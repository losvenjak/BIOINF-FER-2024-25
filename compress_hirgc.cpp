#include <sys/time.h>

#include <algorithm>
#include <cctype>
#include <cstring>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <unordered_map>
#include <vector>
#include <bitset>
#include <map>

using namespace std;

/// Constants
const int MAX_SEQ_LENGTH = 1 << 28;
const int KMER_LENGTH = 20;
const int HASH_TABLE_SIZE = 1 << 20;
const int INITIAL_BUFFER_SIZE = 1024;
const int BITS_PER_BYTE = 8;
const int MAX_DELTA_BITS = 32;

/// Struct for reference and target file names
struct InputFileNames {
  string reference_file;
  string target_file;
};

struct PositionRange {
  int start;
  int length;
};

struct SpecialChar {
  int pos;
  char ch;
};

struct Match {
  int ref_pos;
  int tar_pos;
  int length;
};

vector<char> ref_seq;
vector<char> target_seq;
vector<char> ref_seq_cleaned;
vector<char> target_seq_cleaned;
vector<int> target_seq_encoded;
vector<int> ref_seq_encoded;
unordered_map<uint64_t, vector<int>> kmer_hash_table;
vector<PositionRange> lowercase_ranges;
vector<PositionRange> n_ranges;
vector<SpecialChar> special_chars;
vector<int> line_breaks;
string mismatch_buffer;
string header;
unsigned long timer;
struct timeval timer_start, timer_end;

void showHelpMessage(string reason) {
  /**
   * Displays an error message along with usage instructions.
   * Used when the user provides invalid arguments.
   */
  cout << "Error: " << reason << endl;
  cout << "Usage: ./compress_hirgc -r <reference_file_name> -t "
          "<target_file_name>"
       << endl;
}

void initialize_structures() {
  /**
   * Initialize memory structures for genome sequences and buffers
   */
  ref_seq.reserve(MAX_SEQ_LENGTH);
  target_seq.reserve(MAX_SEQ_LENGTH);
  mismatch_buffer.reserve(INITIAL_BUFFER_SIZE);
  target_seq_encoded.reserve(MAX_SEQ_LENGTH);
  ref_seq_encoded.reserve(MAX_SEQ_LENGTH);
}

void load_sequence(const string& filename, vector<char>& sequence,
                   bool is_target) {
  /**
   * Load genome sequence from FASTA file into memory
   * Store line breaks and header info for target sequence
   */
  ifstream file(filename);
  if (!file) {
    throw runtime_error("Cannot open file: " + filename);
  }

  string line;
  while (getline(file, line)) {
    if (line[0] == '>') {  // Skip header line for reference, store for target
      if (is_target) {
        header = line;
      }
      continue;
    };

    if (line.empty()) {  // Skip empty lines for reference, store for target
      if (is_target) {
        line_breaks.push_back(sequence.size());
      }
      continue;
    }

    for (char c : line) {
      if (c != '\n') {
        if (!is_target) {
          c = toupper(c);
          if (string("ACGT").find(c) == string::npos) {
            continue;
          }
        }
        sequence.push_back(c);
      }
    }

    if (is_target)
      line_breaks.push_back(
          sequence.size());  // Store end of line line breaks for target
  }
}

void build_hash_table() {
  /**
   * Build hash table of k-mers from the reference sequence using rolling hash
   */

  if (ref_seq.size() < KMER_LENGTH) {
    throw runtime_error("Reference sequence too short for k-mer size");
  }

  // Compute first k-mer value
  uint64_t value = 0;
  for (int k = 0; k < KMER_LENGTH; k++) {
    value <<= 2;
    value |= (ref_seq_encoded[k]);
  }
  kmer_hash_table[value].push_back(0);

  // Compute a rolling integer for the rest of the sequence
  const uint64_t mask =
      (1ULL << (2 * KMER_LENGTH)) - 1;  // Mask for removing oldest two bits

  for (int k = KMER_LENGTH; k < ref_seq_encoded.size(); k++) {
    value <<= 2;
    value |= (ref_seq_encoded[k]);
    value &= mask;
    kmer_hash_table[value].push_back(k - KMER_LENGTH + 1);
  }
}

void process_target_sequence() {
  /**
   * Processes the target genome sequence to identify special features :
   * - lowercase regions (often indicate repeats)
   * - n regions (unknown bases)
   * - other special characters
   */
  bool in_lowercase = false;
  bool in_n_region = false;
  int lowercase_start = 0;
  int n_start = 0;
  bool valid;

  for (int i = 0; i < target_seq.size(); ++i) {
    char c = target_seq[i];
    valid = true;

    if (islower(c)) {
      if (!in_lowercase) {
        lowercase_start = i;
        in_lowercase = true;
      }
      c = toupper(c);
    } else if (in_lowercase) {
      lowercase_ranges.push_back({lowercase_start, i - lowercase_start});
      in_lowercase = false;
    }

    if (c == 'N') {
      if (!in_n_region) {
        n_start = i;
        in_n_region = true;
      }
      valid = false;
    } else if (in_n_region) {
      n_ranges.push_back({n_start, i - n_start});
      in_n_region = false;
    }

    if (string("ACGTN").find(toupper(c)) == string::npos) {
      special_chars.push_back({i, c});
      valid = false;
    }

    if (valid) {
      target_seq_cleaned.push_back(c);
    }

    if (i == target_seq.size() - 1) {
      if (in_lowercase) {
        lowercase_ranges.push_back({lowercase_start, i - lowercase_start});
        in_lowercase = false;
      }
      if (in_n_region) {
        n_ranges.push_back({n_start, i - n_start});
        in_n_region = false;
      }
    }
  }
}

void encode_sequence(vector<char>& sequence, vector<int>& encoded_sequence) {
  /**
   * Encode cleaned sequence into a numerical format for compression
   * A=0, C=1, G=2, T=3
   */
  for (int i = 0; i < sequence.size(); ++i) {
    char c = sequence[i];
    switch (c) {
      case 'A':
        encoded_sequence.push_back(0);
        break;
      case 'C':
        encoded_sequence.push_back(1);
        break;
      case 'G':
        encoded_sequence.push_back(2);
        break;
      case 'T':
        encoded_sequence.push_back(3);
        break;
      default:
        continue;
    }
  }
}

void handle_mismatch(int pos, int length) {
  /**
   * Handle mismatched regions between reference and target genomes
   */
  if (pos + length > target_seq.size()) {
    throw out_of_range("Mismatch position out of bounds");
  }
  mismatch_buffer.append(target_seq.begin() + pos,
                         target_seq.begin() + pos + length);
}

template <typename T>
void write_run_length_encoded(ostream& out, const vector<T>& values) {
  /**
   * Compress and write data using run-length encoding (RLE)
   * group consecutive identical values into (value, count) pairs
   */
  if (values.empty()) {
    out << "0";
    return;
  }

  T current = values[0];
  int count = 1;
    
  for (size_t i = 1; i < values.size(); ++i) {
    if (values[i] == current) {
      count++;
    } else {
      out << current << " " << count << " ";
      current = values[i];
      count = 1;
    }
  }
  out << current << " " << count;
}

void write_metadata(const string& output_filename) {
  /**
   * Write all auxiliary metadata needed for decompression
   */
  ofstream out(output_filename);
  if (!out) {
    throw runtime_error("Cannot open output file: " + output_filename);
  }

  out << "HEADER:" << header << "\n\n";
  out << "LINE_BREAKS:";
  
  for (size_t i = 0; i < line_breaks.size(); ++i) {
    if (i != 0) out << ",";
    out << line_breaks[i];
  }
  
  out << "\n\n";
  out << "LOWERCASE_RANGES:" << lowercase_ranges.size() << "\n";
  
  for (const auto& range : lowercase_ranges) {
    out << range.start << " " << range.length << "\n";
  }
  
  out << "\n";
  out << "N_RANGES:" << n_ranges.size() << "\n";
  
  for (const auto& range : n_ranges) {
    out << range.start << " " << range.length << "\n";
  }
  
  out << "\n";
  map<char, int> char_to_code;
  vector<char> unique_chars;
  out << "SPECIAL_CHARS:" << special_chars.size() << "\n";
  
  for (const auto& sc : special_chars) {
    if (char_to_code.find(sc.ch) == char_to_code.end()) {
      char_to_code[sc.ch] = unique_chars.size();
      unique_chars.push_back(sc.ch);
    }
    out << sc.pos << " " << sc.ch << "\n";
  }

  out.close();
}

void find_longest_match(int tar_pos, int& match_ref_pos, int& match_length) {
    /**
     * Finds the longest match between target and reference sequences starting at tar_pos
     * Uses the k-mer hash table to find candidate positions
     */
  match_ref_pos = -1;
  match_length = 0;
  
  if (tar_pos + KMER_LENGTH > target_seq_encoded.size()) {
    return;
  }

  uint64_t target_hash = 0;
  for (int k = 0; k < KMER_LENGTH; k++) {
    target_hash <<= 2;
    target_hash |= target_seq_encoded[tar_pos + k];
  }

  auto it = kmer_hash_table.find(target_hash);
  if (it == kmer_hash_table.end()) {
    return; 
  }

  for (int ref_pos : it->second) {
    int max_possible_length = min(
        (int)ref_seq_encoded.size() - ref_pos,
        (int)target_seq_encoded.size() - tar_pos
    );
        
    int current_length = KMER_LENGTH;
    while (current_length < max_possible_length && 
           ref_seq_encoded[ref_pos + current_length] == target_seq_encoded[tar_pos + current_length]) {
      current_length++;
    }

    if (current_length > match_length) {
      match_length = current_length;
      match_ref_pos = ref_pos;
    }
  }
}

void compress_sequences() {
  // 1. Use kmer_hash_table to find matches
  // 2. Call handle_mismatch() for non-matching regions
  // 3. Generate compressed output
  vector<Match> matches;
  vector<char> mismatches;
  matches.reserve(target_seq_encoded.size() / 100 + 1000);
  mismatches.reserve(10000);
    
  int tar_pos = 0;
  int prev_ref_pos = 0;
  int prev_tar_pos = 0;
  int total_matched = 0;
  int total_mismatched = 0;
    
  while (tar_pos < target_seq_encoded.size()) {
    int match_ref_pos, match_length;
    find_longest_match(tar_pos, match_ref_pos, match_length);

    if (match_length >= KMER_LENGTH) {
      if (!mismatches.empty()) {
        handle_mismatch(prev_tar_pos, mismatches.size());
        total_mismatched += mismatches.size();
        mismatches.clear();
      }
      int delta_ref = match_ref_pos - prev_ref_pos;
      int delta_tar = tar_pos - prev_tar_pos;
      matches.push_back({delta_ref, delta_tar, match_length});
      total_matched += match_length;
      prev_ref_pos = match_ref_pos + match_length;
      prev_tar_pos = tar_pos + match_length;
      tar_pos += match_length;
    } else {
      mismatches.push_back(target_seq[tar_pos]);
      tar_pos++;
    }
  }
    
  if (!mismatches.empty()) {
    handle_mismatch(prev_tar_pos, mismatches.size());
    total_mismatched += mismatches.size();
  }

  string compressed_file = "output_compressed.txt";
  ofstream out(compressed_file);
  if (!out) {
    throw runtime_error("Cannot open output file: " + compressed_file);
  }
    
  out << "MATCHES:" << matches.size() << "\n";
  for (const auto& match : matches) {
    out << match.ref_pos << " " << match.tar_pos << " " << match.length << "\n";
  }
    
  out << "\nMISMATCHES:" << mismatch_buffer.size() << "\n";
  out << mismatch_buffer << "\n";
    
  out.close();

  string metadata_file = "output_metadata.txt";
  write_metadata(metadata_file);
    
  cout << "Total matched bases: " << total_matched << endl;
  cout << "Total mismatched bases: " << total_mismatched << endl;
  cout << "Compression ratio: " 
       << (100.0 * total_matched / (total_matched + total_mismatched)) << "%" << endl;
  cout << "Compressed data written to " << compressed_file << endl;
  cout << "Metadata written to " << metadata_file << endl;
}


void cleanup() {
  /**
   * Cleans up and releases all allocated memory
   */
  ref_seq.clear();
  target_seq.clear();
  kmer_hash_table.clear();
  lowercase_ranges.clear();
  n_ranges.clear();
  special_chars.clear();
  line_breaks.clear();
  mismatch_buffer.clear();
}

int main(int argc, char* argv[]) {
  /**
   * Main function for compressing files using HIRGC algorithm.
   */

  // Start tracking time taken for compression
  gettimeofday(&timer_start, nullptr);

  // Check if passed arguments are valid
  if (argc != 5) {
    showHelpMessage("Invalid number of arguments.");
    return 1;
  }

  if (strcmp(argv[1], "-r") != 0 || strcmp(argv[3], "-t") != 0) {
    showHelpMessage("Invalid arguments.");
    return 1;
  }

  // Assign the reference and target file paths from the command line arguments
  InputFileNames input_file_names;

  input_file_names.reference_file = argv[2];
  input_file_names.target_file = argv[4];

  initialize_structures();

  try {
    load_sequence(input_file_names.reference_file, ref_seq, false);
    load_sequence(input_file_names.target_file, target_seq, true);

    process_target_sequence();

    encode_sequence(target_seq, target_seq_encoded);
    encode_sequence(ref_seq, ref_seq_encoded);

    build_hash_table();
    compress_sequences();

    cout << "Compression completed successfully." << endl;
  } catch (const exception& e) {
    cerr << "Error: " << e.what() << endl;
    cleanup();
    return 1;
  }

  cleanup();

  // Calculate and print the total time taken for compression
  gettimeofday(&timer_end, nullptr);

  timer = 1000000 * (timer_end.tv_sec - timer_start.tv_sec) +
          timer_end.tv_usec - timer_start.tv_usec;
  printf("Total compresion time : %lf ms; \n", timer / 1000.0);

  return 0;
}