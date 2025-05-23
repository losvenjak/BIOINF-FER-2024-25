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
    out.write("\0", 1); // Empty marker
    return;
  }

  T current = values[0];
  int count = 1;
    
  for (size_t i = 1; i < values.size(); ++i) {
    if (values[i] == current) {
      count++;
    } else {
        out.write(reinterpret_cast<const char*>(&current), sizeof(current));
        out.write(reinterpret_cast<const char*>(&count), sizeof(count));
        current = values[i];
        count = 1;
    }
  }
  out.write(reinterpret_cast<const char*>(&current), sizeof(current));
  out.write(reinterpret_cast<const char*>(&count), sizeof(count));
}

void write_metadata(const string& output_filename) {
  /**
   * Write all auxiliary metadata needed for decompression
   */
  ofstream out(output_filename, ios::binary);
  if (!out) {
    throw runtime_error("Cannot open output file: " + output_filename);
  }

  uint32_t header_length = header.size();
  out.write(reinterpret_cast<const char*>(&header_length), sizeof(header_length));
  out.write(header.c_str(), header_length);

  write_run_length_encoded(out, line_breaks);

  uint32_t num_lowercase = lowercase_ranges.size();
  out.write(reinterpret_cast<const char*>(&num_lowercase), sizeof(num_lowercase));
  for (const auto& range : lowercase_ranges) {
    out.write(reinterpret_cast<const char*>(&range.start), sizeof(range.start));
    out.write(reinterpret_cast<const char*>(&range.length), sizeof(range.length));
  }

  uint32_t num_n_ranges = n_ranges.size();
  out.write(reinterpret_cast<const char*>(&num_n_ranges), sizeof(num_n_ranges));
  for (const auto& range : n_ranges) {
    out.write(reinterpret_cast<const char*>(&range.start), sizeof(range.start));
    out.write(reinterpret_cast<const char*>(&range.length), sizeof(range.length));
  }

  map<char, uint8_t> char_to_code;
  vector<char> unique_chars;
  for (const auto& sc : special_chars) {
    if (char_to_code.find(sc.ch) == char_to_code.end()) {
      char_to_code[sc.ch] = unique_chars.size();
      unique_chars.push_back(sc.ch);
    }
  }

  uint8_t num_unique_chars = unique_chars.size();
  out.write(reinterpret_cast<const char*>(&num_unique_chars), sizeof(num_unique_chars));
  out.write(unique_chars.data(), num_unique_chars);

  // delta encoded positions and character codes
  if (!special_chars.empty()) {
    int prev_pos = 0;
    bitset<MAX_DELTA_BITS> delta_bits;
        
    for (const auto& sc : special_chars) {
      int delta = sc.pos - prev_pos;
      prev_pos = sc.pos;
            
      do {
        uint8_t byte = delta & 0x7F;
        delta >>= 7;
        if (delta != 0) byte |= 0x80;
        out.write(reinterpret_cast<const char*>(&byte), sizeof(byte));
      } while (delta != 0);
            
      uint8_t code = char_to_code[sc.ch];
      out.write(reinterpret_cast<const char*>(&code), sizeof(code));
    }
  }

  out.close();
}

void compress_sequences() {
  //
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
    //find_longest_match(tar_pos, match_ref_pos, match_length);

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

  string compressed_file = "output.compressed";
  ofstream out(compressed_file, ios::binary);
  if (!out) {
    throw runtime_error("Cannot open output file: " + compressed_file);
  }
  
  vector<uint8_t> match_buffer;
  match_buffer.reserve(matches.size() * 12);  
  
  uint32_t num_matches = matches.size();
  const uint8_t* num_matches_bytes = reinterpret_cast<const uint8_t*>(&num_matches);
  match_buffer.insert(match_buffer.end(), num_matches_bytes, num_matches_bytes + sizeof(num_matches));
  
  for (const auto& match : matches) {
    auto encode_varint = [&match_buffer](int value) {
      while (value > 0x7F) {
        match_buffer.push_back((value & 0x7F) | 0x80);
        value >>= 7;
      }
      match_buffer.push_back(value & 0x7F);
    };

    encode_varint(match.ref_pos);
    encode_varint(match.tar_pos);
    encode_varint(match.length);
  }
  
  out.write(reinterpret_cast<const char*>(match_buffer.data()), match_buffer.size());

  uint32_t mismatch_length = mismatch_buffer.size();
  out.write(reinterpret_cast<const char*>(&mismatch_length), sizeof(mismatch_length));
  out.write(mismatch_buffer.data(), mismatch_length);

  out.close();

  string metadata_file = "output.meta";
  write_metadata(metadata_file);
  cout << "Metadata written to " << metadata_file << endl;
  cout << "Compressing..." << endl;
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