#include <sys/time.h>
#include <unistd.h>

#include <algorithm>
#include <bitset>
#include <cctype>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <tuple>
#include <unordered_map>
#include <vector>

using namespace std;

/// Constants
const int MAX_SEQ_LENGTH = 1 << 28;
const int KMER_LENGTH = 20;
const int HASH_TABLE_SIZE = 1 << 28;
const int HASH_TABLE_BIT = 28;
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

struct LineLength {
  int length;
  int repeat_count;
};

vector<char> ref_seq;
vector<char> target_seq;
vector<char> ref_seq_cleaned;
vector<char> target_seq_cleaned;
vector<int> target_seq_encoded;
vector<int> ref_seq_encoded;
vector<int> point(HASH_TABLE_SIZE, -1);
vector<int> loc;
vector<PositionRange> lowercase_ranges;
vector<PositionRange> n_ranges;
vector<SpecialChar> special_chars;
vector<LineLength> line_lengths;
vector<int> line_breaks;
string mismatch_buffer;
string header;
unsigned long timer;
struct timeval timer_start, timer_end;
int base_to_index[256];

void show_help_message(string reason) {
  /**
   * Displays an error message along with usage instructions.
   * Used when the user provides invalid arguments.
   */
  cout << "Error: " << reason << endl;
  cout << "Usage: ./compress_hirgc -r <reference_file_name> -t "
          "<target_file_name>"
       << endl;
}

void init_base_index() {
  /**
   * Initialize base to index mapping for encoding
   * A=0, C=1, G=2, T=3
   */
  memset(base_to_index, -1, sizeof(base_to_index));
  base_to_index['A'] = 0;
  base_to_index['C'] = 1;
  base_to_index['G'] = 2;
  base_to_index['T'] = 3;
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
  init_base_index();
}

void load_sequence(const string& filename, vector<char>& sequence,
                   vector<int>& sequence_encoded, bool is_target) {
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
    int length = 0;

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
          switch (c) {
            case 'A':
            case 'C':
            case 'G':
            case 'T':
              break;
            default:
              continue;
          }
        }
        sequence.push_back(c);
        sequence_encoded.push_back(base_to_index[c]);
        length++;
      }
    }

    if (is_target) {
      if (!line_lengths.empty() &&
          line_lengths[line_lengths.size() - 1].length == length) {
        line_lengths[line_lengths.size() - 1].repeat_count++;
      } else {
        LineLength new_line_length;
        new_line_length.length = length;
        new_line_length.repeat_count = 1;
        line_lengths.push_back(new_line_length);
      }
    }
  }

  file.close();
}

void build_hash_table() {
  /**
   * Build hash table of k-mers from the reference sequence using rolling hash
   */
  if (ref_seq_encoded.size() < KMER_LENGTH) {
    throw runtime_error("Reference sequence too short for k-mer size");
  }

  loc.resize(ref_seq_encoded.size(), -1);

  // Compute the hash of the first k-mer
  uint64_t value = 0;
  for (int k = 0; k < KMER_LENGTH; ++k) {
    value <<= 2;
    value += ref_seq_encoded[k];
  }

  int idx = value & (HASH_TABLE_SIZE - 1);
  loc[0] = point[idx];
  point[idx] = 0;

  const uint64_t mask = (1ULL << (2 * KMER_LENGTH)) - 1;

  // Use rolling hash to compute for next k-mers
  for (int i = 1; i <= ref_seq_encoded.size() - KMER_LENGTH; ++i) {
    value <<= 2;
    value += (ref_seq_encoded[i + KMER_LENGTH - 1]);
    value &= mask;

    idx = value & (HASH_TABLE_SIZE - 1);
    loc[i] = point[idx];
    point[idx] = i;
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

    switch (toupper(c)) {
      case 'A':
      case 'C':
      case 'G':
      case 'T':
        break;
      default:
        special_chars.push_back({i, c});
        valid = false;
    }

    if (valid) {
      target_seq_cleaned.push_back(c);
    }

    if (i == target_seq.size() - 1) {
      if (in_lowercase) {
        lowercase_ranges.push_back({lowercase_start, i + 1 - lowercase_start});
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

void write_metadata(const string& output_filename) {
  /**
   * Write all auxiliary metadata needed for decompression as plain text
   */
  ofstream out(output_filename);
  if (!out) {
    throw runtime_error("Cannot open output file: " + output_filename);
  }

  // Write header
  out << header << "\n\n";

  // Write line lengths
  out << line_lengths.size() * 2;
  for (const auto& line_length : line_lengths) {
    out << " " << line_length.length << " " << line_length.repeat_count;
  }
  out << "\n";

  // Lowercase ranges
  int last_position = 0;
  out << lowercase_ranges.size();
  for (const auto& r : lowercase_ranges) {
    out << " " << r.start - last_position << " " << r.length;
    last_position = r.start + r.length;
  }
  out << "\n";

  // N ranges
  last_position = 0;
  out << n_ranges.size();
  for (const auto& r : n_ranges) {
    out << " " << r.start - last_position << " " << r.length;
    last_position = r.start + r.length;
  }
  out << "\n";

  // Special characters
  set<char> seen;
  vector<char> unique_chars;
  last_position = 0;

  out << special_chars.size();
  for (const auto& sc : special_chars) {
    out << " " << sc.pos - last_position;
    last_position = sc.pos + 1;

    // Collect unique characters
    if (seen.insert(sc.ch).second) {
      unique_chars.push_back(sc.ch);
    }
  }

  if (unique_chars.size()) {
    out << " " << unique_chars.size() << " ";
  }

  // Distances from 'A'
  for (char ch : unique_chars) {
    out << toupper(ch) - 'A' << " ";
  }

  // Unique chars pattern
  for (const auto& sc : special_chars) {
    auto it = find(unique_chars.begin(), unique_chars.end(), sc.ch);
    out << distance(unique_chars.begin(), it);
  }

  out << "\n";

  out.close();
}

void find_longest_match(int tar_pos, int& match_ref_pos, int& match_length) {
  /**
   * Finds the longest match between target and reference sequences starting at
   * tar_pos Uses the k-mer hash table to find candidate positions
   */
  match_ref_pos = -1;
  match_length = 0;

  if (tar_pos + KMER_LENGTH > target_seq_encoded.size()) {
    return;
  }

  uint64_t hash = 0;
  for (int i = 0; i < KMER_LENGTH; ++i) {
    hash <<= 2;
    hash += target_seq_encoded[tar_pos + i];
  }

  // Compute the hash value for the current k-mer in the target sequence
  int idx = hash & (HASH_TABLE_SIZE - 1);

  // Find the longest match
  int ht_half_size =
      HASH_TABLE_BIT >> 1;  // hash smo napravili od 20 najnizih bitova k-mera
  int missing_bases_count = KMER_LENGTH - ht_half_size;
  for (int k = point[idx]; k != -1; k = loc[k]) {
    int current_length = 0;

    while (current_length < missing_bases_count &&
           ref_seq_encoded[current_length + k] ==
               target_seq_encoded[tar_pos + current_length])
      current_length++;

    if (current_length == missing_bases_count) {  // našli smo match
      current_length = KMER_LENGTH;
      int max_possible = min((int)ref_seq_encoded.size() - k,
                             (int)target_seq_encoded.size() - tar_pos);
      while (current_length < max_possible &&
             ref_seq_encoded[k + current_length] ==
                 target_seq_encoded[tar_pos + current_length]) {
        current_length++;
      }
    } else {
      current_length = 0;
    }

    if (current_length > match_length) {
      match_length = current_length;
      match_ref_pos = k;
    }
  }
}

void compress_sequences() {
  vector<Match> matches;
  vector<char> mismatches;
  vector<int> encoded_mismatches;
  matches.reserve(target_seq_encoded.size() / 100 + 1000);
  mismatches.reserve(10000);

  int tar_pos = 0;
  int prev_ref_pos = 0;
  int prev_tar_pos = 0;
  int total_matched = 0;
  int total_mismatched = 0;

  string compressed_file = "output.txt";

  ofstream out(compressed_file, ios::app);
  if (!out) {
    throw runtime_error("Cannot open output file: " + compressed_file);
  }

  write_metadata(compressed_file);

  // while (tar_pos < target_seq_encoded.size()) {
  while (tar_pos < target_seq_encoded.size()) {
    int match_ref_pos, match_length;
    find_longest_match(tar_pos, match_ref_pos, match_length);

    if (match_length >= KMER_LENGTH) {
      if (!mismatches.empty()) {
        encode_sequence(mismatches, encoded_mismatches);
        for (size_t i = 0; i < encoded_mismatches.size(); ++i) {
          out << encoded_mismatches[i];
        }
        out << '\n';

        mismatches.clear();
        encoded_mismatches.clear();
      }
      int delta_ref = match_ref_pos - prev_ref_pos;
      int delta_tar = tar_pos - prev_tar_pos;
      matches.push_back({delta_ref, delta_tar, match_length});
      total_matched += match_length;
      prev_ref_pos = match_ref_pos + match_length;
      prev_tar_pos = tar_pos + match_length;
      tar_pos += match_length;
      out << delta_ref << " " << match_length - KMER_LENGTH << '\n';
    } else {
      mismatches.push_back(target_seq[tar_pos]);
      tar_pos++;
    }
  }

  if (!mismatches.empty()) {
    encode_sequence(mismatches, encoded_mismatches);
    for (size_t i = 0; i < encoded_mismatches.size(); ++i) {
      out << encoded_mismatches[i];
    }
  }

  out.close();

  cout << "Total matched bases: " << total_matched << endl;
  cout << "Total mismatched bases: " << total_mismatched << endl;
  cout << "Compression ratio: "
       << (100.0 * (total_matched) / (total_matched + total_mismatched)) << "%"
       << endl;
  cout << "Compressed data written to " << compressed_file << endl;
}

void compress_to_7z(const string& input_file, const string& archive_name) {
  string command = "7z a -t7z " + archive_name + " " + input_file;
  int result = system(command.c_str());
  if (result != 0) {
    throw runtime_error("7z compression failed");
  }
}

void cleanup() {
  /**
   * Cleans up and releases all allocated memory
   */
  ref_seq.clear();
  target_seq.clear();
  lowercase_ranges.clear();
  n_ranges.clear();
  special_chars.clear();
  line_breaks.clear();
  mismatch_buffer.clear();
}

void print_memory_usage() {
  /**
   * Print the current memory usage of the program
   */
  std::ifstream status_file("/proc/self/status");
  std::string line;
  while (std::getline(status_file, line)) {
    if (line.substr(0, 6) == "VmRSS:") {
      std::cout << "Memory usage: " << line.substr(6) << std::endl;
      break;
    }
  }
}

int main(int argc, char* argv[]) {
  /**
   * Main function for compressing files using HIRGC algorithm.
   */

  // Start tracking time taken for compression
  gettimeofday(&timer_start, nullptr);

  // Check if passed arguments are valid
  if (argc != 5) {
    show_help_message("Invalid number of arguments.");
    return 1;
  }

  if (strcmp(argv[1], "-r") != 0 || strcmp(argv[3], "-t") != 0) {
    show_help_message("Invalid arguments.");
    return 1;
  }

  // Assign the reference and target file paths from the command line arguments
  InputFileNames input_file_names;

  input_file_names.reference_file = argv[2];
  input_file_names.target_file = argv[4];

  initialize_structures();

  try {
    load_sequence(input_file_names.reference_file, ref_seq, ref_seq_encoded,
                  false);
    load_sequence(input_file_names.target_file, target_seq, target_seq_encoded,
                  true);

    process_target_sequence();

    // encode_sequence(target_seq, target_seq_encoded);
    // encode_sequence(ref_seq, ref_seq_encoded);

    build_hash_table();
    compress_sequences();

    compress_to_7z("output.txt", "compressed.7z");

    cout << "Compression completed successfully." << endl;
    print_memory_usage();
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