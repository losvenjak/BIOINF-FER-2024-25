#include <sys/time.h>

#include <cctype>
#include <cstring>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <unordered_map>
#include <vector>

using namespace std;

/// Constants
const int MAX_SEQ_LENGTH = 1 << 28;
const int KMER_LENGTH = 20;
const int HASH_TABLE_SIZE = 1 << 20;
const int INITIAL_BUFFER_SIZE = 1024;

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

vector<char> ref_seq;
vector<char> target_seq;
vector<int> target_seq_encoded;
vector<int> ref_seq_encoded;
unordered_map<string, vector<int>> kmer_hash_table;
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
   * Build hash table of k-mers from the reference sequence
   */
  if (ref_seq.size() < KMER_LENGTH) {
    throw runtime_error("Reference sequence too short for k-mer size");
  }

  for (int i = 0; i <= ref_seq.size() - KMER_LENGTH; ++i) {
    string kmer(ref_seq.begin() + i, ref_seq.begin() + i + KMER_LENGTH);
    kmer_hash_table[kmer].push_back(i);
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

  for (int i = 0; i < target_seq.size(); ++i) {
    char c = target_seq[i];

    if (islower(c)) {
      if (!in_lowercase) {
        lowercase_start = i;
        in_lowercase = true;
      }
      target_seq[i] = toupper(c);
    } else if (in_lowercase) {
      lowercase_ranges.push_back({lowercase_start, i - lowercase_start});
      in_lowercase = false;
    }

    if (toupper(c) == 'N') {
      if (!in_n_region) {
        n_start = i;
        in_n_region = true;
      }
      target_seq[i] = '_';  // '_' marks chars to be removed
    } else if (in_n_region) {
      n_ranges.push_back({n_start, i - n_start});
      in_n_region = false;
    }

    if (string("ACGTN").find(toupper(c)) == string::npos) {
      special_chars.push_back({i, c});
      target_seq[i] = '_';  // Any char that is not ACGTN marked to be removed
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

  // Remove all '_' characters from the sequence
  vector<char> cleaned_seq;
  for (char c : target_seq) {
    if (c != '_') {
      cleaned_seq.push_back(c);
    }
  }
  target_seq = move(cleaned_seq);
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

void compress_sequences() {
  //
  // 1. Use kmer_hash_table to find matches
  // 2. Call handle_mismatch() for non-matching regions
  // 3. Generate compressed output
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

    build_hash_table();
    process_target_sequence();

    encode_sequence(target_seq, target_seq_encoded);
    encode_sequence(ref_seq, ref_seq_encoded);

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