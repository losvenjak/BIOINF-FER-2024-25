#include <cstring>
#include <fstream>
#include <iostream>
#include <vector>

using namespace std;

const int MAX_SEQ_LENGTH = 1 << 28;
const int KMER_LENGTH = 20;
const vector<char> decode_into_base = {'A', 'C', 'G', 'T'};

struct InputFileNames {
  string reference_file;
  string compressed_target_file;
};

struct Mismatch {
  vector<int> values;
  int continue_for;
};

vector<char> ref_seq;
vector<char> target_seq;
string header;
vector<int> line_lenghts;
vector<int> lower_case_ranges;
vector<int> n_ranges;
vector<int> special_chars;
vector<int> special_chars_order;
vector<Mismatch> mismatch_data;
int mismatch_offset = 0;

void show_help_message(string reason) {
  /**
   * Displays an error message along with usage instructions.
   * Used when the user provides invalid arguments.
   */
  cout << "Error: " << reason << endl;
  cout << "Usage: ./decompress_hirgc -r <reference_file_name> -t "
          "<target_file_name>"
       << endl;
}

void initialize_structures() {
  /**
   * Initializes all the data structures used in the decompression process.
   */
  ref_seq.reserve(MAX_SEQ_LENGTH);
  target_seq.reserve(MAX_SEQ_LENGTH);
}

void load_and_clean_reference(const string& filename, vector<char>& ref_seq) {
  /**
   * Load and clean reference genome sequence
   * removes any non-ACGT characters, converts to uppercase
   */

  ifstream file(filename);
  if (!file) {
    throw runtime_error("Cannot open file: " + filename);
  }

  string line;
  while (getline(file, line)) {
    if (line[0] == '>') {  // Skip header line
      continue;
    }

    if (line.empty()) {  // Skip empty lines
      continue;
    }

    for (char c : line) {
      c = toupper(c);  // Convert to uppercase
      if (c == 'A' || c == 'C' || c == 'G' || c == 'T') {
        ref_seq.push_back(c);
      }
    }
  }
}

void load_metadata(const string& filename) {
  /**
   * Load metadata from the compressed target file
   * Reads the header, line lengths, lowercase ranges,
   * N ranges, special characters and mismatch offset
   */

  ifstream file(filename, ios::binary);
  if (!file) {
    throw runtime_error("Cannot open file: " + filename);
  }

  string temp;
  int curr;

  // Read the header line
  getline(file, temp);
  header = temp;

  getline(file, temp);  // Skip empty line after header

  // Read the line lengths
  getline(file, temp);
  curr = temp[0] - '0';

  for (size_t i = 1; i < temp.size(); ++i) {
    char c = temp[i];

    if (c == ' ') {
      line_lenghts.push_back(curr);
      curr = 0;
    } else {
      curr = curr * 10 + (c - '0');
    }
  }
  line_lenghts.push_back(curr);

  // Read the lowercase ranges
  getline(file, temp);
  curr = temp[0] - '0';

  for (size_t i = 1; i < temp.size(); ++i) {
    char c = temp[i];

    if (c == ' ') {
      lower_case_ranges.push_back(curr);
      curr = 0;
    } else {
      curr = curr * 10 + (c - '0');
    }
  }
  lower_case_ranges.push_back(curr);

  // Read the N ranges
  getline(file, temp);
  curr = temp[0] - '0';

  for (size_t i = 1; i < temp.size(); ++i) {
    char c = temp[i];

    if (c == ' ') {
      n_ranges.push_back(curr);
      curr = 0;
    } else {
      curr = curr * 10 + (c - '0');
    }
  }
  n_ranges.push_back(curr);

  // Read the special characters
  getline(file, temp);
  int order_list_start = temp.rfind(' ') + 1;
  curr = temp[0] - '0';

  for (size_t i = 1; i < order_list_start; ++i) {
    char c = temp[i];

    if (c == ' ') {
      special_chars.push_back(curr);
      curr = 0;
    } else {
      curr = curr * 10 + (c - '0');
    }
  }
  special_chars.push_back(curr);

  for (int i = order_list_start; i < temp.size(); ++i) {
    char c = temp[i];
    special_chars_order.push_back(c - '0');
  }

  // Check mismatched offset
  getline(file, temp);
  if (temp[0] == '0') {
    int pos = 0;
    for (int i = 2; i < temp.size(); i++) {
      pos = pos * 10 + (temp[i] - '0');
    }
    mismatch_offset = pos;
  }
}

void load_mismatch_data(const string& filename) {
  /**
   * Load mismatch data from the compressed target file
   * Stores the mismatched values
   */
  ifstream file(filename, ios::binary);
  if (!file) {
    throw runtime_error("Cannot open file: " + filename);
  }

  // Skip the header and metadata lines
  string temp;
  for (int i = 0; i < 6; ++i) {
    getline(file, temp);
  }
  // Skip the mismatch offset line if present
  if (mismatch_offset != 0) {
    getline(file, temp);
  }

  // Read mismatch data
  vector<int> values;
  int continue_for;

  while (getline(file, temp)) {
    Mismatch mismatch;
    for (char c : temp) {
      mismatch.values.push_back(c - '0');
    }

    getline(file, temp);
    int pos = 0;
    int index = temp.find(' ') + 1;
    for (int i = index; i < temp.size(); i++) {
      pos = pos * 10 + (temp[i] - '0');
    }

    mismatch.continue_for = pos;

    mismatch_data.push_back(mismatch);
  }
}

void decompress_target_sequence(vector<char>& target_seq) {
  /**
   * Decompress the target sequence from the compressed file
   * Reconstructs the target sequence using the reference sequence
   */
  if (mismatch_offset != 0) {
    for (int i = 0; i < mismatch_offset + KMER_LENGTH; i++) {
      target_seq.push_back(ref_seq[i]);
    }
  }

  for (Mismatch& mismatch : mismatch_data) {
    for (int i = 0; i < mismatch.values.size(); i++) {
      target_seq.push_back(decode_into_base[mismatch.values[i]]);
    }
    for (int i = 0; i < mismatch.continue_for + KMER_LENGTH; i++) {
      target_seq.push_back(ref_seq[target_seq.size()]);
    }
  }
}

void add_special_characters(vector<char>& target_seq) {
  /**
   * Adds special characters to the target sequence
   * based on the special character ranges
   */
  int special_char_num = special_chars[0];
  int unique_special_chars_num = special_chars[special_char_num + 1];
  vector<int> special_chars_positions;
  vector<char> unique_special_chars_decoded;

  if (special_char_num == 0) {
    return;
  }

  // Store positions of special characters
  int pos = 0;
  for (int i = 1; i < special_char_num + 1; i++) {
    pos += special_chars[i];
    special_chars_positions.push_back(pos);
    pos += 1;
  }

  // Decode special characters
  for (int i = special_char_num + 2;
       i < special_char_num + unique_special_chars_num + 2; i++) {
    unique_special_chars_decoded.push_back(special_chars[i] + 'A');
  }

  // Insert special characters into the target sequence
  for (int i = 0; i < special_chars_order.size(); i++) {
    target_seq.insert(target_seq.begin() + special_chars_positions[i],
                      unique_special_chars_decoded[special_chars_order[i]]);
  }
}

void cleanup() {
  /**
   * Cleans up and releases all allocated memory
   */
  ref_seq.clear();
  target_seq.clear();
}

int main(int argc, char* argv[]) {
  /**
   * Main function for compressing files using HIRGC algorithm.
   */

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
  input_file_names.compressed_target_file = argv[4];

  initialize_structures();

  load_and_clean_reference(input_file_names.reference_file, ref_seq);

  load_metadata(input_file_names.compressed_target_file);
  load_mismatch_data(input_file_names.compressed_target_file);

  cout << "Reference Sequence: " << endl;
  for (char c : ref_seq) {
    cout << c;
  }
  cout << endl;

  decompress_target_sequence(target_seq);

  add_special_characters(target_seq);

  cout << "Target Sequence: " << endl;
  for (char c : target_seq) {
    cout << c;
  }
  cout << endl;

  cout << "Decompressing..." << endl;

  cleanup();

  return 0;
}