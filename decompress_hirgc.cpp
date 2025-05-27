#include <cstring>
#include <fstream>
#include <iostream>
#include <vector>

using namespace std;

const int MAX_SEQ_LENGTH = 1 << 28;

/// Struct for reference and target file names
struct InputFileNames {
  string reference_file;
  string compressed_target_file;
};

vector<char> ref_seq;
vector<char> target_seq;

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

void load_and_clean_reference(const string& filename, vector<char>& sequence) {
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
    };

    if (line.empty()) {  // Skip empty lines
      continue;
    }

    for (char c : line) {
      c = toupper(c);  // Convert to uppercase
      if (c == 'A' || c == 'C' || c == 'G' || c == 'T') {
        sequence.push_back(c);
        cout << c;
      }
    }
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

  cout << "Decompressing..." << endl;

  cleanup();

  return 0;
}