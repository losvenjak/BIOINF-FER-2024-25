#include <cstring>
#include <iostream>

using namespace std;

/// Struct for reference and target file names
struct InputFileNames {
  string reference_file;
  string target_file;
};

void showHelpMessage(string reason) {
  /**
   * Displays an error message along with usage instructions.
   * Used when the user provides invalid arguments.
   */
  cout << "Error: " << reason << endl;
  cout << "Usage: ./decompress_hirgc -r <reference_file_name> -t "
          "<target_file_name>"
       << endl;
}

int main(int argc, char* argv[]) {
  /**
   * Main function for compressing files using HIRGC algorithm.
   */

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

  cout << "Reference: " << input_file_names.reference_file << endl;
  cout << "Target: " << input_file_names.target_file << endl;

  cout << "Decompressing..." << endl;

  return 0;
}