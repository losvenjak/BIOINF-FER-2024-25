#include <sys/time.h>

#include <cstring>
#include <fstream>
#include <iostream>
#include <vector>
#include <unistd.h>

using namespace std;

const int MAX_SEQ_LENGTH = 1 << 28;
const int KMER_LENGTH = 20;
const vector<char> decode_into_base = {'A', 'C', 'G', 'T'};

struct InputFileNames {
  string reference_file;
  string compressed_target_file;
};

vector<char> ref_seq;
vector<char> target_seq;
string header;
vector<int> line_lenghts;
vector<int> lower_case_ranges;
vector<int> n_ranges;
vector<int> special_chars;
vector<int> special_chars_order;
int ref_seq_position = 0;
unsigned long timer;
struct timeval timer_start, timer_end;

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
   * N ranges, special characters and initial reference start postion
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
}

void decompress_target_sequence(const string& filename,
                                vector<char>& target_seq) {
  /**
   * Decompress the target sequence using the compressed file info
   * Reconstructs the target sequence using the reference sequence
   * and the mismatch data
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

  while (getline(file, temp)) {
    // If there's no space in the line it's a mismatch
    // the values are copied from compressed file
    if (temp.find(' ') == string::npos) {
      for (int i = 0; i < temp.size(); i++) {
        char c = temp[i];
        target_seq.push_back(decode_into_base[c - '0']);
      }
    } else {  // If there's a space it's a match, read from reference

      // Read the position and length of the match
      int position = 0;
      int counter = 0;
      int curr = 0;
      int mult = 1;

      for (int i = 0; i < temp.size(); i++) {
        if (temp[i] == '-') {
          mult = -1;
          continue;
        }

        if (temp[i] == ' ') {
          position = curr;
          curr = 0;
          continue;
        }
        curr = curr * 10 + mult * (temp[i] - '0');
      }

      counter = curr * mult;

      // Write matched bases from the reference sequence
      ref_seq_position += position;
      for (int i = 0; i < counter + KMER_LENGTH; i++) {
        target_seq.push_back(ref_seq[ref_seq_position]);
        ref_seq_position++;
      }
    }
  }
}

void add_special_characters(vector<char>& target_seq) {
  /**
   * Add special characters to the target sequence
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

void add_n_ranges(vector<char>& target_seq) {
  /**
   * Add N ranges to the target sequence
   * based on the N ranges defined in the metadata
   */
  int n_ranges_num = n_ranges[0];

  if (n_ranges_num == 0) {
    return;
  }

  // For each position-length pair in n_ranges insert 'N's
  int prev_pos = 0;
  for (int i = 1; i < n_ranges_num * 2 + 1; i += 2) {
    int start = n_ranges[i];
    int length = n_ranges[i + 1];

    for (int j = 0; j < length; j++) {
      target_seq.insert(target_seq.begin() + j + start + prev_pos, 'N');
    }
    prev_pos += start + length;
  }
}

void add_lowercase_ranges(vector<char>& target_seq) {
  /**
   * Add lowercase ranges to the target sequence
   * based on the lowercase ranges defined in the metadata
   */
  int lower_case_ranges_num = lower_case_ranges[0];

  if (lower_case_ranges_num == 0) {
    return;
  }

  // For each position-length pair change letters to lowercase
  int prev_pos = 0;
  for (int i = 1; i < lower_case_ranges_num * 2 + 1; i += 2) {
    int start = lower_case_ranges[i];
    int length = lower_case_ranges[i + 1];

    for (int j = 0; j < length; j++) {
      target_seq[j + start + prev_pos] =
          tolower(target_seq[j + start + prev_pos]);
    }

    prev_pos += start + length;
  }
}

void write_reconstructed_sequence_to_file() {
  /**
   * Writes the reconstructed target sequence to a file
   */
  string output_filename = "reconstructed_sequence.txt";
  ofstream out(output_filename);
  if (!out) {
    throw runtime_error("Cannot open output file: " + output_filename);
  }

  out << header << endl;  // Write the header

  out << endl;  // Write an empty line after the header

  // Write out the reconstructed sequence with line breaks
  int line_length_values_num = line_lenghts[0];
  int curr_seq_position = 0;

  for (int i = 1; i < line_length_values_num; i += 2) {
    int lenght = line_lenghts[i];
    int repeat_cnt = line_lenghts[i + 1];

    for (int j = 0; j < repeat_cnt; j++) {
      out.write(&target_seq[curr_seq_position], lenght);
      out.put('\n');
      curr_seq_position += lenght;
    }
  }
  out.close();
}

void decompress_7z(const string& archive_name) {
  string command = "7z x " + archive_name + " -y";
  int result = system(command.c_str());
  if (result != 0) {
    throw runtime_error("7z decompression failed");
  }
}

void cleanup() {
  /**
   * Clean up and release all allocated memory
   */
  ref_seq.clear();
  target_seq.clear();
}

void print_memory_usage() {
  std::ifstream status_file("/proc/self/status");
  std::string line;
  while (std::getline(status_file, line)) {
    if (line.substr(0, 6) == "VmRSS:") {
      std::cout << "Memory usage (Resident Set Size): " << line.substr(6) << std::endl;
      break;
    }
  }
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

  // Start tracking time taken for decompression
  gettimeofday(&timer_start, nullptr);

  // Assign the reference and target file paths from the command line arguments
  InputFileNames input_file_names;

  input_file_names.reference_file = argv[2];
  input_file_names.compressed_target_file = argv[4];

  initialize_structures();

  load_and_clean_reference(input_file_names.reference_file, ref_seq);

  decompress_7z("compressed.7z");

  load_metadata(input_file_names.compressed_target_file);

  decompress_target_sequence(input_file_names.compressed_target_file,
                             target_seq);

  add_special_characters(target_seq);
  add_n_ranges(target_seq);
  add_lowercase_ranges(target_seq);

  write_reconstructed_sequence_to_file();

  cout << "Decompression completed successfully." << endl;
  cout << "Reconstructed sequence written to reconstructed_sequence.txt"
       << endl;
  print_memory_usage();

  // Calculate and print the total time taken for decompression
  gettimeofday(&timer_end, nullptr);

  timer = 1000000 * (timer_end.tv_sec - timer_start.tv_sec) +
          timer_end.tv_usec - timer_start.tv_usec;
  printf("Total decompresion time : %lf ms; \n", timer / 1000.0);

  cleanup();

  return 0;
}