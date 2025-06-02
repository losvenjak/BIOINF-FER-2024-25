# High-speed and high-ration referential genome compression algorithm

The High-speed and High-ratio Referential Genome Compression Algorithm addresses the challenge of compressing genomic data by using a reference-based approach. Instead of storing entire genome sequences, it records only the differences between the target genome and a closely related reference genome, significantly reducing data size. This method achieves a high compression ratio while ensuring fast processing and minimal computational resource usage, making it well-suited for managing the rapidly growing volume of genomic data in research, healthcare, and bioinformatics applications.

This project was created as part of the Bioinformatics course at the Faculty of Electrical Engineering and Computing:
https://www.fer.unizg.hr/predmet/bio1


# Features

- Efficient Data Reduction – Uses referential compression to store only differences, minimizing storage requirements.
- High-speed Processing – Optimized alignment and encoding ensure fast compression and decompression.
- Low Memory Usage – Designed to operate efficiently without excessive computational demands.
- 2-bit Encoding – Compact representation of nucleotide sequences for further space savings.
- Scalability – Suitable for large-scale genomic datasets.


# Parts of the algorithm

- Preprocessing – Removes unnecessary formatting and metadata from genomic files.
- Encoding – Converts nucleotides (A, C, G, T) into a compact 2-bit representation to save space.
- Alignment – Efficiently finds matching regions between the target and reference genome.
- Compression – Stores only the differences while maintaining the ability to reconstruct the original sequence.
- Decompression – Reconstructs the original genome by applying stored variations to the reference genome.

## Installation and Running Instructions

# Install dependencies
    make init
        
# Compile
    make compress_hirgc
    make decompress_hirgc

# Compress
    ./compress_hirgc -r <reference_file_name> -t <target_file_name>

# Decompress
    ./decompress_hirgc -r <reference_file_name> -t <target_file_name>

# Run example
    follow previous steps for installing dependencies and compiling
    
    compress using command
    ./compress_hirgc -r ref.fna -t tar.fna

    decompress using command
    ./compress_hirgc -r ref.fna -t output.txt
