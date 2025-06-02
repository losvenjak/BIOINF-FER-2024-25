CC = g++
CFLAG = -O3 -w -Wall -std=c++0x

init:
	@echo "Updating and installing dependencies..."
	@sudo apt-get update
	@sudo apt-get install -y p7zip-full

compress_hirgc: compress_hirgc.cpp
	@$(CC) compress_hirgc.cpp -o compress_hirgc $(CFLAG)
	@echo "Compiled successfully"

decompress_hirgc: decompress_hirgc.cpp
	@$(CC) decompress_hirgc.cpp -o decompress_hirgc $(CFLAG)
	@echo "Compiled successfully"