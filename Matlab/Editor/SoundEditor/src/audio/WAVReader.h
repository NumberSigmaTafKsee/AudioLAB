#pragma once

#include <cstdint>
#include <cstdlib>

#include <string>
#include <vector>
#include <fstream>

struct WAV_HEADER {
	//RIFF
	int RIFF_chunk_ID;
	int RIFF_chunk_size;
	int RIFF_format;
	//FMT
	int FMT_chunk_ID;
	int FMT_chunk_size;
	short FMT_audio_format;
	short FMT_num_channels;
	int FMT_sample_rate;
	int FMT_byte_rate;
	short FMT_block_align;
	short FMT_bits_per_sample;
	//DATA
	int DATA_chunk_ID;
	int DATA_chunk_size;
};

struct WAV_FILE {
	WAV_HEADER* header;
	std::vector<short>* data;
};

class WAVReader {
public:
	WAVReader();
	~WAVReader();
public:
	WAV_FILE* ReadFileData(std::string url);
	WAV_FILE* GetLastSound();
private:
	WAV_FILE* sound;
};