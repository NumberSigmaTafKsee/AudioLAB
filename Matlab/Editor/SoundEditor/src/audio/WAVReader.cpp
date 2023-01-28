#include "WAVReader.h"

#include <iostream>

#include <wx/filesys.h>

WAVReader::WAVReader() {
	sound = nullptr;
}

WAVReader::~WAVReader() {
	delete sound;
}

WAV_FILE* WAVReader::GetLastSound() {
	return sound;
}

WAV_FILE* WAVReader::ReadFileData(std::string url) {
	WAV_FILE* temp = new WAV_FILE();

	FILE* file = fopen(url.c_str(), "r+");

	if (file != NULL) {
		// Read the header
		temp->header = new WAV_HEADER();

		size_t bytesRead = fread(temp->header, 1, 44, file);

		if (bytesRead > 0){
			// Read the data
			temp->data = new std::vector<short>(temp->header->DATA_chunk_size/2);

			bytesRead = fread(temp->data->data(), 2, temp->header->DATA_chunk_size, file);
		}

		fclose(file);

		sound = temp;
	}else {
		std::cout << "Problem opening the file...\n";
	}

	return temp;
}

