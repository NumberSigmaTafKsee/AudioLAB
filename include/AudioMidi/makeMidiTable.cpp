#include <cmath>
#include "MusicFunctions.hpp"
#include <cstdio>
#include <cstdlib>

int main()
{
	FILE * f = fopen("midi_table.h","w");
	fprintf(f,"double midi_table[128] = {");
	for(size_t i = 0; i < 127; i++)
		fprintf(f,"%f,",midi_to_freq(i));
	fprintf(f,"};\n");
	fclose(f);
}
