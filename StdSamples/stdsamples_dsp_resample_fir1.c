#include "src.h"
// 44.1kHz * 640 / 147 = 192kHz
FIR_Filter* filter = src_generate(640, 147); 
...
float* input_buffer[IN_SIZE]; // Incoming stream
float* output_buffer[IN_SIZE * 5]; // Outcoming stream
int output_length = src_filt(filter, input_buffer, input_length, output_buffer);