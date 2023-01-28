#include "src.h"
int interpolation = 2; // x 2
int decimation = 1;    // / 1
float cutoff_frequency = 0.5; // half the sampling frequency
int num_taps = 24 * interpolation;
float* coefficients = src_generate_fir_coeffs(num_taps, cutoff_frequency / interpolation);
FIR_Filter* filter_left = src_generate_fir_filter(coefficients, num_taps, interpolation, decimation);
FIR_Filter* filter_right = src_generate_fir_filter(coefficients, num_taps, interpolation, decimation);
free(coefficients);
...
float* input_buffer_left[IN_SIZE]; // Incoming stream
float* output_buffer_left[IN_SIZE * interpolation / decimation + 1]; // Outcoming stream
int output_length = src_filt(filter_left, input_buffer, input_length, output_buffer);