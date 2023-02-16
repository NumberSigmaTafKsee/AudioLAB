/*	Copyright (C) 2018, 2020 Harris M. Snyder

	This file is part of guitardsp.

	guitardsp is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	guitardsp is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with guitardsp.  If not, see <https://www.gnu.org/licenses/>.
*/
#ifndef CONVOLUTION_H
#define CONVOLUTION_H


typedef struct convolution
{
	char * data;
	unsigned int periodsz;
	unsigned int n;
} convolution_t;

// returns 0 on success, -1 on error. 
// -1 means a critical error (you need to abort). 
// errMsg might be set even on success (a warning). You can print this or ignore it.
int convolution_construct(convolution_t * self, char ** errMsg, unsigned int * rate, unsigned int period_sz, const char * IR_filename, unsigned int IR_max_size_truncate);

void convolution_destruct(convolution_t * self);

// returns a pointer to the input buffer for this convolution object. Write samples to that buffer before calling convolution_apply
float * convolution_getInputPtr(convolution_t * self);


#include <string.h>
#include <arm_neon.h>

static inline void convolution_apply(convolution_t * self, float * output)
{
	// this is a pretty straightforward convolution implementation, though it's perhaps a bit obfuscated by the explicit SIMD intrinsics. 
	// a previous non-SIMD version did not give good enough performance for inaudible latency on a raspberry pi 3b.

	float * buffer = (float*) (self->data + sizeof(float32x4_t) * self->n);
	float32x4_t* kernel_reverse = (float32x4_t*)self->data;
	
	float32x4_t out, data_block;
	for(int i=0; i<self->periodsz; i+=4){
		out = vmovq_n_f32(0.0f);
		// After this loop, we have computed 4 output samples for the price of one.
		for(int k=0; k<self->n; k++)
		{
			data_block = vld1q_f32(&buffer[i+k]);
			out = vfmaq_f32(out,data_block,kernel_reverse[k]);
		}
		vst1q_f32(output+i, out);

	}
	memmove(buffer, buffer + self->periodsz, (self->n - 1)*sizeof(float));
}

#endif
