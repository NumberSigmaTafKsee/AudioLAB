#pragma once 

#include "DSP/resampler.cc"   
#include "DSP/resampler-table.cc"
#include "DSP/resampler.h"

class SimpleResampler {
 private:
    Resampler r_up, r_down;
    int m_fact;
 public:
    SimpleResampler(): r_up(), r_down(), m_fact() {}
    void setup(int sampleRate, unsigned int fact);
    void up(int count, DspFloatType *input, DspFloatType *output);
    void down(int count, DspFloatType *input, DspFloatType *output);
};

void SimpleResampler::setup(int sampleRate, unsigned int fact)
{
	m_fact = fact;
	const int qual = 16; // resulting in a total delay of 2*qual (0.7ms @44100)
	// upsampler
	r_up.setup(sampleRate, sampleRate*fact, 1, qual);
	// k == inpsize() == 2 * qual
	// pre-fill with k-1 zeros
	r_up.inp_count = r_up.inpsize() - 1;
	r_up.out_count = 1;
	r_up.inp_data = r_up.out_data = 0;
	r_up.process();
	// downsampler
	r_down.setup(sampleRate*fact, sampleRate, 1, qual);
	// k == inpsize() == 2 * qual * fact
	// pre-fill with k-1 zeros
	r_down.inp_count = r_down.inpsize() - 1;
	r_down.out_count = 1;
	r_down.inp_data = r_down.out_data = 0;
	r_down.process();
}

void SimpleResampler::up(int count, DspFloatType *input, DspFloatType *output)
{
	float temp[count],out[m_fact*count];
	for(size_t i = 0; i < count; i++) temp[i] = input[i];
	r_up.inp_count = count;
	r_up.inp_data = temp;
	r_up.out_count = count * m_fact;
	r_up.out_data = out;
	r_up.process();
	assert(r_up.inp_count == 0);
	assert(r_up.out_count == 0);
	for(size_t i = 0; i < m_fact*count; i++) output[i] = out[i];
}

void SimpleResampler::down(int count, DspFloatType *input, DspFloatType *output)
{
	float temp[m_fact*count],out[count];
	for(size_t i = 0; i < m_fact*count; i++) temp[i] = input[i];	
	r_down.inp_count = count * m_fact;
	r_down.inp_data = temp;
	r_down.out_count = count +1;// == trick to drain input
	r_down.out_data = out;
	r_down.process();
	assert(r_down.inp_count == 0);
	assert(r_down.out_count == 1);
	for(size_t i = 0; i < count; i++) output[i] = out[i];
}