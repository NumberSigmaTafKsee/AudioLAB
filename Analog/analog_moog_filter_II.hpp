/*
This class is an adaption of the Moog VCF variation 2, posted here:
https://www.musicdsp.org/en/latest/Filters/26-moog-vcf-variation-2.html

Comes from a Stilson/Smith CCRMA paper

added highpass and bandpass options
*/

#pragma once
#include "GenericSoundObject.hpp"


namespace Analog::Filters::Moog::MoogFilterII
{
	template<typename DSP>
	struct MoogFilterII : public GSSoundProcessor<DSP>
	{	
		enum {
			LOWPASS,
			HIGHPASS,
			BANDPASS
		};
		MoogFilterII() : FilterProcessor() {
			init(sampleRate);
			set(1000.0,0.5);			
		}

		void init(DSP sampleRate);
		void set(DSP cutoff, DSP resonance);
		DSP processSample(DSP in);
		
		void  ProcessSIMD(size_t n, DSP * inputs, DSP * outputs);
		void  ProcessBlock(size_t n, DSP * inputs, DSP * outputs) {
			ProcessSIMD(n,inputs,outputs);
		}
		void ProcessInplace(size_t n, DSP * p) {
			ProcessSIMD(n,p,p);
		}

		enum
        {
            PORT_CUTOFF,
            PORT_RESONANCE,			
        };
        void setPort(int port, DSP v)
        {
            switch (port)
            {
            case PORT_CUTOFF:
                set(v,resonance);
                break;
            case PORT_RESONANCE:
                set(cutoff,v);
                break;		
            }
        }
		DSP Tick(DSP I, DSP A=1, DSP X=1, DSP Y=1) {
			return A*processSample(I);
		}

		int passMode = LOWPASS;
		DSP sampleRate;
		DSP cutoff, resonance;
		DSP out1, out2, out3, out4;
		DSP in1, in2, in3, in4;
	};


	void MoogFilterII::init(DSP sampleRate)
	{
		this->sampleRate = sampleRate;
		out1 = 0.0f;
		out2 = 0.0f;
		out3 = 0.0f;
		out4 = 0.0f;
		in1 = 0.0f;
		in2 = 0.0f;
		in3 = 0.0f;
		in4 = 0.0f;
	}

	// Set coefficients given frequency & resonance [0.0...1.0]
	void MoogFilterII::set(DSP cutoff, DSP resonance)
	{
		cutoff = cutoff / (sampleRate/2.0f);
		resonance = resonance * 4.0f;
		this->cutoff = cutoff;
		this->resonance = resonance;
	}


	// Filter (in [-1.0...+1.0])
	DSP MoogFilterII::processSample(DSP input)
	{
		Undenormal denormals;
		DSP in = input;
		DSP f = cutoff * 1.16;
		DSP fb = resonance * (1.0 - 0.15 * f * f);
		in -= out4 * fb;
		in *= 0.35013 * (f * f) * (f * f);
		out1 = in + 0.3 * in1 + (1 - f) * out1; // Pole 1
		in1 = in;
		out2 = out1 + 0.3 * in2 + (1 - f) * out2;  // Pole 2
		in2 = out1;
		out3 = out2 + 0.3 * in3 + (1 - f) * out3;  // Pole 3
		in3 = out2;
		out4 = out3 + 0.3 * in4 + (1 - f) * out4;  // Pole 4
		in4 = out3;

		switch (passMode)
		{
		case LOWPASS:
			return out4;
			break;
		case HIGHPASS:
			return input - out4 - out1;
			break;
		case BANDPASS:
			return out4 - out1;
			break;
		default:
			return out4;
		}

		// Lowpass = out4
		// Highpass = input - out4 - out1
		// Bandpass = out4 - out1
	}

	void MoogFilterII::ProcessSIMD(size_t numSamples, DSP * inputs, DSP * outputs)
	{
		Undenormal denormals;
		#pragma omp simd aligned(inputs,outputs)
		for (int s = 0; s < numSamples; s++)
		{
			DSP in = inputs[i];
			DSP f = cutoff * 1.16;
			DSP fb = resonance * (1.0 - 0.15 * f * f);
			in -= out4 * fb;
			in *= 0.35013 * (f * f) * (f * f);
			out1 = in + 0.3 * in1 + (1 - f) * out1; // Pole 1
			in1 = in;
			out2 = out1 + 0.3 * in2 + (1 - f) * out2;  // Pole 2
			in2 = out1;
			out3 = out2 + 0.3 * in3 + (1 - f) * out3;  // Pole 3
			in3 = out2;
			out4 = out3 + 0.3 * in4 + (1 - f) * out4;  // Pole 4
			in4 = out3;

			switch (passMode)
			{
			case LOWPASS:
				outputs[i] = out4;
				break;
			case HIGHPASS:
				outputs[i] = in - out4 - out1;
				break;
			case BANDPASS:
				outputs[i] = out4 - out1;
				break;
			default:
				outputs[i] = out4;
			}

		// Lowpass = out4
		// Highpass = input - out4 - out1
		// Bandpass = out4 - out1
		}
	}
}
#undef LOWPASS
#undef HIGHPASS
#undef BANDPASS
