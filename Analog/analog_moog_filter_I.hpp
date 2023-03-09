/*
The following code was taken from the following post:
https://www.musicdsp.org/en/latest/Filters/24-moog-vcf.html

A C++ class implementation of the filter was posted as a comment by moc.erehwon@tsaot.

That class was revised and added to for this implementation:
 - several variable names were changed for better readability
 - The processBlock function was added for easy use in JUCE projects.
 - Highpass and Bandpass filtering have been added as options
 - A saturation effect was applied. 
   (adapted from the saturation function in the StilsonMoogFilter class)
*/


#pragma once
#include "GenericSoundObject.hpp"

namespace Analog::Filters::Moog::MoogFilterI
{
	template<typename DSP>
	struct MoogFilterI : public GSSoundProcessor<DSP>
	{
		enum {
			LOWPASS,
			HIGHPASS,
			BANDPASS
		};
		MoogFilterI() : GSSoundProcessor<DSP>() {
			init(sampleRate);
			setCutoff(1000.0);
			setResonance(0.01);
		}
		
		void init(DSP sampleRate);
		
		void setCutoff(DSP cutoff);
		void setResonance(DSP resonance);
		
		// saturationAmount ranges [0..1] and is a dry/wet ratio
		// 0.0 will effectively turn off the saturation
		void setSaturation(DSP saturationAmount);

		void ProcessSIMD(size_t n, DSP * input, DSP * output);
		void ProcessBlock(size_t n, DSP * input, DSP * output) {
			ProcessSIMD(n,input,output);
		}
		void ProcessBlock(size_t n, DSP * input) {
			ProcessSIMD(n,input,input);
		}
	
		void calculateCoefficients();
		DSP process(DSP x);
		DSP saturate(DSP input);

		DSP Tick(DSP I, DSP A=1, DSP X=1, DSP Y=1) {
			return A*process(I);
		}

		enum
        {
            PORT_CUTOFF,
            PORT_RESONANCE,			
			PORT_SATURATION,
        };
        void setPort(int port, DSP v)
        {
            switch (port)
            {
            case PORT_CUTOFF:
                setCutoff(v);
                break;
            case PORT_RESONANCE:
                setResonance(v);
                break;		
			case PORT_SATURATION:
				setSaturation(v);
				break;
            }
        }
		int passMode = LOWPASS;
		DSP cutoff;
		DSP resonance;
		DSP saturationAmount;
		
		DSP sampleRate;
		DSP out1, out2, out3, out4;
		DSP in1, in2, in3, in4;

		// coefficients determined by cutoff and resonance
		DSP r, p, k;

		const DSP saturationLimit = 0.95;
	};

	
	void MoogFilterI::init(DSP sampleRate)
	{
		this->sampleRate = sampleRate;

		out1 = out2 = out3 = out4 = 0.0f;
		in1 = in2 = in3 = in4 = 0.0f;
		calculateCoefficients();
	};

	void MoogFilterI::setCutoff(DSP cutoff)
	{
		this->cutoff = cutoff;
		calculateCoefficients();
	}

	void MoogFilterI::setResonance(DSP resonance)
	{
		if (resonance >= 0.0 && resonance <= 1.0)
		{
			this->resonance = resonance;
			calculateCoefficients();
		}
	}

	void MoogFilterI::setSaturation(DSP saturationAmount)
	{
		if (saturationAmount >= 0.0 && saturationAmount <= 1.0)
		{
			this->saturationAmount = saturationAmount;
		}
	}

	void MoogFilterI::ProcessBlock(size_t n, DSP * inputs, DSP * outputs)
	{
		for (int s = 0; s < n; s++)
		{
			outputs[s] =  process(inputs[s]);
		}
	}

	void MoogFilterI::calculateCoefficients()
	{
		DSP f = (cutoff + cutoff) / this->sampleRate;
		p = f * (1.8f - 0.8f * f);
		k = p + p - 1.f;

		DSP t = (1.f - p) * 1.386249f;
		DSP t2 = 12.f + t * t;
		r = resonance * (t2 + 6.f * t) / (t2 - 6.f * t);
	};

	DSP MoogFilterI::process(DSP input)
	{
		Undenormal denormals;
		// process input
		input = saturate(input);
		input -= r * out4;

		//Four cascaded onepole filters (bilinear transform)
		out1 = input * p + in1 * p - k * out1;
		out2 =  out1 * p + in2 * p - k * out2;
		out3 =  out2 * p + in1 * p - k * out3;
		out4 =  out3 * p + in4 * p - k * out4;

		in1 = input;
		in2 = out1;
		in3 = out2;
		in4 = out3;

		//Clipper band limited sigmoid
		out4 -= (out4 * out4 * out4) / 6.f;

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
	}
	void ProcessSIMD(size_t n, DSP * in, DSP * output) {
		Undenormal denormals;
		#pragma omp simd aligned(input,output)
		for(size_t i = 0; i < n; i++)
		{
			// process input
			DSP input = saturate(in[i]);
			input -= r * out4;

			//Four cascaded onepole filters (bilinear transform)
			out1 = input * p + in1 * p - k * out1;
			out2 =  out1 * p + in2 * p - k * out2;
			out3 =  out2 * p + in1 * p - k * out3;
			out4 =  out3 * p + in4 * p - k * out4;

			in1 = input;
			in2 = out1;
			in3 = out2;
			in4 = out3;

			//Clipper band limited sigmoid
			out4 -= (out4 * out4 * out4) / 6.f;

			switch (passMode)
			{
			case LOWPASS:
				output[i] =  out4;
				break;
			case HIGHPASS:
				output[i] = input - out4 - out1;
				break;
			case BANDPASS:
				output[i] = out4 - out1;
				break;
			default:
				output[i] = out4;
			}
		}
	}
	
	DSP MoogFilterI::saturate(DSP input)
	{
		DSP drySignal = input;
		input *= 1.5f;

		DSP x1 = fabsf(input + saturationLimit);
		DSP x2 = fabsf(input - saturationLimit);
		DSP wetSignal = (0.5f * (x1 - x2));
		
		return (wetSignal * saturationAmount) + (drySignal * (1.0f - saturationAmount));
	}
}

