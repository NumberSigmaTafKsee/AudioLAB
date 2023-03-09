#pragma once


#include "MoogFilters.hpp"
#include "GenericSoundObject.hpp"

namespace Liquid
{
	template<typename DSP>
	struct LiquidMoog : public GSSoundProcessor<DSP>
	{
		Analog::MoogFilters::MoogLadderFilter moog;
		DSP cutoff,Q;
		DSP low,high,band,notch;
		enum {
			LOWPASS,
			HIGHPASS,
			BANDPASS,
			NOTCH,
			APF,
			UBS,
			PEAK,
			SHELF,        
		};

		int type  = LOWPASS;

		LiquidMoog(Analog::MoogFilters::MoogModelType type = Analog::MoogFilters::MoogModelType::FINN_MOOG) :        
		GSSoundProcessor<DSP>(),
		moog(type,sampleRate)
		{
			cutoff = Q = low = high = band = notch = 0;
		}

		void setCutoff(DSP c) {
			cutoff = MusicFunctions::cv2freq(c);
			moog.SetCutoff(cutoff);
		}
		void setResonance(DSP q) {
			Q = q;
			moog.SetResonance(q);
		}
		void setType(Analog::MoogFilters::MoogModelType type) {
			moog.setType(type);
		}
		DSP Tick(DSP I, DSP A=1, DSP X=0, DSP Y=0) {
			Undenormal denormal;
			DSP out;

			low = moog.Tick(I,A,X,Y);

			DSP x_scale = 1.5*(1.0-Q);
			if(x_scale == 0.0f) x_scale = 0.001f;
			DSP x_f     = 2*std::sin(2 * M_PI * cutoff/sampleRate);        

			high  = x_scale*I - low - x_scale*band;
			band  = x_f * high + band;
			notch = high + low;
			/*
			ubp = 2 * scale * bp;              
			shelf = input + 2*K*scale*bp;        
			apf   = xn - 4*scale*bp;
			peak  = lp - hp;
			*/
			switch(type)
			{
				case LOWPASS: out = low; break;
				case HIGHPASS: out = high; break;
				case BANDPASS: out = band; break;
				case NOTCH: out = notch; break;
			}
			out *= A;
			return out;
		}        
		void ProcessSIMD(size_t n, DSP * input, DSP * output)
		{
		   Undenormal denormal;
		   #pragma omp simd aligned(input,output)
		   for(size_t i = 0; i < n; i++)
		   {
				DSP out;
				DSP I = input[i];
				low = moog.Tick(I);

				DSP x_scale = 1.5*(1.0-Q);
				if(x_scale == 0.0f) x_scale = 0.001f;
				DSP x_f     = 2*std::sin(2 * M_PI * cutoff/sampleRate);        

				high  = x_scale*I - low - x_scale*band;
				band  = x_f * high + band;
				notch = high + low;
				/*
				ubp = 2 * scale * bp;              
				shelf = input + 2*K*scale*bp;        
				apf   = xn - 4*scale*bp;
				peak  = lp - hp;
				*/
				switch(type)
				{
					case LOWPASS: out = low; break;
					case HIGHPASS: out = high; break;
					case BANDPASS: out = band; break;
					case NOTCH: out = notch; break;
				}
				output[i] = out;
		}
		void ProcessBlock(size_t n, DSP * in, DSP * out)
		{
			ProcessSIMD(n,in,out);
		}
		void ProcessInplace(size_t n, DSP * samples)
		{
			ProcessSIMD(n,samples,samples);
		}
    };
}
