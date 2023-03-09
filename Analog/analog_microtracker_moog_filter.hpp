#pragma once

#include <algorithm>
#include <cmath>
#include "Undenormal.hpp"
#include "GenericSoundObject.hpp"

namespace Analog::Filters::Moog
{
	///////////////////////////////////////////////////////////////////////////////////////////
	// Stolen from microtracker
	///////////////////////////////////////////////////////////////////////////////////////////
	template<typename DSP>
	class MicrotrackerMoog : public GSSoundProcessor<DSP>
	{

	public:

		MicrotrackerMoog(DSP sr) : GSSoundProcessor<DSP>(),sampleRate(sr)
		{
			p0 = p1 = p2 = p3 = p32 = p33 = p34 = 0.0;
			SetCutoff(1000.0f);
			SetResonance(0.10f);
		}

		virtual ~MicrotrackerMoog() {}

		DSP fast_tanh(DSP x) 
		{
			DSP x2 = x * x;
			return x * (27.0 + x2) / (27.0 + 9.0 * x2);
		}

		void ProcessSIMD(size_t n, DSP * samples, DSP * output)
		{
			Undenormal denormal;
			DSP k = resonance * 4;
			#pragma omp simd aligned(samples,output)
			for (uint32_t s = 0; s < n; ++s)
			{
				// Coefficients optimized using differential evolution
				// to make feedback gain 4.0 correspond closely to the
				// border of instability, for all values of omega.
				DSP out = p3 * 0.360891 + p32 * 0.417290 + p33 * 0.177896 + p34 * 0.0439725;

				p34 = p33;
				p33 = p32;
				p32 = p3;

				p0 += (fast_tanh(samples[s] - k * out) - fast_tanh(p0)) * cutoff;
				p1 += (fast_tanh(p0) - fast_tanh(p1)) * cutoff;
				p2 += (fast_tanh(p1) - fast_tanh(p2)) * cutoff;
				p3 += (fast_tanh(p2) - fast_tanh(p3)) * cutoff;

				output[s] = out;
			}
		}
	
		void ProcessBlock(size_t n, DSP * in, DSP * out)
		{
			ProcessSIMD(n,in,out);
		}
		void ProcessInplace(size_t n, DSP * samples)
		{
			ProcessSIMD(n,samples,samples);
		}

		
		DSP Tick(DSP I, DSP A=1, DSP X=1, DSP Y=1) {
			Undenormal denormal;
			DSP k = resonance * 4;
			DSP out = p3 * 0.360891 + p32 * 0.417290 + p33 * 0.177896 + p34 * 0.0439725;

			p34 = p33;
			p33 = p32;
			p32 = p3;

			p0 += (fast_tanh(I - k * out) - fast_tanh(p0)) * cutoff;
			p1 += (fast_tanh(p0) - fast_tanh(p1)) * cutoff;
			p2 += (fast_tanh(p1) - fast_tanh(p2)) * cutoff;
			p3 += (fast_tanh(p2) - fast_tanh(p3)) * cutoff;

			return A * out;

		}

		void SetResonance(DSP r)
		{
			resonance = r;
		}

		void SetCutoff(DSP c)
		{
			cutoff = c * 2 * M_PI / sampleRate;
			cutoff = std::min(cutoff, (DSP)1);
		}

		DSP GetResonance() { return resonance; }
		DSP GetCutoff() { return cutoff; }

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
                SetCutoff(v);
                break;
            case PORT_RESONANCE:
                SetResonance(v);
                break;			
            }
        }
	private:

		DSP p0;
		DSP p1;
		DSP p2;
		DSP p3;
		DSP p32;
		DSP p33;
		DSP p34;
		DSP cutoff,resonance,sampleRate;
	};
}
