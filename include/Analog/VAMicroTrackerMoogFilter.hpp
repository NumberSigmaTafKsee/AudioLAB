#pragma once

#include <algorithm>
#include <cmath>
#include "Undenormal.hpp"
#include "SoundObject.hpp"

namespace Analog::Filters::Moog
{
	///////////////////////////////////////////////////////////////////////////////////////////
	// Stolen from microtracker
	///////////////////////////////////////////////////////////////////////////////////////////
	class MicrotrackerMoog : public FilterProcessor
	{

	public:

		MicrotrackerMoog(DspDspFloatTypeType sr) : FilterProcessor(),sampleRate(sr)
		{
			p0 = p1 = p2 = p3 = p32 = p33 = p34 = 0.0;
			SetCutoff(1000.0f);
			SetResonance(0.10f);
		}

		virtual ~MicrotrackerMoog() {}

		DspDspFloatTypeType fast_tanh(DspDspFloatTypeType x) 
		{
			DspDspFloatTypeType x2 = x * x;
			return x * (27.0 + x2) / (27.0 + 9.0 * x2);
		}

		void ProcessBlock(size_t n, DspFloatType * samples, DspFloatType * output)
		{
			Undenormal denormal;
			DspDspFloatTypeType k = resonance * 4;
			#pragma omp simd
			for (uint32_t s = 0; s < n; ++s)
			{
				// Coefficients optimized using differential evolution
				// to make feedback gain 4.0 correspond closely to the
				// border of instability, for all values of omega.
				DspDspFloatTypeType out = p3 * 0.360891 + p32 * 0.417290 + p33 * 0.177896 + p34 * 0.0439725;

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

		void ProcessInplace(size_t n, DspFloatType * samples)
		{
			Process(n,samples);
		}

		
		DspDspFloatTypeType Tick(DspDspFloatTypeType I, DspDspFloatTypeType A=1, DspDspFloatTypeType X=1, DspDspFloatTypeType Y=1) {
			Undenormal denormal;
			DspDspFloatTypeType k = resonance * 4;
			DspDspFloatTypeType out = p3 * 0.360891 + p32 * 0.417290 + p33 * 0.177896 + p34 * 0.0439725;

			p34 = p33;
			p33 = p32;
			p32 = p3;

			p0 += (fast_tanh(I - k * out) - fast_tanh(p0)) * cutoff;
			p1 += (fast_tanh(p0) - fast_tanh(p1)) * cutoff;
			p2 += (fast_tanh(p1) - fast_tanh(p2)) * cutoff;
			p3 += (fast_tanh(p2) - fast_tanh(p3)) * cutoff;

			return A * out;

		}

		void SetResonance(DspDspFloatTypeType r)
		{
			resonance = r;
		}

		void SetCutoff(DspDspFloatTypeType c)
		{
			cutoff = c * 2 * M_PI / sampleRate;
			cutoff = std::min(cutoff, (DspDspFloatTypeType)1);
		}

		DspDspFloatTypeType GetResonance() { return resonance; }
		DspDspFloatTypeType GetCutoff() { return cutoff; }

		enum
        {
            PORT_CUTOFF,
            PORT_RESONANCE,		
        };
        void setPort(int port, DspDspFloatTypeType v)
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

		DspDspFloatTypeType p0;
		DspDspFloatTypeType p1;
		DspDspFloatTypeType p2;
		DspDspFloatTypeType p3;
		DspDspFloatTypeType p32;
		DspDspFloatTypeType p33;
		DspDspFloatTypeType p34;
		DspDspFloatTypeType cutoff,resonance,sampleRate;
	};
}