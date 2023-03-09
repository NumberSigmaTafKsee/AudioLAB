#pragma once

namespace Analog::Moog
{
	template<typename DSP>
	class MicrotrackerMoog : public LadderFilterBase
	{

	public:

		MicrotrackerMoog(DspFloatType sampleRate) : LadderFilterBase(sampleRate)
		{
			p0 = p1 = p2 = p3 = p32 = p33 = p34 = 0.0;
			SetCutoff(1000.0f);
			SetResonance(0.10f);
		}

		virtual ~MicrotrackerMoog() {}

		void Process(size_t n, DspFloatType * samples, DspFloatType * output)
		{
			Undenormal denormal;
			DspFloatType k = resonance * 4;
			for (uint32_t s = 0; s < n; ++s)
			{
				// Coefficients optimized using differential evolution
				// to make feedback gain 4.0 correspond closely to the
				// border of instability, for all values of omega.
				DspFloatType out = p3 * 0.360891 + p32 * 0.417290 + p33 * 0.177896 + p34 * 0.0439725;

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

		void Process(size_t n, DspFloatType * samples)
		{
			Undenormal denormal;
			DspFloatType k = resonance * 4;
			for (uint32_t s = 0; s < n; ++s)
			{
				// Coefficients optimized using differential evolution
				// to make feedback gain 4.0 correspond closely to the
				// border of instability, for all values of omega.
				DspFloatType out = p3 * 0.360891 + p32 * 0.417290 + p33 * 0.177896 + p34 * 0.0439725;

				p34 = p33;
				p33 = p32;
				p32 = p3;

				p0 += (fast_tanh(samples[s] - k * out) - fast_tanh(p0)) * cutoff;
				p1 += (fast_tanh(p0) - fast_tanh(p1)) * cutoff;
				p2 += (fast_tanh(p1) - fast_tanh(p2)) * cutoff;
				p3 += (fast_tanh(p2) - fast_tanh(p3)) * cutoff;

				samples[s] = out;
			}
		}

		void Process(DspFloatType * samples, DspFloatType * modulation, uint32_t n)
		{
			Undenormal denormal;
			DspFloatType k = resonance * 4;
			for (uint32_t s = 0; s < n; ++s)
			{
				SetCutoff(modulation[s]);
				// Coefficients optimized using differential evolution
				// to make feedback gain 4.0 correspond closely to the
				// border of instability, for all values of omega.
				DspFloatType out = p3 * 0.360891 + p32 * 0.417290 + p33 * 0.177896 + p34 * 0.0439725;

				p34 = p33;
				p33 = p32;
				p32 = p3;

				p0 += (fast_tanh(samples[s] - k * out) - fast_tanh(p0)) * cutoff;
				p1 += (fast_tanh(p0) - fast_tanh(p1)) * cutoff;
				p2 += (fast_tanh(p1) - fast_tanh(p2)) * cutoff;
				p3 += (fast_tanh(p2) - fast_tanh(p3)) * cutoff;

				samples[s] = out;
			}
		}

		DspFloatType Tick(DspFloatType input) {
			DspFloatType r = 0.0;
			Process(1,&input,&r);
			return r;

		}

		virtual void SetResonance(DspFloatType r) override
		{
			resonance = r;
		}

		virtual void SetCutoff(DspFloatType c) override
		{
			cutoff = c * 2 * MOOG_PI / sampleRate;
			cutoff = moog_min(cutoff, 1);
		}

	private:

		DspFloatType p0;
		DspFloatType p1;
		DspFloatType p2;
		DspFloatType p3;
		DspFloatType p32;
		DspFloatType p33;
		DspFloatType p34;
	};
}
