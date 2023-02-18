#pragma once

namespace Analog::Moog
{
	///////////////////////////////////////////////////////////////////////////////////////////
	// Finnish Vampire
	///////////////////////////////////////////////////////////////////////////////////////////
	class HuovilainenMoog : public LadderFilterBase
	{
	public:

		HuovilainenMoog(DspFloatType sampleRate) : LadderFilterBase(sampleRate), thermal(0.000025)
		{
			memset(stage, 0, sizeof(stage));
			memset(delay, 0, sizeof(delay));
			memset(stageTanh, 0, sizeof(stageTanh));
			SetCutoff(1000.0f);
			SetResonance(0.10f);
		}

		virtual ~HuovilainenMoog()
		{

		}

		void ProcessBlock(size_t n, DspFloatType * _input, DspFloatType * _output)
		{
			Undenormal denormal;
			#pragma omp simd aligned(_input,_output)
			for (uint32_t s = 0; s < n; ++s)
			{
				// Oversample
				for (int j = 0; j < 2; j++)
				{
					DspFloatType input = _input[s] - resQuad * delay[5];
					delay[0] = stage[0] = delay[0] + tune * (tanh(input * thermal) - stageTanh[0]);
					for (int k = 1; k < 4; k++)
					{
						input = stage[k-1];
						stage[k] = delay[k] + tune * ((stageTanh[k-1] = tanh(input * thermal)) - (k != 3 ? stageTanh[k] : tanh(delay[k] * thermal)));
						delay[k] = stage[k];
					}
					// 0.5 sample delay for phase compensation
					delay[5] = (stage[3] + delay[4]) * 0.5;
					delay[4] = stage[3];
				}
				_output[s] = delay[5];
			}
		}
		DspFloatType Tick(DspFloatType _input) {		
			// Oversample
			#pragma omp simd
			for (int j = 0; j < 2; j++)
			{
				DspFloatType input = _input - resQuad * delay[5];
				delay[0] = stage[0] = delay[0] + tune * (tanh(input * thermal) - stageTanh[0]);
				for (int k = 1; k < 4; k++)
				{
					input = stage[k-1];
					stage[k] = delay[k] + tune * ((stageTanh[k-1] = tanh(input * thermal)) - (k != 3 ? stageTanh[k] : tanh(delay[k] * thermal)));
					delay[k] = stage[k];
				}
				// 0.5 sample delay for phase compensation
				delay[5] = (stage[3] + delay[4]) * 0.5;
				delay[4] = stage[3];
			}		
			return delay[5];
		}

		virtual void SetResonance(DspFloatType r) override
		{
			resonance = r;
			resQuad = 4.0 * resonance * acr;
		}

		virtual void SetCutoff(DspFloatType c) override
		{
			cutoff = c;

			DspFloatType fc =  cutoff / sampleRate;
			DspFloatType f  =  fc * 0.5; // oversampled
			DspFloatType fc2 = fc * fc;
			DspFloatType fc3 = fc * fc * fc;

			DspFloatType fcr = 1.8730 * fc3 + 0.4955 * fc2 - 0.6490 * fc + 0.9988;
			acr = -3.9364 * fc2 + 1.8409 * fc + 0.9968;

			tune = (1.0 - exp(-((2 * MOOG_PI) * f * fcr))) / thermal;

			SetResonance(resonance);
		}

	private:

		DspFloatType stage[4];
		DspFloatType stageTanh[3];
		DspFloatType delay[6];

		DspFloatType thermal;
		DspFloatType tune;
		DspFloatType acr;
		DspFloatType resQuad;

	};
}
