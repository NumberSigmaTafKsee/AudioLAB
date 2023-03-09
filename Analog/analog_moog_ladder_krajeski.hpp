#pragma once

namespace Analog::Moog
{
	template<typename DSP>
	class KrajeskiMoog : public LadderFilterBase
	{

	public:

		KrajeskiMoog(DspFloatType sampleRate) : LadderFilterBase(sampleRate)
		{
			memset(state, 0, sizeof(state));
			memset(delay, 0, sizeof(delay));

			drive = 1.0;
			gComp = 1.0;

			SetCutoff(1000.0f);
			SetResonance(0.1f);
		}

		virtual ~KrajeskiMoog() { }

		void Process(size_t n, DspFloatType * samples, DspFloatType * output)
		{
			Undenormal denormal;
			for (uint32_t s = 0; s < n; ++s)
			{
				state[0] = tanh(drive * (samples[s] - 4 * gRes * (state[4] - gComp * samples[s])));

				for(int i = 0; i < 4; i++)
				{
					state[i+1] = g * (0.3 / 1.3 * state[i] + 1 / 1.3 * delay[i] - state[i + 1]) + state[i + 1];
					delay[i] = state[i];
				}
				output[s] = state[4];
			}
		}

		void Process(size_t n, DspFloatType * samples)
		{
			Undenormal denormal;
			for (uint32_t s = 0; s < n; ++s)
			{
				state[0] = tanh(drive * (samples[s] - 4 * gRes * (state[4] - gComp * samples[s])));

				for(int i = 0; i < 4; i++)
				{
					state[i+1] = g * (0.3 / 1.3 * state[i] + 1 / 1.3 * delay[i] - state[i + 1]) + state[i + 1];
					delay[i] = state[i];
				}
				samples[s] = state[4];
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
			gRes = resonance * (1.0029 + 0.0526 * wc - 0.926 * std::pow(wc, 2) + 0.0218 * std::pow(wc, 3));
		}

		virtual void SetCutoff(DspFloatType c) override
		{
			cutoff = c;
			wc = 2 * MOOG_PI * cutoff / sampleRate;
			g = 0.9892 * wc - 0.4342 * std::pow(wc, 2) + 0.1381 * std::pow(wc, 3) - 0.0202 * std::pow(wc, 4);
		}

	private:

		DspFloatType state[5];
		DspFloatType delay[5];
		DspFloatType wc; // The angular frequency of the cutoff.
		DspFloatType g; // A derived parameter for the cutoff frequency
		DspFloatType gRes; // A similar derived parameter for resonance.
		DspFloatType gComp; // Compensation factor.
		DspFloatType drive; // A parameter that controls intensity of nonlinearities.

	};
}
