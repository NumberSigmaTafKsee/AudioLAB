#pragma once

namespace Analog::Moog
{
	///////////////////////////////////////////////////////////////////////////////////////////
	// Runge-Kutta Simulator
	///////////////////////////////////////////////////////////////////////////////////////////
	template<typename DSP>
	class RKSimulationMoog : public LadderFilterBase
	{
	public:

		RKSimulationMoog(DspFloatType sampleRate) : LadderFilterBase(sampleRate)
		{
			memset(state, 0, sizeof(state));

			saturation = 3.0;
			saturationInv = 1.0 / saturation;

			oversampleFactor = 1;

			stepSize = 1.0 / (oversampleFactor * sampleRate);

			SetCutoff(1000.f);
			SetResonance(1.0f);
		}

		virtual ~RKSimulationMoog()
		{
		}

		void Process(size_t n, DspFloatType * samples, DspFloatType * output)
		{
			Undenormal denormal;
			for (uint32_t s = 0; s < n; ++s)
			{
				for (int j = 0; j < oversampleFactor; j++)
				{
					rungekutteSolver(samples[s], state);
				}

				output[s] = state[3];
			}
		}
		void Process(size_t n, DspFloatType * samples)
		{
			Undenormal denormal;
			for (uint32_t s = 0; s < n; ++s)
			{
				for (int j = 0; j < oversampleFactor; j++)
				{
					rungekutteSolver(samples[s], state);
				}
				samples[s] = state[3];
			}
		}
		DspFloatType Tick(DspFloatType input) {
			DspFloatType r;
			Process(1,&input,&r);
			return r;
		}

		virtual void SetResonance(DspFloatType r) override
		{
			// 0 to 10
			resonance = r;
		}

		virtual void SetCutoff(DspFloatType c) override
		{
			cutoff = (2.0 * MOOG_PI * c);
		}

	private:

		void calculateDerivatives(DspFloatType input, DspFloatType * dstate, DspFloatType * state)
		{
			DspFloatType satstate0 = clip(state[0], saturation, saturationInv);
			DspFloatType satstate1 = clip(state[1], saturation, saturationInv);
			DspFloatType satstate2 = clip(state[2], saturation, saturationInv);

			dstate[0] = cutoff * (clip(input - resonance * state[3], saturation, saturationInv) - satstate0);
			dstate[1] = cutoff * (satstate0 - satstate1);
			dstate[2] = cutoff * (satstate1 - satstate2);
			dstate[3] = cutoff * (satstate2 - clip(state[3], saturation, saturationInv));
		}

		void rungekutteSolver(DspFloatType input, DspFloatType * state)
		{
			int i;
			DspFloatType deriv1[4], deriv2[4], deriv3[4], deriv4[4], tempState[4];

			calculateDerivatives(input, deriv1, state);

			for (i = 0; i < 4; i++)
				tempState[i] = state[i] + 0.5 * stepSize * deriv1[i];

			calculateDerivatives(input, deriv2, tempState);

			for (i = 0; i < 4; i++)
				tempState[i] = state[i] + 0.5 * stepSize * deriv2[i];

			calculateDerivatives(input, deriv3, tempState);

			for (i = 0; i < 4; i++)
				tempState[i] = state[i] + stepSize * deriv3[i];

			calculateDerivatives(input, deriv4, tempState);

			for (i = 0; i < 4; i++)
				state[i] += (1.0 / 6.0) * stepSize * (deriv1[i] + 2.0 * deriv2[i] + 2.0 * deriv3[i] + deriv4[i]);
		}

		DspFloatType state[4];
		DspFloatType saturation, saturationInv;
		int oversampleFactor;
		DspFloatType stepSize;

	};
}
