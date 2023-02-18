#pragma once

namespace Analog::Moog
{
	///////////////////////////////////////////////////////////////////////////////////////////
	// The new and improved moog filter
	///////////////////////////////////////////////////////////////////////////////////////////

	class ImprovedMoog : public LadderFilterBase
	{
	public:
		const DspFloatType VT = 0.312;

		ImprovedMoog(DspFloatType sampleRate) : LadderFilterBase(sampleRate)
		{
			memset(V, 0, sizeof(V));
			memset(dV, 0, sizeof(dV));
			memset(tV, 0, sizeof(tV));

			drive = 1.0f;

			SetCutoff(1000.0f); // normalized cutoff frequency
			SetResonance(0.1f); // [0, 4]
		}

		virtual ~ImprovedMoog() { }

		void ProcessBlock(size_t n,DspFloatType * samples, DspFloatType * output)
		{
			Undenormal denormal;
			DspFloatType dV0, dV1, dV2, dV3;
			#pragma omp simd aligned(samples,output)
			for (uint32_t i = 0; i < n; i++)
			{
				dV0 = -g * (std::tanh((drive * samples[i] + resonance * V[3]) / (2.0 * VT)) + tV[0]);
				V[0] += (dV0 + dV[0]) / (2.0 * sampleRate);
				dV[0] = dV0;
				tV[0] = tanh(V[0] / (2.0 * VT));

				dV1 = g * (tV[0] - tV[1]);
				V[1] += (dV1 + dV[1]) / (2.0 * sampleRate);
				dV[1] = dV1;
				tV[1] = std::tanh(V[1] / (2.0 * VT));

				dV2 = g * (tV[1] - tV[2]);
				V[2] += (dV2 + dV[2]) / (2.0 * sampleRate);
				dV[2] = dV2;
				tV[2] = std::tanh(V[2] / (2.0 * VT));

				dV3 = g * (tV[2] - tV[3]);
				V[3] += (dV3 + dV[3]) / (2.0 * sampleRate);
				dV[3] = dV3;
				tV[3] = std::tanh(V[3] / (2.0 * VT));

				output[i] = V[3];
			}
		}	
		
		virtual void SetResonance(DspFloatType r) override
		{
			resonance = r;
		}

		virtual void SetCutoff(DspFloatType c) override
		{
			cutoff = c;
			x = (MOOG_PI * cutoff) / sampleRate;
			g = 4.0 * MOOG_PI * VT * cutoff * (1.0 - x) / (1.0 + x);
		}

		DspFloatType Tick(DspFloatType input) {
			Undenormal denormal;
			DspFloatType dV0, dV1, dV2, dV3;
			dV0 = -g * (std::tanh((drive * input + resonance * V[3]) / (2.0 * VT)) + tV[0]);
			V[0] += (dV0 + dV[0]) / (2.0 * sampleRate);
			dV[0] = dV0;
			tV[0] = tanh(V[0] / (2.0 * VT));

			dV1 = g * (tV[0] - tV[1]);
			V[1] += (dV1 + dV[1]) / (2.0 * sampleRate);
			dV[1] = dV1;
			tV[1] = std::tanh(V[1] / (2.0 * VT));

			dV2 = g * (tV[1] - tV[2]);
			V[2] += (dV2 + dV[2]) / (2.0 * sampleRate);
			dV[2] = dV2;
			tV[2] = std::tanh(V[2] / (2.0 * VT));

			dV3 = g * (tV[2] - tV[3]);
			V[3] += (dV3 + dV[3]) / (2.0 * sampleRate);
			dV[3] = dV3;
			tV[3] = std::tanh(V[3] / (2.0 * VT));
		
			return V[3];
		}
	
	private:

		DspFloatType V[4];
		DspFloatType dV[4];
		DspFloatType tV[4];

		DspFloatType x;
		DspFloatType g;
		DspFloatType drive;
	};
}
