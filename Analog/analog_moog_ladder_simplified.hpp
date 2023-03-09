#pragma once

namespace Analog::Moog
{
	template<typename DSP>
	class SimplifiedMoog : public LadderFilterBase
	{
	public:

		SimplifiedMoog(DspFloatType sampleRate) : LadderFilterBase(sampleRate)
		{
			// To keep the overall level approximately constant, comp should be set
			// to 0.5 resulting in a 6 dB passband gain decrease at the maximum resonance
			// (compared to a 12 dB decrease in the original Moog model
			gainCompensation = 0.5;

			memset(stage, 0, sizeof(stage));
			memset(stageZ1, 0, sizeof(stageZ1));
			memset(stageTanh, 0, sizeof(stageTanh));

			SetCutoff(1000.0f);
			SetResonance(0.10f);
		}

		virtual ~SimplifiedMoog()
		{

		}

		// This system is nonlinear so we are probably going to create a signal with components that exceed nyquist.
		// To prevent aliasing distortion, we need to oversample this processing chunk. Where do these extra samples
		// come from? Todo! We can use polynomial interpolation to generate the extra samples, but this is expensive.
		// The cheap solution is to zero-stuff the incoming sample buffer.
		// With resampling, numSamples should be 2x the frame size of the existing sample rate.
		// The output of this filter needs to be run through a decimator to return to the original samplerate.
		void Process(size_t n, DspFloatType * samples, DspFloatType * _output)
		{
			Undenormal denormal;
			// Processing still happens at sample rate...
			for (size_t s = 0; s < n; ++s)
			{
				for (int stageIdx = 0; stageIdx < 4; ++stageIdx)
				{
					if (stageIdx)
					{
						input = stage[stageIdx-1];
						stageTanh[stageIdx-1] = tanh(input);
						stage[stageIdx] = (h * stageZ1[stageIdx] + h0 * stageTanh[stageIdx-1]) + (1.0 - g) * (stageIdx != 3 ? stageTanh[stageIdx] : tanh(stageZ1[stageIdx]));
					}
					else
					{
						input = samples[s] - ((4.0 * resonance) * (output - gainCompensation * samples[s]));
						stage[stageIdx] = (h * tanh(input) + h0 * stageZ1[stageIdx]) + (1.0 - g) * stageTanh[stageIdx];
					}

					stageZ1[stageIdx] = stage[stageIdx];
				}

				output = stage[3];
				SNAP_TO_ZERO(output);
				_output[s] = output;
			}
		}

		void Process(size_t n, DspFloatType * samples)
		{
			Undenormal denormal;
			// Processing still happens at sample rate...
			for (int s = 0; s < n; ++s)
			{
				for (int stageIdx = 0; stageIdx < 4; ++stageIdx)
				{
					if (stageIdx)
					{
						input = stage[stageIdx-1];
						stageTanh[stageIdx-1] = tanh(input);
						stage[stageIdx] = (h * stageZ1[stageIdx] + h0 * stageTanh[stageIdx-1]) + (1.0 - g) * (stageIdx != 3 ? stageTanh[stageIdx] : tanh(stageZ1[stageIdx]));
					}
					else
					{
						input = samples[s] - ((4.0 * resonance) * (output - gainCompensation * samples[s]));
						stage[stageIdx] = (h * tanh(input) + h0 * stageZ1[stageIdx]) + (1.0 - g) * stageTanh[stageIdx];
					}

					stageZ1[stageIdx] = stage[stageIdx];
				}

				output = stage[3];
				SNAP_TO_ZERO(output);
				samples[s] = output;
			}
		}

		DspFloatType Tick(DspFloatType input) {
			DspFloatType r;
			Process(1,&input,&r);
			return r;
		}
		virtual void SetResonance(DspFloatType r) override
		{
			resonance = r;
		}

		virtual void SetCutoff(DspFloatType c) override
		{
			cutoff = c;

			// Not being oversampled at the moment... * 2 when functional
			DspFloatType fs2 = sampleRate;

			// Normalized cutoff [0, 1] in radians: ((2*pi) * cutoff / samplerate)
			g = (2 * MOOG_PI) * cutoff / fs2; // feedback coefficient at fs*2 because of doublesampling
			g *= MOOG_PI / 1.3; // correction factor that allows _cutoff to be supplied Hertz

			// FIR part with gain g
			h = g / 1.3;
			h0 = g * 0.3 / 1.3;
		}

	private:

		DspFloatType output;
		DspFloatType lastStage;

		DspFloatType stage[4];
		DspFloatType stageZ1[4];
		DspFloatType stageTanh[3];

		DspFloatType input;
		DspFloatType h;
		DspFloatType h0;
		DspFloatType g;

		DspFloatType gainCompensation;
	};
}
