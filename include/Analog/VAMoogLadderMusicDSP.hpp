#pragma once

namespace Analog::Moog
{
    ///////////////////////////////////////////////////////////////////////////////////////////
    // MusicDSP kaka
    ///////////////////////////////////////////////////////////////////////////////////////////
    class MusicDSPMoog : public LadderFilterBase
    {

    public:

        MusicDSPMoog(DspFloatType sampleRate) : LadderFilterBase(sampleRate)
        {
            memset(stage, 0, sizeof(stage));
            memset(delay, 0, sizeof(delay));
            SetCutoff(1000.0f);
            SetResonance(0.10f);
        }
        virtual ~MusicDSPMoog()
        {

        }

        void ProcessBlock(size_t n, DspFloatType * samples, DspFloatType * output)
        {
            Undenormal denormal;
            #pragma omp simd aligned(samples,output)
            for (uint32_t s = 0; s < n; ++s)
            {
                DspFloatType x = samples[s] - resonance * stage[3];

                // Four cascaded one-pole filters (bilinear transform)
                stage[0] = x * p + delay[0]  * p - k * stage[0];
                stage[1] = stage[0] * p + delay[1] * p - k * stage[1];
                stage[2] = stage[1] * p + delay[2] * p - k * stage[2];
                stage[3] = stage[2] * p + delay[3] * p - k * stage[3];

                // Clipping band-limited sigmoid
                stage[3] -= (stage[3] * stage[3] * stage[3]) / 6.0;

                delay[0] = x;
                delay[1] = stage[0];
                delay[2] = stage[1];
                delay[3] = stage[2];

                output[s] = stage[3];
            }
        }

        DspFloatType Tick(DspFloatType input) {
            Undenormal denormal;
			DspFloatType x = input - resonance * stage[3];

			// Four cascaded one-pole filters (bilinear transform)
			stage[0] = x * p + delay[0]  * p - k * stage[0];
			stage[1] = stage[0] * p + delay[1] * p - k * stage[1];
			stage[2] = stage[1] * p + delay[2] * p - k * stage[2];
			stage[3] = stage[2] * p + delay[3] * p - k * stage[3];

			// Clipping band-limited sigmoid
			stage[3] -= (stage[3] * stage[3] * stage[3]) / 6.0;

			delay[0] = x;
			delay[1] = stage[0];
			delay[2] = stage[1];
			delay[3] = stage[2];

			return stage[3];
		}

        virtual void SetResonance(DspFloatType r) override
        {
            resonance = r * (t2 + 6.0 * t1) / (t2 - 6.0 * t1);
        }

        virtual void SetCutoff(DspFloatType c) override
        {
            cutoff = 2.0 * c / sampleRate;

            p = cutoff * (1.8 - 0.8 * cutoff);
            k = 2.0 * std::sin(cutoff * MOOG_PI * 0.5) - 1.0;
            t1 = (1.0 - p) * 1.386249;
            t2 = 12.0 + t1 * t1;

            SetResonance(resonance);
        }

    private:

        DspFloatType stage[4];
        DspFloatType delay[4];

        DspFloatType p;
        DspFloatType k;
        DspFloatType t1;
        DspFloatType t2;

    };

}
