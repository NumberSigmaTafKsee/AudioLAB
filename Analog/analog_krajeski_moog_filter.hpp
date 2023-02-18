#pragma once

#include <cmath>
#include <cstdint>
#include "Undenormal.hpp"
#include "GenericSoundObject.hpp"

namespace Analog::Filters::Moog
{
    ///////////////////////////////////////////////////////////////////////////////////////////
    // The Krajtski 
    ///////////////////////////////////////////////////////////////////////////////////////////
    template<typename DSP>
    struct KrajeskiMoog : public GSSoundProcessor<DSP>
    {
        KrajeskiMoog(DSP sr) : GSSoundProcessor<DSP>(),sampleRate(sr)
        {
            memset(state, 0, sizeof(state));
            memset(delay, 0, sizeof(delay));

            drive = 1.0;
            gComp = 1.0;

            SetCutoff(1000.0f);
            SetResonance(0.1f);
        }

        virtual ~KrajeskiMoog() { }

        void Process(size_t n, DSP * samples, DSP * output)
        {
            Undenormal denormal;
            #pragma omp simd
            for (uint32_t s = 0; s < n; ++s)
            {
                state[0] = std::tanh(drive * (samples[s] - 4 * gRes * (state[4] - gComp * samples[s])));

                for(int i = 0; i < 4; i++)
                {
                    state[i+1] = g * (0.3 / 1.3 * state[i] + 1 / 1.3 * delay[i] - state[i + 1]) + state[i + 1];
                    delay[i] = state[i];
                }
                output[s] = state[4];
            }
        }
        void ProcessSIMD(size_t n, DSP * samples, DSP * output)
        {
            #pragma omp simd aligned(samples,output)
            for(size_t s = 0; s < n; s++)
            {
                state[0] = std::tanh(drive * (samples[s] - 4 * gRes * (state[4] - gComp * samples[s])));

                for(int i = 0; i < 4; i++)
                {
                    state[i+1] = g * (0.3 / 1.3 * state[i] + 1 / 1.3 * delay[i] - state[i + 1]) + state[i + 1];
                    delay[i] = state[i];
                }
                output[s] = state[4];
            }
        }
        void ProcessBlock(size_t n, DSP * samples, DSP * output) {
			ProcessSIMD(n,samples,output);
		}
        void ProcessInplace(size_t n, DSP * samples)
        {
            ProcessSIMD(n,samples,samples);
        }
        
        DSP Tick(DSP I, DSP A=1, DSP X=1, DSP Y=1) {
            Undenormal denormal;
            DSP c = GetCutoff();
            DSP r = GetResonance();
            SetCutoff(c * fabs(X));
            SetResonance(r * fabs(Y));
            state[0] = std::tanh(drive * (I - 4 * gRes * (state[4] - gComp * I)));
            for(int i = 0; i < 4; i++)
            {
                state[i+1] = g * (0.3 / 1.3 * state[i] + 1 / 1.3 * delay[i] - state[i + 1]) + state[i + 1];
                delay[i] = state[i];
            }
            SetCutoff(c);
            SetResonance(r);
            return A * state[4];
        }
        void SetResonance(DSP r)
        {
            resonance = r;
            gRes = resonance * (1.0029 + 0.0526 * wc - 0.926 * std::pow(wc, 2) + 0.0218 * std::pow(wc, 3));
        }
        void SetCutoff(DSP c)
        {
            cutoff = c;
            wc = 2 * M_PI * cutoff / sampleRate;
            g = 0.9892 * wc - 0.4342 * std::pow(wc, 2) + 0.1381 * std::pow(wc, 3) - 0.0202 * std::pow(wc, 4);
        }

        void setDrive(DSP d) {
            drive = d;
        }

        DSP GetResonance() { return resonance; }
        DSP GetCutoff() { return cutoff; }
        
        enum
        {
            PORT_CUTOFF,
            PORT_RESONANCE,
            PORT_DRIVE,
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
            case PORT_DRIVE:
                setDrive(v);
                break;
            }
        }
        DSP state[5];
        DSP delay[5];
        DSP wc; // The angular frequency of the cutoff.
        DSP g; // A derived parameter for the cutoff frequency
        DSP gRes; // A similar derived parameter for resonance.
        DSP gComp; // Compensation factor.
        DSP drive; // A parameter that controls intensity of nonlinearities.
        DSP sampleRate;
        DSP resonance;
        DSP cutoff;        
    };
}
