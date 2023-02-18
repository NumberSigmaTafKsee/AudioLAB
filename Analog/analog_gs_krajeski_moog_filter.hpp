#pragma once

#include <cmath>
#include <cstdint>
#include "Undenormal.hpp"
#include "GenericSoundObject.hpp"

namespace Analog::Filters::Moog
{
    ///////////////////////////////////////////////////////////////////////////////////////////
    // The Krajtski 5
    ///////////////////////////////////////////////////////////////////////////////////////////
    template<typename T>
    struct KrajeskiMoog : public GSSoundProcessor<T>
    {
        KrajeskiMoog(T sr) : GSSoundProcessor<T>(),sampleRate(sr)
        {
            memset(state, 0, sizeof(state));
            memset(delay, 0, sizeof(delay));

            drive = 1.0;
            gComp = 1.0;

            SetCutoff(1000.0f);
            SetResonance(0.1f);
        }

        virtual ~KrajeskiMoog() { }

        void ProcessSIMD(size_t n, T * samples, T * output)
        {
            Undenormal denormal;
            #pragma omp simd aligned(samples,output)
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
        void ProcessBlock(size_t n, T * in, T * out) {
			ProcessSIMD(n,in,out);
		}
		void ProcessInplace(size_t n, T * out) {
			ProcessSIMD(n,out);
		}
        T Tick(T I, T A=1, T X=1, T Y=1) {
            Undenormal denormal;
            T c = GetCutoff();
            T r = GetResonance();
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
        void SetResonance(T r)
        {
            resonance = r;
            gRes = resonance * (1.0029 + 0.0526 * wc - 0.926 * std::pow(wc, 2) + 0.0218 * std::pow(wc, 3));
        }
        void SetCutoff(T c)
        {
            cutoff = c;
            wc = 2 * M_PI * cutoff / sampleRate;
            g = 0.9892 * wc - 0.4342 * std::pow(wc, 2) + 0.1381 * std::pow(wc, 3) - 0.0202 * std::pow(wc, 4);
        }

        void setDrive(T d) {
            drive = d;
        }

        T GetResonance() { return resonance; }
        T GetCutoff() { return cutoff; }
        
        enum
        {
            PORT_CUTOFF,
            PORT_RESONANCE,
            PORT_DRIVE,
        };
        void setPort(int port, T v)
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
        T state[5];
        T delay[5];
        T wc; // The angular frequency of the cutoff.
        T g; // A derived parameter for the cutoff frequency
        T gRes; // A similar derived parameter for resonance.
        T gComp; // Compensation factor.
        T drive; // A parameter that controls intensity of nonlinearities.
        T sampleRate;
        T resonance;
        T cutoff;        
    };
}
