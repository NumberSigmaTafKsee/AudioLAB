#pragma once

#include <cmath>
#include "GenericSoundObject.hpp"
#include "Undenormal.hpp"


namespace Analog::Filters::Moog
{
    ///////////////////////////////////////////////////////////////////////////////////////////
    // The new and improved moog filter
    ///////////////////////////////////////////////////////////////////////////////////////////
    
    template<typename DSP>
    struct ImprovedMoog : public GSSoundProcessor<DSP>
    {    
        DSP VT = 0.312;

        ImprovedMoog(DSP sr)         
        : GSSoundProcessor<DSP>(), sampleRate(sr)
        {
            memset(V, 0, sizeof(V));
            memset(dV, 0, sizeof(dV));
            memset(tV, 0, sizeof(tV));

            drive = 1.0f;

            SetCutoff(1000.0f); // normalized cutoff frequency
            SetResonance(0.1f); // [0, 4]
        }

        virtual ~ImprovedMoog() { }

        
        void ProcessBlock(uint32_t n, DSP * samples, DSP * output, DSP * A=NULL,DSP * X=NULL,DSP * Y=NULL)
        {
            Undenormal denormal;
            DSP dV0, dV1, dV2, dV3;
            DSP c,r;
            c = cutoff;
            r = resonance;
            #pragma omp simd
            for (uint32_t i = 0; i < n; i++)
            {
                
                if(X != NULL) {
                    SetCutoff(c * fabs(X[i]));
                }
                if(Y != NULL) {
                    SetResonance(r * fabs(Y[i]));
                }
                

                dV0 = -g * (std::tanh((drive * samples[i] + resonance * V[3]) / (2.0 * VT)) + tV[0]);
                V[0] += (dV0 + dV[0]) / (2.0 * sampleRate);
                dV[0] = dV0;
                tV[0] = std::tanh(V[0] / (2.0 * VT));

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

                DSP amp = 1.0;
                if(A != NULL) amp = A[i];
                output[i] = amp*V[3];
            }
            SetCutoff(c);
            SetResonance(r);
        }

        void ProcessInplace(size_t n, DSP * i)
        {
            ProcessBlock(n,i,i);
        }

        void SetResonance(DSP r)
        {
            resonance = r;
        }
        void SetCutoff(DSP c)
        {
            cutoff = c;
            x = (M_PI * cutoff) / sampleRate;
            g = 4.0 * M_PI * VT * cutoff * (1.0 - x) / (1.0 + x);
        }
        DSP GetResonance() { return resonance; }
        DSP GetCutoff() { return cutoff; }

        enum
        {
            PORT_CUTOFF,
            PORT_RESONANCE,
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
            default: printf("No port %d\n",port);
            }
        }
        
        DSP Tick(DSP I, DSP A=1, DSP X=1, DSP Y=1) {
            Undenormal denormal;
            DSP dV0, dV1, dV2, dV3;
            dV0 = -g * (std::tanh((drive * I + resonance * V[3]) / (2.0 * VT)) + tV[0]);
            V[0] += (dV0 + dV[0]) / (2.0 * sampleRate);
            dV[0] = dV0;
            tV[0] = std::tanh(V[0] / (2.0 * VT));

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

            return A*V[3];
        }
    
        DSP V[4];
        DSP dV[4];
        DSP tV[4];

        DSP x;
        DSP g;
        DSP drive;

        DSP cutoff,resonance,sampleRate;
    };
}