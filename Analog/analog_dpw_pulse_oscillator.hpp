#pragma once

#include "GenericSoundObject.hpp"
#include "FX/OnePole.hpp"
#include "FX/Filters.h"


namespace Analog::Oscillators
{
    template<typename DSP>
    struct DPWPulse : public GSSoundProcessor<DSP>
    {
        DSP freq,fs,inc;
        DSP phase,lastPhase;
        DSP lastValueA,lastValueB,position;
        DSP positionA,positionB;
        DSP scaleFactor;
        DSP invSampleRate=1.0/44100.0;
        
        DPWPulse(DSP sampleRate=44100)
        {
            freq = 440.0f;
            fs   = sampleRate;
            invSampleRate = 1.0/fs;
            inc  = freq/fs;
            lastValueA = lastValueB = phase = lastPhase = position = 0.0f;
            positionA = 0.5f;
            positionB = 0.5f;
            scaleFactor = 0.5f * sampleRate /(4.0f * freq);    
            phase = 0.5;
        }
        void setFrequency(DSP f) {
            freq = f;
            inc  = f/fs;
            scaleFactor = 0.5f * fs /(4.0f * freq);    
        }
        void setDuty(DSP d) {
            phase = clamp(d,0.01f,0.99f);
        }
        enum {
            PORT_FREQ,
            PORT_DUTY,
        };
        void setPort(int port, DSP v) {
            switch(port) {
                case PORT_FREQ: setFrequency(v); break;
                case PORT_DUTY: setDuty(v); break;
                default: printf("No port %d\n", port);
            }
        }
        DSP Tick(DSP I=1, DSP A=1, DSP X=0, DSP Y=0) {
                            
            positionB += phase - lastPhase;
            lastPhase = phase;

            positionA = std::fmod(positionA, 1.0f);
            positionB = std::fmod(positionB, 1.0f);

            DSP valueA = positionA * 2.0f - 1.0f;
            DSP valueB = positionB * 2.0f - 1.0f;
            valueA = valueA * valueA;
            valueB = valueB * valueB;
            DSP out = ((valueA - lastValueA) -(valueB - lastValueB)) * scaleFactor;
            lastValueA = valueA;
            lastValueB = valueB;

            positionA += freq * invSampleRate;
            positionB += freq * invSampleRate;

            return out;        
        }
        void ProcessBlock(size_t n, DSP * in, DSP * out) {
            #pragma omp simd
            for(size_t i = 0; i < n; i++)
            {
                positionB += phase - lastPhase;
                lastPhase = phase;

                positionA = std::fmod(positionA, 1.0f);
                positionB = std::fmod(positionB, 1.0f);

                DSP valueA = positionA * 2.0f - 1.0f;
                DSP valueB = positionB * 2.0f - 1.0f;
                valueA = valueA * valueA;
                valueB = valueB * valueB;
                DSP output = ((valueA - lastValueA) -(valueB - lastValueB)) * scaleFactor;
                lastValueA = valueA;
                lastValueB = valueB;

                positionA += freq * invSampleRate;
                positionB += freq * invSampleRate;

                out[i] = output;        
                if(in) out[i] *= in;
            }
        }
        void ProcessInplace(size_t n, DSP * in) {
            ProcessBlock(n,nullptr,in);
        }
    };
}