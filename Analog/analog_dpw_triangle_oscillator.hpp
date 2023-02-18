#pragma once

#include "GenericSoundObject.hpp"
#include "FX/OnePole.hpp"
#include "FX/Filters.h"


namespace Analog::Oscillators
{
    template<typename DSP>
    struct DPWTriangle : public GSSoundProcessor<DSP>
    {
        DSP freq,fs,inc;
        DSP phase,lastPhase;
        DSP lastValue,position;    
        DSP scaleFactor;
        DSP invSampleRate=1.0/44100.0;        

        DPWTriangle(DSP sampleRate=44100) : GSSoundProcessor<DSP>()
        {
            freq = 440.0f;
            fs   = sampleRate;
            inc  = freq/fs;
            invSampleRate = 1.0/fs;
            lastValue = phase = lastPhase = position = 0.0f;
            position = 0.0f;        
            scaleFactor =  sampleRate / (2.0f * freq);
            phase = 0.5;
        }
        void setFrequency(DSP f) {
            freq = f;
            inc  = f/fs;
            scaleFactor =  fs / (2.0f * freq);
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
        DSP Tick(DSP I=1, DSP A=1, DSP X=1, DSP Y=1)
        {        
            position += phase - lastPhase;
            lastPhase = phase;
            position = fmod(position, 1.0f);                
            DSP out = std::abs(position - 0.5) * 4 - 1;                
            position += freq * invSampleRate;        
            return out;
        }
        void ProcessSIMD(size_t n, DSP * in, DSP * out) {
            #pragma omp simd aligned(in,out)
            for(size_t i = 0; i < n; i++)
            {
                position += phase - lastPhase;
                lastPhase = phase;
                position = std::fmod(position, 1.0f);                
                DSP output = std::abs(position - 0.5) * 4 - 1;                
                position += freq * invSampleRate;        
                out[i] = output;
                if(in) out[i] *= in[i];
            }
        }
        void ProcessBlock(size_t n, DSP * in, DSP * out) {
			ProcessSIMD(n,in,out);
		}
        void ProcessInplace(size_t n, DSP * in) {
            ProcessSIMD(n,nullptr,in);
        }
    };
}
