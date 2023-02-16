#pragma once

#include "GenericSoundObject.hpp"
#include "FX/OnePole.hpp"
#include "FX/Filters.h"


namespace Analog::Oscillators
{

    /////////////////////////////////////////////////////////////////////
    // Differential Parablic Wave
    /////////////////////////////////////////////////////////////////////
    tempalte<typename DSP
    struct DPWSaw : public GSSoundProcessor<DSP>
    {
        DSP freq,fs,inc;
        DSP phase,lastPhase;
        DSP lastValue,position;
        DSP scaleFactor;
        

        DPWSaw(DSP sampleRate=44100) : GSSoundProcessor<DSP>()
        {
            freq = 440.0f;
            fs   = sampleRate;
            inc  = freq/fs;
            lastValue = phase = lastPhase = position = 0.0f;
            scaleFactor = sampleRate / (4.0f * freq);
        }
        void setFrequency(DSP f) {
            freq = f;
            inc  = f/fs;
            scaleFactor = fs / (4.0f * freq);
        }
        enum {
            PORT_FREQ,
        };
        void setPort(int port, DSP v) {
            switch(port) {
                case PORT_FREQ: setFrequency(v); break;
            }
        }
        DSP Tick(DSP I=1, DSP A=1, DSP X=1, DSP Y=1)
        {                                    
            position += phase - lastPhase;
            lastPhase = phase;

            position = fmod(position, 1.0f);

            DSP value = position * 2 - 1;
            value = value * value;
            
            DSP out = scaleFactor * (value - lastValue);
            lastValue = value;

            phase = fmod(phase + inc,1.0f);
            return out;
        }   
        void ProcessBlock(size_t n, DSP * in, DSP * out) {
            #pragma omp simd
            for(size_t i = 0; i < n; i++)
            {
                position += phase - lastPhase;
                lastPhase = phase;
                position = std::fmod(position, 1.0f);
                DSP value = position * 2 - 1;
                value = value * value;                
                DSP output = scaleFactor * (value - lastValue);
                lastValue = value;
                phase = std::fmod(phase + inc,1.0f);
                out[i] = output;
                if(in) out[i] *= in[i];
            }
        }
        void ProcessInplace(size_t n, DSP * in) {
            ProcessBlock(n,nullptr,in);
        }
    };
}