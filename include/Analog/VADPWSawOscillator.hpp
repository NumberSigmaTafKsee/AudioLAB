#pragma once

#include "SoundObject.hpp"
#include "FX/OnePole.hpp"
#include "FX/Filters.h"


namespace Analog::Oscillators
{

    /////////////////////////////////////////////////////////////////////
    // Differential Parablic Wave
    /////////////////////////////////////////////////////////////////////
    struct DPWSaw : public OscillatorProcessor
    {
        DspFloatType freq,fs,inc;
        DspFloatType phase,lastPhase;
        DspFloatType lastValue,position;
        DspFloatType scaleFactor;
        

        DPWSaw(DspFloatType sampleRate=44100) : OscillatorProcessor()
        {
            freq = 440.0f;
            fs   = sampleRate;
            inc  = freq/fs;
            lastValue = phase = lastPhase = position = 0.0f;
            scaleFactor = sampleRate / (4.0f * freq);
        }
        void setFrequency(DspFloatType f) {
            freq = f;
            inc  = f/fs;
            scaleFactor = fs / (4.0f * freq);
        }
        enum {
            PORT_FREQ,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_FREQ: setFrequency(v); break;
            }
        }
        DspFloatType Tick(DspFloatType I=1, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1)
        {                                    
            position += phase - lastPhase;
            lastPhase = phase;
            position = std::fmod(position, 1.0f);
            DspFloatType value = position * 2 - 1;
            value = value * value;            
            DspFloatType out = scaleFactor * (value - lastValue);
            lastValue = value;
            phase = std::fmod(phase + inc,1.0f);
            return out;
        }   
        void ProcessSIMD(size_t n, DspFloatType * in, DspFloatType * output) {
            #pragma omp simd aligned(in,output)
            for(size_t i = 0; i < n; i++) {
                position += phase - lastPhase;
                lastPhase = phase;
                position = std::fmod(position, 1.0f);
                DspFloatType value = position * 2 - 1;
                value = value * value;                
                DspFloatType out = scaleFactor * (value - lastValue);
                lastValue = value;
                phase = std::fmod(phase + inc,1.0f);
                output[i] = out;
            }
        }
        void ProcessBlock(size_t n, DspFloatType * input, DspFloatType * output) {
            ProcessSIMD(n,input,output);
        }
            
        void ProcessInplace(size_t n, DspFloatType * input) {
            ProcessBlock(n,nullptr,input);
        }
    };
}
