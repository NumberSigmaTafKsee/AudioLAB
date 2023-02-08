#pragma once

#include "SoundObject.hpp"
#include "FX/OnePole.hpp"
#include "FX/Filters.h"


namespace Analog::Oscillators
{
    struct DPWPulse
    {
        DspFloatType freq,fs,inc;
        DspFloatType phase,lastPhase;
        DspFloatType lastValueA,lastValueB,position;
        DspFloatType positionA,positionB;
        DspFloatType scaleFactor;
        DspFloatType invSampleRate=1.0/44100.0;
        
        DPWPulse(DspFloatType sampleRate=44100)
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
        void setFrequency(DspFloatType f) {
            freq = f;
            inc  = f/fs;
            scaleFactor = 0.5f * fs /(4.0f * freq);    
        }
        void setDuty(DspFloatType d) {
            phase = clamp(d,0.01f,0.99f);
        }
        enum {
            PORT_FREQ,
            PORT_DUTY,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_FREQ: setFrequency(v); break;
                case PORT_DUTY: setDuty(v); break;
                default: printf("No port %d\n", port);
            }
        }
        DspFloatType Tick(DspFloatType I=1, DspFloatType A=1, DspFloatType X=0, DspFloatType Y=0) {
                            
            positionB += phase - lastPhase;
            lastPhase = phase;

            positionA = fmod(positionA, 1.0f);
            positionB = fmod(positionB, 1.0f);

            DspFloatType valueA = positionA * 2.0f - 1.0f;
            DspFloatType valueB = positionB * 2.0f - 1.0f;
            valueA = valueA * valueA;
            valueB = valueB * valueB;
            DspFloatType out = ((valueA - lastValueA) -(valueB - lastValueB)) * scaleFactor;
            lastValueA = valueA;
            lastValueB = valueB;

            positionA += freq * invSampleRate;
            positionB += freq * invSampleRate;

            return out;        
        }
        void ProcessSIMD(size_t n, DspFloatType * in, DspFloatType * output) {
            #pragma omp simd
            for(size_t i = 0; i < n; i++) {
                positionB += phase - lastPhase;
                lastPhase = phase;

                positionA = fmod(positionA, 1.0f);
                positionB = fmod(positionB, 1.0f);

                DspFloatType valueA = positionA * 2.0f - 1.0f;
                DspFloatType valueB = positionB * 2.0f - 1.0f;
                valueA = valueA * valueA;
                valueB = valueB * valueB;
                DspFloatType out = ((valueA - lastValueA) -(valueB - lastValueB)) * scaleFactor;
                lastValueA = valueA;
                lastValueB = valueB;

                positionA += freq * invSampleRate;
                positionB += freq * invSampleRate;
                output[i] = out;
            }
        }
    };
}