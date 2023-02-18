#pragma once

#include <cmath>
#include "SoundObject.hpp"

namespace Analog::Oscillators::DPW
{
    struct DPWSaw : public OscillatorProcessor
    {
        DspFloatType freq,fs,inc;
        DspFloatType phase,lastPhase;
        DspFloatType lastValue,position;
        DspFloatType scaleFactor;
        
        DPWSaw(DspFloatType sampleRate=44100.0) : OscillatorProcessor()
        {
            freq = 440.0f;
            fs   = sampleRate;
            inc  = freq/fs;
            lastValue = phase = lastPhase = position = 0.0f;
            scaleFactor = fs / (4.0f * freq);
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
            if(port == PORT_FREQ) setFrequency(v);
            else printf("No port %d\n", port);
        }
        DspFloatType Tick(DspFloatType I=1, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1)
        {                                    
            position += phase - lastPhase;
            lastPhase = phase;

            position = fmod(position, 1.0f);

            DspFloatType value = position * 2 - 1;
            value = value * value;
            
            DspFloatType out = scaleFactor * (value - lastValue);
            lastValue = value;

            phase = std::fmod(phase + inc,1.0f);
            return A*out;
        }

        void ProcessSIMD(size_t n, DspFloatType * in, DspFloatType * output)  {
            #pragma omp simd aligned(in,output)
            for(size_t i = 0; i < n; i++)
            {
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

    struct DPWPulse : public OscillatorProcessor
    {
        DspFloatType freq,fs,inc;
        DspFloatType phase,lastPhase;
        DspFloatType lastValueA,lastValueB,position;
        DspFloatType positionA,positionB;
        DspFloatType scaleFactor;
        DspFloatType invSampleRate = 1.0/44100.0;

        DPWPulse(DspFloatType sampleRate=44100.0) : OscillatorProcessor()
        {
            freq = 440.0f;
            fs   = sampleRate;
            inc  = freq/fs;
            invSampleRate = 1.0/fs;
            lastValueA = lastValueB = phase = lastPhase = position = 0.0f;
            positionA = 0.5f;
            positionB = 0.5f;
            scaleFactor = 0.5f * fs /(4.0f * freq);    
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
        DspFloatType Tick(DspFloatType I=1, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {
                            
            positionB += phase - lastPhase;
            lastPhase = phase;

            positionA = std::fmod(positionA, 1.0f);
            positionB = std::fmod(positionB, 1.0f);

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
        void ProcessSIMD(size_t n, DspFloatType * in, DspFloatType * output)  {
            #pragma omp simd aligned(in,output)
            for(size_t i = 0; i < n; i++)
            {
                positionB += phase - lastPhase;
                lastPhase = phase;

                positionA = std::fmod(positionA, 1.0f);
                positionB = std::fmod(positionB, 1.0f);

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
        void ProcessBlock(size_t n, DspFloatType * input, DspFloatType * output) {
            ProcessSIMD(n,input,output);
        }
            
        void ProcessInplace(size_t n, DspFloatType * input) {
            ProcessBlock(n,nullptr,input);
        }
    };

    struct DPWTriangle : public OscillatorProcessor
    {
        DspFloatType freq,fs,inc;
        DspFloatType phase,lastPhase;
        DspFloatType lastValue,position;    
        DspFloatType scaleFactor;
        DspFloatType invSampleRate = 1.0/44100.0;

        DPWTriangle(DspFloatType sampleRate=44100.0) : OscillatorProcessor()
        {
            freq = 440.0f;
            fs   = sampleRate;
            invSampleRate=1.0/fs;
            inc  = freq/fs;
            lastValue = phase = lastPhase = position = 0.0f;
            position = 0.0f;        
            scaleFactor =  fs / (2.0f * freq);
            phase = 0.5;
        }
        void setFrequency(DspFloatType f) {
            freq = f;
            inc  = f/fs;
            scaleFactor =  fs / (2.0f * freq);
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
        DspFloatType Tick(DspFloatType I=1, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1)
        {        
            position += phase - lastPhase;
            lastPhase = phase;
            position = std::fmod(position, 1.0f);                
            DspFloatType out = std::abs(position - 0.5) * 4 - 1;                
            position += freq * invSampleRate;        
            return A*out;
        }
        void ProcessSIMD(size_t n, DspFloatType * input, DspFloatType * output) {
            #pragma omp simd aligned(input,output)
            for(size_t i = 0; i < n; i++) {
                position += phase - lastPhase;
                lastPhase = phase;
                position = std::fmod(position, 1.0f);                
                DspFloatType out = std::abs(position - 0.5) * 4 - 1;                
                position += freq * invSampleRate;        
                output[i] =  out;
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
