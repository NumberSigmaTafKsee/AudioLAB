#pragma once

#include <cmath>
#include "GenericSoundObject.hpp"

namespace Analog::Oscillators::DPW
{
    template<typename DSP>
    struct DPWSaw : public GSSoundProcessor<DSP>
    {
        DSP freq,fs,inc;
        DSP phase,lastPhase;
        DSP lastValue,position;
        DSP scaleFactor;
        
        DPWSaw(DSP sampleRate=44100.0) : GSSoundProcessor<DSP>()
        {
            freq = 440.0f;
            fs   = sampleRate;
            inc  = freq/fs;
            lastValue = phase = lastPhase = position = 0.0f;
            scaleFactor = fs / (4.0f * freq);
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
            if(port == PORT_FREQ) setFrequency(v);
            else printf("No port %d\n", port);
        }
        DSP Tick(DSP I=1, DSP A=1, DSP X=1, DSP Y=1)
        {                                    
            position += phase - lastPhase;
            lastPhase = phase;

            position = std::fmod(position, 1.0f);

            DSP value = position * 2 - 1;
            value = value * value;
            
            DSP out = scaleFactor * (value - lastValue);
            lastValue = value;

            phase = std::fmod(phase + inc,1.0f);
            return A*out;
        }   
        void ProcessSIMD(size_t n, DSP * in, DSP * out) {
            #pragma omp simd aligned(in,out)
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
        void ProcessBlock(size_t n, DSP * in, DSP * out)
        {
			ProcessSIMD(n,in,out);
		}
        void ProcessInplace(size_t n, DSP * in) {
            ProcessSIMD(n,nullptr,in);
        }
    };

    template<typename DSP>
    struct DPWPulse : public GSSoundProcessor<DSP>
    {
        DSP freq,fs,inc;
        DSP phase,lastPhase;
        DSP lastValueA,lastValueB,position;
        DSP positionA,positionB;
        DSP scaleFactor;
        DSP invSampleRate = 1.0/44100.0;

        DPWPulse(DSP sampleRate=44100.0) : GSSoundProcessor<DSP>()
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
        DSP Tick(DSP I=1, DSP A=1, DSP X=1, DSP Y=1) {
                            
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

            return 2*out;        
        }
        void ProcessSIMD(size_t n, DSP * in, DSP * out) {
            #pragma omp simd aligned(in,out)
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
        void ProcessBlock(size_t n, DSP * in, DSP * out) {
			ProcessSIMD(n,in,out);
		}			
        void ProcessInplace(size_t n, DSP * in) {
            ProcessSIMD(n,nullptr,in);
        }
    };

    template<typename DSP>
    struct DPWTriangle : public GSSoundProcessor<DSP>
    {
        DSP freq,fs,inc;
        DSP phase,lastPhase;
        DSP lastValue,position;    
        DSP scaleFactor;
        DSP invSampleRate = 1.0/44100.0;

        DPWTriangle(DSP sampleRate=44100.0) : GSSoundProcessor<DSP>()
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
            position = std::fmod(position, 1.0f);                
            DSP out = std::abs(position - 0.5) * 4 - 1;                
            position += freq * invSampleRate;        
            return A*out;
        }
        void ProcessSIMD(size_t n, DSP * in, DSP * out) {
            #pragma omp simd
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
