#pragma once
#include <cmath>
#include "ClipFunctions.hpp"

namespace Analog::Filters::SVF
{
    struct StateVariableFilter : public FilterProcessor
    {
        /*
        cutoff = cutoff freq in Hz
        fs = sampling frequency //(e.g. 44100Hz)
        f = 2 sin (pi * cutoff / fs) //[approximately]
        q = resonance/bandwidth [0 < q <= 1]  most res: q=1, less: q=0
        low = lowpass output
        high = highpass output
        band = bandpass output
        notch = notch output
        */
        DspFloatType cutoff,scale,fs,low,high,band,notch;
            
        enum {
            LP,
            HP,
            BP,
            NOTCH
        };
        int type = LP;
        
        StateVariableFilter() {
            scale = 0.5;
            cutoff= 1000;
            fs    = 44100.0;
            low=high=band=notch=0;
        }
        StateVariableFilter(DspFloatType Fc, DspFloatType Fs, DspFloatType Q) {
            scale = Q;
            cutoff= Fc;
            fs    = Fs;
            low=high=band=notch=0;
        }
        void setCutoff(DspFloatType F) { cutoff = F; }
        void setResonance(DspFloatType R) { scale = 1.25*(1.0-R); }

        enum {
            PORT_CUTOFF,
            PORT_RESONANCE,
            PORT_LP,
            PORT_HP,
            PORT_BP,
            PORT_NOTCH
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_CUTOFF: setCutoff(v); break;
                case PORT_RESONANCE: setResonance(v); break;
                case PORT_LP: type = LP; break;
                case PORT_HP: type = HP; break;
                case PORT_BP: type = BP; break;
                case PORT_NOTCH: type = NOTCH; break;               
            }
        }
        DspFloatType Tick(DspFloatType I, DspFloatType A = 1, DspFloatType X=0, DspFloatType Y=0)
        {
            Undenormal denormal;
            DspFloatType f     = 2 * std::sin(2 * M_PI * cutoff/fs);        
            //--beginloop
            //I = Distortion::tanhify(I);
            low = low + f * band;
            high = scale * I - low - scale*band;
            band = f * high + band;
            notch = high + low;
            DspFloatType out;
            switch(type) {
                case LP: out = low; break;
                case HP: out = high; break;
                case BP: out = band; break;
                case NOTCH: out = notch; break;
            }
            return out;
        }
        void ProcessSIMD(size_t n, DspFloatType * input, DspFloatType * output) {
            Undenormal denormal;
            #pragma omp simd aligned(input,output)
            for(size_t i = 0; i < n; i++) {
                DspFloatType f = 2 * std::sin(2 * M_PI * cutoff/fs);        
                DspFloatType I = input[i];   
                //--beginloop
                //I = Distortion::tanhify(I);
                low = low + f * band;
                high = scale * I - low - scale*band;
                band = f * high + band;
                notch = high + low;
                DspFloatType out;
                switch(type) {
                    case LP: out = low; break;
                    case HP: out = high; break;
                    case BP: out = band; break;
                    case NOTCH: out = notch; break;
                }
                output[i] = out;
            }
        }
        void ProcessBlock(size_t n, DspFloatType * input, DspFloatType * output) {
            ProcessSIMD(n,input,output);
        }
        void ProcessInplace(size_t n, DspFloatType * input) {
            ProcessSIMD(n,input,input);
        }
    };
}
