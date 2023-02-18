#pragma once
#include <cmath>
#include "FX/ClipFunctions.hpp"

namespace Analog::Filters::SVF
{
    struct ChamberlinSVF : public FilterProcessor
    {
        /*
        //Input/Output
        I - input sample
        L - lowpass output sample
        B - bandpass output sample
        H - highpass output sample
        N - notch output sample
        F1 - Frequency control parameter
        Q1 - Q control parameter
        D1 - delay associated with bandpass output
        D2 - delay associated with low-pass output
        */
        DspFloatType x,L,B,H,N,F1,Q1,D1,D2;
        DspFloatType Fc,Fs,R;
        enum {
            LP,
            HP,
            BP,
            NOTCH,
        };
        int Type = LP;
        ChamberlinSVF(DspFloatType sr, DspFloatType fc, DspFloatType q) : FilterProcessor() {        
            Fc = fc;
            Fs = sr;
            R  = q;
            Q1 = 1.5*q + 0.5;
            F1 = 2 * std::sin( M_PI * Fc / Fs );
            x=L=B=H=N=F1=Q1=D1=D2 = 0;
        }
        void setCutoff(DspFloatType fc) { Fc = fc; F1 = 2 * std::sin( M_PI * Fc / Fs );}
        void setResonance(DspFloatType r) { R = r; Q1 = 1.0-r; }

        enum {
            PORT_CUTOFF,
            PORT_RESONANCE
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_CUTOFF: setCutoff(v); break;
                case PORT_RESONANCE: setResonance(v); break;
            }
        }
        DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=0, DspFloatType Y=0)
        {    
            Undenormal denormal;
            x = I;
            // algorithm
            // loop            
            L = D2 + F1 * D1;
            H = I - L - Q1*D1;
            B = F1 * H + D1;
            N = H + L;

            // store delays
            D1 = B;
            D2 = L;

            // outputs
            //L,H,B,N
            return L;
        }
        void ProcessSIMD(size_t n, DspFloatType * in, DspFloatType * output) {
            Undenormal denormal;
            #pragma omp simd aligned(in,output)
            for(size_t i = 0; i < n; i++) {
                x = in[i];
                // algorithm
                // loop
                L = D2 + F1 * D1;
                H = x - L - Q1*D1;
                B = F1 * H + D1;
                N = H + L;

                // store delays
                D1 = B;
                D2 = L;

                // outputs
                //L,H,B,N
                DspFloatType out = L;
                switch(Type) {
                    case LP: 
                    default: out = L; break;
                    case HP: out = H; break;
                    case BP: out = B; break;
                    case NOTCH: out = N; break;
                }
                output[i] = out;
			}
        }
        void ProcessBlock(size_t n, DspFloatType * input, DspFloatType * output) {
            ProcessSIMD(n,input,output);
        }
        void InplaceProcess(size_t n, DspFloatType * input) {
            ProcessSIMD(n,input,input);
        }
    };
}
