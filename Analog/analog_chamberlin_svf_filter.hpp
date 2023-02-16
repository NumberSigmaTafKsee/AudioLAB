#pragma once
#include <cmath>
#include "FX/ClipFunctions.hpp"
#include "GenericSoundProcessor.hpp"

namespace Analog::Filters::SVF
{
    template<typename DSP>
    struct ChamberlinSVF : public GSSoundProcessor<DSP>
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
        DSP x,L,B,H,N,F1,Q1,D1,D2;
        DSP Fc,Fs,R;

        ChamberlinSVF(DSP sr, DSP fc, DSP q) : GSSoundProcessor<DSP> {        
            Fc = fc;
            Fs = sr;
            R  = q;
            Q1 = 1.5*q + 0.5;
            F1 = 2 * std::sin( M_PI * Fc / Fs );            
            x=L=B=H=N=F1=Q1=D1=D2 = 0;
            setCutoff(fc);
            setResonance(q);
        }
        void setCutoff(DSP fc) { Fc = fc; F1 = 2 * std::sin( M_PI * Fc / Fs );}
        void setResonance(DSP r) { R = r; Q1 = 1.0-r; }

        enum {
            PORT_CUTOFF,
            PORT_RESONANCE
        };
        void setPort(int port, DSP v) {
            switch(port) {
                case PORT_CUTOFF: setCutoff(v); break;
                case PORT_RESONANCE: setResonance(v); break;
            }
        }
        DSP Tick(DSP I, DSP A=1, DSP X=0, DSP Y=0)
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
            return std::tanh(A*L);
        }
        void ProcessSIMD(size_t n, DSP * in, DSP * out) {
            Undenormal denormal;            
            #pragma omp simd
            for(size_t i = 0; i < n; i++)
            {            
                DSP I = in[i];
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
                out[i] = std::tanh(A*L);
            }
        }
        void ProcessBlock(size_t n, DSP * in, DSP * out) {
            ProcessSIMD(n,in,out);
        }
        void ProcessInplace(size_t n, DSP * out) {
            ProcessSIMD(n,out);
        }
    };
}