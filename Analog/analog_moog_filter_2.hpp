#pragma once

#include <cmath>
#include "GenericSoundObject.hpp"

namespace Analog::MoogFilters::MoogFilter2
{
	template<typename DSP>
    struct MoogFilter2 : public GSSoundProcessor<DSP>
    {
        // Moog 24 dB/oct resonant lowpass VCF
        // References: CSound source code, Stilson/Smith CCRMA paper.
        // Modified by paul.kellett@maxim.abel.co.uk July 2000

        DSP f, p, q;             //filter coefficients
        DSP b0, b1, b2, b3, b4;  //filter buffers (beware denormals!)
        DSP t1, t2;              //temporary buffers
        DSP fs,fc,res;

        // Set coefficients given frequency & resonance [0.0...1.0]
        MoogFilter2(DSP sr, DSP cutoff, DSP r) : GSSoundProcessor<DSP>()
        {
            fs = sr;
            fc = cutoff/sr;
            res = r;
            calc();
            b0=b1=b2=b3=b4=0;
        }
        void calc()
        {
            q = 1.0f - fc;
            p = fc + 0.8f * fc * q;
            f = p + p - 1.0f;
            q = res * (1.0f + 0.5f * q * (1.0f - q + 5.6f * q * q));
        }
        void setCutoff(DSP f) { fc = f/fs; }
        void setResonance(DSP r) { res = r; }
        enum
        {
            PORT_CUTOFF,
            PORT_RESONANCE,			
        };
        void setPort(int port, DSP v)
        {
            switch (port)
            {
            case PORT_CUTOFF:
                setCutoff(v);
                break;
            case PORT_RESONANCE:
                setResonance(v);
                break;		
            }
        }
        DSP Tick(DSP I, DSP A=1, DSP X = 0, DSP Y=0)
        {
            Undenormal denormals;
            calc();
            DSP in = I - q*b4;       
            t1 = b1; //std::tanh(b1);  
            b1 = (in + b0) * p - b1 * f;
            t2 = b2; //std::tanh(b2);  
            b2 = (b1 + t1) * p - b2 * f;
            t1 = b3; //std::tanh(b3); 
            b3 = (b2 + t2) * p - b3 * f;
            b4 = (b3 + t1) * p - b4 * f;
            b4 = b4 - b4 * b4 * b4 * 0.166667f;
            b0 = in;
            return b4;
        }
        void ProcessSIMD(size_t n, DSP * I, DSP * out)
        {
			Undenormal denormals;
            calc();
            #pragma omp simd aligned(in,out)
            for(size_t i = 0; i < n; i++)
            {
				DSP in = I[i] - q*b4;       
				t1 = b1; //std::tanh(b1);  
				b1 = (in + b0) * p - b1 * f;
				t2 = b2; //std::tanh(b2);  
				b2 = (b1 + t1) * p - b2 * f;
				t1 = b3; //std::tanh(b3); 
				b3 = (b2 + t2) * p - b3 * f;
				b4 = (b3 + t1) * p - b4 * f;
				b4 = b4 - b4 * b4 * b4 * 0.166667f;
				b0 = in;
				out[i] = b4;
			}
		}
		void ProcessBlock(size_t n, DSP * in, DSP * out) {
			ProcessSIMD(n,in,out);
		}
		void ProcessInplace(size_t n, DSP * buffer) {
			ProcessSIMD(n,buffer,buffer);
		}
    };
}
