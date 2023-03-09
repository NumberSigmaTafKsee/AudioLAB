#pragma once

#include <cmath>
#include "GenericSoundObject.hpp"

namespace Analog::MoogFilters::MoogFilter3
{
	template<typename DSP>
    struct MoogVCF : public GSSoundProcessor<DSP>
    {
        //Init
        DSP fc;
        DSP fs;
        DSP res;
        DSP out1,out2,out3,out4;
        DSP in1,in2,in3,in4;
        
        MoogVCF(DSP sr, DSP Fc, DSP R) : GSSoundProcessor<DSP>()
        {
            fs = sr;
            fc = Fc/sr;
            res= R;
            out1=out2=out3=out4=0;
            in1=in2=in3=in4=0;
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
        DSP Tick(DSP I, DSP A=1, DSP X=0, DSP Y=0) {
            DSP f = fc * 1.16;
            DSP fb = res * (1.0 - 0.15 * f * f);
            DSP input = I;
            input -= out4 * fb;
            input *= 0.35013 * (f*f)*(f*f);
            out1 = input + 0.3 * in1 + (1 - f) * out1; // Pole 1
            in1  = input;
            out2 = out1 + 0.3 * in2 + (1 - f) * out2;  // Pole 2
            in2  = out1;
            out3 = out2 + 0.3 * in3 + (1 - f) * out3;  // Pole 3
            in3  = out2;
            out4 = out3 + 0.3 * in4 + (1 - f) * out4;  // Pole 4
            in4  = out3;
            return out4;
        }
        void ProcessSIMD(size_t n, DSP * in, DSP *out) {
			Undenormal denormals;
			#pragma omp simd aligned(in,out)
			for(size_t i = 0; i < n i++)
			{
				DSP f = fc * 1.16;
				DSP fb = res * (1.0 - 0.15 * f * f);
				DSP input = in[i];
				input -= out4 * fb;
				input *= 0.35013 * (f*f)*(f*f);
				out1 = input + 0.3 * in1 + (1 - f) * out1; // Pole 1
				in1  = input;
				out2 = out1 + 0.3 * in2 + (1 - f) * out2;  // Pole 2
				in2  = out1;
				out3 = out2 + 0.3 * in3 + (1 - f) * out3;  // Pole 3
				in3  = out2;
				out4 = out3 + 0.3 * in4 + (1 - f) * out4;  // Pole 4
				in4  = out3;
				out[i] = out4;
			}
		}
        void ProcessBlock(size_t n, DSP * in, DSP  out) {
			ProcessSIMD(n,in,out);
		}
		void ProcessInplace(size_t n, DSP * buffer) {
			ProcessSIMD(n,buffer,buffer);
		}
    };
}
