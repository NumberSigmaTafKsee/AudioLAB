#pragma once

#include <cmath>
#include "Undenormal.hpp"
#include "GenericSoundObject.hpp"

namespace Analog::Filters::Moog::MoogFilter
{
    // Classes definition
    template<typename DSP>
    struct MoogFilter : public GSSoundProcessor<DSP>
    {
        DSP s0, s1, s2, s3;
        DSP zi;
        DSP xx, y0, y1, y2, y3;
        DSP  fc,q,sr;
        DSP dbGain;
        
        MoogFilter(DSP Fs) : GSSoundProcessor<DSP>()
        {
            sr = Fs;
            fc = 440.0f/sr;
            q  = 0.5;
            s0 = s1 = s2 = s3 = zi = xx = y0 = y1 = y2 =y3 = 0.0;
            dbGain = 0;
        }
        void setCutoff(DSP f) {
            fc = f/sr;
        }
        void setResonance(DSP r) {
            q = r;
        }
        void setGain(DSP db) {
            dbGain = db;
        }
        enum
        {
            PORT_CUTOFF,
            PORT_RESONANCE,
			PORT_DBGAIN,
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
            case PORT_DBGAIN:
                setGain(v);
                break;
            }
        }
        DSP Tick(DSP sample, DSP A=1, DSP X=1, DSP Y=1)
        {
            Undenormal denormals;
            
            DSP f = tan(M_PI * fc/sr);
            DSP r = (DSP)(40.0/9.0) * q;
            
            DSP volume = std::pow(10.0, (dbGain/20.0));

            const DSP input=sample;
            
            // input with half delay, for non-linearities
            const DSP ih = (DSP)0.5 * (input + zi);

            // evaluate the non-linear gains
            // optimization of function calls 
            const DSP in0 = ih - r * s3;
            const DSP a0 = in0 * in0;
            const DSP t0 = f * (((a0 + 105)*a0 + 945) / ((15*a0 + 420)*a0 + 945));

            const DSP a1 = s0 * s0;
            const DSP t1 = f * (((a1 + 105)*a1 + 945) / ((15*a1 + 420)*a1 + 945));

            const DSP a2 = s1 * s1;
            const DSP t2 = f * (((a2 + 105)*a2 + 945) / ((15*a2 + 420)*a2 + 945));

            const DSP a3 = s2 * s2;
            const DSP t3 = f * (((a3 + 105)*a3 + 945) / ((15*a3 + 420)*a3 + 945));

            const DSP a4 = s3 * s3;
            const DSP t4 = f * (((a4 + 105)*a4 + 945) / ((15*a4 + 420)*a4 + 945));

            // This formula gives the main result
            DSP y3 = (s3*(1+t3) + s2*t3)*(1+t2);
            y3 = (y3 + t2*t3*s1)*(1+t1);
            y3 = (y3 + t1*t2*t3*(s0+t0*zi));
            y3 = y3 / ((1+t1)*(1+t2)*(1+t3)*(1+t4) + r*t0*t1*t2*t3);

            const DSP xx = t0 * (zi - r*y3);
            const DSP y0 = t1 * (s0 + xx) / (1+t1);			
            const DSP y1 = t2 * (s1 + y0) / (1+t2);
            const DSP y2 = t3 * (s2 + y1) / (1+t3);

            // update state
            s0 += 2 * (xx - y0);
            s1 += 2 * (y0 - y1);
            s2 += 2 * (y1 - y2);
            s3 += 2 * (y2 - t4*y3);

            zi = input;

            return y3*volume;                
        }
		void ProcessSIMD(size_t n, DSP * in, DSP * out)
		{
			Undenormal denormals;
			#pragma omp simd aligned(in,out)
			for(size_t i = 0; i < n; i++)
			{
				const DSP input = in[i];
				DSP f = std::tan(M_PI * fc/sr);
				DSP r = (DSP)(40.0/9.0) * q;
				
				DSP volume = std::pow(10.0, (dbGain/20.0));
				
				// input with half delay, for non-linearities
				const DSP ih = (DSP)0.5 * (input + zi);

				// evaluate the non-linear gains
				// optimization of function calls 
				const DSP in0 = ih - r * s3;
				const DSP a0 = in0 * in0;
				const DSP t0 = f * (((a0 + 105)*a0 + 945) / ((15*a0 + 420)*a0 + 945));

				const DSP a1 = s0 * s0;
				const DSP t1 = f * (((a1 + 105)*a1 + 945) / ((15*a1 + 420)*a1 + 945));

				const DSP a2 = s1 * s1;
				const DSP t2 = f * (((a2 + 105)*a2 + 945) / ((15*a2 + 420)*a2 + 945));

				const DSP a3 = s2 * s2;
				const DSP t3 = f * (((a3 + 105)*a3 + 945) / ((15*a3 + 420)*a3 + 945));

				const DSP a4 = s3 * s3;
				const DSP t4 = f * (((a4 + 105)*a4 + 945) / ((15*a4 + 420)*a4 + 945));

				// This formula gives the main result
				DSP y3 = (s3*(1+t3) + s2*t3)*(1+t2);
				y3 = (y3 + t2*t3*s1)*(1+t1);
				y3 = (y3 + t1*t2*t3*(s0+t0*zi));
				y3 = y3 / ((1+t1)*(1+t2)*(1+t3)*(1+t4) + r*t0*t1*t2*t3);

				const DSP xx = t0 * (zi - r*y3);
				const DSP y0 = t1 * (s0 + xx) / (1+t1);			
				const DSP y1 = t2 * (s1 + y0) / (1+t2);
				const DSP y2 = t3 * (s2 + y1) / (1+t3);

				// update state
				s0 += 2 * (xx - y0);
				s1 += 2 * (y0 - y1);
				s2 += 2 * (y1 - y2);
				s3 += 2 * (y2 - t4*y3);

				zi = input;
				out[i] = y3*volume;
			}
		}
		void ProcessBlock(size_t n, DSP * in, DSP * out) {
			ProcessSIMD(n,in,out);
		}
		void ProcessInplace(size_t n, DSP * in) {
			ProcessSIMD(n,in,in);
		}
        void reset()
        {  	zi = 0;
            s0 = s1 = s2 = s3 = 0;
        }
        
    };
}
