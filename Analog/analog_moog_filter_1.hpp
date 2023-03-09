#pragma once 
#include <cmath>
#include "GenericSoundObject.hpp"
#include "FX/ClipFunctions.hpp"
#include "Undenormal.hpp"

namespace Analog::MoogFilters::MoogFilter1
{
	template<typename DSP>
    struct MoogFilter1 : public GSSoundProcessor<DSP>
    {
    //Init
    // cutoff = cutoff freq in Hz
    //fs = sampling frequency //(e.g. 44100Hz)
    //res = resonance [0 - 1] //(minimum - maximum)

        DSP f,fs,k,p,scale,r,y1,y2,y3,y4,oldx,oldy1,oldy2,oldy3;
        DSP cutoff,Q;
        DSP x;

        MoogFilter1(DSP sampleRate, DSP cutoff, DSP resonance) : GSSoundProcessor<DSP>() {
                    
            coefficients(sampleRate,cutoff,resonance);
            x=y1=y2=y3=y4=oldx=oldy1=oldy2=oldy3=0;
        }

        void coefficients(DSP sampleRate,DSP frequency, DSP resonance) 
        {
            fs = sampleRate;
            cutoff = frequency;
            Q = resonance;

            f = 2 * cutoff / fs; //[0 - 1]
            k = 3.6*f - 1.6*f*f -1; //(Empirical tuning)
            p = (k+1)*0.5;

            // resonance sucks 
            scale = std::exp((1-p)*1.386249);
            r = resonance*scale;        
            //DSP t=(1.f-p)*1.386249f;
            //DSP t2=12.f+t*t;
            //r = Q*(t2+6.f*t)/(t2-6.f*t);
        }
        void setCutoff(DSP c) {        
            c = clamp(c,0,fs/2);
            coefficients(fs,c,Q);
        }
        void setResonance(DSP res) {
            res = clamp(res,0,1);
            coefficients(fs,cutoff,res);
        }
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
        DSP Tick(DSP input, DSP A=1, DSP X=0, DSP Y=0)
        {
            Undenormal denormal;
            DSP c = cutoff;
            DSP res = Q;
            coefficients(fs,c + 0.5*X*c,Q + 0.5*Y*Q);
            //Loop
            //--Inverted feed back for corner peaking
            x = input - r*y4;                
            
            //Four cascaded onepole filters (bilinear transform)
            y1=x*p + oldx*p - k*y1;        
            y2=y1*p+oldy1*p - k*y2;        
            y3=y2*p+oldy2*p - k*y3;        
            y4=y3*p+oldy3*p - k*y4;        

            coefficients(fs,c,res);

            //Clipper band limited sigmoid
            y4 = y4 - (y4*y4*y4)/6;        
            oldx  = x;
            oldy1 = y1;
            oldy2 = y2;
            oldy3 = y3;

            return A*y4;
        }
        void ProcessSIMD(size_t n, DSP * in, DSP * out)
        {
			Undenormal denormal;
			#pragma omp simd aligned(in,out)
			for(size_t i = 0; i < n; i++)
			{
				DSP c = cutoff;
				DSP res = Q;
				
				//Loop
				//--Inverted feed back for corner peaking
				x = in[i] - r*y4;                
				
				//Four cascaded onepole filters (bilinear transform)
				y1=x*p + oldx*p - k*y1;        
				y2=y1*p+oldy1*p - k*y2;        
				y3=y2*p+oldy2*p - k*y3;        
				y4=y3*p+oldy3*p - k*y4;        
				
				//Clipper band limited sigmoid
				y4 = y4 - (y4*y4*y4)/6;        
				oldx  = x;
				oldy1 = y1;
				oldy2 = y2;
				oldy3 = y3;

				out[i] = y4;
			}
		}
		void ProcessBlock(size_t n, DSP * in, DSP *out) {
			ProcessSIMD(n,in,out);
		}
		void ProcessInplace(size_t n, DSP * buffer) {
			ProcessSIMD(n,buffer,buffer);
		}
    };
}
