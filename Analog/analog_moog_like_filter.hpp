#pragma once 

#include <cmath>
#include "Undenormal.hpp"

namespace Analog::Filters::MoogLike
{
	template<typename DSP>
    struct MoogLike : public GSSoundProcessor<DSP>
    {
        
        enum {
            LOWPASS,
            HIGHPASS
        };
        int type = LOWPASS;
        DSP coef[9];
        DSP d[4];
        DSP omega; //peak freq        
        DSP  fs,res;
        DSP  _in,_out;

        // calculating coefficients:

        DSP k,p,q,g,a;
        DSP a0,a1,a2,a3,a4;
        
        MoogLike(DSP Fc, DSP Q, DSP G, DSP Fs) : GSSoundProcessor<DSP>()
        {
            omega = Fc;
            q  = Q;
            fs = Fs;
            g  = G;
            k =p=q=a=a0=a1=a2=a3=a4=0;
        }

        void SetCoefficients(DSP Fc, DSP R)
        {
            omega = Fc;
            q     = R;
            k=(4.0*g-3.0)/(g+1.0);
            p=1.0-0.25*k;
            p*=p;
            
            if(type == LOWPASS) {
                // LP:
                a=1.0/(std::tan(0.5*omega)*(1.0+p));
                p=1.0+a;
                q=1.0-a;

                a0= 1.0/(k+p*p*p*p);
                a1= 4.0*(k+p*p*p*q);
                a2= 6.0*(k+p*p*q*q);
                a3= 4.0*(k+p*q*q*q);
                a4= (k+q*q*q*q);
                p = a0*(k+1.0);

                coef[0]=p;
                coef[1]=4.0*p;
                coef[2]=6.0*p;
                coef[3]=4.0*p;
                coef[4]=p;
                coef[5]=-a1*a0;
                coef[6]=-a2*a0;
                coef[7]=-a3*a0;
                coef[8]=-a4*a0;
            }
            else {
                // or HP:
                a=std::tan(0.5*omega)/(1.0+p);
                p=a+1.0;
                q=a-1.0;

                a0=1.0/(p*p*p*p+k);
                a1=4.0*(p*p*p*q-k);
                a2=6.0*(p*p*q*q+k);
                a3=4.0*(p*q*q*q-k);
                a4=    (q*q*q*q+k);
                p=a0*(k+1.0);

                coef[0]=p;
                coef[1]=-4.0*p;
                coef[2]=6.0*p;
                coef[3]=-4.0*p;
                coef[4]=p;
                coef[5]=-a1*a0;
                coef[6]=-a2*a0;
                coef[7]=-a3*a0;
                coef[8]=-a4*a0;
            }
        }
        DSP Tick(DSP I, DSP A = 1, DSP X=0, DSP Y=0)
        {
            Undenormal denormal;
                    
            _in = I;
            // per sample:
            _out=coef[0]*_in+d[0];
            d[0]=coef[1]*_in+coef[5]*_out+d[1];
            d[1]=coef[2]*_in+coef[6]*_out+d[2];
            d[2]=coef[3]*_in+coef[7]*_out+d[3];
            d[3]=coef[4]*_in+coef[8]*_out;
            return _out;
        }
        void ProcessSIMD(size_t n, DSP * input, DSP * output) {
			Undenormal denormal;
			#pragma omp simd aligned(input,output)
			for(size_t i = 0; i < n; i++) {
				_in = input[i];				
				_out=coef[0]*_in+d[0];
				d[0]=coef[1]*_in+coef[5]*_out+d[1];
				d[1]=coef[2]*_in+coef[6]*_out+d[2];
				d[2]=coef[3]*_in+coef[7]*_out+d[3];
				d[3]=coef[4]*_in+coef[8]*_out;
				output[i] = _out;
			}
		}
		void ProcessBlock(size_t n, DSP * input, DSP * output) {
			ProcessSIMD(n,input,output);
		}
		void ProcessBlock(size_t n, DSP * input) {
			ProcessSIMD(n,input,input);
		}
    };
}
