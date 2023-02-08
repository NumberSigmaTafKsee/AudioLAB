#pragma once 

#include <cmath>
#include "Undenormal.hpp"

namespace Analog::Filters::MoogLike
{
    struct MoogLike : public FilterProcessor
    {
        
        enum {
            LOWPASS,
            HIGHPASS
        };
        int type = LOWPASS;
        DspFloatType coef[9];
        DspFloatType d[4];
        DspFloatType omega; //peak freq        
        DspFloatType  fs,res;
        DspFloatType  _in,_out;

        // calculating coefficients:

        DspFloatType k,p,q,g,a;
        DspFloatType a0,a1,a2,a3,a4;
        
        MoogLike(DspFloatType Fc, DspFloatType Q, DspFloatType G, DspFloatType Fs) : FilterProcessor()
        {
            omega = Fc;
            q  = Q;
            fs = Fs;
            g  = G;
            k =p=q=a=a0=a1=a2=a3=a4=0;
        }

        void SetCoefficients(DspFloatType Fc, DspFloatType R)
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
        DspFloatType Tick(DspFloatType I, DspFloatType A = 1, DspFloatType X=0, DspFloatType Y=0)
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
        void ProcessSIMD(size_t n, DspFloatType * in, DspFloatType * out) {
            Undenormal denormal;
            #pragma omp simd
            for(size_t i = 0; i < n; i++)  {
                _in = in[i];
                // per sample:
                _out=coef[0]*_in+d[0];
                d[0]=coef[1]*_in+coef[5]*_out+d[1];
                d[1]=coef[2]*_in+coef[6]*_out+d[2];
                d[2]=coef[3]*_in+coef[7]*_out+d[3];
                d[3]=coef[4]*_in+coef[8]*_out;
                out[i] = _out;
            }
        }
    };
}