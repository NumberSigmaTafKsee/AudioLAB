#pragma once

#include "SoundObject.hpp"

namespace Analog::Filters::AnalogSVF
{
    struct AnalogSVF : public FilterProcessor
    {    
        DspFloatType fc,fs,q,K;
        DspFloatType lp,hp,bp,ubp,shelf,notch,apf,peak;
        DspFloatType z1,z2;
        DspFloatType minC = -1;
        DspFloatType maxC = 1;
        DspFloatType gain = 1.0;

        enum {
            LP,
            HP,
            BP,
            UBP,
            SHELF,
            NOTCH,
            APF,
            PEAK,
        };
        int type = LP;

        AnalogSVF(DspFloatType Fs, DspFloatType Fc, DspFloatType Q) : FilterProcessor() {
            init(Fs,Fc,Q);
        }
        AnalogSVF() = default;
        void init(DspFloatType Fs, DspFloatType Fc, DspFloatType Q)
        {
            fc = Fc;
            fs = Fs;
            q = Q;
            if(q == 0) q = 0.01;            
            K = 1;
            reset();
        }
        void reset() {
			z1 = z2 = 0.0;
			lp = bp = hp = ubp = peak = shelf = notch = apf = 0;
		}
        void setCutoff(DspFloatType f) {
            //if(f < 30 || f >= fs/2) return;                  
            if(f < 0) return;
            if(f >= 0.5) return;
            fc = f;
        }
        void setQ(DspFloatType Q) {   
            //if(Q < 0.5) return;
            if(Q <= 0) return;
            if(Q >= 1) return;
            q = Q;
        }
        enum {
            PORT_CUTOFF,
            PORT_Q,
            PORT_TYPE,
            PORT_GAIN,
            PORT_MINC,
            PORT_MAXC,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_CUTOFF: setCutoff(v); break;
                case PORT_Q: setQ(v); break;
                case PORT_TYPE: type = (int)v; break;
                case PORT_GAIN: gain = pow(10.0,v/20.0); break;
                case PORT_MINC: minC = v; break;
                case PORT_MAXC: maxC = v; break;
                default: printf("No port %d\n",port);
            }
        }
        DspFloatType Tick(DspFloatType I, DspFloatType A = 1, DspFloatType X = 1, DspFloatType Y = 1)
        {         
			/*           
            DspFloatType wd = 2*M_PI*fc;
            DspFloatType T  = 1/fs;
            DspFloatType temp = wd*T/2;
            DspFloatType wa;

            if(fabs(temp) != M_PI/2)
                wa = (2/T)*tan(temp);
            else
                wa = (2/T)*tan(temp*0.995);
            
            wa *= fabs(X);

            DspFloatType g  = wa*T/2;
            DspFloatType xn = A*I;
            DspFloatType qt = fabs(Y)*q;
            if(qt < 0.5) qt = 0.5;
            DspFloatType  R  = 1.0/(2*qt);
                        
            if(xn < minC) xn = minC;
            if(xn > maxC) xn = maxC;
			*/
			
			DspFloatType wd = 2*M_PI*fc;
			DspFloatType T  = 1/fs;
			DspFloatType wa = 2/T*tan(wd/2);
			DspFloatType g  = wa*T/2;
			DspFloatType xn = A*I;
			DspFloatType  R  = q;
			
            hp = (xn - (2*R+g)*z1 - z2) / (1.0 + 2*R*g + g*g);
            bp = g*hp + z1;
            lp = g*bp + z2;                    
            ubp = 2 * R * bp;            
            shelf = xn + 2*K*R*bp;
            notch = xn - 2*R*bp;
            apf   = xn - 4*R*bp;
            peak  = lp - hp;

            // delay feedback
            z1 = g*hp + bp;
            z2 = g*bp + lp;
                    
            DspFloatType out;
            switch(type) {
                case LP: out = lp; break;
                case HP: out = hp; break;
                case BP: out = bp; break;
                case UBP: out = ubp; break;
                case SHELF: out = shelf; break;
                case NOTCH: out = notch; break;
                case APF: out = apf; break;
                case PEAK: out = peak; break;
            }
            return tanh(out);
        }
        
        
        void ProcessSIMD(size_t n, DspFloatType * input, DspFloatType * output)
        {
            #pragma omp simd aligned(input,output)
            for(size_t i = 0; i < n; i++)
            {
                DspFloatType wd = 2*M_PI*fc;
                DspFloatType T  = 1/fs;
                DspFloatType temp = wd*T/2;
                DspFloatType wa;

                if(fabs(temp) != M_PI/2)
                    wa = (2/T)*tan(temp);
                else
                    wa = (2/T)*tan(temp*0.995);
                                
                DspFloatType I  = input[i];
                DspFloatType g  = wa*T/2;
                DspFloatType xn = I;
                DspFloatType qt = q;
                if(qt < 0.5) qt = 0.5;
                DspFloatType  R  = 1.0/(2*qt);
                
                
                if(xn < minC) xn = minC;
                if(xn > maxC) xn = maxC;

                hp = (xn - (2*R+g)*z1 - z2) / (1.0 + 2*R*g + g*g);
                bp = g*hp + z1;
                lp = g*bp + z2;        
                // not sure these work right yet
                ubp = 2 * R * bp;
                // K is Gain
                shelf = xn + 2*K*R*bp;
                notch = xn - 2*R*bp;
                apf   = xn - 4*R*bp;
                peak  = lp - hp;

                // delay feedback
                z1 = g*hp + bp;
                z2 = g*bp + lp;
                        
                
                DspFloatType out;
                switch(type) {
                    case LP: out = lp; break;
                    case HP: out = hp; break;
                    case BP: out = bp; break;
                    case UBP: out = ubp; break;
                    case SHELF: out = shelf; break;
                    case NOTCH: out = notch; break;
                    case APF: out = apf; break;
                    case PEAK: out = peak; break;                
                }   
                output[i] = out;
            }
        }
        void ProcessBlock(size_t n, DspFloatType * input, DspFloatType * output) {
            ProcessSIMD(n,input,output);
        }
        
        void ProcessInplace(size_t n, DspFloatType * input) {
            ProcessBlock(n,input,input);
        }

    };

    
}
