// DspFilters
#pragma once
#include "IIRDspFilters.hpp"
#include "DspFilters/Butterworth.h"

namespace Filters
{    
    // todo: get analog prototype and apply Radius/Q
    struct ButterworthLowPassFilterQ : public FilterBase
    {
        Dsp::Butterworth::LowPass<32> prototype;        
        
        ButterworthLowPassFilterQ(size_t Order, DspFloatType Fc, DspFloatType Fs)
        {            
            order = Order;
            sr    = Fs;
            R     = 1.0/10.0;
            setCutoff(Fc);            
        }
        void setCutoff(DspFloatType Fc) {     
            if(Fc < 0 || Fc >= sr/2.0) return;                         
            fc = Fc;                        
            prototype.setup(order,sr,fc);            
            int total = prototype.getNumStages();
            BiquadSOS sos;
            sos.resize(total);      
            std::vector<Dsp::PoleZeroPair> pz = prototype.getPoleZeros();     
            double beta = tan(M_PI*fc/sr) ;
            double k    = 1/beta;
            for(size_t i = 0; i < total; i++) {

                std::complex<double> p1 = pz[i].poles.first;
                std::complex<double> p2 = pz[i].poles.second;
                std::complex<double> z1 = pz[i].zeros.first;
                std::complex<double> z2 = pz[i].zeros.second;

                double b0 = abs(z1*z2);
                double b1 = abs(-z1-z2);
                double b2 = 1.0;
                double a0 = abs(p1*p2);
                double a1 = abs(-p1-p2);
                double a2 = 1.0;
                
                b0 /= a0;
                b1 /= a0;
                b2 /= a0;
                a1 /= a0;
                a2 /= a0;
                a0  = 1.0;

                sos[i].z[0] = b0;
                sos[i].z[1] = k*b1;
                sos[i].z[2] = k*b2;
                sos[i].p[0] = k*a1;
                sos[i].p[1] = k*a2;
                sos[i].p[2] = 0;
            }
            setCoefficients(sos);
        }
        void setResonance(DspFloatType r) {
            
        }
        void setQ(DspFloatType q) {
            
        }
        enum {
            PORT_ORDER,
            PORT_CUTOFF
        };
        void setPort(int port, DspFloatType v) {
            switch(port)
            {
                case PORT_ORDER: order = (int)v; break;
                case PORT_CUTOFF: fc = v; break;
            }
            setCutoff(fc);
        }       
    };
}