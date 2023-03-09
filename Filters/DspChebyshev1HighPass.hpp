#pragma once

#include "IIRDspFilters.hpp"
#include "DspFilters/ChebyshevI.h"

namespace Filters
{
    struct ChebyshevIHighPassFilter : public FilterBase
    {
        Dsp::ChebyshevI::HighPass<32> prototype;
        

        ChebyshevIHighPassFilter(size_t Order, DspFloatType ripdb, DspFloatType Fc, DspFloatType Fs)
        {            
            order = Order;
            sr    = Fs;
            ripple = ripdb;
            setCutoff(Fc);            
        }
        void setCutoff(DspFloatType Fc) {
            if(Fc < 0 || Fc >= sr/2.0) return;                         
            fc = Fc;                        
            prototype.setup(order,sr,fc,ripple);            
            int total = prototype.getNumStages();
            BiquadSOS sos;
            sos.resize(total);            
            for(size_t i = 0; i < total; i++) {
                Dsp::Cascade::Stage s = prototype[i];
                sos[i].z[0] = s.m_b0;
                sos[i].z[1] = s.m_b1;
                sos[i].z[2] = s.m_b2;
                sos[i].p[0] = s.m_a1;
                sos[i].p[1] = s.m_a2;
                sos[i].p[2] = 0;
            }
            setCoefficients(sos);
        }         
        void setQ(DspFloatType q) {
            // not used yet
        }
        void setRipple(DspFloatType r) {
            ripple = r;
            setCutoff(fc);
        }
        enum {
            PORT_ORDER,
            PORT_CUTOFF,            
            PORT_RIPPLE,            
        };
        void setPort(int port, DspFloatType v) {
            switch(port)
            {
                case PORT_ORDER: order = (int)v; break;
                case PORT_CUTOFF: fc = v; break;
                case PORT_RIPPLE: ripple = v; break;                
            }
            setCutoff(fc);
        }        
    };
}    