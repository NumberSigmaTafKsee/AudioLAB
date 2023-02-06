#pragma once
#include "IIRDspFilters.hpp"
#include "DspFilters/ChebyshevII.h"
namespace Filters
{
    struct ChebyshevIILowShelfFilter : public FilterBase
    {
        Dsp::ChebyshevII::LowShelf<32> prototype;
        
        ChebyshevIILowShelfFilter(size_t Order, DspFloatType G, DspFloatType stopdb, DspFloatType Fc, DspFloatType Fs)
        {            
            order = Order;
            sr    = Fs;
            stop=stopdb;            
            g  = G;
            setCutoff(Fc);            
        }
        void setCutoff(DspFloatType Fc) {
            if(Fc < 0 || Fc >= sr/2.0) return;                         
            fc = Fc;                        
            prototype.setup(order,sr,fc,g,stop);            
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
        void setStopband(DspFloatType s) {
            stop = s;
            setCutoff(fc);
        }
        enum {
            PORT_ORDER,
            PORT_CUTOFF,
            PORT_BANDWIDTH,
            PORT_STOPBAND,
            PORT_GAIN,
        };
        void setPort(int port, DspFloatType v) {
            switch(port)
            {
                case PORT_ORDER: order = (int)v; break;
                case PORT_CUTOFF: fc = v; break;
                case PORT_STOPBAND: stop = v; break;
                case PORT_BANDWIDTH: bw = v; break;
                case PORT_GAIN: g = v; break;
            }
            setCutoff(fc);
        }
    };        
}