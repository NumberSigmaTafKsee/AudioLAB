// DspFilters
#pragma once
#include "IIRDspFilters.hpp"
#include "DspFilters/Bessel.h"

namespace Filters
{
    // bessel filter doesn't preserve group delay in digital
    struct BesselBandStopFilter : public FilterBase
    {
        Dsp::Bessel::BandStop<32> prototype;                
        
        BesselBandStopFilter(size_t Order, DspFloatType BW, DspFloatType Fc, DspFloatType Fs)
        {            
            order = Order;
            sr    = Fs;
            bw    = BW;
            setCutoff(Fc);            
        }
        void setCutoff(DspFloatType Fc) {            
            if(Fc < 0 || Fc >= sr/2.0) {
                std::cout << "Cutoff " << Fc << " out of range.\n";
                return;
            }
            fc = Fc;                                    
            BiquadSOS sos;            
            prototype.setup(order,sr,fc,bw);
            int total = prototype.getNumStages();
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
        void setBandWidth(DspFloatType BW) {
            bw = BW;
            setCutoff(fc);
        }
        enum {
            PORT_ORDER,
            PORT_CUTOFF,
            PORT_BANDWIDTH,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) 
            {
                case PORT_ORDER: order = (int)v; break;
                case PORT_CUTOFF: fc = v; break;
                case PORT_BANDWIDTH: bw = v; break;
            }
            setCutoff(fc);
        }        
    };
}
