#pragma once

#define GAMMA_H_INC_ALL
#include "Gamma/Gamma.h"
#include "Gamma/Voices.h"


namespace Gamma::FX
{
    struct Pan : public StereoFunctionProcessorPlugin<gam::Pan<DspFloatType>>
    {
        Pan() : StereoFunctionProcessorPlugin<gam::Pan<DspFloatType>>()
        {

        }
        enum {
            PORT_PAN,
            PORT_POS,
            PORT_POSL,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_PAN: this->pan(v); break;
                case PORT_POSL: this->posL(v); break;
            }
        }
        void setPort2(int port, DspFloatType a, DspFloatType b) {
            switch(port) {
                case PORT_POS: this->pos(a,b);
            }
        }
        DspFloatType Tick(DspFloatType iL, DspFloatType iR, DspFloatType &L, DspFloatType &R, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {
            (*this)(iL,iR,L,R);
            L *= A;
            R *= A;
            return 0.5*(L+R);
        }
        void ProcessBlock(size_t n, float ** in, float ** out) {
            for(size_t i = 0; i < n; i++) 
            {
                DspFloatType L = out[0][i];
                DspFloatType R = out[1][i];
                Tick(in[0][i],in[1][i],L,R);
                out[0][i] = L;
                out[1][i] = R;
            }                
        }
    };
}