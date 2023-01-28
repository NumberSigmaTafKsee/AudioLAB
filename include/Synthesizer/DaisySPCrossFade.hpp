#pragma once

#include "DaisySP.hpp"

namespace DaisySP::Dynamics
{   
    struct CrossFade : public Parameter2ProcessorPlugin<daisysp::CrossFade>
    {
        enum
        {
            CROSSFADE_LIN,
            CROSSFADE_CPOW,
            CROSSFADE_LOG,
            CROSSFADE_EXP,
            CROSSFADE_LAST,
        };
        CrossFade(int curve) : Parameter2ProcessorPlugin<daisysp::CrossFade>()
        {
            this->Init(curve);
        }
        enum {
            PORT_POS,
            PORT_CURVE,            
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_POS: this->SetPos(v); break;
                case PORT_CURVE: this->SetCurve((int)v); break;
            }
        }
        DspFloatType Tick(DspFloatType X, DspFloatType Y) {
            float x = X;
            float y = Y;
            return this->Process(x,y);
        }
    };
}