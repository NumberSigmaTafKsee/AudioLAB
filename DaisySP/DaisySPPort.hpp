#pragma once

#include "DaisySP.hpp"
#include "Utility/port.h"

namespace DaisySP::Util
{
    struct Port : public FunctionProcessorPlugin<daisysp::Port>
    {
        Port(float htime,DspFloatType sampleRate=44100.0) : FunctionProcessorPlugin<daisysp::Port>()
        {
            this->Init(sampleRate,htime);
        }
        enum {
            PORT_HTIME,            
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_HTIME: this->SetHtime(v); break;
            }
        }
        DspFloatType getPort(int port)
        {
            if(port != PORT_HTIME) return 0;
            return this->GetHtime();
        }
        DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {
            return A*this->Process(I);
        }
    };
}