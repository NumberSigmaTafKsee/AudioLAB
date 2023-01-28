#pragma once

#include "DaisySP.hpp"

namespace DaisySP::PhysicalModels
{
    template<int num_voices=4>
    struct PolyPluck : public GeneratorProcessorPlugin<daisysp::PolyPluck<num_voices>>
    {
        float trig=0;
        float note=0;
        PolyPluck(DspFloatType sampleRate=44100.0) : GeneratorProcessorPlugin<daisysp::PolyPluck<num_voices>>()
        {
            this->Init(sampleRate);
        }
    
        enum {
            PORT_DECAY,        
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_DECAY: this->SetDecay(v); break;            
            }
        }    
        DspFloatType Tick(DspFloatType I=1, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {
            return A*this->Process(trig,note);
        }
    };
}