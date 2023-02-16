#pragma once

#include "DaisySP.hpp"

namespace DaisySP::PhysicalModels
{
    template<int batch_size>
    struct Resonator : public FilterProcessorPlugin<daisysp::ResonatorSvf<batch_size>>
    {
        enum {
            LP,
            BP,
            BPN,
            HP
        };
        int type = LP;
        DspFloatType f,q,g;

        Resonator() : FilterProcessorPlugin<daisysp::ResonatorSvf<batch_size>>()
        {
            this->Init();
            f = 1000.0;
            q = 0.5;
            g = 1.0;
        }
        enum {
            PORT_TYPE,
            PORT_FREQ,
            PORT_Q,
            PORT_GAIN
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_TYPE: type = v; break;
                case PORT_FREQ: f = v; break;
                case PORT_Q: q = v; break;
                case PORT_GAIN: g = v; break;
            }
        }
        DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {
            float out = I;
            switch(type)
            {
                case  HP:  this->Process<daisysp::ResonatorSvf<batch_size>::FilterMode::HIGH_PASS,true>(&f,&q,&g,I,&out); break;
                case  BP:  this->Process<daisysp::ResonatorSvf<batch_size>::FilterMode::BAND_PASS,true>(&f,&q,&g,I,&out); break;
                case  BPN: this->Process<daisysp::ResonatorSvf<batch_size>::FilterMode::BAND_PASS_NORMALIZED,true>(&f,&q,&g,I,&out); break;
                default:   this->Process<daisysp::ResonatorSvf<batch_size>::FilterMode::LOW_PASS,true>(&f,&q,&g,I,&out); break;
            }
            return A*out;
        }
    };
}