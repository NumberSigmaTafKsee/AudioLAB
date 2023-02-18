#pragma once

#include "DaisySP.hpp"
#include <functional>

namespace DaisySP::Control
{
    struct Line : public GeneratorProcessorPlugin<daisysp::Line>
    {
        uint8_t finished = 0;
        std::function<void (Line *)> callback;// = [](Line * p){};
        Line(DspFloatType sampleRate=44100.0) : GeneratorProcessorPlugin<daisysp::Line>()
        {
            this->Init(sampleRate);
            callback = [](Line * p){};
        }
        enum {
            PORT_RESET,
        };
        void setPort(int port, DspFloatType v) {
            if(port == PORT_RESET) finished = 0;
        }
        DspFloatType Tick(DspFloatType I=1, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {
            DspFloatType out = A*this->Process(&finished);
            if(finished == 1) callback(this);
            return out;
        }
    };
}