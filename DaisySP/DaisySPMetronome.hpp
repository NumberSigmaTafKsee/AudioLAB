#pragma once

#include "DaisySP.hpp"
#include "Utility/metro.h"

namespace DaisySP::Util
{
    struct Metronome : public GeneratorProcessorPlugin<daisysp::Metro>
    {
        std::function<void (Metronome*)> callback;// = [](Metronome* ptr){};
        Metronome(DspFloatType freq, DspFloatType sampleRate=44100.0) : GeneratorProcessorPlugin<daisysp::Metro>()
        {
            this->Init(freq,sampleRate);
            callback = [](Metronome* ptr){};
        }
        enum {
            PORT_RESET,
            PORT_FREQ,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_RESET: this->Reset(); break;
                case PORT_FREQ: this->SetFreq(v); break;
            }
        }
        DspFloatType Tick(DspFloatType I=1, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {
            if(this->Process()) callback(this);
            return A*I;
        }
    };
}