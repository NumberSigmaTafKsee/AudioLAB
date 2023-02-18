#pragma once

#include "DaisySP.hpp"

namespace DaisySP::Synthesis
{
    struct BLOsc : public GeneratorProcessorPlugin<daisysp::BlOsc>
    {
        BLOsc(DspFloatType sampleRate=44100.0) : GeneratorProcessorPlugin<daisysp::BlOsc>()
        {
            this->Init(sampleRate);
        }
        enum {
            PORT_WAVEFORM,
            PORT_FREQ,
            PORT_AMP,
            PORT_PW,
            PORT_RESET
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_WAVEFORM: this->SetWaveform((uint8_t)v); break;
                case PORT_FREQ: this->SetFreq(v); break;
                case PORT_AMP: this->SetAmp(v); break;
                case PORT_PW: this->SetPw(v); break;
            }
        }
        DspFloatType Tick(DspFloatType I=1, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {            
            return A*this->Process();
        }        
    };
}