#pragma once

#include "DaisySP.hpp"

namespace DaisySP::Synthesis
{
    struct Oscillator : public GeneratorProcessorPlugin<daisysp::Oscillator>
    {
        Oscillator(DspFloatType sampleRate=44100.0) : GeneratorProcessorPlugin<daisysp::Oscillator>()
        {
            this->Init(sampleRate);
        }
        enum {
            PORT_FREQ,
            PORT_AMP,
            PORT_WAVEFORM,
            PORT_PW,
            PORT_EOR,
            PORT_EOC,
            PORT_RISING,
            PORT_FALLING,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_FREQ: this->SetFreq(v); break;
                case PORT_AMP: this->SetAmp(v); break;
                case PORT_WAVEFORM: this->SetWaveform((uint8_t)v); break;
                case PORT_PW: this->SetPw(v); break;
            }
        }
        DspFloatType getPort(int port)
        {
            switch(port) {
                case PORT_EOR: return this->IsEOR(); 
                case PORT_EOC: return this->IsEOC();
                case PORT_RISING: return this->IsRising();
                case PORT_FALLING: return this->IsFalling();
            }
            return 0;
        }
        
        DspFloatType Tick(DspFloatType I=1, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {
            return A*this->Process();
        }
    };
}