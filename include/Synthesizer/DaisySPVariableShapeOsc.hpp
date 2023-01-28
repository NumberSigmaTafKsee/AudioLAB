#pragma once

#include "DaisySP.hpp"

namespace DaisySP::Synthesis
{
    struct VariableShapeOsc : public GeneratorProcessorPlugin<daisysp::VariableShapeOscillator>
    {
        VariableShapeOsc(DspFloatType sampleRate=44100.0) : GeneratorProcessorPlugin<daisysp::VariableShapeOscillator>()
        {            
            this->Init(sampleRate);
        }
        enum {
            PORT_FREQ,
            PORT_PW,
            PORT_WAVESHAPE,
            PORT_SYNC,
            PORT_SYNCFREQ,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_FREQ: this->SetFreq(v); break;
                case PORT_PW: this->SetPW(v); break;
                case PORT_WAVESHAPE: this->SetWaveshape(v); break;
                case PORT_SYNC: this->SetSync(v); break;
                case PORT_SYNCFREQ: this->SetSyncFreq(v); break;
            }
        }
        DspFloatType Tick(DspFloatType I=1, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {
            return A*this->Process();
        }
    };
}