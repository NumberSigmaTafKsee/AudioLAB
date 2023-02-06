#pragma once

#include "DaisySP.hpp"

namespace DaisySP::Synthesis
{
    template<int num_harmonics=16>
    struct HarmonicOsc : public GeneratorProcessorPlugin<daisysp::HarmonicOscillator<num_harmonics>>
    {
        int amp_idx=0;
        HarmonicOsc(DspFloatType freq, DspFloatType * amps, DspFloatType sampleRate=44100.0) : GeneratorProcessorPlugin<daisysp::HarmonicOscillator<num_harmonics>>()
        {
            this->Init(sampleRate);
            this->SetFreq(freq);
            this->SetAmplitudes(amps);
        }
        enum {
            PORT_FREQ,
            PORT_FIRSTIDX,            
            PORT_AMP,        
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_FREQ: this->SetFreq(v); break;
                case PORT_FIRSTIDX: this->SetFirstHarmIdx(v); break;
            }            
        }
        void setPort2(int port, DspFloatType x, DspFloatType y) {
            switch(port) {
                case PORT_AMP: this->SetSingleAmp(x,y);
            }
        }
        DspFloatType Tick(DspFloatType I=1, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {
            return A*this->Process();
        }
    };
}