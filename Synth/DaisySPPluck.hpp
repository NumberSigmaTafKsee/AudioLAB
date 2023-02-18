#pragma once

#include "DaisySP.hpp"

namespace DaisySP::PhysicalModels
{
    struct Pluck : public GeneratorProcessorPlugin<daisysp::Pluck>
    {
        float * memory;
        float trig = 0;
        Pluck(size_t n, int mode = daisysp::PLUCK_MODE_RECURSIVE,DspFloatType sampleRate=44100.0) {
            memory = new float[n];
            this->Init(sampleRate,memory,n,mode);
        }
        ~Pluck() {
            if(memory) delete [] memory;
        }
        enum {
            PORT_AMP,
            PORT_FREQ,
            PORT_DECAY,
            PORT_DAMP,
            PORT_MODE, 
            PORT_TRIG,           
        };
        void setPort(int port, DspFloatType v) {
            switch(port)
            {
                case PORT_AMP: this->SetAmp(v); break;
                case PORT_FREQ: this->SetFreq(v); break;
                case PORT_DECAY: this->SetDecay(v); break;
                case PORT_DAMP: this->SetDamp(v); break;
                case PORT_MODE: this->SetMode(v); break;
                case PORT_TRIG: trig = v; break;
            }
        }
        DspFloatType getPort(int port) {
            switch(port) {
                case PORT_TRIG: return trig;
            }
            return 0.0;
        }
        DspFloatType Tick(DspFloatType I=1, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {
            return A*this->Process(trig);
        }
    };
}