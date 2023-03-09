#pragma once

#include "DaisySP.hpp"

namespace DaisySP::Filters
{
    struct SVF : public FilterProcessorPlugin<daisysp::Svf>
    {
        enum {
            LP,
            HP,
            BP,
            NOTCH,
            PEAK,
        };
        int type = LP;

        SVF(DspFloatType sampleRate=44100.0) :FilterProcessorPlugin<daisysp::Svf>()
        {
            this->Init(sampleRate);
        }
        enum {
            PORT_CUTOFF,
            PORT_RES,
            PORT_DRIVE,
            PORT_TYPE,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_CUTOFF: this->SetFreq(v); break;
                case PORT_RES: this->SetRes(v); break;
                case PORT_DRIVE: this->SetDrive(v); break;
                case PORT_TYPE: this->type = (int)v % (PEAK+1); break;
            }
        }
        DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1)
        {
            DspFloatType out;
            this->Process(I);
            switch(type) {
                case LP: out = this->Low(); break;
                case HP: out = this->High(); break;
                case BP: out = this->Band(); break;
                case NOTCH: out = this->Notch(); break;
                case PEAK: out = this->Peak(); break;
            }
            return A*out;
        }
    };
}