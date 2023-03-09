#pragma once

#include "DaisySP.hpp"
#include "Utility/looper.h"

namespace DaisySP::Util
{
    struct Looper : public MonoFXProcessorPlugin<daisysp::Looper>
    {
        std::vector<float> memory;
        size_t  size;
        enum 
        {
            NORMAL,
            ONETIME_DUB,
            REPLACE,
            FRIPPERTRONICS,
        };

        Looper(size_t n) : MonoFXProcessorPlugin<daisysp::Looper>()
        {
            memory.resize(n);            
            this->Init(memory.data(),size);
        }
        ~Looper() {
            
        }

        enum {
            PORT_CLEAR,
            PORT_TRIGRECORD,
            PORT_INCMODE,
            PORT_SETMODE,
            PORT_TOGGLE_REVERSE,
            PORT_SET_REVERSE,
            PORT_TOGGLE_HALFSPEED,
            PORT_SET_HALFSPEED,
            PORT_RECORDING_QUEUED,
            PORT_RECORDING,
            PORT_GETMODE,
            PORT_GET_REVERSE,
            PORT_GET_HALFSPEED,
            PORT_IS_NEAR_BEGINNING,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_CLEAR: this->Clear(); break;
                case PORT_TRIGRECORD: this->TrigRecord(); break;
                case PORT_INCMODE: this->IncrementMode(); break;
                case PORT_SETMODE: this->SetMode((Mode)v); break;
                case PORT_TOGGLE_REVERSE: this->ToggleReverse(); break;
                case PORT_SET_REVERSE: this->SetReverse(v); break;
                case PORT_TOGGLE_HALFSPEED: this->ToggleHalfSpeed(); break;
                case PORT_SET_HALFSPEED: this->SetHalfSpeed(v); break;
            }
        }
        DspFloatType getPort(int port) {
            switch(port) {
                case PORT_RECORDING: return this->Recording();
                case PORT_RECORDING_QUEUED: return this->RecordingQueued();
                case PORT_GETMODE: return (DspFloatType)this->GetMode();
                case PORT_GET_REVERSE: return this->GetReverse();
                case PORT_GET_HALFSPEED: return this->GetHalfSpeed();
                case PORT_IS_NEAR_BEGINNING: return this->IsNearBeginning(); 
            }
            return 0;
        }
        DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {
            DspFloatType out = A*this->Process(I);            
            return out;
        }
         void ProcessBlock(size_t n, DspFloatType * in, DspFloatType * out) {
            std::vector<float> _in(n),_out(n);
            for(size_t i = 0; i < n; i++) _in[i] = in[i];
            for(size_t i = 0; i < n; i++) _out[i] = Tick(_in[i]);
            for(size_t i = 0; i < n; i++) out[i] = _out[i];
        }
    };
}