#pragma once
#include "StkHeaders.hpp"

namespace Stk::FX
{
    // todo: ports
    struct FreeVerb : public StereoFXProcessorPlugin<stk::FreeVerb>
    {
        FreeVerb() : StereoFXProcessorPlugin<stk::FreeVerb>()
        {

        }
        enum {
            PORT_MIX,
            PORT_ROOMSIZE,
            PORT_DAMPING,
            PORT_WIDTH,
            PORT_MODE,
            PORT_CLEAR
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_MIX: this->setEffectMix(v); break;
                case PORT_ROOMSIZE: this->setRoomSize(v); break;
                case PORT_DAMPING: this->setDamping(v); break;
                case PORT_WIDTH: this->setWidth(v); break;
                case PORT_MODE: this->setMode((bool)v); break;
                case PORT_CLEAR: this->clear(); break;
            }
        }
        DspFloatType Tick(DspFloatType iL, DspFloatType iR, DspFloatType &L, DspFloatType &R, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1)
        {
            L = this->tick(iL,iR,0);
            R = this->tick(iL,iR,1);
            return (L+R)*0.5;
        }
        void ProcessBlock(size_t n, DspFloatType ** in, DspFloatType ** out) {
            for(size_t i = 0; i < n; i++) {
                DspFloatType L,R;
                Tick(in[0][i],in[1][i],L,R);
                out[0][i] = L;
                out[1][i] = R;
            }
        }
    };
}