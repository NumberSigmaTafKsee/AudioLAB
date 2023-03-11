#pragma once
#include "StkHeaders.hpp"

namespace Stk::FX
{
    struct NRev : public MonoFXProcessorPlugin<stk::NRev>
    {
        NRev() : MonoFXProcessorPlugin<stk::NRev>()
        {

        }
        enum {
            PORT_CLEAR,
            PORT_T60,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_CLEAR: this->clear(); break;
                case PORT_T60: this->setT60(v); break;
            }
        }
        DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1)
        {
            return this->tick(I);
        }
        void ProcessBlock(size_t n, DspFloatType * input, DspFloatType * output)
        {
            for(size_t i = 0; i < n; i++) output[i] = Tick(input[i]);
        }
    };
    struct StereoNRev : public StereoFXProcessorPlugin<stk::NRev>
    {
        StereoNRev() : StereoFXProcessorPlugin<stk::NRev>()
        {

        }
        enum {
            PORT_CLEAR,
            PORT_T60,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_CLEAR: this->clear(); break;
                case PORT_T60: this->setT60(v); break;
            }
        }
        DspFloatType Tick(DspFloatType iL, DspFloatType iR, DspFloatType &L, DspFloatType &R, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1)
        {
            L = A*this->tick(iL,0);
            R = A*this->tick(iR,1);
            return 0.5*(L+R);
        }
        void ProcessBlock(size_t n, DspFloatType ** in, DspFloatType ** out) {
            for(size_t i = 0; i < n; i++) {                
                Tick(in[0][i],in[i][1],out[i][0],out[i][1]);                
            }
        }
    };
}