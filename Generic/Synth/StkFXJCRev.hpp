#pragma once
#include "StkHeaders.hpp"

namespace Stk::FX
{
    struct JCRev : public MonoFXProcessorPlugin<stk::JCRev>
    {
        JCRev() : MonoFXProcessorPlugin<stk::JCRev>()
        {

        }
        enum {
            PORT_T60,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_T60: this->setT60(v); break;
            }
        }
        DspFloatType Tick(DspFloatType I=0, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1)
        {
            return this->tick(I);
        }
        void ProcessBlock(size_t n, DspFloatType * input, DspFloatType * output) {
            for(size_t i = 0; i < n; i++) output[i] = Tick(input[i]);
        }
    };

    struct StereoJCRev : public StereoFXProcessorPlugin<stk::JCRev>
    {
        StereoJCRev() : StereoFXProcessorPlugin<stk::JCRev>()
        {

        }
        enum {
            PORT_T60,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_T60: this->setT60(v); break;
            }
        }
        DspFloatType Tick(DspFloatType iL, DspFloatType iR, DspFloatType &L, DspFloatType &R, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1)
        {
            L = A*this->tick(iL,0);
            R = A*this->tick(iR,1);
            return 0.5*(L+R);
        }
        void ProcessBlock(size_t n, DspFloatType ** input, DspFloatType ** output) {
            for(size_t i = 0; i < n; i++) 
            {                   
                Tick(input[i][0],input[i][1],output[i][0],output[i][1]);
            }                
        }
    };
}