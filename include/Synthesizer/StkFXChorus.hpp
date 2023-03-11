#pragma once

#include "StkHeaders.hpp"

namespace Stk::FX
{
    struct Chorus : public MonoFXProcessorPlugin<stk::Chorus>
    {
        Chorus() : MonoFXProcessorPlugin<stk::Chorus>()
        {

        }
        enum {
            PORT_CLEAR,
            PORT_MODDEPTH,
            PORT_MODFREQ,
            PORT_LASTOUT,
        };
        DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1)
        {
            return A*this->tick(I);            
        }        
        void ProcessBlock(size_t n, DspFloatType * in, DspFloatType * out) {
            for(size_t i = 0; i < n; i++) out[i] = Tick(in[i]);
        }
    };
    struct StereoChorus : public StereoFXProcessorPlugin<stk::Chorus>
    {
        StereoChorus() : StereoFXProcessorPlugin<stk::Chorus>()
        {

        }
        enum {
            PORT_CLEAR,
            PORT_MODDEPTH,
            PORT_MODFREQ,
            PORT_LASTOUT,
        };   
        DspFloatType Tick(DspFloatType iL, DspFloatType iR, DspFloatType &L, DspFloatType &R, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1)
        {
            L = A*this->tick(iL,0);
            R = A*this->tick(iR,1);
            return 0.5*(L+R);
        }
        void ProcessBlock(size_t n, DspFloatType ** in, DspFloatType ** out) {
            for(size_t i = 0; i < n; i++) {
                DspFloatType L = out[0][i];
                DspFloatType R = out[1][i];
                Tick(in[0][i],in[1][i],L,R);
                out[0][i] = L;
                out[1][i] = R;
            }
        }
    };
}