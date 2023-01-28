#pragma once

#define GAMMA_H_INC_ALL
#include "Gamma/Gamma.h"
#include "Gamma/Voices.h"


namespace Gamma::FX
{
    struct Chorus : public MonoFXProcessorPlugin<gam::Chorus<DspFloatType>>
    {
        Chorus() : MonoFXProcessorPlugin<gam::Chorus<DspFloatType>>()
        {

        }
        enum {
            PORT_MODULATE,
            PORT_MAXDELAY,
            PORT_DELAY,
            PORT_FBK,
            PORT_FFD,
            PORT_FREQ,
            PORT_DEPTH,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_MODULATE: this->modulate(); break;
                case PORT_MAXDELAY: this->maxDelay(v); break;
                case PORT_DELAY: this->delay(v); break;
                case PORT_FBK: this->fbk(v); break;
                case PORT_FFD: this->ffd(v); break;
                case PORT_FREQ: this->freq(v); break;
                case PORT_DEPTH: this->depth(v); break;
            }
        }
        DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {
            return A*(*this)(I);
        }
        void ProcessBlock(size_t n, float * in, float * out) {
            for(size_t i = 0; i < n; i++) out[i] = Tick(in[i]);
        }
    };

    struct StereoChorus : public StereoFXProcessorPlugin<gam::Chorus<DspFloatType>>
    {
        StereoChorus() : StereoFXProcessorPlugin<gam::Chorus<DspFloatType>>()
        {

        }
        enum {
            PORT_MODULATE,
            PORT_MAXDELAY,
            PORT_DELAY,
            PORT_FBK,
            PORT_FFD,
            PORT_FREQ,
            PORT_DEPTH,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_MODULATE: this->modulate(); break;
                case PORT_MAXDELAY: this->maxDelay(v); break;
                case PORT_DELAY: this->delay(v); break;
                case PORT_FBK: this->fbk(v); break;
                case PORT_FFD: this->ffd(v); break;
                case PORT_FREQ: this->freq(v); break;
                case PORT_DEPTH: this->depth(v); break;
            }
        }
        DspFloatType Tick(DspFloatType iL, DspFloatType iR, DspFloatType &L, DspFloatType &R, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {
            (*this)(iL,iR,L,R);
            L *= A;
            R *= A;
            return 0.5*(L+R);
        }
        void ProcessBlock(size_t n, float ** in, float ** out) {
            for(size_t i = 0; i < n; i++) 
                Tick(in[i][0],in[i][1],out[0][i],out[1][i]);
        }
    };
}