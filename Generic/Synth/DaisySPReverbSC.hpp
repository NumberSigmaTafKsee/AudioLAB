#pragma once

#include "DaisySP.hpp"

namespace DaisySP::FX
{
    struct ReverbSC : public StereoFXProcessorPlugin<daisysp::ReverbSc>
    {
        ReverbSC(DspFloatType sampleRate=44100.0) : StereoFXProcessorPlugin<daisysp::ReverbSc>()
        {
            this->Init(sampleRate);
        }
        enum {
            PORT_FEEDBACK,
            PORT_LPFREQ,
        };
        void setPort(int port, DspFloatType v) {
            float x = v;
            switch(port) {
                case PORT_FEEDBACK: this->SetFeedback(x); break;
                case PORT_LPFREQ: this->SetLpFreq(x); break;
            }
        }
        DspFloatType Tick(DspFloatType iL, DspFloatType iR, DspFloatType & oL, DspFloatType &oR, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {
            float L,R;
            L = iL;
            R = iR;
            float o1,o2;
            this->Process(L,R,&o1,&o2);
            oL = L;
            oR = R;
            return 0.5*(oL+oR);
        }
        void ProcessBlock(size_t n, DspFloatType ** in, DspFloatType ** out) {
            for(size_t i = 0; i < n; i++) 
            {
                float iL = in[i][0];
                float iR = in[i][1];
                DspFloatType oR,oL;
                Tick(iL,iR,oL,oR);
                out[i][0] = oL;
                out[i][1] = oR;
            }
        }
    };
}