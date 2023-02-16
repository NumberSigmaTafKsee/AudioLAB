#pragma once
#include "StkHeaders.hpp"

namespace Stk::FX
{
    struct TapDelay : public MonoFXProcessorPlugin<stk::TapDelay>
    {        
        size_t numTaps = 0;

        TapDelay(std::vector<unsigned long> taps = std::vector<unsigned long>(1,0),
                 unsigned long maxDelay=4095) 
        : MonoFXProcessorPlugin<stk::TapDelay>()
        {
            numTaps = taps.size();
            this->setTapDelays(taps);
            this->setMaximumDelay(maxDelay);
        }
        enum {
            PORT_MAXDELAY,            
        };
        void setPort(int port, DspFloatType v) {
            if(port == PORT_MAXDELAY) this->setMaximumDelay(v);
        }
        // if callback is set it is called with the StkFrames of each tap
        // returns the average of all taps
        DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1)
        {
            stk::StkFrames outputs(numTaps);
            outputs = this->tick(I,outputs);            
            DspFloatType o = 0;
            for(size_t i = 0; i < numTaps; i++) o += outputs[i];
            return o/(DspFloatType)numTaps;
        }
        void ProcessBlock(size_t n, DspFloatType * input, DspFloatType * output) {
            for(size_t i = 0; i < n; i++) output[i] = Tick(input[i]);
        }
    };

    /*
    struct StereoTapDelay : public StereoFXProcessorPlugin<stk::TapDelay>
    {        
        size_t numTaps = 0;

        StereoTapDelay(std::vector<unsigned long> taps = std::vector<unsigned long>(1,0),
                 unsigned long maxDelay=4095)  : StereoFXProcessorPlugin<stk::TapDelay>()
        {
            numTaps = taps.size();
            this->setTapDelays(taps);
            this->setMaximumDelay(maxDelay);
        }
        enum {
            PORT_MAXDELAY,            
        };
        void setPort(int port, DspFloatType v) {
            if(port == PORT_MAXDELAY) this->setMaximumDelay(v);
        }
        // if callback is set it is called with the StkFrames of each tap
        // returns the average of all taps
        DspFloatType Tick(DspFloatType iL, DspFloatType iR, DspFloatType &L, DspFloatType & R, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1)
        {
            stk::StkFrames outputs(numTaps);
            outputs = this->tick(I,outputs);            
            DspFloatType o = 0;
            L = 0;
            R = 0;
            for(size_t i = 0; i < numTaps; i++) 
            { 
                o += outputs[i];
                if(i % 2 == 0) R += outputs[i];
                else L += outputs[i];
            }
            return o/(DspFloatType)numTaps;
        }
        void ProcessBlock(size_t n, DspFloatType ** input, DspFloatType ** output) {
            for(size_t i = 0; i < n; i++) 
            {
                DspFloatType L,R;
                Tick(input[0][i],input[1][i],L,R);
                output[0][i] = L;
                output[1][i] = R;
            }
        }
    };
    */
}