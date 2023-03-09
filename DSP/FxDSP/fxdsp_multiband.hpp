#pragma once

#include <cstdlib>
#include "fxdsp.hpp"
#include "fxdsp_rbj.hpp"
#include "fxdp_linkwitzreily.hpp"

// Sqrt(2)/2
#define FILT_Q (0.70710681186548)

namespace FXDSP
{
    /*******************************************************************************
    MultibandFilter */

    template<typename T>
    struct MultibandFilter
    {
        LRFilter<T>*   LPA;
        LRFilter<T>*   HPA;
        LRFilter<T>*   LPB;
        LRFilter<T>*   HPB;
        RBJFilter<T>*  APF;
        T       lowCutoff;
        T       highCutoff;
        T       sampleRate;

        MultibandFilter(T lowCutoff, T highCutoff, T sampleRate)
        {            
            this->lowCutoff = lowCutoff;
            this->highCutoff = highCutoff;
            this->sampleRate = sampleRate;
            LPA = new LRFilter<T>(LOWPASS, lowCutoff, FILT_Q, sampleRate);
            HPA = new LRFilter<T>(HIGHPASS, lowCutoff, FILT_Q, sampleRate);
            LPB = new LRFilter<T>(LOWPASS, highCutoff, FILT_Q, sampleRate);
            HPB = new LRFilter<T>(HIGHPASS, highCutoff, FILT_Q, sampleRate);
            APF = new RBJFilter<T>(ALLPASS, sampleRate/2.0, sampleRate);
            APF->setQ(0.5);
        }
        ~MultibandFilter() {
            if(LPA) delete LPA;
            if(HPA) delete HPA;
            if(LPB) delete LPB;
            if(HPB) delete HPB;
            if(APF) delete APF;
        }
        void flush() {
            LPA->flush();
            HPA->flush();
            LPB->flush();
            HPB->flush();
            APF->flush();
        }
        void setLowCutoff(T low)
        {
            lowCutoff = low;
            LPA->setParams(LOWPASS,lowCutoff,FILT_Q);
            HPA->setParams(HIGHPASS,lowCutoff,FILT_Q);
        }
        void setHighCutoff(T hi)
        {
            highCutoff = hi;
            LPA->setParams(LOWPASS,highCutoff,FILT_Q);
            HPA->setParams(HIGHPASS,highCutoff,FILT_Q);
        }
        void update(T low, T high)
        {
            setLowCutoff(low);
            setHighCutoff(high);
        }

        void ProcessBlock(size_t n, T * inBuffer, T * lowOut, T * midOut, T* highOut)
        {
            T tempLow[n];
            T tempHi[n];

            LPA->ProcessBlock(n,inBuffer,tempLow);
            HPA->ProcessBlock(n,inBuffer,tempHigh);
            APF->ProcessBlock(n,tempLow,lowOut);
            LPB->ProcessBlock(n,tempHi,midOut);
            HPB->ProcessBlock(n,tempHi,highOut);
        }
    };
}