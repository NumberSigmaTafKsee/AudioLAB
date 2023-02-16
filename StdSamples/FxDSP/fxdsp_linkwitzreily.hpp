#pragma once

#include "fxdsp.hpp"
#include "fxdsp_rbj.hpp"

namespace FXDSP
{
    /* LRFilter ***************************************************************/
    template<typename T>
    struct LRFilter
    {
        RBJFilter<T>*  filterA;
        RBJFilter<T>*  filterB;
        Filter_t    type;
        T       cutoff;
        T       Q;
        T       sampleRate;

        LRFilter(Filter_t type, T cutoff, T Q, T sampleRate)
        {
            
            this->type = type;
            this->cutoff = cutoff;
            this->Q = Q;
            this->sampleRate = sampleRate;
            filterA = new RBJFilter<T>(type, cutoff, sampleRate);
            filterB = new RBJFilter<T>(type, cutoff, sampleRate);
            filterA->setQ(Q);
            filterB->setQ(Q);            
        }
        ~LRFilter()
        {
            if(filterA) delete filterA;
            if(filterB) delete filterB;
        }
        void flush() {
            filterA->flush();
            filterB->flush();
        }
        void setParams(Filter_t type, T cutoff, T Q)
        {
            this->type = type;
            this->cutoff = cutoff;
            this->Q = Q;
            filterA->setParams(type, cutoff, Q);
            filterB->setParams(type, cutoff, Q);
        }
        void ProcessBlock(size_t n, T * inBuffer, T * outBuffer)
        {
            T tempBuffer[n_samples];
            filterA->ProcessBlock(n, inBuffer, tempBuffer);
            filterB->ProcessBlock(n,tempBuffer, outBuffer);            
        }
        T Tick(T I, T A=1, T X=1, T Y=1)
        {
            return 0;
        }
    };
}