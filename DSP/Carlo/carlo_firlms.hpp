#pragma once

namespace Casino::IPP
{
    // warning - do not know if this works yet
    struct FIRLMS32F
    {
        Ipp32f * pDlyLine;
        Ipp8u * pBuffer;
        int  size,nTaps;
        Ipp32f *pTaps;
        IppsFIRLMSState_32f * pState;
    
        FIRLMS32F(size_t n, Ipp32f * taps) {            
            IppStatus status = ippsFIRLMSGetStateSize_32f(n,0,&size);
            checkStatus(status);
            nTaps = n;
            pTaps = ippsMalloc_32f(n);
            pDlyLine = ippsMalloc_32f(2*n);
            pBuffer = ippsMalloc_8u(size);            
            status = ippsFIRLMSInit_32f(&pState,pTaps,nTaps,pDlyLine,0,pBuffer);
            checkStatus(status);
        }
        ~FIRLMS32F() {            
            if(pBuffer) ippsFree(pBuffer);
            if(pDlyLine) ippFree(pDlyLine);            
            if(pTaps) ippsFree(pTaps);
        }        
        void Execute(const Ipp32f * pSrc, Ipp32f * pRef, Ipp32f* pDst, int len, Ipp32f mu) {
            IppStatus status = ippsFIRLMS_32f(pSrc,pRef,pDst,len,mu,pState);
            checkStatus(status);
        }
    };

    void filter(FIRLMS32F& filter,const float* pSrc, float* pRef, float* pDst, int numIters, float mu) {
        filter.Execute(pSrc,pRef,pDst,numIters,mu);
    }
}