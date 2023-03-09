#pragma once

namespace Casino::IPP
{
        template<typename T>
    void FIRMRGetSize(int tapsLen, int upFactor, int downFactor, int * pSpecSize, int *pBufSize)
    {
        IppDataType dType = GetDataType<T>();
        IppStatus status = ippsFIRMRGetSize(tapsLen,upFactor,downFactor,dType,pSpecSize,pBufSize);
        checkStatus(status);
    }
    template<typename T1>
    void FIRMRInit(const T1* pTaps, int tapsLen, int upFactor, int upPhase, int downFact, int downPhase, void* pSpec)
    {
        assert(1==0);
    }
    template<>
    void FIRMRInit<float>(const float* pTaps, int tapsLen, int upFactor, int upPhase,int downFactor, int downPhase, void * pSpec)
    {        
        IppStatus status = ippsFIRMRInit_32f(pTaps,tapsLen,upFactor,upPhase,downFactor,downPhase,(IppsFIRSpec_32f*)pSpec);
        checkStatus(status);
    }
    template<>
    void FIRMRInit<Ipp64f>(const Ipp64f* pTaps, int tapsLen, int upFactor, int upPhase, int downFactor, int downPhase, void * pSpec)
    {        
        IppStatus status = ippsFIRMRInit_64f(pTaps,tapsLen,upFactor,upPhase,downFactor,downPhase,(IppsFIRSpec_64f*)pSpec);
        checkStatus(status);
    }

    template<typename T>
    void FIRMR_(const T * pSrc, T * pDst, int numIters, void* pSpec, const T* pDlySrc, T* pDlyDst, Ipp8u * pBuffer)
    {
        assert(1==0);
    }
    template<>
    void FIRMR_<float>(const float * pSrc, float * pDst, int numIters, void * pSpec, const float* pDlySrc, float* pDlyDst, Ipp8u * pBuffer)
    {
        IppStatus status = ippsFIRMR_32f(pSrc,pDst,numIters,(IppsFIRSpec_32f*)pSpec,pDlySrc,pDlyDst,pBuffer);
        checkStatus(status);
    }
    template<>
    void FIRMR_<Ipp64f>(const Ipp64f * pSrc, Ipp64f * pDst, int numIters, void * pSpec, const Ipp64f* pDlySrc, Ipp64f* pDlyDst, Ipp8u * pBuffer)
    {
        IppStatus status = ippsFIRMR_64f(pSrc,pDst,numIters,(IppsFIRSpec_64f*)pSpec,pDlySrc,pDlyDst,pBuffer);
        checkStatus(status);
    }


    template<typename T>
    struct FIRMR
    {
        Ipp8u * pSpec;
        Ipp8u * pBuffer;
        size_t  size;
        T     * pTaps;
        FIRMR(size_t n, int up, int down, int upPhase=0,int downPhase=0) {
            int specSize,bufferSize;
            FIRMRGetSize<T>(n,up,down,&specSize,&bufferSize);            
            pSpec = Malloc<Ipp8u>(specSize);
            pBuffer= Malloc<Ipp8u>(bufferSize);
            //pDlySrc= ippsMalloc_32f(((n + up -1) / up);
            //pDlyDst= ippsMalloc_32f(((n + down -1) / down);
            pTaps   = Malloc<T>(n);
            FIRMRInit<T>(pTaps,n,up,upPhase,down,downPhase,(void*)pSpec);            
        }
        ~FIRMR() {
            if(pSpec) Free(pSpec);
            if(pBuffer) Free(pBuffer);
            //if(pDlySrc) ippFree(pDlySrc);
            //if(pDlyDst) ippFree(pDlyDst);
            if(pTaps) Free(pTaps);
        }
        void Execute(const T* pSrc, T* pDst, int numIters) {
            FIRMR_<T>(pSrc,pDst,numIters,(void*)pSpec,NULL,NULL,pBuffer);            
        }
    };    

    template<typename T>
    void filter(FIRMR<T>& filter,const T* pSrc, T* pDst, int numIters) {
        filter.Execute(pSrc,pDst,numIters);
    }
    
}