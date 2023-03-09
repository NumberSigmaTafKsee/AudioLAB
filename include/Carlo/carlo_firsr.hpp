#pragma once

namespace Casino::IPP
{
    template<typename T>
    void FIRSRGetSize(int tapsLen, int * pSpecSize, int *pBufSize)
    {
        IppDataType dType = GetDataType<T>();
        IppStatus status = ippsFIRSRGetSize(tapsLen,dType,pSpecSize,pBufSize);
        checkStatus(status);
    }
    template<typename T1>
    void FIRSRInit(const T1* pTaps, int tapsLen, IppAlgType algType, void* pSpec)
    {
        assert(1==0);
    }
    template<>
    void FIRSRInit<float>(const float* pTaps, int tapsLen, IppAlgType algType, void * pSpec)
    {
        IppStatus status = ippsFIRSRInit_32f(pTaps,tapsLen, algType,(IppsFIRSpec_32f*)pSpec);
        checkStatus(status);
    }
    template<>
    void FIRSRInit<Ipp64f>(const Ipp64f* pTaps, int tapsLen, IppAlgType algType, void * pSpec)
    {
        IppStatus status = ippsFIRSRInit_64f(pTaps,tapsLen, algType,(IppsFIRSpec_64f*)pSpec);
        checkStatus(status);
    }

    template<typename T>
    void FIRSR_(const T * pSrc, T * pDst, int numIters, void * pSpect, const T * pDlySrc, T * pDlyDst, Ipp8u* pBuf)
    {
        assert(1==0);
    }
    template<>
    void FIRSR_<float>(const float * pSrc, float * pDst, int numIters, void * pSpect, const float * pDlySrc, float * pDlyDst, Ipp8u* pBuf)
    {
        IppStatus status = ippsFIRSR_32f(pSrc,pDst,numIters,(IppsFIRSpec_32f*)pSpect,pDlySrc,pDlyDst,pBuf);
        checkStatus(status);
    }
    template<>
    void FIRSR_<Ipp64f>(const Ipp64f * pSrc, Ipp64f * pDst, int numIters, void * pSpect, const Ipp64f * pDlySrc, Ipp64f * pDlyDst, Ipp8u* pBuf)
    {
        IppStatus status = ippsFIRSR_64f(pSrc,pDst,numIters,(IppsFIRSpec_64f*)pSpect,pDlySrc,pDlyDst,pBuf);
        checkStatus(status);
    }

    template<typename T>
    struct FIRSR
    {
        Ipp8u * pSpec;
        Ipp8u * pBuffer;
        size_t  size,nTaps;
        T *pTaps;

        enum {
            LP,
            HP,
            BP,
            BS
        };
        int type = LP;

        FIRSR() = default;
        FIRSR(size_t n, T* taps, IppAlgType algType = ippAlgFFT) {
            int specSize,bufferSize;
            FIRSRGetSize<T>(n,&specSize,&bufferSize);            
            pSpec = Malloc<Ipp8u>(specSize);
            pBuffer= Malloc<Ipp8u>(bufferSize);
            //pDlySrc= ippsMalloc_32f(((n + up -1) / up);
            //pDlyDst= ippsMalloc_32f(((n + down -1) / down);
            nTaps   = n;
            pTaps   = Malloc<T>(n);
            memcpy(pTaps,taps,n*sizeof(T));
            FIRSRInit<T>(pTaps,n,algType,(void*)pSpec);            
        }
        FIRSR(int type, size_t n, T fc1, T fc2=0, IppAlgType algType = ippAlgFFT) {
            int specSize,bufferSize;
            FIRSRGetSize<T>(n,&specSize,&bufferSize);            
            pSpec = Malloc<Ipp8u>(specSize);
            pBuffer= Malloc<Ipp8u>(bufferSize);
            //pDlySrc= ippsMalloc_32f(((n + up -1) / up);
            //pDlyDst= ippsMalloc_32f(((n + down -1) / down);
            pTaps   = Malloc<T>(n);
            makeFilter(type,fc1,fc2);
            FIRSRInit<T>(pTaps,n,algType,(void*)pSpec);            
        }
        ~FIRSR() {
            if(pSpec) Free(pSpec);
            if(pBuffer) Free(pBuffer);
            //if(pDlySrc) ippFree(pDlySrc);
            //if(pDlyDst) ippFree(pDlyDst);
            if(pTaps) Free(pTaps);
        }

        // bandpass/bandstop - fc1 = lower, fc2 = upper
        void makeFilter(int filterType, T fc1, T fc2=0) {
            type = filterType;
            switch(type) {
                case LP: makeLowPass(fc1);  break;
                case HP: makeHighPass(fc1);  break;
                case BP: makeBandPass(fc1,fc2);  break;
                case BS: makeBandStop(fc1,fc2);  break;
            }
        }
        void makeLowPass(T fc) {            
            for(size_t i = 0; i < nTaps; i++)
            {   
                float x= fc*(T)(i-nTaps/2);            
                pTaps[i] = fc * sin(M_PI*x)/(M_PI*x);
            }
            WinHamming(pTaps,pTaps,nTaps);
        }
        void makeHighPass(T fc) {            
            for(size_t i = 0; i < nTaps; i++)
            {   
                float x= fc*(float)(i-nTaps/2);            
                pTaps[i] = (sin(M_PI*x)/(M_PI*x)) - fc * sin(M_PI*x)/(M_PI*x);
            }            
            WinHamming(pTaps,pTaps,nTaps);
        }
        void makeBandPass(T fc1, T fc2) {
            for(size_t i = 0; i < nTaps; i++)
            {   
                float xl= fc1*(T)(i-nTaps/2);            
                float xu= fc2*(T)(i-nTaps/2);            
                pTaps[i] = (xu * sin(M_PI*xu)/(M_PI*xu)) - (xl * sin(M_PI*xl)/(M_PI*xl)) ;
            }
            WinHamming(pTaps,pTaps,nTaps);
        }        
        void makeBandStop(T fc1, T fc2) {
            for(size_t i = 0; i < nTaps; i++)
            {   
                float xl= fc1*(T)(i-nTaps/2);            
                float xu= fc2*(T)(i-nTaps/2);            
                pTaps[i] = (xu * sin(M_PI*xu)/(M_PI*xu)) - (xl * sin(M_PI*xl)/(M_PI*xl)) ;
            }
            WinHamming(pTaps,pTaps,nTaps);
        }    
        
        void Execute(const T* pSrc, T* pDst, int numIters) {
            FIRSR_<T>(pSrc,pDst,numIters,(void*)pSpec,NULL,NULL,pBuffer);            
        }
    };

    struct FIRLowPass
    {
        Ipp8u * pSpec;
        Ipp8u * pBuffer;
        size_t  size,nTaps;
        Ipp64f *pTaps;

        FIRLowPass(size_t n, Ipp64f freq, IppAlgType algType = ippAlgFFT)
        {        
            int specSize,bufferSize;
            FIRSRGetSize<Ipp64f>(n,&specSize,&bufferSize);            
            pSpec = Malloc<Ipp8u>(specSize);
            pBuffer= Malloc<Ipp8u>(bufferSize);
            nTaps   = n;
            pTaps   = Malloc<Ipp64f>(n);            
            genLowPass(n,freq);
            FIRSRInit<Ipp64f>(pTaps,n,algType,(void*)pSpec);       
        }
        ~FIRLowPass() {
            if(pSpec) Free(pSpec);
            if(pBuffer) Free(pBuffer);
            if(pTaps) Free(pTaps);
        }
        void genLowPass(size_t n, Ipp64f freq)
        {
            Ipp64f taps[n];            
            IppStatus status = ippsFIRGenLowpass_64f(freq,taps,size,ippWinHamming,ippTrue,pBuffer);
            checkStatus(status);            
        }

        void Execute(const double* pSrc, double* pDst, int numIters) {
            FIRSR_<double>(pSrc,pDst,numIters,(void*)pSpec,NULL,NULL,pBuffer);            
        }    
    };

    struct FIRHighPass
    {
        Ipp8u * pSpec;
        Ipp8u * pBuffer;
        size_t  size,nTaps;
        Ipp64f *pTaps;

        FIRHighPass(size_t n, Ipp64f freq, IppAlgType algType = ippAlgFFT)
        {        
            int specSize,bufferSize;
            FIRSRGetSize<Ipp64f>(n,&specSize,&bufferSize);            
            pSpec = Malloc<Ipp8u>(specSize);
            pBuffer= Malloc<Ipp8u>(bufferSize);
            nTaps   = n;
            pTaps   = Malloc<Ipp64f>(n);            
            genHighPass(n,freq);
            FIRSRInit<Ipp64f>(pTaps,n,algType,(void*)pSpec);       
        }
        ~FIRHighPass() {
            if(pSpec) Free(pSpec);
            if(pBuffer) Free(pBuffer);
            if(pTaps) Free(pTaps);
        }
        void genHighPass(size_t n, Ipp64f freq)
        {
            Ipp64f taps[n];            
            IppStatus status = ippsFIRGenHighpass_64f(freq,taps,size,ippWinHamming,ippTrue,pBuffer);
            checkStatus(status);            
        }

        void Execute(const double* pSrc, double* pDst, int numIters) {
            FIRSR_<double>(pSrc,pDst,numIters,(void*)pSpec,NULL,NULL,pBuffer);            
        }    
    };


    struct FIRBandPass
    {
        Ipp8u * pSpec;
        Ipp8u * pBuffer;
        size_t  size,nTaps;
        Ipp64f *pTaps;

        FIRBandPass(size_t n, Ipp64f lowFreq, Ipp64f highFreq, IppAlgType algType = ippAlgFFT)
        {        
            int specSize,bufferSize;
            FIRSRGetSize<Ipp64f>(n,&specSize,&bufferSize);            
            pSpec = Malloc<Ipp8u>(specSize);
            pBuffer= Malloc<Ipp8u>(bufferSize);
            nTaps   = n;
            pTaps   = Malloc<Ipp64f>(n);            
            genBandPass(n,lowFreq, highFreq);
            FIRSRInit<Ipp64f>(pTaps,n,algType,(void*)pSpec);       
        }
        ~FIRBandPass() {
            if(pSpec) Free(pSpec);
            if(pBuffer) Free(pBuffer);
            if(pTaps) Free(pTaps);
        }
        void genBandPass(size_t n, Ipp64f lowFreq, Ipp64f highFreq)
        {
            Ipp64f taps[n];            
            IppStatus status = ippsFIRGenBandpass_64f(lowFreq,highFreq,taps,size,ippWinHamming,ippTrue,pBuffer);
            checkStatus(status);            
        }

        void Execute(const double* pSrc, double* pDst, int numIters) {
            FIRSR_<double>(pSrc,pDst,numIters,(void*)pSpec,NULL,NULL,pBuffer);            
        }    
    };


    struct FIRBandStop
    {
        Ipp8u * pSpec;
        Ipp8u * pBuffer;
        size_t  size,nTaps;
        Ipp64f *pTaps;

        FIRBandStop(size_t n, Ipp64f lowFreq, Ipp64f highFreq, IppAlgType algType = ippAlgFFT)
        {        
            int specSize,bufferSize;
            FIRSRGetSize<Ipp64f>(n,&specSize,&bufferSize);            
            pSpec = Malloc<Ipp8u>(specSize);
            pBuffer= Malloc<Ipp8u>(bufferSize);
            nTaps   = n;
            pTaps   = Malloc<Ipp64f>(n);            
            genBandStop(n,lowFreq, highFreq);
            FIRSRInit<Ipp64f>(pTaps,n,algType,(void*)pSpec);       
        }
        ~FIRBandStop() {
            if(pSpec) Free(pSpec);
            if(pBuffer) Free(pBuffer);
            if(pTaps) Free(pTaps);
        }
        void genBandStop(size_t n, Ipp64f lowFreq, Ipp64f highFreq)
        {
            Ipp64f taps[n];            
            IppStatus status = ippsFIRGenBandstop_64f(lowFreq,highFreq,taps,size,ippWinHamming,ippTrue,pBuffer);
            checkStatus(status);            
        }
        void Execute(const double* pSrc, double* pDst, int numIters) {
            FIRSR_<double>(pSrc,pDst,numIters,(void*)pSpec,NULL,NULL,pBuffer);            
        }    
    };


    template<typename T>
    void filter(FIRSR<T>& filter,const T* pSrc, T* pDst, int numIters) {
        filter.Execute(pSrc,pDst,numIters);
    }
    void filter(FIRLowPass& filter,const double* pSrc, double* pDst, int numIters) {
        filter.Execute(pSrc,pDst,numIters);
    }
    void filter(FIRHighPass& filter,const double* pSrc, double* pDst, int numIters) {
        filter.Execute(pSrc,pDst,numIters);
    }
    void filter(FIRBandPass& filter,const double* pSrc, double* pDst, int numIters) {
        filter.Execute(pSrc,pDst,numIters);
    }
    void filter(FIRBandStop& filter,const double* pSrc, double* pDst, int numIters) {
        filter.Execute(pSrc,pDst,numIters);
    }
}