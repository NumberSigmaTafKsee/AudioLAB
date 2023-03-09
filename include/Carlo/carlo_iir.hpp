#pragma once

namespace Casino::IPP
{
    template<typename T>
    void IIRInit(void ** ppState, const T* pTaps, int order, const T* pDlyLine, Ipp8u * pBuf)
    {
        assert(1==0);
    }
    template<>
    void IIRInit<float>(void ** ppState, const float* pTaps, int order, const float* pDlyLine, Ipp8u * pBuf)
    {
        IppStatus status = ippsIIRInit_32f((IppsIIRState_32f**)ppState,pTaps,order,pDlyLine,pBuf);
        checkStatus(status);
    }
    template<>
    void IIRInit<double>(void ** ppState, const double* pTaps, int order, const double* pDlyLine, Ipp8u * pBuf)
    {
        IppStatus status = ippsIIRInit_64f((IppsIIRState_64f**)ppState,pTaps,order,pDlyLine,pBuf);
        checkStatus(status);
    }
    template<typename T>
    void IIRInit_BiQuad(void ** ppState, const T* pTaps, int numBq, const T* pDlyLine, Ipp8u * pBuf)
    {
        assert(1==0);
    }
    template<>
    void IIRInit_BiQuad<float>(void ** ppState, const float* pTaps, int numBq, const float* pDlyLine, Ipp8u * pBuf)
    {
        IppStatus status = ippsIIRInit_BiQuad_32f((IppsIIRState_32f**)ppState,pTaps,numBq,pDlyLine,pBuf);
        checkStatus(status);
    }
    template<>
    void IIRInit_BiQuad<double>(void ** ppState, const double* pTaps, int numBq, const double* pDlyLine, Ipp8u * pBuf)
    {
        IppStatus status = ippsIIRInit_BiQuad_64f((IppsIIRState_64f**)ppState,pTaps,numBq,pDlyLine,pBuf);
        checkStatus(status);
    }
    template<typename T>
    void IIRGetStateSize(int order, int * pBufferSize)
    {
        assert(1==0);
    }
    template<>
    void IIRGetStateSize<float>(int order, int * pBufferSize)
    {
        IppStatus status = ippsIIRGetStateSize_32f(order,pBufferSize);
        checkStatus(status);
    }
    template<>
    void IIRGetStateSize<double>(int order, int * pBufferSize)
    {
        IppStatus status = ippsIIRGetStateSize_64f(order,pBufferSize);
        checkStatus(status);
    }
    template<typename T>
    void IIRGetStateSize_BiQuad(int order, int * pBufferSize)
    {
        assert(1==0);
    }
    template<>
    void IIRGetStateSize_BiQuad<float>(int order, int * pBufferSize)
    {
        IppStatus status = ippsIIRGetStateSize_BiQuad_32f(order,pBufferSize);
        checkStatus(status);
    }
    template<>
    void IIRGetStateSize_BiQuad<double>(int order, int * pBufferSize)
    {
        IppStatus status = ippsIIRGetStateSize_BiQuad_64f(order,pBufferSize);
        checkStatus(status);
    }
    template<typename T>
    void IIRGetDlyLine(void * pState, T * pDlyLine)
    {
        assert(1==0);
    }
    template<>
    void IIRGetDlyLine<float>(void * pState, float * pDlyLine)
    {
        IppStatus status = ippsIIRGetDlyLine_32f((IppsIIRState_32f*)pState,pDlyLine);
        checkStatus(status);
    }
    template<>
    void IIRGetDlyLine<double>(void * pState, double * pDlyLine)
    {
        IppStatus status = ippsIIRGetDlyLine_64f((IppsIIRState_64f*)pState,pDlyLine);
        checkStatus(status);
    }
    template<typename T>
    void IIRSetDlyLine(void * pState, T * pDlyLine)
    {
        assert(1==0);
    }
    template<>
    void IIRSetDlyLine<float>(void * pState, float * pDlyLine)
    {
        IppStatus status = ippsIIRSetDlyLine_32f((IppsIIRState_32f*)pState,pDlyLine);
        checkStatus(status);
    }
    template<>
    void IIRSetDlyLine<double>(void * pState, double * pDlyLine)
    {
        IppStatus status = ippsIIRSetDlyLine_64f((IppsIIRState_64f*)pState,pDlyLine);
        checkStatus(status);
    }
    template<typename T>
    void IIR_(const T* pSrc, T * pDst, int len, void * pState)
    {
        assert(1==0);
    }
    template<>
    void IIR_<float>(const float* pSrc, float * pDst, int len, void * pState)
    {
        IppStatus status = ippsIIR_32f(pSrc,pDst,len,(IppsIIRState_32f*)pState);
        checkStatus(status);
    }
    template<>
    void IIR_<double>(const double* pSrc, double * pDst, int len, void * pState)
    {
        IppStatus status = ippsIIR_64f(pSrc,pDst,len,(IppsIIRState_64f*)pState);
        checkStatus(status);
    }    
    void IIRGenGetBufferSize(int order, int * pBufferSize)
    {
        IppStatus status = ippsIIRGenGetBufferSize(order,pBufferSize);
        checkStatus(status);
    }
    
    void IIRGenLowpassButterworth(double freq, double ripple, int order, double * pTaps, Ipp8u * pBuffer)
    {
        IppStatus status = ippsIIRGenLowpass_64f(freq,ripple,order,pTaps,ippButterworth,pBuffer);
        checkStatus(status);
    }
    void IIRGenHighpassButterworth(double freq, double ripple, int order, double * pTaps, Ipp8u * pBuffer)
    {
        IppStatus status = ippsIIRGenHighpass_64f(freq,ripple,order,pTaps,ippButterworth,pBuffer);
        checkStatus(status);
    }
    void IIRGenLowpassChebyshev1(double freq, double ripple, int order, double * pTaps, Ipp8u * pBuffer)
    {
        IppStatus status = ippsIIRGenLowpass_64f(freq,ripple,order,pTaps,ippChebyshev1,pBuffer);
        checkStatus(status);
    }
    void IIRGenHighpassChebyshev1(double freq, double ripple, int order, double * pTaps, Ipp8u * pBuffer)
    {
        IppStatus status = ippsIIRGenHighpass_64f(freq,ripple,order,pTaps,ippChebyshev1,pBuffer);
        checkStatus(status);
    }


    
    template<typename T>
    struct IIR
    {
        Ipp8u * pBuffer;
        size_t  len;                
        void  * pState;
        T*      dlyLine;
        
        IIR() {
            pBuffer = nullptr;
            len = 0;
            pState = nullptr;
            dlyLine = nullptr;
        }
        // B0,B1..Border,A0,A1..Aorder = 2*(order+1)
        IIR(size_t n, int order, const T * taps) {
            initCoefficients(n,order,taps);
        }
        ~IIR() {
            if(pBuffer) Free(pBuffer);      
            if(dlyLine) Free(dlyLine);          
        }        
        void initCoefficients(int n, int order, const T* taps) {
            if(pBuffer) Free(pBuffer);
            if(dlyLine) Free(dlyLine);
            int bufferSize;
            len = n;
            // have to save it so it doesn't pop
            IIRGetDlyLine<T>(pState,dlyLine);            
            IIRGetStateSize<T>(n,&bufferSize);            
            pBuffer = Malloc<Ipp8u>(bufferSize);            
            IIRInit<T>(&pState,taps,order,NULL,pBuffer);            
            IIRSetDlyLine<T>(pState,dlyLine);            
        }
        void Execute(const T* pSrc, T* pDst)
        {
            IIR_<T>(pSrc,pDst,len,pState);            
        }
    };

    template<typename T>
    struct IIRBiquad
    {        
        Ipp8u * pBuffer;
        size_t  len,numb;
        void  * pState;
        T*      dlyLine;
        
        IIRBiquad() {
            pBuffer = nullptr;
            len = 0;
            numb=0;
            pState = nullptr;
            dlyLine = nullptr;
        }
        // B0,B1..Border,A0,A1..Aorder = 2*(order+1)
        IIRBiquad(size_t n, int numBiquads, const T* taps) {                                    
            initCoefficients(n,numBiquads,taps);
        }
        ~IIRBiquad() {
            if(pBuffer) Free(pBuffer);        
            if(dlyLine) Free(dlyLine);
        }        
        void initCoefficients(size_t n, int numBiquads, const T* taps) {
            
            int bufferSize;
            
            len = n;            
            if(dlyLine) IIRGetDlyLine<T>(pState,dlyLine);            
            else dlyLine = Malloc<T>(2*numBiquads);            
            if( numb != numBiquads)
            {                
                if(pBuffer) Free(pBuffer);                
                IIRGetStateSize<T>(numBiquads,&bufferSize);                                                    
                pBuffer = Malloc<Ipp8u>(bufferSize);                            
                numb = numBiquads;
            }         
            IIRInit_BiQuad<T>(&pState,taps,numBiquads,dlyLine,pBuffer);                                    
            IIRSetDlyLine<T>(pState,dlyLine);            
        }
        void Execute(const T* pSrc, T* pDst)
        {            
            IIR_<T>(pSrc,pDst,len,pState);                              
        }
    };


    struct IIRButterworthLowpassFilter
    {
        IIR<double> * iir;
        double      * taps;
        Ipp8u       * buffer;
        IIRButterworthLowpassFilter(size_t n, double freq, int order)
        {
            int bufferSize;
            IIRGetStateSize<double>(n,&bufferSize);   
            buffer = Malloc<Ipp8u>(bufferSize);
            taps   = Malloc<double>(n);
            IIRGenLowpassButterworth(freq,0,order,taps,buffer);
            iir = new IIR<double>(n,order,taps);
            assert(iir != nullptr);
        }
        ~IIRButterworthLowpassFilter() {
            if(iir) delete iir;
            if(taps) Free(taps);
            if(buffer) Free(buffer);
        }
        void Execute(const double* pSrc, double* pDst)
        {
            iir->Execute(pSrc,pDst);
        }
    };

    struct IIRButterworthHighpassFilter
    {
        IIR<double> * iir;
        double      * taps;
        Ipp8u       * buffer;
        IIRButterworthHighpassFilter(size_t n, double freq, int order)
        {
            int bufferSize;
            IIRGetStateSize<double>(n,&bufferSize);   
            buffer = Malloc<Ipp8u>(bufferSize);
            taps   = Malloc<double>(n);
            IIRGenHighpassButterworth(freq,0,order,taps,buffer);
            iir = new IIR<double>(n,order,taps);
            assert(iir != nullptr);
        }
        ~IIRButterworthHighpassFilter() {
            if(iir) delete iir;
            if(taps) Free(taps);
            if(buffer) Free(buffer);
        }
        void Execute(const double* pSrc, double* pDst)
        {
            iir->Execute(pSrc,pDst);
        }
    };

    struct IIRChebyshevLowpassFilter
    {
        IIR<double> * iir;
        double      * taps;
        Ipp8u       * buffer;


        IIRChebyshevLowpassFilter(int order,size_t n, double freq, double ripple)
        {
            int bufferSize;
            IIRGetStateSize<double>(n,&bufferSize);   
            buffer = Malloc<Ipp8u>(bufferSize);
            taps   = Malloc<double>(n);
            IIRGenLowpassChebyshev1(freq,ripple,order,taps,buffer);
            iir = new IIR<double>(n,order,taps);
            assert(iir != nullptr);
        }
        ~IIRChebyshevLowpassFilter() {
            if(iir) delete iir;
            if(taps) Free(taps);
            if(buffer) Free(buffer);
        }
        void Execute(const double* pSrc, double* pDst)
        {
            iir->Execute(pSrc,pDst);
        }
    };

    struct IIRChebyshevHighpassFilter
    {
        IIR<double> * iir;
        double      * taps;
        Ipp8u       * buffer;


        IIRChebyshevHighpassFilter(int order,size_t n, double freq, double ripple)
        {
            int bufferSize;
            IIRGetStateSize<double>(n,&bufferSize);   
            buffer = Malloc<Ipp8u>(bufferSize);
            taps   = Malloc<double>(n);
            IIRGenHighpassChebyshev1(freq,ripple,order,taps,buffer);
            iir = new IIR<double>(n,order,taps);
            assert(iir != nullptr);
        }
        ~IIRChebyshevHighpassFilter() {
            if(iir) delete iir;
            if(taps) Free(taps);
            if(buffer) Free(buffer);
        }
        void Execute(const double* pSrc, double* pDst)
        {
            iir->Execute(pSrc,pDst);
        }
    };

    template<typename T>
    void filter(IIR<T> & filter, const T* pSrc, T* pDst)
    {
        filter.Execute(pSrc,pDst);
    }
    template<typename T>
    void filter(IIRBiquad<T> & filter, const T* pSrc, T* pDst)
    {
        filter.Execute(pSrc,pDst);
    }

    void filter(IIRButterworthLowpassFilter & filter, const double* pSrc, double* pDst)
    {
        filter.Execute(pSrc,pDst);
    }
    void filter(IIRButterworthHighpassFilter & filter, const double* pSrc, double* pDst)
    {
        filter.Execute(pSrc,pDst);
    }
    void filter(IIRChebyshevLowpassFilter & filter, const double* pSrc, double* pDst)
    {
        filter.Execute(pSrc,pDst);
    }
    void filter(IIRChebyshevHighpassFilter & filter, const double* pSrc, double* pDst)
    {
        filter.Execute(pSrc,pDst);
    }
   
}