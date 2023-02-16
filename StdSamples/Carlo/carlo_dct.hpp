#pragma once

namespace Casino::IPP
{
    template<typename T>
    void DCTFwdInit(void ** ppDCTSpec, int len, IppHintAlgorithm hint, Ipp8u* pSpec, Ipp8u* pSpecBuffer)
    {
        assert(1==0);
    }
    template<>        
    void DCTFwdInit<Ipp32f>(void ** ppDCTSpec, int len, IppHintAlgorithm hint, Ipp8u* pSpec, Ipp8u* pSpecBuffer)
    {
        IppStatus status = ippsDCTFwdInit_32f((IppsDCTFwdSpec_32f**)ppDCTSpec,len,hint,pSpec,pSpecBuffer);
        checkStatus(status);
    }
    template<>        
    void DCTFwdInit<Ipp64f>(void ** ppDCTSpec, int len, IppHintAlgorithm hint, Ipp8u* pSpec, Ipp8u* pSpecBuffer)
    {
        IppStatus status = ippsDCTFwdInit_64f((IppsDCTFwdSpec_64f**)ppDCTSpec,len,hint,pSpec,pSpecBuffer);
        checkStatus(status);
    }
    template<typename T>
    void DCTInvInit(void ** ppDCTSpec, int len, IppHintAlgorithm hint, Ipp8u* pSpec, Ipp8u* pSpecBuffer)
    {
        assert(1==0);
    }
    template<>        
    void DCTInvInit<Ipp32f>(void ** ppDCTSpec, int len, IppHintAlgorithm hint, Ipp8u* pSpec, Ipp8u* pSpecBuffer)
    {
        IppStatus status = ippsDCTInvInit_32f((IppsDCTInvSpec_32f**)ppDCTSpec,len,hint,pSpec,pSpecBuffer);
        checkStatus(status);
    }
    template<>        
    void DCTInvInit<Ipp64f>(void ** ppDCTSpec, int len, IppHintAlgorithm hint, Ipp8u* pSpec, Ipp8u* pSpecBuffer)
    {
        IppStatus status = ippsDCTInvInit_64f((IppsDCTInvSpec_64f**)ppDCTSpec,len,hint,pSpec,pSpecBuffer);
        checkStatus(status);
    }
    template<typename T>
    void DCTFwdGetSize(int len, IppHintAlgorithm hint, int * pSpectSize, int * pSpecBufferSize, int * pBuffersize)
    {
        assert(1==0);
    }
    template<>
    void DCTFwdGetSize<Ipp32f>(int len, IppHintAlgorithm hint, int * pSpectSize, int * pSpecBufferSize, int * pBufferSize)
    {
        IppStatus status = ippsDCTFwdGetSize_32f(len,hint,pSpectSize,pSpecBufferSize,pBufferSize);
        checkStatus(status);
    }
    template<>
    void DCTFwdGetSize<Ipp64f>(int len, IppHintAlgorithm hint, int * pSpectSize, int * pSpecBufferSize, int * pBufferSize)
    {
        IppStatus status = ippsDCTFwdGetSize_64f(len,hint,pSpectSize,pSpecBufferSize,pBufferSize);
        checkStatus(status);
    }
    template<typename T>
    void DCTInvGetSize(int len, IppHintAlgorithm hint, int * pSpectSize, int * pSpecBufferSize, int * pBufferSize)
    {
        assert(1==0);
    }
    template<>
    void DCTInvGetSize<Ipp32f>(int len, IppHintAlgorithm hint, int * pSpectSize, int * pSpecBufferSize, int * pBufferSize)
    {
        IppStatus status = ippsDCTInvGetSize_32f(len,hint,pSpectSize,pSpecBufferSize,pBufferSize);
        checkStatus(status);
    }
    template<>
    void DCTInvGetSize<Ipp64f>(int len, IppHintAlgorithm hint, int * pSpectSize, int * pSpecBufferSize, int * pBufferSize)
    {
        IppStatus status = ippsDCTInvGetSize_64f(len,hint,pSpectSize,pSpecBufferSize,pBufferSize);
        checkStatus(status);
    }
    template<typename T>
    void DCTFwd(const T* pSrc, T * pDst, void *, Ipp8u* pBuffer)
    {
        assert(1==0);
    }
    template<>
    void DCTFwd<Ipp32f>(const Ipp32f* pSrc, Ipp32f * pDst, void * pDCTSpec, Ipp8u* pBuffer)    
    {
        IppStatus status = ippsDCTFwd_32f(pSrc,pDst,(const IppsDCTFwdSpec_32f*)pDCTSpec,pBuffer);
        checkStatus(status);
    }
    template<>
    void DCTFwd<Ipp64f>(const Ipp64f* pSrc, Ipp64f * pDst, void * pDCTSpec, Ipp8u* pBuffer)    
    {
        IppStatus status = ippsDCTFwd_64f(pSrc,pDst,(const IppsDCTFwdSpec_64f*)pDCTSpec,pBuffer);
        checkStatus(status);
    }
    template<typename T>
    void DCTInv(const T* pSrc, T * pDst, void *, Ipp8u* pBuffer)
    {
        assert(1==0);
    }
    template<>
    void DCTInv<Ipp32f>(const Ipp32f* pSrc, Ipp32f * pDst, void * pDCTSpec, Ipp8u* pBuffer)    
    {
        IppStatus status = ippsDCTInv_32f(pSrc,pDst,(const IppsDCTInvSpec_32f*)pDCTSpec,pBuffer);
        checkStatus(status);
    }
    template<>
    void DCTInv<Ipp64f>(const Ipp64f* pSrc, Ipp64f * pDst, void * pDCTSpec, Ipp8u* pBuffer)    
    {
        IppStatus status = ippsDCTInv_64f(pSrc,pDst,(const IppsDCTInvSpec_64f*)pDCTSpec,pBuffer);
        checkStatus(status);
    }




    template<typename T>
    struct DCT
    {
        Ipp8u * pBuffer[2];        
        Ipp8u * pSpec[2];
        Ipp8u * pSpecBuffer[2];
        void * forward;
        void *inverse;
                        
        DCT(size_t n) {
            int size,specbuffer,spec;            
            int flag = IPP_FFT_DIV_FWD_BY_N;

            DCTFwdGetSize<T>(n,ippAlgHintNone,&spec,&specbuffer,&size);            
            pSpec[0]       = Malloc<Ipp8u>(spec);
            pSpecBuffer[0] = specbuffer > 0? Malloc<Ipp8u>(specbuffer) : NULL;
            pBuffer[0]     = Malloc<Ipp8u>(size);

            int order = n;
            DCTFwdInit<T>(&forward,n,ippAlgHintNone,pSpec[0],pSpecBuffer[0]);            
            

            DCTInvGetSize<T>(n,ippAlgHintNone,&spec,&specbuffer,&size);            
            pSpec[1]       = Malloc<Ipp8u>(spec);
            pSpecBuffer[1] = specbuffer > 0? Malloc<Ipp8u>(specbuffer) : NULL;
            pBuffer[1]     = Malloc<Ipp8u>(size);

            DCTInvInit<T>(&inverse,n,ippAlgHintNone,pSpec[1],pSpecBuffer[1]);                        
        }
        ~DCT() {
            if(forward) Free(forward);
            if(inverse) Free(inverse);
            if(pSpecBuffer[0]) Free(pSpecBuffer[0]);
            if(pBuffer[0]) Free(pBuffer[0]);
            if(pSpecBuffer[1]) Free(pSpecBuffer[1]);
            if(pBuffer[1]) Free(pBuffer[1]);        
        }        
        void Forward(const T* pSrc, T* pDst)
        {                   
            DCTFwd<T>(pSrc,pDst,forward,pSpecBuffer[0]);            
        }
        void Inverse(const T* pSrc, T* pDst)
        {                   
            DCTInv<T>(pSrc,pDst,inverse,pSpecBuffer[1]);            
        }
    };
    
    void dct(size_t n, const float * pSrc, float * pDst)
    {
        DCT<float> d(n);
        d.Forward(pSrc,pDst);
    }
    void dct(DCT<float>& d, const float * pSrc, float * pDst)
    {        
        d.Forward(pSrc,pDst);
    }
    void idct(size_t n, const float * pSrc, float * pDst)
    {
        DCT<float> d(n);
        d.Inverse(pSrc,pDst);
    }
    void idct(DCT<float> &d, const float * pSrc, float * pDst)
    {        
        d.Inverse(pSrc,pDst);
    }
        
    void dct(size_t n, const double * pSrc, double * pDst)
    {
        DCT<double> d(n);
        d.Forward(pSrc,pDst);
    }
    void dct(DCT<double>& d, const double * pSrc, double * pDst)
    {        
        d.Forward(pSrc,pDst);
    }
    void idct(size_t n, const double * pSrc, double * pDst)
    {
        DCT<double> d(n);
        d.Inverse(pSrc,pDst);
    }
    void idct(DCT<double> &d, const double * pSrc, double * pDst)
    {        
        d.Inverse(pSrc,pDst);
    }      
 
}