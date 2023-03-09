#pragma once

namespace Casino::IPP
{
      struct HilbertTransform32
    {   
        Ipp8u * pSpec;
        Ipp8u * pBuffer;
        size_t  size;        

        HilbertTransform32(size_t n) {
            int specSize,bufferSize;
            size = n;
            IppStatus status = ippsHilbertGetSize_32f32fc(n,ippAlgHintNone,&specSize,&bufferSize);
            checkStatus(status);
            pSpec = Malloc<Ipp8u>(specSize);
            pBuffer= Malloc<Ipp8u>(bufferSize);
            status = ippsHilbertInit_32f32fc(n,ippAlgHintNone,(IppsHilbertSpec*)pSpec,pBuffer);
            checkStatus(status);
        }
        ~HilbertTransform32() {
            if(pSpec) Free(pSpec);
            if(pBuffer) Free(pBuffer);
        }
        void Execute(const Ipp32f * pSrc, Ipp32fc * pDst) {
            IppStatus status = ippsHilbert_32f32fc(pSrc,pDst,(IppsHilbertSpec*)pSpec,pBuffer);
            checkStatus(status);            
        }
    };

    struct HilbertTransform64 
    {    
        Ipp8u * pSpec;
        Ipp8u * pBuffer;
        size_t  size;        
        HilbertTransform64(size_t n) {
            int specSize,bufferSize;
            size = n;
            IppStatus status = ippsHilbertGetSize_64f64fc(n,ippAlgHintNone,&specSize,&bufferSize);
            checkStatus(status);
            pSpec = Malloc<Ipp8u>(specSize);
            pBuffer= Malloc<Ipp8u>(bufferSize);
            status = ippsHilbertInit_64f64fc(n,ippAlgHintNone,(IppsHilbertSpec*)pSpec,pBuffer);
            checkStatus(status);
        }
        ~HilbertTransform64() {
            if(pSpec) Free(pSpec);
            if(pBuffer) Free(pBuffer);
        }
        void Execute(const Ipp64f * pSrc, Ipp64fc * pDst) {
            IppStatus status = ippsHilbert_64f64fc(pSrc,pDst,(IppsHilbertSpec*)pSpec,pBuffer);
            checkStatus(status);            
        }
    };

    void hilbert(size_t n, const Ipp32f* pSrc, Ipp32fc * pDst)
    {
        HilbertTransform32 h(n);
        h.Execute(pSrc,pDst);
    }
    void hilbert(HilbertTransform32& h, const Ipp32f* pSrc, Ipp32fc * pDst)
    {
        h.Execute(pSrc,pDst);
    }
    void hilbert(size_t n, const Ipp64f* pSrc, Ipp64fc * pDst)
    {
        HilbertTransform64 h(n);
        h.Execute(pSrc,pDst);
    }
    void hilbert(HilbertTransform64& h, const Ipp64f* pSrc, Ipp64fc * pDst)
    {
        h.Execute(pSrc,pDst);
    }
}