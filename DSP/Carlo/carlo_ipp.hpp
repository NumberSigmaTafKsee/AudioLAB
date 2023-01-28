#pragma once


#include <complex>
#include <algorithm>
#include <vector>
#include <random>
#include <functional>
#include <memory>
#include <chrono>
#include <cmath>
#include <cassert>


#include <ippcore.h>
#include <ipps.h>

namespace Casino::IPP
{    
    enum Algorithm
    {
        ALG_AUTO = ippAlgAuto,
        ALG_DIRECT= ippAlgDirect,
        ALG_FFT = ippAlgFFT,
    };

    enum Norm
    {
        NORM_NONE = ippsNormNone,
        NORM_NORMA = ippsNormA,
        NORM_NORMB = ippsNormB,
    };    
    

    inline const char* GetStatusString(IppStatus status) {
        return ippGetStatusString(status);
    }
    inline void checkStatus(IppStatus status) {
        if(status != ippStsNoErr) {
            std::string err = GetStatusString(status);
            std::cout << "Status Error: " << err << std::endl;
            throw std::runtime_error(err);
        }
    }
    void Init() {
        IppStatus status = ippInit();
        checkStatus(status);
    }

    template<typename T>
    T * Malloc(size_t n) {
        throw std::runtime_error("Called the abstract Malloc\n");
    }
    template<>
    Ipp8u* Malloc<Ipp8u>(size_t n) {
        return ippsMalloc_8u(n);
    }
    template<>
    Ipp32f* Malloc<Ipp32f>(size_t n) {
        return ippsMalloc_32f(n);
    }
    template<>
    Ipp32fc* Malloc<Ipp32fc>(size_t n) {
        return ippsMalloc_32fc(n);
    }
    template<>
    Ipp64f* Malloc<Ipp64f>(size_t n) {
        return ippsMalloc_64f(n);
    }
    template<>
    Ipp64fc* Malloc<Ipp64fc>(size_t n) {
        return ippsMalloc_64fc(n);
    }

    template<typename T>
    void Free(T * ptr) {
        ippFree(ptr);
    }
    
    // Intel IPP
    template<typename T>
    inline void Move(const T* pSrc, T* pDst, int len) {
        std::memcpy(pDst,pSrc,len*sizeof(T));
    }
    template<>
    inline void Move<Ipp8u>(const Ipp8u* pSrc, Ipp8u* pDst, int len) {
        IppStatus status = ippsMove_8u(pSrc,pDst,len);
        checkStatus(status);
    }
    template<>
    inline void Move<Ipp16s>(const Ipp16s* pSrc, Ipp16s* pDst, int len) {
        IppStatus status = ippsMove_16s(pSrc,pDst,len);
        checkStatus(status);
    }
    template<>
    inline void Move<Ipp32f>(const Ipp32f* pSrc, Ipp32f* pDst, int len) {
        IppStatus status = ippsMove_32f(pSrc,pDst,len);
        checkStatus(status);
    }
    template<>
    inline void Move<Ipp32fc>(const Ipp32fc* pSrc, Ipp32fc* pDst, int len) {
        IppStatus status = ippsMove_32fc(pSrc,pDst,len);
        checkStatus(status);
    }    
    template<>
    inline void Move<Ipp64f>(const Ipp64f* pSrc, Ipp64f* pDst, int len) {
        IppStatus status = ippsMove_64f(pSrc,pDst,len);
        checkStatus(status);
    }
    template<>
    inline void Move<Ipp64fc>(const Ipp64fc* pSrc, Ipp64fc* pDst, int len) {
        IppStatus status = ippsMove_64fc(pSrc,pDst,len);
        checkStatus(status);
    }    
    template<>
    inline void Move<Ipp32s>(const Ipp32s* pSrc, Ipp32s* pDst, int len) {
        IppStatus status = ippsMove_32s(pSrc,pDst,len);
        checkStatus(status);
    }    
    template<>
    inline void Move<Ipp64s>(const Ipp64s* pSrc, Ipp64s* pDst, int len) {
        IppStatus status = ippsMove_64s(pSrc,pDst,len);
        checkStatus(status);
    }    
    template<>
    inline void Move<Ipp16sc>(const Ipp16sc* pSrc, Ipp16sc* pDst, int len) {
        IppStatus status = ippsMove_16sc(pSrc,pDst,len);
        checkStatus(status);
    }
    template<>
    inline void Move<Ipp32sc>(const Ipp32sc* pSrc, Ipp32sc* pDst, int len) {
        IppStatus status = ippsMove_32sc(pSrc,pDst,len);
        checkStatus(status);
    }    
    template<>
    inline void Move<Ipp64sc>(const Ipp64sc* pSrc, Ipp64sc* pDst, int len) {
        IppStatus status = ippsMove_64sc(pSrc,pDst,len);
        checkStatus(status);
    }

    // i dont know if there is a differen with copy and move
    template<typename T>
    inline void Copy(const T* pSrc, T* pDst, int len) {
        //std::memcpy(pDst,pSrc,len*sizeof(T));        
        std::copy(&pSrc[0],&pSrc[len-1],&pDst[0]);
    }

    template<>
    inline void Copy<Ipp8u>(const Ipp8u* pSrc, Ipp8u* pDst, int len) {
        IppStatus status = ippsCopy_8u(pSrc,pDst,len);
        checkStatus(status);
    }

    template<>
    inline void Copy<Ipp16s>(const Ipp16s* pSrc, Ipp16s* pDst, int len) {
        IppStatus status = ippsCopy_16s(pSrc,pDst,len);
        checkStatus(status);
    }

    template<>
    inline void Copy<Ipp32s>(const Ipp32s* pSrc, Ipp32s* pDst, int len) {
        IppStatus status = ippsCopy_32s(pSrc,pDst,len);
        checkStatus(status);
    }

    template<>
    inline void Copy<Ipp64s>(const Ipp64s* pSrc, Ipp64s* pDst, int len) {
        IppStatus status = ippsCopy_64s(pSrc,pDst,len);
        checkStatus(status);
    }

    template<>
    inline void Copy<Ipp16sc>(const Ipp16sc* pSrc, Ipp16sc* pDst, int len) {
        IppStatus status = ippsCopy_16sc(pSrc,pDst,len);
        checkStatus(status);
    }

    template<>
    inline void Copy<Ipp32sc>(const Ipp32sc* pSrc, Ipp32sc* pDst, int len) {
        IppStatus status = ippsCopy_32sc(pSrc,pDst,len);
        checkStatus(status);
    }   

    template<>
    inline void Copy<Ipp64sc>(const Ipp64sc* pSrc, Ipp64sc* pDst, int len) {
        IppStatus status = ippsCopy_64sc(pSrc,pDst,len);
        checkStatus(status);
    }

    template<>
    inline void Copy<Ipp32f>(const Ipp32f* pSrc, Ipp32f* pDst, int len) {
        IppStatus status = ippsCopy_32f(pSrc,pDst,len);
        checkStatus(status);
    }

    template<>
    inline void Copy<Ipp32fc>(const Ipp32fc* pSrc, Ipp32fc* pDst, int len) {
        IppStatus status = ippsCopy_32fc(pSrc,pDst,len);
        checkStatus(status);
    }

    template<typename T>
    inline void Set(const T val, T* pDst, int len) {
        //for(size_t i = 0; i < len; i++) pDst[i] = val;
        std::fill(pDst,val,len);
    }
    template<>
    inline void Set<Ipp32f>(const Ipp32f val, Ipp32f* pDst, int len) {
        IppStatus status = ippsSet_32f(val,pDst,len);
        checkStatus(status);
    }

    template<>
    inline void Set<Ipp32fc>(const Ipp32fc val, Ipp32fc* pDst, int len) {
        IppStatus status = ippsSet_32fc(val,pDst,len);
        checkStatus(status);
    }    
    template<>
    inline void Set<Ipp64f>(const Ipp64f val, Ipp64f* pDst, int len) {
        IppStatus status = ippsSet_64f(val,pDst,len);
        checkStatus(status);
    }

    template<>
    inline void Set<Ipp64fc>(const Ipp64fc val, Ipp64fc* pDst, int len) {
        IppStatus status = ippsSet_64fc(val,pDst,len);
        checkStatus(status);
    }    

    template<>
    inline void Set<Ipp8u>(const Ipp8u val, Ipp8u* pDst, int len) {
        IppStatus status = ippsSet_8u(val,pDst,len);
        checkStatus(status);
    }

    template<>
    inline void Set<Ipp16s>(const Ipp16s val, Ipp16s* pDst, int len) {
        IppStatus status = ippsSet_16s(val,pDst,len);
        checkStatus(status);
    }

    template<>
    inline void Set<Ipp32s>(const Ipp32s val, Ipp32s* pDst, int len) {
        IppStatus status = ippsSet_32s(val,pDst,len);
        checkStatus(status);
    }

    template<>
    inline void Set<Ipp64s>(const Ipp64s val, Ipp64s* pDst, int len) {
        IppStatus status = ippsSet_64s(val,pDst,len);
        checkStatus(status);
    }

    template<>
    inline void Set<Ipp16sc>(const Ipp16sc val, Ipp16sc* pDst, int len) {
        IppStatus status = ippsSet_16sc(val,pDst,len);
        checkStatus(status);
    }

    template<>
    inline void Set<Ipp32sc>(const Ipp32sc val, Ipp32sc* pDst, int len) {
        IppStatus status = ippsSet_32sc(val,pDst,len);
        checkStatus(status);
    }

    template<>
    inline void Set<Ipp64sc>(const Ipp64sc val, Ipp64sc* pDst, int len) {
        IppStatus status = ippsSet_64sc(val,pDst,len);
        checkStatus(status);
    }    

    template<typename T>
    inline void Zero(T* pDst, int len) {
        std::fill(&pDst[0],&pDst[len-1],len);
    }
    template<>
    inline void Zero<Ipp8u>(Ipp8u* pDst, int len) {
        IppStatus status = ippsZero_8u(pDst,len);
        checkStatus(status);
    }

    template<>
    inline void Zero<Ipp16s>(Ipp16s* pDst, int len) {
        IppStatus status = ippsZero_16s(pDst,len);
        checkStatus(status);
    }

    template<>
    inline void Zero<Ipp32s>(Ipp32s* pDst, int len) {
        IppStatus status = ippsZero_32s(pDst,len);
        checkStatus(status);
    }

    template<>
    inline void Zero<Ipp64s>(Ipp64s* pDst, int len) {
        IppStatus status = ippsZero_64s(pDst,len);
        checkStatus(status);
    }

    template<>
    inline void Zero<Ipp16sc>(Ipp16sc* pDst, int len) {
        IppStatus status = ippsZero_16sc(pDst,len);
        checkStatus(status);
    }

    template<>
    inline void Zero<Ipp32sc>(Ipp32sc* pDst, int len) {
        IppStatus status = ippsZero_32sc(pDst,len);
        checkStatus(status);
    } 

    template<>
    inline void Zero<Ipp64sc>(Ipp64sc* pDst, int len) {
        IppStatus status = ippsZero_64sc(pDst,len);
        checkStatus(status);
    }

    template<>
    inline void Zero<Ipp32f>(Ipp32f* pDst, int len) {
        IppStatus status = ippsZero_32f(pDst,len);
        checkStatus(status);
    }

    template<>
    inline void Zero<Ipp32fc>(Ipp32fc* pDst, int len) {
        IppStatus status = ippsZero_32fc(pDst,len);
        checkStatus(status);
    }

    template<>
    inline void Zero<Ipp64f>(Ipp64f* pDst, int len) {
        IppStatus status = ippsZero_64f(pDst,len);
        checkStatus(status);
    }

    template<>
    inline void Zero<Ipp64fc>(Ipp64fc* pDst, int len) {
        IppStatus status = ippsZero_64fc(pDst,len);
        checkStatus(status);
    }


    inline void SetNumThreads(int num) {
        IppStatus status = ippSetNumThreads(num);
        checkStatus(status);
    }

    inline void SetFlushToZero(int value, unsigned int *uMask = NULL) {
        IppStatus status = ippSetFlushToZero(value,uMask);
        checkStatus(status);
    }
    inline int GetNumThreads() {
        int num = 1;
        IppStatus status = ippGetNumThreads(&num);
        checkStatus(status);
        return num;
    }
    template<typename T>
    T Sinc(const T x)
    {
        if(std::abs((double)x) < 1e-6)
            return 1.0;
        return std::sin(M_PI*(double)x) / (M_PI*(double)x);
    }        
    template<typename T>
    void Sinc(const T * src, T * dst, size_t n) {
        for(size_t i = 0; i < n; i++)
            dst[i] = Sinc(src[i]); 
    }
    
    void RealToComplex(const Ipp32f * pSrcR, const Ipp32f * pSrcD, Ipp32fc* pDst, int len) {
        IppStatus status = ippsRealToCplx_32f(pSrcR,pSrcD,pDst,len);
        checkStatus(status);
    }
    void ComplexToReal(const Ipp32fc * src, Ipp32f * real, Ipp32f * imag, int len) {
        IppStatus status = ippsCplxToReal_32fc(src,real,imag,len);
        checkStatus(status);
    }
    void Magnitude(const Ipp32f * pSrcR, const Ipp32f * pSrcD, Ipp32f* pDst, int len) {
        IppStatus status = ippsMagnitude_32f(pSrcR,pSrcD,pDst,len);
        checkStatus(status);
    }
    void Phase(const Ipp32f * pSrcR, const Ipp32f * pSrcD, Ipp32f* pDst, int len) {
        IppStatus status = ippsPhase_32f(pSrcR,pSrcD,pDst,len);
        checkStatus(status);
    }
    void CartToPolar(const Ipp32f * pSrcR, const Ipp32f * pSrcI, Ipp32f* pDstMag, Ipp32f* pDstPhase, int len) {
        IppStatus status = ippsCartToPolar_32f(pSrcR,pSrcI,pDstMag,pDstPhase,len);
        checkStatus(status);
    }
    void PolarToCart(const Ipp32f* pMag, const Ipp32f* pPhase, Ipp32f* real, Ipp32f* imag, int len) {
        IppStatus status = ippsPolarToCart_32f(pMag,pPhase,real,imag,len);
        checkStatus(status);
    }
    void Magnitude(const Ipp32fc * pSrc, Ipp32f* pDst, int len) {
        IppStatus status = ippsMagnitude_32fc(pSrc,pDst,len);
        checkStatus(status);
    }
    void Phase(const Ipp32fc * pSrc, Ipp32f* pDst, int len) {
        IppStatus status = ippsPhase_32fc(pSrc,pDst,len);
        checkStatus(status);
    }
    void CartToPolar(const Ipp32fc * src, Ipp32f* pDstMag, Ipp32f* pDstPhase, int len) {
        IppStatus status = ippsCartToPolar_32fc(src,pDstMag,pDstPhase,len);
        checkStatus(status);
    }
    void PolarToCart(const Ipp32f* pMag, const Ipp32f* pPhase, Ipp32fc * dst, int len) {
        IppStatus status = ippsPolarToCart_32fc(pMag,pPhase,dst,len);
        checkStatus(status);
    }
    void Conj(const Ipp32fc * pSrc, Ipp32fc* pDst, int len) {
        IppStatus status = ippsConj_32fc(pSrc,pDst,len);
        checkStatus(status);
    }
    void ConjFlip(const Ipp32fc * pSrc, Ipp32fc* pDst, int len) {
        IppStatus status = ippsConjFlip_32fc(pSrc,pDst,len);
        checkStatus(status);
    }
    void PowerSpectrum(const Ipp32fc * pSrc, Ipp32f* pDst, int len) {
        IppStatus status = ippsPowerSpectr_32fc(pSrc,pDst,len);
        checkStatus(status);
    }
    void Real(const Ipp32fc * pSrc, Ipp32f* pDst, int len) {
        IppStatus status = ippsReal_32fc(pSrc,pDst,len);
        checkStatus(status);
    }
    void Imag(const Ipp32fc * pSrc, Ipp32f* pDst, int len) {
        IppStatus status = ippsImag_32fc(pSrc,pDst,len);
        checkStatus(status);
    }
    void Threshold(const Ipp32f * pSrc, Ipp32f * pDst, size_t len, Ipp32f level, IppCmpOp op = ippCmpGreater) {
        IppStatus status = ippsThreshold_32f(pSrc,pDst,level,len, op);
        checkStatus(status);
    }
    void Threshold(const Ipp32fc * pSrc, Ipp32fc * pDst, size_t len, Ipp32f level, IppCmpOp op = ippCmpGreater) {
        IppStatus status = ippsThreshold_32fc(pSrc,pDst,level,len, op);
        checkStatus(status);
    }
    void WinBartlett(const Ipp32f* src, Ipp32f * dst, int len)
    {
        IppStatus status = ippsWinBartlett_32f(src,dst,len);
        checkStatus(status);
    }
    void WinBartlett(const Ipp32fc* src, Ipp32fc * dst, int len)
    {
        IppStatus status = ippsWinBartlett_32fc(src,dst,len);
        checkStatus(status);
    }
    void WinBartlett(Ipp32f * dst, int len)
    {
        IppStatus status = ippsWinBartlett_32f_I(dst,len);
        checkStatus(status);
    }
    void WinBartlett(Ipp32fc * dst, int len)
    {
        IppStatus status = ippsWinBartlett_32fc_I(dst,len);
        checkStatus(status);
    }

    void WinBlackman(const Ipp32f* src, Ipp32f * dst, int len, Ipp32f alpha)
    {
        IppStatus status = ippsWinBlackman_32f(src,dst,len,alpha);
        checkStatus(status);
    }
    void WinBlackman(const Ipp32fc* src, Ipp32fc * dst, int len, Ipp32f alpha)
    {
        IppStatus status = ippsWinBlackman_32fc(src,dst,len,alpha);
        checkStatus(status);
    }
    void WinBlackman(Ipp32f * dst, int len, Ipp32f alpha)
    {
        IppStatus status = ippsWinBlackman_32f_I(dst,len,alpha);
        checkStatus(status);
    }
    void WinBlackman(Ipp32fc * dst, int len, Ipp32f alpha)
    {
        IppStatus status = ippsWinBlackman_32fc_I(dst,len,alpha);
        checkStatus(status);
    }
    void WinKaiser(const Ipp32f* src, Ipp32f * dst, int len, Ipp32f alpha)
    {
        IppStatus status = ippsWinKaiser_32f(src,dst,len,alpha);
        checkStatus(status);
    }
    void WinKaiser(const Ipp32fc* src, Ipp32fc * dst, int len, Ipp32f alpha)
    {
        IppStatus status = ippsWinKaiser_32fc(src,dst,len,alpha);
        checkStatus(status);
    }
    void WinKaiser(Ipp32f * dst, int len, Ipp32f alpha)
    {
        IppStatus status = ippsWinKaiser_32f_I(dst,len,alpha);
        checkStatus(status);
    }
    void WinKaiser(Ipp32fc * dst, int len, Ipp32f alpha)
    {
        IppStatus status = ippsWinKaiser_32fc_I(dst,len,alpha);
        checkStatus(status);
    }
    void WinHamming(const Ipp32f* src, Ipp32f * dst, int len)
    {
        IppStatus status = ippsWinHamming_32f(src,dst,len);
        checkStatus(status);
    }
    void WinHamming(const Ipp32fc* src, Ipp32fc * dst, int len)
    {
        IppStatus status = ippsWinHamming_32fc(src,dst,len);
        checkStatus(status);
    }
    void WinHamming(Ipp32f * dst, int len)
    {
        IppStatus status = ippsWinHamming_32f_I(dst,len);
        checkStatus(status);
    }
    void WinHamming(Ipp32fc * dst, int len)
    {
        IppStatus status = ippsWinHamming_32fc_I(dst,len);
        checkStatus(status);
    }

    void WinHann(const Ipp32f* src, Ipp32f * dst, int len)
    {
        IppStatus status = ippsWinHann_32f(src,dst,len);
        checkStatus(status);
    }
    void WinHann(const Ipp32fc* src, Ipp32fc * dst, int len)
    {
        IppStatus status = ippsWinHann_32fc(src,dst,len);
        checkStatus(status);
    }
    void WinHann(Ipp32f * dst, int len)
    {
        IppStatus status = ippsWinHann_32f_I(dst,len);
        checkStatus(status);
    }
    void WinHann(Ipp32fc * dst, int len)
    {
        IppStatus status = ippsWinHann_32fc_I(dst,len);
        checkStatus(status);
    }
   
    void RealToComplex(const Ipp64f * pSrcR, const Ipp64f * pSrcD, Ipp64fc* pDst, int len) {
        IppStatus status = ippsRealToCplx_64f(pSrcR,pSrcD,pDst,len);
        checkStatus(status);
    }
    void ComplexToReal(const Ipp64fc * src, Ipp64f * real, Ipp64f * imag, int len) {
        IppStatus status = ippsCplxToReal_64fc(src,real,imag,len);
        checkStatus(status);
    }
    void Magnitude(const Ipp64f * pSrcR, const Ipp64f * pSrcD, Ipp64f* pDst, int len) {
        IppStatus status = ippsMagnitude_64f(pSrcR,pSrcD,pDst,len);
        checkStatus(status);
    }
    void Phase(const Ipp64f * pSrcR, const Ipp64f * pSrcD, Ipp64f* pDst, int len) {
        IppStatus status = ippsPhase_64f(pSrcR,pSrcD,pDst,len);
        checkStatus(status);
    }
    void CartToPolar(const Ipp64f * pSrcR, const Ipp64f * pSrcI, Ipp64f* pDstMag, Ipp64f* pDstPhase, int len) {
        IppStatus status = ippsCartToPolar_64f(pSrcR,pSrcI,pDstMag,pDstPhase,len);
        checkStatus(status);
    }
    void PolarToCart(const Ipp64f* pMag, const Ipp64f* pPhase, Ipp64f* real, Ipp64f* imag, int len) {
        IppStatus status = ippsPolarToCart_64f(pMag,pPhase,real,imag,len);
        checkStatus(status);
    }
    void Magnitude(const Ipp64fc * pSrc, Ipp64f* pDst, int len) {
        IppStatus status = ippsMagnitude_64fc(pSrc,pDst,len);
        checkStatus(status);
    }
    void Phase(const Ipp64fc * pSrc, Ipp64f* pDst, int len) {
        IppStatus status = ippsPhase_64fc(pSrc,pDst,len);
        checkStatus(status);
    }
    void CartToPolar(const Ipp64fc * src, Ipp64f* pDstMag, Ipp64f* pDstPhase, int len) {
        IppStatus status = ippsCartToPolar_64fc(src,pDstMag,pDstPhase,len);
        checkStatus(status);
    }
    void PolarToCart(const Ipp64f* pMag, const Ipp64f* pPhase, Ipp64fc * dst, int len) {
        IppStatus status = ippsPolarToCart_64fc(pMag,pPhase,dst,len);
        checkStatus(status);
    }
    void Conj(const Ipp64fc * pSrc, Ipp64fc* pDst, int len) {
        IppStatus status = ippsConj_64fc(pSrc,pDst,len);
        checkStatus(status);
    }
    void ConjFlip(const Ipp64fc * pSrc, Ipp64fc* pDst, int len) {
        IppStatus status = ippsConjFlip_64fc(pSrc,pDst,len);
        checkStatus(status);
    }
    void PowerSpectrum(const Ipp64fc * pSrc, Ipp64f* pDst, int len) {
        IppStatus status = ippsPowerSpectr_64fc(pSrc,pDst,len);
        checkStatus(status);
    }
    void Real(const Ipp64fc * pSrc, Ipp64f* pDst, int len) {
        IppStatus status = ippsReal_64fc(pSrc,pDst,len);
        checkStatus(status);
    }
    void Imag(const Ipp64fc * pSrc, Ipp64f* pDst, int len) {
        IppStatus status = ippsImag_64fc(pSrc,pDst,len);
        checkStatus(status);
    }
    void Threshold(const Ipp64f * pSrc, Ipp64f * pDst, size_t len, Ipp64f level, IppCmpOp op = ippCmpGreater) {
        IppStatus status = ippsThreshold_64f(pSrc,pDst,level,len, op);
        checkStatus(status);
    }
    void Threshold(const Ipp64fc * pSrc, Ipp64fc * pDst, size_t len, Ipp64f level, IppCmpOp op = ippCmpGreater) {
        IppStatus status = ippsThreshold_64fc(pSrc,pDst,level,len, op);
        checkStatus(status);
    }
    void WinBartlett(const Ipp64f* src, Ipp64f * dst, int len)
    {
        IppStatus status = ippsWinBartlett_64f(src,dst,len);
        checkStatus(status);
    }
    void WinBartlett(const Ipp64fc* src, Ipp64fc * dst, int len)
    {
        IppStatus status = ippsWinBartlett_64fc(src,dst,len);
        checkStatus(status);
    }
    void WinBartlett(Ipp64f * dst, int len)
    {
        IppStatus status = ippsWinBartlett_64f_I(dst,len);
        checkStatus(status);
    }
    void WinBartlett(Ipp64fc * dst, int len)
    {
        IppStatus status = ippsWinBartlett_64fc_I(dst,len);
        checkStatus(status);
    }

    void WinBlackman(const Ipp64f* src, Ipp64f * dst, int len, Ipp64f alpha)
    {
        IppStatus status = ippsWinBlackman_64f(src,dst,len,alpha);
        checkStatus(status);
    }
    void WinBlackman(const Ipp64fc* src, Ipp64fc * dst, int len, Ipp64f alpha)
    {
        IppStatus status = ippsWinBlackman_64fc(src,dst,len,alpha);
        checkStatus(status);
    }
    void WinBlackman(Ipp64f * dst, int len, Ipp64f alpha)
    {
        IppStatus status = ippsWinBlackman_64f_I(dst,len,alpha);
        checkStatus(status);
    }
    void WinBlackman(Ipp64fc * dst, int len, Ipp64f alpha)
    {
        IppStatus status = ippsWinBlackman_64fc_I(dst,len,alpha);
        checkStatus(status);
    }
    void WinKaiser(const Ipp64f* src, Ipp64f * dst, int len, Ipp64f alpha)
    {
        IppStatus status = ippsWinKaiser_64f(src,dst,len,alpha);
        checkStatus(status);
    }
    void WinKaiser(const Ipp64fc* src, Ipp64fc * dst, int len, Ipp64f alpha)
    {
        IppStatus status = ippsWinKaiser_64fc(src,dst,len,alpha);
        checkStatus(status);
    }
    void WinKaiser(Ipp64f * dst, int len, Ipp64f alpha)
    {
        IppStatus status = ippsWinKaiser_64f_I(dst,len,alpha);
        checkStatus(status);
    }
    void WinKaiser(Ipp64fc * dst, int len, Ipp64f alpha)
    {
        IppStatus status = ippsWinKaiser_64fc_I(dst,len,alpha);
        checkStatus(status);
    }
    void WinHamming(const Ipp64f* src, Ipp64f * dst, int len)
    {
        IppStatus status = ippsWinHamming_64f(src,dst,len);
        checkStatus(status);
    }
    void WinHamming(const Ipp64fc* src, Ipp64fc * dst, int len)
    {
        IppStatus status = ippsWinHamming_64fc(src,dst,len);
        checkStatus(status);
    }
    void WinHamming(Ipp64f * dst, int len)
    {
        IppStatus status = ippsWinHamming_64f_I(dst,len);
        checkStatus(status);
    }
    void WinHamming(Ipp64fc * dst, int len)
    {
        IppStatus status = ippsWinHamming_64fc_I(dst,len);
        checkStatus(status);
    }

    void WinHann(const Ipp64f* src, Ipp64f * dst, int len)
    {
        IppStatus status = ippsWinHann_64f(src,dst,len);
        checkStatus(status);
    }
    void WinHann(const Ipp64fc* src, Ipp64fc * dst, int len)
    {
        IppStatus status = ippsWinHann_64fc(src,dst,len);
        checkStatus(status);
    }
    void WinHann(Ipp64f * dst, int len)
    {
        IppStatus status = ippsWinHann_64f_I(dst,len);
        checkStatus(status);
    }
    void WinHann(Ipp64fc * dst, int len)
    {
        IppStatus status = ippsWinHann_64fc_I(dst,len);
        checkStatus(status);
    }
    template<typename T>
    void Sum(const T * pSrc, int len, T* pDst, int hint = 0) {
        throw std::runtime_error("Called abstract Sum");
    }

    template<>
    void Sum<Ipp32f>(const Ipp32f * pSrc, int len, Ipp32f* pDst, int hint) {
        IppStatus status = ippsSum_32f(pSrc,len,pDst,(IppHintAlgorithm)hint);
        checkStatus(status);
    }
    template<>
    void Sum<Ipp64f>(const Ipp64f * pSrc, int len, Ipp64f* pDst, int hint) {
        IppStatus status = ippsSum_64f(pSrc,len,pDst);
        checkStatus(status);
    }

    template<typename T>
    void AddC(const T* src, T val, T* dst, int len)
    {
        throw std::runtime_error("foobar");    
    }

    template<>
    void AddC<Ipp32f>(const Ipp32f * src, Ipp32f val, Ipp32f* dst, int len) {
        IppStatus status = ippsAddC_32f(src,val,dst,len);
        checkStatus(status);
    }
    template<>
    void AddC<Ipp64f>(const Ipp64f * src, Ipp64f val, Ipp64f* dst, int len) {
        IppStatus status = ippsAddC_64f(src,val,dst,len);
        checkStatus(status);
    }
    template<>
    void AddC<Ipp32fc>(const Ipp32fc * src, Ipp32fc val, Ipp32fc* dst, int len) {
        IppStatus status = ippsAddC_32fc(src,val,dst,len);
        checkStatus(status);
    }
    template<>
    void AddC<Ipp64fc>(const Ipp64fc * src, Ipp64fc val, Ipp64fc* dst, int len) {
        IppStatus status = ippsAddC_64fc(src,val,dst,len);
        checkStatus(status);
    }

    template<typename T>
    void AddC(T val, T* dst, int len)
    {
        throw std::runtime_error("foobar");    
    }
    template<>
    void AddC<Ipp32f>(Ipp32f val, Ipp32f* dst, int len) {
        IppStatus status = ippsAddC_32f_I(val,dst,len);
        checkStatus(status);
    }
    template<>
    void AddC<Ipp32fc>(Ipp32fc val, Ipp32fc* dst, int len) {
        IppStatus status = ippsAddC_32fc_I(val,dst,len);
        checkStatus(status);
    }
    template<>
    void AddC<Ipp64f>(Ipp64f val, Ipp64f* dst, int len) {
        IppStatus status = ippsAddC_64f_I(val,dst,len);
        checkStatus(status);
    }
    template<>
    void AddC<Ipp64fc>(Ipp64fc val, Ipp64fc* dst, int len) {
        IppStatus status = ippsAddC_64fc_I(val,dst,len);
        checkStatus(status);
    }

    template<typename T>
    void Add(const T * a, T *b, T* dst, int len) {
        throw std::runtime_error("foobar");
    }
    template<>
    void Add<Ipp32f>(const Ipp32f * a, Ipp32f *b, Ipp32f* dst, int len) {
        IppStatus status = ippsAdd_32f(a,b,dst,len);
        checkStatus(status);
    }
    template<>
    void Add<Ipp32fc>(const Ipp32fc * a, Ipp32fc *b, Ipp32fc* dst, int len) {
        IppStatus status = ippsAdd_32fc(a,b,dst,len);
        checkStatus(status);
    }
    template<>
    void Add<Ipp64f>(const Ipp64f * a, Ipp64f *b, Ipp64f* dst, int len) {
        IppStatus status = ippsAdd_64f(a,b,dst,len);
        checkStatus(status);
    }
    template<>
    void Add<Ipp64fc>(const Ipp64fc * a, Ipp64fc *b, Ipp64fc* dst, int len) {
        IppStatus status = ippsAdd_64fc(a,b,dst,len);
        checkStatus(status);
    }

    template<typename T>
    void Add(const T * pSrc, T * pDst, int len) {
        throw std::runtime_error("foobar");
    }
    template<>
    void Add<Ipp32f>(const Ipp32f * src, Ipp32f* dst, int len) {
        IppStatus status = ippsAdd_32f_I(src,dst,len);
        checkStatus(status);
    }
    template<>
    void Add<Ipp32fc>(const Ipp32fc * src, Ipp32fc* dst, int len) {
        IppStatus status = ippsAdd_32fc_I(src,dst,len);
        checkStatus(status);
    }
    template<>
    void Add<Ipp64f>(const Ipp64f * src, Ipp64f* dst, int len) {
        IppStatus status = ippsAdd_64f_I(src,dst,len);
        checkStatus(status);
    }
    template<>
    void Add<Ipp64fc>(const Ipp64fc * src, Ipp64fc* dst, int len) {
        IppStatus status = ippsAdd_64fc_I(src,dst,len);
        checkStatus(status);
    }

    template<typename T>
    void SubC(const T* src, T val, T* dst, int len) {
        throw std::runtime_error("foobar");
    }
    template<>
    void SubC<Ipp32f>(const Ipp32f * src, Ipp32f val, Ipp32f* dst, int len) {
        IppStatus status = ippsSubC_32f(src,val,dst,len);
        checkStatus(status);
    }
    template<>
    void SubC<Ipp32fc>(const Ipp32fc * src, Ipp32fc val, Ipp32fc* dst, int len) {
        IppStatus status = ippsSubC_32fc(src,val,dst,len);
        checkStatus(status);
    }
    template<>
    void SubC<Ipp64f>(const Ipp64f * src, Ipp64f val, Ipp64f* dst, int len) {
        IppStatus status = ippsSubC_64f(src,val,dst,len);
        checkStatus(status);
    }
    template<>
    void SubC<Ipp64fc>(const Ipp64fc * src, Ipp64fc val, Ipp64fc* dst, int len) {
        IppStatus status = ippsSubC_64fc(src,val,dst,len);
        checkStatus(status);
    }

    template<typename T>
    void SubC(T val, T* dst, int len) {
        throw std::runtime_error("foobar");
    }
    template<>
    void SubC<Ipp32f>(Ipp32f val, Ipp32f* dst, int len) {
        IppStatus status = ippsSubC_32f_I(val,dst,len);
        checkStatus(status);
    }
    template<>
    void SubC<Ipp32fc>(Ipp32fc val, Ipp32fc* dst, int len) {
        IppStatus status = ippsSubC_32fc_I(val,dst,len);
        checkStatus(status);
    }
    template<>
    void SubC<Ipp64f>(Ipp64f val, Ipp64f* dst, int len) {
        IppStatus status = ippsSubC_64f_I(val,dst,len);
        checkStatus(status);
    }
    template<>
    void SubC<Ipp64fc>(Ipp64fc val, Ipp64fc* dst, int len) {
        IppStatus status = ippsSubC_64fc_I(val,dst,len);
        checkStatus(status);
    }

    template<typename T>
    void Sub(const T * a, T *b, T* dst, int len) {
        assert(1==0);
    }

    template<>
    void Sub<Ipp32f>(const Ipp32f * a, Ipp32f *b, Ipp32f* dst, int len) {
        IppStatus status = ippsSub_32f(a,b,dst,len);
        checkStatus(status);
    }
    template<>
    void Sub<Ipp32fc>(const Ipp32fc * a, Ipp32fc *b, Ipp32fc* dst, int len) {
        IppStatus status = ippsSub_32fc(a,b,dst,len);
        checkStatus(status);
    }
    template<>
    void Sub<Ipp64f>(const Ipp64f * a, Ipp64f *b, Ipp64f* dst, int len) {
        IppStatus status = ippsSub_64f(a,b,dst,len);
        checkStatus(status);
    }
    template<>
    void Sub<Ipp64fc>(const Ipp64fc * a, Ipp64fc *b, Ipp64fc* dst, int len) {
        IppStatus status = ippsSub_64fc(a,b,dst,len);
        checkStatus(status);
    }

    template<typename T>
    void Sub(const T * a, T* dst, int len) {
        assert(1==0);
    }
    template<>
    void Sub<Ipp32f>(const Ipp32f * a, Ipp32f* dst, int len) {
        IppStatus status = ippsSub_32f_I(a,dst,len);
        checkStatus(status);
    }
    template<>
    void Sub<Ipp32fc>(const Ipp32fc * a, Ipp32fc* dst, int len) {
        IppStatus status = ippsSub_32fc_I(a,dst,len);
        checkStatus(status);
    }
    template<>
    void Sub<Ipp64f>(const Ipp64f * a, Ipp64f* dst, int len) {
        IppStatus status = ippsSub_64f_I(a,dst,len);
        checkStatus(status);
    }
    template<>
    void Sub<Ipp64fc>(const Ipp64fc * a, Ipp64fc* dst, int len) {
        IppStatus status = ippsSub_64fc_I(a,dst,len);
        checkStatus(status);
    }
    
    template<typename T>
    void SubCRev(const T * src, T val, T* dst, int len) {
        assert(1==0);
    }

    template<>
    void SubCRev<Ipp32f>(const Ipp32f * src, Ipp32f val, Ipp32f* dst, int len) {
        IppStatus status = ippsSubCRev_32f(src,val,dst,len);
        checkStatus(status);
    }
    template<>
    void SubCRev<Ipp32fc>(const Ipp32fc * src, Ipp32fc val, Ipp32fc* dst, int len) {
        IppStatus status = ippsSubCRev_32fc(src,val,dst,len);
        checkStatus(status);
    }
    template<>
    void SubCRev<Ipp64f>(const Ipp64f * src, Ipp64f val, Ipp64f* dst, int len) {
        IppStatus status = ippsSubCRev_64f(src,val,dst,len);
        checkStatus(status);
    }
    template<>
    void SubCRev<Ipp64fc>(const Ipp64fc * src, Ipp64fc val, Ipp64fc* dst, int len) {
        IppStatus status = ippsSubCRev_64fc(src,val,dst,len);
        checkStatus(status);
    }

    template<typename T>
    void SubRev(const T * a, T b, T* dst, int len) {
        assert(1==0);
    }

    template<>
    void SubRev<Ipp32f>(const Ipp32f * a, Ipp32f b, Ipp32f* dst, int len) {
        IppStatus status = ippsSubCRev_32f(a,b,dst,len);
        checkStatus(status);
    }
    template<>
    void SubRev<Ipp32fc>(const Ipp32fc * a, Ipp32fc b, Ipp32fc* dst, int len) {
        IppStatus status = ippsSubCRev_32fc(a,b,dst,len);
        checkStatus(status);
    }
    template<>
    void SubRev<Ipp64f>(const Ipp64f * a, Ipp64f b, Ipp64f* dst, int len) {
        IppStatus status = ippsSubCRev_64f(a,b,dst,len);
        checkStatus(status);
    }
    template<>
    void SubRev<Ipp64fc>(const Ipp64fc * a, Ipp64fc b, Ipp64fc* dst, int len) {
        IppStatus status = ippsSubCRev_64fc(a,b,dst,len);
        checkStatus(status);
    }

    template<typename T>
    void MulC(const T* src, T val, T* dst, int len) {
        assert(1==0);
    }
    template<>
    void MulC<Ipp32f>(const Ipp32f * src, Ipp32f val, Ipp32f* dst, int len) {
        IppStatus status = ippsMulC_32f(src,val,dst,len);
        checkStatus(status);
    }
    template<>
    void MulC<Ipp32fc>(const Ipp32fc * src, Ipp32fc val, Ipp32fc* dst, int len) {
        IppStatus status = ippsMulC_32fc(src,val,dst,len);
        checkStatus(status);
    }
    template<>
    void MulC<Ipp64f>(const Ipp64f * src, Ipp64f val, Ipp64f* dst, int len) {
        IppStatus status = ippsMulC_64f(src,val,dst,len);
        checkStatus(status);
    }
    template<>
    void MulC<Ipp64fc>(const Ipp64fc * src, Ipp64fc val, Ipp64fc* dst, int len) {
        IppStatus status = ippsMulC_64fc(src,val,dst,len);
        checkStatus(status);
    }

    template<typename T>
    void MulC(T val, T* dst, int len) {
        assert(1==0);
    }
    template<>
    void MulC<Ipp32f>(Ipp32f val, Ipp32f* dst, int len) {
        IppStatus status = ippsMulC_32f_I(val,dst,len);
        checkStatus(status);
    }
    template<>
    void MulC<Ipp32fc>(Ipp32fc val, Ipp32fc* dst, int len) {
        IppStatus status = ippsMulC_32fc_I(val,dst,len);
        checkStatus(status);
    }
    template<>
    void MulC<Ipp64f>(Ipp64f val, Ipp64f* dst, int len) {
        IppStatus status = ippsMulC_64f_I(val,dst,len);
        checkStatus(status);
    }
    template<>
    void MulC<Ipp64fc>(Ipp64fc val, Ipp64fc* dst, int len) {
        IppStatus status = ippsMulC_64fc_I(val,dst,len);
        checkStatus(status);
    }

    template<typename T>
    void Mul(const T* a, T*b, T* dst, int len) {
        assert(1==0);
    }
    template<>
    void Mul<Ipp32f>(const Ipp32f * a, Ipp32f *b, Ipp32f* dst, int len) {
        IppStatus status = ippsMul_32f(a,b,dst,len);
        checkStatus(status);
    }
    template<>
    void Mul<Ipp32fc>(const Ipp32fc * a, Ipp32fc *b, Ipp32fc* dst, int len) {
        IppStatus status = ippsMul_32fc(a,b,dst,len);
        checkStatus(status);
    }
    template<>
    void Mul<Ipp64f>(const Ipp64f * a, Ipp64f *b, Ipp64f* dst, int len) {
        IppStatus status = ippsMul_64f(a,b,dst,len);
        checkStatus(status);
    }
    template<>
    void Mul<Ipp64fc>(const Ipp64fc * a, Ipp64fc *b, Ipp64fc* dst, int len) {
        IppStatus status = ippsMul_64fc(a,b,dst,len);
        checkStatus(status);
    }
    
    template<typename T>
    void Mul(const T* a, T* dst, int len) {
        assert(1==0);
    }
    template<>
    void Mul<Ipp32f>(const Ipp32f * a, Ipp32f* dst, int len) {
        IppStatus status = ippsMul_32f_I(a,dst,len);
        checkStatus(status);
    }
    template<>
    void Mul<Ipp32fc>(const Ipp32fc * a, Ipp32fc* dst, int len) {
        IppStatus status = ippsMul_32fc_I(a,dst,len);
        checkStatus(status);
    }    
    template<>
    void Mul<Ipp64f>(const Ipp64f * a, Ipp64f* dst, int len) {
        IppStatus status = ippsMul_64f_I(a,dst,len);
        checkStatus(status);
    }
    template<>
    void Mul<Ipp64fc>(const Ipp64fc * a, Ipp64fc* dst, int len) {
        IppStatus status = ippsMul_64fc_I(a,dst,len);
        checkStatus(status);
    }

    template<typename T>
    void DivC(const T* src, T val, T* dst, int len)
    {
        assert(1==0);
    }

    template<>
    void DivC<Ipp32f>(const Ipp32f * src, Ipp32f val, Ipp32f* dst, int len) {
        IppStatus status = ippsDivC_32f(src,val,dst,len);
        checkStatus(status);
    }
    template<>
    void DivC<Ipp32fc>(const Ipp32fc * src, Ipp32fc val, Ipp32fc* dst, int len) {
        IppStatus status = ippsDivC_32fc(src,val,dst,len);
        checkStatus(status);
    }
    template<>
    void DivC<Ipp64f>(const Ipp64f * src, Ipp64f val, Ipp64f* dst, int len) {
        IppStatus status = ippsDivC_64f(src,val,dst,len);
        checkStatus(status);
    }
    template<>
    void DivC<Ipp64fc>(const Ipp64fc * src, Ipp64fc val, Ipp64fc* dst, int len) {
        IppStatus status = ippsDivC_64fc(src,val,dst,len);
        checkStatus(status);
    }

    template<typename T>
    void DivC(T val, T* dst, int len) {
        assert(1==0);
    }

    template<>
    void DivC<Ipp32f>(Ipp32f val, Ipp32f* dst, int len) {
        IppStatus status = ippsDivC_32f_I(val,dst,len);
        checkStatus(status);
    }
    template<>
    void DivC<Ipp32fc>(Ipp32fc val, Ipp32fc* dst, int len) {
        IppStatus status = ippsDivC_32fc_I(val,dst,len);
        checkStatus(status);
    }
    template<>
    void DivC<Ipp64f>(Ipp64f val, Ipp64f* dst, int len) {
        IppStatus status = ippsDivC_64f_I(val,dst,len);
        checkStatus(status);
    }
    template<>
    void DivC<Ipp64fc>(Ipp64fc val, Ipp64fc* dst, int len) {
        IppStatus status = ippsDivC_64fc_I(val,dst,len);
        checkStatus(status);
    }

    template<typename T>
    void Div(const T* a, T *b, T* dst, int len) {
        assert(1==0);
    }
    template<>
    void Div<Ipp32f>(const Ipp32f * a, Ipp32f *b, Ipp32f* dst, int len) {
        IppStatus status = ippsDiv_32f(a,b,dst,len);
        checkStatus(status);
    }
    template<>
    void Div<Ipp32fc>(const Ipp32fc * a, Ipp32fc *b, Ipp32fc* dst, int len) {
        IppStatus status = ippsDiv_32fc(a,b,dst,len);
        checkStatus(status);
    }
    template<>
    void Div<Ipp64f>(const Ipp64f * a, Ipp64f *b, Ipp64f* dst, int len) {
        IppStatus status = ippsDiv_64f(a,b,dst,len);
        checkStatus(status);
    }
    template<>
    void Div<Ipp64fc>(const Ipp64fc * a, Ipp64fc *b, Ipp64fc* dst, int len) {
        IppStatus status = ippsDiv_64fc(a,b,dst,len);
        checkStatus(status);
    }

    template<typename T>
    void Div(const T* a, T* dst, int len) {
        assert(1==0);
    }
    template<>
    void Div<Ipp32f>(const Ipp32f * a, Ipp32f* dst, int len) {
        IppStatus status = ippsDiv_32f_I(a,dst,len);
        checkStatus(status);
    }
    template<>
    void Div<Ipp32fc>(const Ipp32fc * a, Ipp32fc* dst, int len) {
        IppStatus status = ippsDiv_32fc_I(a,dst,len);
        checkStatus(status);
    }
    template<>
    void Div<Ipp64f>(const Ipp64f * a, Ipp64f* dst, int len) {
        IppStatus status = ippsDiv_64f_I(a,dst,len);
        checkStatus(status);
    }
    template<>
    void Div<Ipp64fc>(const Ipp64fc * a, Ipp64fc* dst, int len) {
        IppStatus status = ippsDiv_64fc_I(a,dst,len);
        checkStatus(status);
    }

    template<typename T>
    void DivCRev(const T* src, T val, T* dst, int len) {
        assert(1==0);
    }

    template<>
    void DivCRev<Ipp32f>(const Ipp32f * src, Ipp32f val, Ipp32f* dst, int len) {
        IppStatus status = ippsDivCRev_32f(src,val,dst,len);
        checkStatus(status);
    }    
     
    template<typename T>
    void DivCRev(T val, T* dst, int len) {
        assert(1==0);
    }

    template<>
    void DivCRev<Ipp32f>(Ipp32f val, Ipp32f* dst, int len) {
        IppStatus status = ippsDivCRev_32f_I(val,dst,len);
        checkStatus(status);
    }    
    

    
    void AddProductC(const Ipp32f * array, int len, const Ipp32f val, Ipp32f * dst) {            
        IppStatus status = ippsAddProductC_32f(array,val,dst,len);
        checkStatus(status);        
    }
    
    template<typename T>
    void AddProduct(const T* array1, int len, const T* array2, T* dst) {
        assert(1==0);
    }
    template<>
    void AddProduct<Ipp32f>(const Ipp32f * array1, int len, const Ipp32f *array2, Ipp32f * dst) {            
        IppStatus status = ippsAddProduct_32f(array1,array2,dst,len);
        checkStatus(status);        
    }
    template<>
    void AddProduct<Ipp64f>(const Ipp64f * array1, int len, const Ipp64f *array2, Ipp64f * dst) {            
        IppStatus status = ippsAddProduct_64f(array1,array2,dst,len);
        checkStatus(status);        
    }
    template<>
    void AddProduct<Ipp32fc>(const Ipp32fc * array1, int len, const Ipp32fc *array2, Ipp32fc * dst) {            
        IppStatus status = ippsAddProduct_32fc(array1,array2,dst,len);
        checkStatus(status);        
    }
    template<>
    void AddProduct<Ipp64fc>(const Ipp64fc * array1, int len, const Ipp64fc *array2, Ipp64fc * dst) {            
        IppStatus status = ippsAddProduct_64fc(array1,array2,dst,len);
        checkStatus(status);        
    }

    template<typename T>
    void Abs(const T* src, T* dst, int len) {
        assert(1==0);
    }
    template<>
    void Abs<Ipp32f>(const Ipp32f * src, Ipp32f * dst, int len) {
        IppStatus status = ippsAbs_32f(src,dst,len);
        checkStatus(status);        
    }
    template<>
    void Abs<Ipp64f>(const Ipp64f * src, Ipp64f * dst, int len) {
        IppStatus status = ippsAbs_64f(src,dst,len);
        checkStatus(status);        
    }
    template<typename T>
    void Abs(T* dst, int len) {
        assert(1==0);
    }
    template<>
    void Abs<Ipp32f>( Ipp32f * dst, int len) {
        IppStatus status = ippsAbs_32f_I(dst,len);
        checkStatus(status);        
    }
    template<>
    void Abs<Ipp64f>(Ipp64f * dst, int len) {
        IppStatus status = ippsAbs_64f_I(dst,len);
        checkStatus(status);        
    }


    template<typename T>
    void Sqr(const T* src, T* dst, int len) {
        assert(1==0);
    }
    template<>
    void Sqr<Ipp32f>(const Ipp32f * src, Ipp32f * dst, int len) {
        IppStatus status = ippsSqr_32f(src,dst,len);
        checkStatus(status);        
    }
    template<>
    void Sqr<Ipp32fc>(const Ipp32fc * src, Ipp32fc * dst, int len) {
        IppStatus status = ippsSqr_32fc(src,dst,len);
        checkStatus(status);        
    }
    template<>
    void Sqr<Ipp64f>(const Ipp64f * src, Ipp64f * dst, int len) {
        IppStatus status = ippsSqr_64f(src,dst,len);
        checkStatus(status);        
    }
    template<>
    void Sqr<Ipp64fc>(const Ipp64fc * src, Ipp64fc * dst, int len) {
        IppStatus status = ippsSqr_64fc(src,dst,len);
        checkStatus(status);        
    }

    template<typename T>
    void Sqr(T* dst, int len) {
        assert(1==0);
    }
    template<>
    void Sqr<Ipp32f>(Ipp32f * dst, int len) {
        IppStatus status = ippsSqr_32f_I(dst,len);
        checkStatus(status);        
    }
    template<>
    void Sqr<Ipp32fc>(Ipp32fc * dst, int len) {
        IppStatus status = ippsSqr_32fc_I(dst,len);
        checkStatus(status);        
    }
    template<>
    void Sqr<Ipp64f>(Ipp64f * dst, int len) {
        IppStatus status = ippsSqr_64f_I(dst,len);
        checkStatus(status);        
    }
    template<>
    void Sqr<Ipp64fc>(Ipp64fc * dst, int len) {
        IppStatus status = ippsSqr_64fc_I(dst,len);
        checkStatus(status);        
    }

    template<typename T>
    void Sqrt(const T * src, T * dst, int len) {
        assert(1==0);
    }
    template<>
    void Sqrt<Ipp32f>(const Ipp32f * src, Ipp32f * dst, int len) {
        IppStatus status = ippsSqrt_32f(src,dst,len);
        checkStatus(status);        
    }
    template<>
    void Sqrt<Ipp32fc>(const Ipp32fc * src, Ipp32fc * dst, int len) {
        IppStatus status = ippsSqrt_32fc(src,dst,len);
        checkStatus(status);        
    }
    template<>
    void Sqrt<Ipp64f>(const Ipp64f * src, Ipp64f * dst, int len) {
        IppStatus status = ippsSqrt_64f(src,dst,len);
        checkStatus(status);        
    }
    template<>
    void Sqrt<Ipp64fc>(const Ipp64fc * src, Ipp64fc * dst, int len) {
        IppStatus status = ippsSqrt_64fc(src,dst,len);
        checkStatus(status);        
    }

    template<typename T>
    void Sqrt(T * dst, int len) {
        assert(1==0);
    }
    template<>
    void Sqrt<Ipp32f>(Ipp32f * dst, int len) {
        IppStatus status = ippsSqrt_32f_I(dst,len);
        checkStatus(status);        
    }
    template<>
    void Sqrt<Ipp32fc>(Ipp32fc * dst, int len) {
        IppStatus status = ippsSqrt_32fc_I(dst,len);
        checkStatus(status);        
    }
    template<>
    void Sqrt<Ipp64f>(Ipp64f * dst, int len) {
        IppStatus status = ippsSqrt_64f_I(dst,len);
        checkStatus(status);        
    }
    template<>
    void Sqrt<Ipp64fc>(Ipp64fc * dst, int len) {
        IppStatus status = ippsSqrt_64fc_I(dst,len);
        checkStatus(status);        
    }

    void Cubrt(const Ipp32f * src, Ipp32f * dst, int len) {
        IppStatus status = ippsCubrt_32f(src,dst,len);
        checkStatus(status);        
    }

    template<typename T>
    void Exp(const T * src, T * dst, int len) {
        assert(1==0);
    }

    template<>
    void Exp<Ipp32f>(const Ipp32f * src, Ipp32f * dst, int len) {
        IppStatus status = ippsExp_32f(src,dst,len);
        checkStatus(status);        
    }
    template<>
    void Exp<Ipp64f>(const Ipp64f * src, Ipp64f * dst, int len) {
        IppStatus status = ippsExp_64f(src,dst,len);
        checkStatus(status);        
    }
    
    template<typename T>
    void Exp( T * dst, int len) {
        assert(1==0);
    }
    template<>
    void Exp<Ipp32f>(Ipp32f * dst, int len) {
        IppStatus status = ippsExp_32f_I(dst,len);
        checkStatus(status);        
    }
    template<>
    void Exp<Ipp64f>(Ipp64f * dst, int len) {
        IppStatus status = ippsExp_64f_I(dst,len);
        checkStatus(status);        
    }
    
    template<typename T>
    void Ln(const T* src, T* dst, int len) {
        assert(1==0);
    }
    template<>
    void Ln<Ipp32f>(const Ipp32f * src, Ipp32f * dst, int len) {
        IppStatus status = ippsLn_32f(src,dst,len);
        checkStatus(status);        
    }
    template<>
    void Ln<Ipp64f>(const Ipp64f * src, Ipp64f * dst, int len) {
        IppStatus status = ippsLn_64f(src,dst,len);
        checkStatus(status);        
    }
    template<typename T>
    void Ln(T* dst, int len) {
        assert(1==0);
    }
    template<>
    void Ln<Ipp32f>(Ipp32f * dst, int len) {
        IppStatus status = ippsLn_32f_I(dst,len);
        checkStatus(status);        
    }
    template<>
    void Ln<Ipp64f>(Ipp64f * dst, int len) {
        IppStatus status = ippsLn_64f_I(dst,len);
        checkStatus(status);        
    }

    template<typename T>
    void SumLn(const T* src, T* dst, int len) {
        assert(1==0);
    }
    template<>
    void SumLn<Ipp32f>(const Ipp32f * src, Ipp32f * dst, int len) {
        IppStatus status = ippsLn_32f(src,dst,len);
        checkStatus(status);        
    }
    template<>
    void SumLn<Ipp64f>(const Ipp64f * src, Ipp64f * dst, int len) {
        IppStatus status = ippsLn_64f(src,dst,len);
        checkStatus(status);        
    }

    template<typename T>
    void Arctan(const T* src, T* dst, int len) {
        assert(1==0);
    }
    template<>
    void Arctan<Ipp32f>(const Ipp32f * src, Ipp32f * dst, int len) {
        IppStatus status = ippsArctan_32f(src,dst,len);
        checkStatus(status);        
    }
    template<>
    void Arctan<Ipp64f>(const Ipp64f * src, Ipp64f * dst, int len) {
        IppStatus status = ippsArctan_64f(src,dst,len);
        checkStatus(status);        
    }
    template<typename T>
    void Arctan(T* dst, int len) {
        assert(1==0);
    }
    template<>
    void Arctan<Ipp32f>(Ipp32f * dst, int len) {
        IppStatus status = ippsArctan_32f_I(dst,len);
        checkStatus(status);        
    }
    template<>
    void Arctan<Ipp64f>(Ipp64f * dst, int len) {
        IppStatus status = ippsArctan_64f_I(dst,len);
        checkStatus(status);        
    }

    template<typename T>
    T Max(const T* src,int len) {
        assert(1==0);
    }

    template<>
    Ipp32f Max<Ipp32f>(const Ipp32f * src,int len) {
        Ipp32f max=0;
        IppStatus status = ippsMax_32f(src,len,&max);
        checkStatus(status);        
        return max;
    }
    template<>
    Ipp64f Max<Ipp64f>(const Ipp64f * src,int len) {
        Ipp64f max=0;
        IppStatus status = ippsMax_64f(src,len,&max);
        checkStatus(status);        
        return max;
    }
    
    template<typename T>
    T Min(const T* src,int len) {
        assert(1==0);
    }

    template<>
    Ipp32f Min<Ipp32f>(const Ipp32f * src,int len) {
        Ipp32f min=0;
        IppStatus status = ippsMin_32f(src,len,&min);
        checkStatus(status);        
        return min;
    }
    template<>
    Ipp64f Min<Ipp64f>(const Ipp64f * src,int len) {
        Ipp64f min=0;
        IppStatus status = ippsMin_64f(src,len,&min);
        checkStatus(status);        
        return min;
    }

    template<typename T>
    int MaxIndex(const T* src, T * max, int len) {
        assert(1==0);
    }
    template<>
    int MaxIndex<Ipp32f>(const Ipp32f * src,Ipp32f * max, int len) {
        int index=-1;
        IppStatus status = ippsMaxIndx_32f(src,len,max,&index);
        checkStatus(status);        
        return index;
    }
    template<>
    int MaxIndex<Ipp64f>(const Ipp64f * src,Ipp64f * max, int len) {
        int index=-1;
        IppStatus status = ippsMaxIndx_64f(src,len,max,&index);
        checkStatus(status);        
        return index;
    }

    template<typename T>
    int MinIndex(const T * src, int len, T * min) {
        assert(1==0);
    }

    template<>
    int MinIndex<Ipp32f>(const Ipp32f * src, int len, Ipp32f * min) {
        int index=-1;
        IppStatus status = ippsMinIndx_32f(src,len,min,&index);
        checkStatus(status);        
        return index;
    }
    template<>
    int MinIndex<Ipp64f>(const Ipp64f * src, int len, Ipp64f * min) {
        int index=-1;
        IppStatus status = ippsMinIndx_64f(src,len,min,&index);
        checkStatus(status);        
        return index;
    }

    template<typename T>
    T MaxAbs(const T* src, int len) {
        assert(1==0);
    }

    template<>
    Ipp32f MaxAbs<Ipp32f>(const Ipp32f * src, int len) {
        Ipp32f max=0;
        IppStatus status = ippsMaxAbs_32f(src,len,&max);
        checkStatus(status);        
        return max;
    }
    template<>
    Ipp64f MaxAbs<Ipp64f>(const Ipp64f * src, int len) {
        Ipp64f max=0;
        IppStatus status = ippsMaxAbs_64f(src,len,&max);
        checkStatus(status);        
        return max;
    }

    template<typename T>
    T MinAbs(const T * src, int len) {
        assert(1==0);
    }

    template<>
    Ipp32f MinAbs<Ipp32f>(const Ipp32f * src, int len) {
        Ipp32f min=0;
        IppStatus status = ippsMinAbs_32f(src,len,&min);
        checkStatus(status);        
        return min;
    }
    template<>
    Ipp64f MinAbs<Ipp64f>(const Ipp64f * src, int len) {
        Ipp64f min=0;
        IppStatus status = ippsMinAbs_64f(src,len,&min);
        checkStatus(status);        
        return min;
    }

    template<typename T>
    int MaxAbsIndex(const T* src, T* max, int len) {
        assert(1==0);
    }

    template<>
    int MaxAbsIndex<Ipp32f>(const Ipp32f * src, Ipp32f * max, int len) {
        int index=-1;
        IppStatus status = ippsMaxAbsIndx_32f(src,len,max, &index);
        checkStatus(status);        
        return index;
    }
    template<>
    int MaxAbsIndex<Ipp64f>(const Ipp64f * src, Ipp64f * max, int len) {
        int index=-1;
        IppStatus status = ippsMaxAbsIndx_64f(src,len,max, &index);
        checkStatus(status);        
        return index;
    }

    template<typename T>
    int MinAbsIndex(const T* src, T* min, int len) {
        assert(1==0);
    }
    template<>
    int MinAbsIndex<Ipp32f>(const Ipp32f * src, Ipp32f * min, int len) {
        int index=-1;
        IppStatus status = ippsMinAbsIndx_32f(src,len,min,&index);
        checkStatus(status);        
        return index;
    }
    
    template<typename T>
    void MinMax(const T* src, int len, T* min, T* max) {    
        assert(1==0);
    }
    template<>
    void MinMax<Ipp32f>(const Ipp32f * src, int len, Ipp32f * min, Ipp32f * max) {    
        IppStatus status = ippsMinMax_32f(src,len,min,max);
        checkStatus(status);            
    }

    template<typename T>
    void MinMaxIndex(const T* src, int len, T *min, int * min_index, T *max, int * max_index) {
        assert(1==0);
    }
        
    template<>
    void MinMaxIndex<Ipp32f>(const Ipp32f * src, int len, Ipp32f *min, int * min_index, Ipp32f *max, int * max_index) {
        IppStatus status = ippsMinMaxIndx_32f(src,len,min,min_index,max,max_index);
        checkStatus(status);            
    }


    void ReplaceNAN(Ipp32f * p, Ipp32f v, int len) {
        IppStatus status = ippsReplaceNAN_32f_I(p,len,v);
        checkStatus(status);            
    }
    
    template<typename T>
    T Mean(const T * src, int len) {
        assert(1==0);
    }

    template<>
    Ipp32f Mean<Ipp32f>(const Ipp32f * src, int len)
    {
        Ipp32f mean;
        IppStatus status = ippsMean_32f(src,len,&mean,ippAlgHintNone);
        checkStatus(status);            
        return mean;
    }
    template<>
    Ipp64f Mean<Ipp64f>(const Ipp64f * src, int len)
    {
        Ipp64f mean;
        IppStatus status = ippsMean_64f(src,len,&mean);
        checkStatus(status);            
        return mean;
    }

    template<typename T>
    T StdDev(const T* src, int len)
    {
        assert(1==0);
    }

    template<>
    Ipp32f StdDev<Ipp32f>(const Ipp32f * src, int len)
    {
        Ipp32f mean;
        IppStatus status = ippsStdDev_32f(src,len,&mean,ippAlgHintNone);
        checkStatus(status);            
        return mean;
    }
    template<>
    Ipp64f StdDev<Ipp64f>(const Ipp64f * src, int len)
    {
        Ipp64f mean;
        IppStatus status = ippsStdDev_64f(src,len,&mean);
        checkStatus(status);            
        return mean;
    }

    template<typename T>
    void MeanStdDev(const T* src, int len, T* mean, T* dev)
    {
        assert(1==0);
    }
    template<>
    void MeanStdDev<Ipp32f>(const Ipp32f * src, int len, Ipp32f * mean, Ipp32f * dev)
    {        
        IppStatus status = ippsMeanStdDev_32f(src,len,mean,dev,ippAlgHintNone);
        checkStatus(status);                    
    }
    template<>
    void MeanStdDev<Ipp64f>(const Ipp64f * src, int len, Ipp64f * mean, Ipp64f * dev)
    {        
        IppStatus status = ippsMeanStdDev_64f(src,len,mean,dev);
        checkStatus(status);                    
    }

    template<typename T>
    T NormInf(const T* src, int len) {
        assert(1==0);
    }
    template<>
    Ipp32f NormInf<Ipp32f>(const Ipp32f * src, int len) {
        Ipp32f norm;
        IppStatus status = ippsNorm_Inf_32f(src,len,&norm);
        checkStatus(status);                    
        return norm;
    }
    template<>
    Ipp64f NormInf<Ipp64f>(const Ipp64f * src, int len) {
        Ipp64f norm;
        IppStatus status = ippsNorm_Inf_64f(src,len,&norm);
        checkStatus(status);                    
        return norm;
    }

    template<typename T>
    T NormL1(const T* src, int len) {
        assert(1==0);
    }

    template<>
    Ipp32f NormL1<Ipp32f>(const Ipp32f * src, int len) {
        Ipp32f norm;
        IppStatus status = ippsNorm_L1_32f(src,len,&norm);
        checkStatus(status);                    
        return norm;
    }
    template<>
    Ipp64f NormL1<Ipp64f>(const Ipp64f * src, int len) {
        Ipp64f norm;
        IppStatus status = ippsNorm_L1_64f(src,len,&norm);
        checkStatus(status);                    
        return norm;
    }

    template<typename T>
    T NormL2(const T* src, int len) {
        assert(1==0);
    }

    
    template<>
    Ipp32f NormL2<Ipp32f>(const Ipp32f * src, int len) {
        Ipp32f norm;
        IppStatus status = ippsNorm_L2_32f(src,len,&norm);
        checkStatus(status);                    
        return norm;
    }
    template<>
    Ipp64f NormL2<Ipp64f>(const Ipp64f * src, int len) {
        Ipp64f norm;
        IppStatus status = ippsNorm_L2_64f(src,len,&norm);
        checkStatus(status);                    
        return norm;
    }

    template<typename T>
    T NormDiffInf(const T* src1,const T* src2, int len) {
        assert(1==0);
    }

    template<>
    Ipp32f NormDiffInf<Ipp32f>(const Ipp32f * src1,const Ipp32f * src2, int len) {
        Ipp32f norm;
        IppStatus status = ippsNormDiff_Inf_32f(src1,src2,len,&norm);
        checkStatus(status);                    
        return norm;
    }
    template<>
    Ipp64f NormDiffInf<Ipp64f>(const Ipp64f * src1,const Ipp64f * src2, int len) {
        Ipp64f norm;
        IppStatus status = ippsNormDiff_Inf_64f(src1,src2,len,&norm);
        checkStatus(status);                    
        return norm;
    }

    template<typename T>
    T NormDiffL1(const T* src1, const T* src2, int len) {
        assert(1==0);
    }
    template<>
    Ipp32f NormDiffL1<Ipp32f>(const Ipp32f * src1, const Ipp32f * src2, int len) {
        Ipp32f norm;
        IppStatus status = ippsNormDiff_L1_32f(src1,src2,len,&norm);
        checkStatus(status);                    
        return norm;
    }
    template<>
    Ipp64f NormDiffL1<Ipp64f>(const Ipp64f * src1, const Ipp64f * src2, int len) {
        Ipp64f norm;
        IppStatus status = ippsNormDiff_L1_64f(src1,src2,len,&norm);
        checkStatus(status);                    
        return norm;
    }

    template<typename T>
    T NormDiffL2(const T* src1, const T* src2, int len) {
        assert(1==0);
    }
    template<>
    Ipp32f NormDiffL2<Ipp32f>(const Ipp32f * src1, const Ipp32f * src2, int len) {
        Ipp32f norm;
        IppStatus status = ippsNormDiff_L2_32f(src1,src2,len,&norm);
        checkStatus(status);                    
        return norm;
    }
    template<>
    Ipp64f NormDiffL2<Ipp64f>(const Ipp64f * src1, const Ipp64f * src2, int len) {
        Ipp64f norm;
        IppStatus status = ippsNormDiff_L2_64f(src1,src2,len,&norm);
        checkStatus(status);                    
        return norm;
    }

    template<typename T>
    T DotProduct(const T * src1, const T * src2, int len) 
    {
        assert(1==0);
    }

    template<>
    Ipp32f DotProduct<Ipp32f>(const Ipp32f * src1, const Ipp32f * src2, int len) {
        Ipp32f dp;
        IppStatus status = ippsDotProd_32f(src1,src2,len,&dp);
        checkStatus(status);     
        return dp;
    }
    template<>
    Ipp64f DotProduct<Ipp64f>(const Ipp64f * src1, const Ipp64f * src2, int len) {
        Ipp64f dp;
        IppStatus status = ippsDotProd_64f(src1,src2,len,&dp);
        checkStatus(status);     
        return dp;
    }

    void MaxEvery(const Ipp32f * src1, const Ipp32f * src2, Ipp32f * dst, int len) {        
        IppStatus status = ippsMaxEvery_32f(src1,src2,dst,len);
        checkStatus(status);             
    }
    void MinEvery(const Ipp32f * src1, const Ipp32f * src2, Ipp32f * dst, int len) {        
        IppStatus status = ippsMinEvery_32f(src1,src2,dst,len);
        checkStatus(status);             
    }
 
    Ipp32f ZeroCrossing(const Ipp32f * src1, int len, IppsZCType zcType = ippZCR) {
        Ipp32f zc;
        IppStatus status = ippsZeroCrossing_32f(src1,len,&zc,zcType);
        checkStatus(status);     
        return zc;        
    }
    /*
    Ipp64f ZeroCrossing(const Ipp64f * src1, int len, IppsZCType zcType = ippZCR) {
        Ipp64f zc;
        IppStatus status = ippsZeroCrossing_64f(src1,len,&zc,zcType);
        checkStatus(status);     
        return zc;
    }
    */
    template<typename T>
    T MSE(const T * thetaAbs, const T * thetaModel, int len) {
        T r = 0;
        for(size_t i = 0; i <  len; i++)
            r += pow(thetaAbs[i] - thetaModel[i],2.0);
        return (1.0f / (T)len) * r;
    }
    template<typename T>
    T RMSE(const T * thetaAbs, const T * thetaModel, int len) {
        T r = MSE(thetaAbs,thetaModel,len);
        return sqrt(r);
    }
    template<typename T>
    T MeanSquare(const T * x, int len) {
        T r = 0;
        for(size_t i = 0; i < len; i++) r += pow(x[i],2.0);
        return (1.0f / (T)len) * r;
    }

    template<typename T>
    T AmpToDb(T a) {
        return pow(10.0,a/20.0);
    }
    template<typename T>
    T DbToAmp(T a) {
        return 20.0*log10(a);
    }
    void Tone(Ipp32f * array, int len, Ipp32f mag, Ipp32f freq, Ipp32f * phase) {
        IppStatus status = ippsTone_32f(array,len,mag,freq,phase,ippAlgHintNone);
        checkStatus(status);
    }
    void Triangle(Ipp32f * array, int len, Ipp32f m, Ipp32f f, Ipp32f a, Ipp32f * p) {
        IppStatus status = ippsTriangle_32f(array,len,m,f,a,p);
        checkStatus(status);
    }
    void Sort(Ipp32f * array, int len, int dir=1) {
        IppStatus status;
        if(dir >= 0)
            status = ippsSortAscend_32f_I(array,len);
        else
            status = ippsSortDescend_32f_I(array,len);
        checkStatus(status);
    }
    
    template<typename T>
    void UpSample(const T * pSrc, int srcLen, T* pDst, int *dstLen, int fact, int * pPhase)
    {
        throw std::runtime_error("");
    }

    template<>
    void UpSample<Ipp32f>(const Ipp32f * pSrc, int srcLen, Ipp32f* pDst, int *dstLen, int fact, int * pPhase) {
        IppStatus status = ippsSampleUp_32f(pSrc,srcLen,pDst,dstLen,fact,pPhase);
        checkStatus(status);
    }
    template<>
    void UpSample<Ipp32fc>(const Ipp32fc * pSrc, int srcLen, Ipp32fc* pDst, int *dstLen, int fact, int * pPhase) {
        IppStatus status = ippsSampleUp_32fc(pSrc,srcLen,pDst,dstLen,fact,pPhase);
        checkStatus(status);
    }
    template<>
    void UpSample<Ipp64f>(const Ipp64f * pSrc, int srcLen, Ipp64f* pDst, int *dstLen, int fact, int * pPhase) {
        IppStatus status = ippsSampleUp_64f(pSrc,srcLen,pDst,dstLen,fact,pPhase);
        checkStatus(status);
    }
    template<>
    void UpSample<Ipp64fc>(const Ipp64fc * pSrc, int srcLen, Ipp64fc* pDst, int *dstLen, int fact, int * pPhase) {
        IppStatus status = ippsSampleUp_64fc(pSrc,srcLen,pDst,dstLen,fact,pPhase);
        checkStatus(status);
    }

    template<typename T>
    void DownSample(const T * pSrc, int srcLen, T* pDst, int *dstLen, int fact, int * pPhase) {
        throw std::runtime_error("");
    }

    template<>
    void DownSample<Ipp32f>(const Ipp32f * pSrc, int srcLen, Ipp32f* pDst, int *dstLen, int fact, int * pPhase) {
        IppStatus status = ippsSampleDown_32f(pSrc,srcLen,pDst,dstLen,fact,pPhase);
        checkStatus(status);
    }
    template<>
    void DownSample<Ipp32fc>(const Ipp32fc * pSrc, int srcLen, Ipp32fc* pDst, int *dstLen, int fact, int * pPhase) {
        IppStatus status = ippsSampleDown_32fc(pSrc,srcLen,pDst,dstLen,fact,pPhase);
        checkStatus(status);
    }
    template<>
    void DownSample<Ipp64f>(const Ipp64f * pSrc, int srcLen, Ipp64f* pDst, int *dstLen, int fact, int * pPhase) {
        IppStatus status = ippsSampleDown_64f(pSrc,srcLen,pDst,dstLen,fact,pPhase);
        checkStatus(status);
    }
    template<>
    void DownSample<Ipp64fc>(const Ipp64fc * pSrc, int srcLen, Ipp64fc* pDst, int *dstLen, int fact, int * pPhase) {
        IppStatus status = ippsSampleDown_64fc(pSrc,srcLen,pDst,dstLen,fact,pPhase);
        checkStatus(status);
    }

    template<typename T>
    T kaiserBeta(T As)
    {
        if(As > (T)50.0)
            return (T)0.1102 * (As - (T)8.7);
        else if(As >= (T)21.)
            return (T)0.5842 * std::pow(As - (T)21.0, (T)0.4) + (T)0.07886 * (As - (T)21.0);
        else
            return (T)0.;
    }

    
    void Goertzal(int len, const Ipp32fc * src, Ipp32fc * dst, Ipp32f freq) {
        IppStatus status = ippsGoertz_32fc(src,len,dst,freq);
    }
    
    void Goertzal(int len, const Ipp32f * src, Ipp32fc * dst, Ipp32f freq) {
        IppStatus status = ippsGoertz_32f(src,len,dst,freq);
    }

    void DoubleToFloat(size_t len, double * src, Ipp32f * dst) {
        for(size_t i = 0; i < len; i++) dst[i] = (Ipp32f)src[i];
    }
    void FloatToDouble(size_t len, Ipp32f * src, double * dst) {
        for(size_t i = 0; i < len; i++) dst[i] = (double)src[i];
    }

    template<typename T1, typename T2>
    void ConjPack(const T1* pSrc, T2* pDst, int dstLen)
    {
        assert(1==0);
    }
    template<>
    void ConjPack<Ipp32f,Ipp32fc>(const Ipp32f* pSrc, Ipp32fc * pDst, int dstLen)
    {
        IppStatus status = ippsConjPack_32fc(pSrc,pDst,dstLen);
        checkStatus(status);
    }
    template<>
    void ConjPack<Ipp64f,Ipp64fc>(const Ipp64f* pSrc, Ipp64fc * pDst, int dstLen)
    {
        IppStatus status = ippsConjPack_64fc(pSrc,pDst,dstLen);
        checkStatus(status);
    }

    template<typename T1, typename T2>
    void ConjPerm(const T1* pSrc, T2* pDst, int dstLen)
    {
        assert(1==0);
    }
    template<>
    void ConjPerm<Ipp32f,Ipp32fc>(const Ipp32f* pSrc, Ipp32fc * pDst, int dstLen)
    {
        IppStatus status = ippsConjPerm_32fc(pSrc,pDst,dstLen);
        checkStatus(status);
    }
    template<>
    void ConjPerm<Ipp64f,Ipp64fc>(const Ipp64f* pSrc, Ipp64fc * pDst, int dstLen)
    {
        IppStatus status = ippsConjPerm_64fc(pSrc,pDst,dstLen);
        checkStatus(status);
    }

    template<typename T1, typename T2>
    void ConjCcs(const T1* pSrc, T2* pDst, int dstLen)
    {
        assert(1==0);
    }
    template<>
    void ConjCcs<Ipp32f,Ipp32fc>(const Ipp32f* pSrc, Ipp32fc * pDst, int dstLen)
    {
        IppStatus status = ippsConjCcs_32fc(pSrc,pDst,dstLen);
        checkStatus(status);
    }
    template<>
    void ConjCcs<Ipp64f,Ipp64fc>(const Ipp64f* pSrc, Ipp64fc * pDst, int dstLen)
    {
        IppStatus status = ippsConjCcs_64fc(pSrc,pDst,dstLen);
        checkStatus(status);
    }


    template<typename T>
    void MulPack(const T* pSrc1, const T * pSrc2, T* pDst, int dstLen)
    {
        assert(1==0);
    }
    template<>
    void MulPack<Ipp32f>(const Ipp32f* pSrc1, const Ipp32f * pSrc2, Ipp32f * pDst, int dstLen)
    {
        IppStatus status = ippsMulPack_32f(pSrc1,pSrc2,pDst,dstLen);
        checkStatus(status);
    }
    template<>
    void MulPack<Ipp64f>(const Ipp64f* pSrc1, const Ipp64f * pSrc2, Ipp64f * pDst, int dstLen)
    {
        IppStatus status = ippsMulPack_64f(pSrc1,pSrc2,pDst,dstLen);
        checkStatus(status);
    }

    template<typename T>
    void MulPerm(const T* pSrc1, const T* pSrc2, T* pDst, int dstLen)
    {
        assert(1==0);
    }
    template<>
    void MulPerm<Ipp32f>(const Ipp32f* pSrc1, const Ipp32f* pSrc2, Ipp32f * pDst, int dstLen)
    {
        IppStatus status = ippsMulPerm_32f(pSrc1,pSrc2,pDst,dstLen);
        checkStatus(status);
    }
    template<>
    void MulPerm<Ipp64f>(const Ipp64f* pSrc1, const Ipp64f* pSrc2, Ipp64f * pDst, int dstLen)
    {
        IppStatus status = ippsMulPerm_64f(pSrc1,pSrc2,pDst,dstLen);
        checkStatus(status);
    }

    template<typename T>
    void MulPackConj(const T* pSrc1, T* pDst, int dstLen)
    {
        assert(1==0);
    }
    template<>
    void MulPackConj<Ipp32f>(const Ipp32f* pSrc1,Ipp32f * pDst, int dstLen)
    {
        IppStatus status = ippsMulPackConj_32f_I(pSrc1,pDst,dstLen);
        checkStatus(status);
    }
    template<>
    void MulPackConj<Ipp64f>(const Ipp64f* pSrc1,Ipp64f * pDst, int dstLen)
    {
        IppStatus status = ippsMulPackConj_64f_I(pSrc1,pDst,dstLen);
        checkStatus(status);
    }

    template<typename T>
    IppDataType GetDataType()
    {
        IppDataType dataType;
        if(typeid(T) == typeid(Ipp32f)) dataType = ipp32f;
        else if(typeid(T) == typeid(Ipp64f)) dataType = ipp64f;
        else if(typeid(T) == typeid(Ipp32fc)) dataType = ipp32fc;
        else if(typeid(T) == typeid(Ipp64fc)) dataType = ipp64fc;
        else throw std::runtime_error("Type not supported yet");
        return dataType;
    }

}