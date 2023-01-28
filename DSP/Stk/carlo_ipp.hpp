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

#include "ipp.hpp"

// this is just std::vector
#include "Samples.hpp"
#include "SamplesDSP.hpp"

// this is MKL stuff
//#include "mkl_array.hpp"
//#include "mkl_vector.hpp"
//#include "mkl_matrix.hpp"

// we can do SIMD kernels with these
/*
#include "vector2d.hpp"
#include "vector4d.hpp"
#include "vector4f.hpp"
#include "vector8f.hpp"
*/
// Eigen
//#include "TinyEigen.hpp"
//#include "SampleVector.h"

namespace Casino::IPP
{
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
    // Intel IPP
    inline void Move(const Ipp8u* pSrc, Ipp8u* pDst, int len) {
        IppStatus status = ippsMove_8u(pSrc,pDst,len);
        checkStatus(status);
    }
    inline void Move(const Ipp16s* pSrc, Ipp16s* pDst, int len) {
        IppStatus status = ippsMove_16s(pSrc,pDst,len);
        checkStatus(status);
    }
    inline void Move(const Ipp32s* pSrc, Ipp32s* pDst, int len) {
        IppStatus status = ippsMove_32s(pSrc,pDst,len);
        checkStatus(status);
    }
    inline void Move(const Ipp32f* pSrc, Ipp32f* pDst, int len) {
        IppStatus status = ippsMove_32f(pSrc,pDst,len);
        checkStatus(status);
    }
    inline void Move(const Ipp64s* pSrc, Ipp64s* pDst, int len) {
        IppStatus status = ippsMove_64s(pSrc,pDst,len);
        checkStatus(status);
    }
    inline void Move(const Ipp64f* pSrc, Ipp64f* pDst, int len) {
        IppStatus status = ippsMove_64f(pSrc,pDst,len);
        checkStatus(status);
    }
    inline void Move(const Ipp16sc* pSrc, Ipp16sc* pDst, int len) {
        IppStatus status = ippsMove_16sc(pSrc,pDst,len);
        checkStatus(status);
    }
    inline void Move(const Ipp32sc* pSrc, Ipp32sc* pDst, int len) {
        IppStatus status = ippsMove_32sc(pSrc,pDst,len);
        checkStatus(status);
    }
    inline void Move(const Ipp32fc* pSrc, Ipp32fc* pDst, int len) {
        IppStatus status = ippsMove_32fc(pSrc,pDst,len);
        checkStatus(status);
    }
    inline void Move(const Ipp64sc* pSrc, Ipp64sc* pDst, int len) {
        IppStatus status = ippsMove_64sc(pSrc,pDst,len);
        checkStatus(status);
    }
    inline void Move(const Ipp64fc* pSrc, Ipp64fc* pDst, int len) {
        IppStatus status = ippsMove_64fc(pSrc,pDst,len);
        checkStatus(status);
    }
        inline void Copy(const Ipp8u* pSrc, Ipp8u* pDst, int len) {
        IppStatus status = ippsCopy_8u(pSrc,pDst,len);
        checkStatus(status);
    }
    inline void Copy(const Ipp16s* pSrc, Ipp16s* pDst, int len) {
        IppStatus status = ippsCopy_16s(pSrc,pDst,len);
        checkStatus(status);
    }
    inline void Copy(const Ipp32s* pSrc, Ipp32s* pDst, int len) {
        IppStatus status = ippsCopy_32s(pSrc,pDst,len);
        checkStatus(status);
    }
    inline void Copy(const Ipp32f* pSrc, Ipp32f* pDst, int len) {
        IppStatus status = ippsCopy_32f(pSrc,pDst,len);
        checkStatus(status);
    }
    inline void Copy(const Ipp64s* pSrc, Ipp64s* pDst, int len) {
        IppStatus status = ippsCopy_64s(pSrc,pDst,len);
        checkStatus(status);
    }
    inline void Copy(const Ipp64f* pSrc, Ipp64f* pDst, int len) {
        IppStatus status = ippsCopy_64f(pSrc,pDst,len);
        checkStatus(status);
    }
    inline void Copy(const Ipp16sc* pSrc, Ipp16sc* pDst, int len) {
        IppStatus status = ippsCopy_16sc(pSrc,pDst,len);
        checkStatus(status);
    }
    inline void Copy(const Ipp32sc* pSrc, Ipp32sc* pDst, int len) {
        IppStatus status = ippsCopy_32sc(pSrc,pDst,len);
        checkStatus(status);
    }
    inline void Copy(const Ipp32fc* pSrc, Ipp32fc* pDst, int len) {
        IppStatus status = ippsCopy_32fc(pSrc,pDst,len);
        checkStatus(status);
    }
    inline void Copy(const Ipp64sc* pSrc, Ipp64sc* pDst, int len) {
        IppStatus status = ippsCopy_64sc(pSrc,pDst,len);
        checkStatus(status);
    }
    inline void Copy(const Ipp64fc* pSrc, Ipp64fc* pDst, int len) {
        IppStatus status = ippsCopy_64fc(pSrc,pDst,len);
        checkStatus(status);
    }
    inline void Set(const Ipp8u val, Ipp8u* pDst, int len) {
        IppStatus status = ippsSet_8u(val,pDst,len);
        checkStatus(status);
    }
    inline void Set(const Ipp16s val, Ipp16s* pDst, int len) {
        IppStatus status = ippsSet_16s(val,pDst,len);
        checkStatus(status);
    }
    inline void Set(const Ipp32s val, Ipp32s* pDst, int len) {
        IppStatus status = ippsSet_32s(val,pDst,len);
        checkStatus(status);
    }
    inline void Set(const Ipp32f val, Ipp32f* pDst, int len) {
        IppStatus status = ippsSet_32f(val,pDst,len);
        checkStatus(status);
    }
    inline void Set(const Ipp64s val, Ipp64s* pDst, int len) {
        IppStatus status = ippsSet_64s(val,pDst,len);
        checkStatus(status);
    }
    inline void Set(const Ipp64f val, Ipp64f* pDst, int len) {
        IppStatus status = ippsSet_64f(val,pDst,len);
        checkStatus(status);
    }
    inline void Set(const Ipp16sc val, Ipp16sc* pDst, int len) {
        IppStatus status = ippsSet_16sc(val,pDst,len);
        checkStatus(status);
    }
    inline void Set(const Ipp32sc val, Ipp32sc* pDst, int len) {
        IppStatus status = ippsSet_32sc(val,pDst,len);
        checkStatus(status);
    }
    inline void Set(const Ipp32fc val, Ipp32fc* pDst, int len) {
        IppStatus status = ippsSet_32fc(val,pDst,len);
        checkStatus(status);
    }
    inline void Set(const Ipp64sc val, Ipp64sc* pDst, int len) {
        IppStatus status = ippsSet_64sc(val,pDst,len);
        checkStatus(status);
    }
    inline void Set(const Ipp64fc val, Ipp64fc* pDst, int len) {
        IppStatus status = ippsSet_64fc(val,pDst,len);
        checkStatus(status);
    }
    inline void Zero(Ipp8u* pDst, int len) {
        IppStatus status = ippsZero_8u(pDst,len);
        checkStatus(status);
    }
    inline void Zero(Ipp16s* pDst, int len) {
        IppStatus status = ippsZero_16s(pDst,len);
        checkStatus(status);
    }
    inline void Zero(Ipp32s* pDst, int len) {
        IppStatus status = ippsZero_32s(pDst,len);
        checkStatus(status);
    }
    inline void Zero(Ipp32f* pDst, int len) {
        IppStatus status = ippsZero_32f(pDst,len);
        checkStatus(status);
    }
    inline void Zero(Ipp64s* pDst, int len) {
        IppStatus status = ippsZero_64s(pDst,len);
        checkStatus(status);
    }
    inline void Zero(Ipp64f* pDst, int len) {
        IppStatus status = ippsZero_64f(pDst,len);
        checkStatus(status);
    }
    inline void Zero(Ipp16sc* pDst, int len) {
        IppStatus status = ippsZero_16sc(pDst,len);
        checkStatus(status);
    }
    inline void Zero(Ipp32sc* pDst, int len) {
        IppStatus status = ippsZero_32sc(pDst,len);
        checkStatus(status);
    }
    inline void Zero(Ipp32fc* pDst, int len) {
        IppStatus status = ippsZero_32fc(pDst,len);
        checkStatus(status);
    }
    inline void Zero(Ipp64sc* pDst, int len) {
        IppStatus status = ippsZero_64sc(pDst,len);
        checkStatus(status);
    }
    inline void Zero(Ipp64fc* pDst, int len) {
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
    double Sinc(const double x)
    {
        if(std::abs(x) < 1e-6)
            return 1.0;
        return std::sin(M_PI*x) / (M_PI*x);
    }
    void Sinc(const Ipp32f * src, Ipp32f * dst, size_t n) {
        for(size_t i = 0; i < n; i++)
            dst[i] = Sinc(src[i]); 
    }
    void Sinc(const Ipp64f * src, Ipp64f * dst, size_t n) {
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
    void Sum(const Ipp32f * pSrc, int len, Ipp32f* pDst, int hint = 0) {
        IppStatus status = ippsSum_32f(pSrc,len,pDst,(IppHintAlgorithm)hint);
        checkStatus(status);
    }

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

    void AddC(const Ipp32f * src, Ipp32f val, Ipp32f* dst, int len) {
        IppStatus status = ippsAddC_32f(src,val,dst,len);
        checkStatus(status);
    }
    void Add(const Ipp32f * a, Ipp32f *b, Ipp32f* dst, int len) {
        IppStatus status = ippsAdd_32f(a,b,dst,len);
        checkStatus(status);
    }
    void SubC(const Ipp32f * src, Ipp32f val, Ipp32f* dst, int len) {
        IppStatus status = ippsSubC_32f(src,val,dst,len);
        checkStatus(status);
    }
    void Sub(const Ipp32f * a, Ipp32f *b, Ipp32f* dst, int len) {
        IppStatus status = ippsSub_32f(a,b,dst,len);
        checkStatus(status);
    }
    void SubCRev(const Ipp32f * src, Ipp32f val, Ipp32f* dst, int len) {
        IppStatus status = ippsSubCRev_32f(src,val,dst,len);
        checkStatus(status);
    }
    void SubRev(const Ipp32f * a, Ipp32f b, Ipp32f* dst, int len) {
        IppStatus status = ippsSubCRev_32f(a,b,dst,len);
        checkStatus(status);
    }
    void MulC(const Ipp32f * src, Ipp32f val, Ipp32f* dst, int len) {
        IppStatus status = ippsMulC_32f(src,val,dst,len);
        checkStatus(status);
    }
    void Mul(const Ipp32f * a, Ipp32f *b, Ipp32f* dst, int len) {
        IppStatus status = ippsMul_32f(a,b,dst,len);
        checkStatus(status);
    }
    void DivC(const Ipp32f * src, Ipp32f val, Ipp32f* dst, int len) {
        IppStatus status = ippsDivC_32f(src,val,dst,len);
        checkStatus(status);
    }
    void Div(const Ipp32f * a, Ipp32f *b, Ipp32f* dst, int len) {
        IppStatus status = ippsDiv_32f(a,b,dst,len);
        checkStatus(status);
    }
    void DivCRev(const Ipp32f * src, Ipp32f val, Ipp32f* dst, int len) {
        IppStatus status = ippsDivCRev_32f(src,val,dst,len);
        checkStatus(status);
    }    
    void AddProductC(const Ipp32f * array, int len, const Ipp32f val, Ipp32f * dst) {            
        IppStatus status = ippsAddProductC_32f(array,val,dst,len);
        checkStatus(status);        
    }
    void AddProduct(const Ipp32f * array1, int len, const Ipp32f *array2, Ipp32f * dst) {            
        IppStatus status = ippsAddProduct_32f(array1,array2,dst,len);
        checkStatus(status);        
    }
    void Abs(const Ipp32f * src, Ipp32f * dst, int len) {
        IppStatus status = ippsAbs_32f(src,dst,len);
        checkStatus(status);        
    }
    void Sqr(const Ipp32f * src, Ipp32f * dst, int len) {
        IppStatus status = ippsSqr_32f(src,dst,len);
        checkStatus(status);        
    }
    void Sqrt(const Ipp32f * src, Ipp32f * dst, int len) {
        IppStatus status = ippsSqrt_32f(src,dst,len);
        checkStatus(status);        
    }
    void Cubrt(const Ipp32f * src, Ipp32f * dst, int len) {
        IppStatus status = ippsCubrt_32f(src,dst,len);
        checkStatus(status);        
    }
    void Exp(const Ipp32f * src, Ipp32f * dst, int len) {
        IppStatus status = ippsExp_32f(src,dst,len);
        checkStatus(status);        
    }
    void Ln(const Ipp32f * src, Ipp32f * dst, int len) {
        IppStatus status = ippsLn_32f(src,dst,len);
        checkStatus(status);        
    }
    void SumLn(const Ipp32f * src, Ipp32f * dst, int len) {
        IppStatus status = ippsLn_32f(src,dst,len);
        checkStatus(status);        
    }
    void Arctan(const Ipp32f * src, Ipp32f * dst, int len) {
        IppStatus status = ippsArctan_32f(src,dst,len);
        checkStatus(status);        
    }
    Ipp32f Max(const Ipp32f * src,int len) {
        Ipp32f max=0;
        IppStatus status = ippsMax_32f(src,len,&max);
        checkStatus(status);        
        return max;
    }
    Ipp32f Min(const Ipp32f * src,int len) {
        Ipp32f min=0;
        IppStatus status = ippsMin_32f(src,len,&min);
        checkStatus(status);        
        return min;
    }
    int MaxIndex(const Ipp32f * src,Ipp32f * max, int len) {
        int index=-1;
        IppStatus status = ippsMaxIndx_32f(src,len,max,&index);
        checkStatus(status);        
        return index;
    }
    int MinIndex(const Ipp32f * src, int len, Ipp32f * min) {
        int index=-1;
        IppStatus status = ippsMinIndx_32f(src,len,min,&index);
        checkStatus(status);        
        return index;
    }
    Ipp32f MaxAbs(const Ipp32f * src, int len) {
        Ipp32f max=0;
        IppStatus status = ippsMaxAbs_32f(src,len,&max);
        checkStatus(status);        
        return max;
    }
    Ipp32f MinAbs(const Ipp32f * src, int len) {
        Ipp32f min=0;
        IppStatus status = ippsMinAbs_32f(src,len,&min);
        checkStatus(status);        
        return min;
    }
    int MaxAbsIndex(const Ipp32f * src, Ipp32f * max, int len) {
        int index=-1;
        IppStatus status = ippsMaxAbsIndx_32f(src,len,max, &index);
        checkStatus(status);        
        return index;
    }
    int MinAbsIndex(const Ipp32f * src, Ipp32f * min, int len) {
        int index=-1;
        IppStatus status = ippsMinAbsIndx_32f(src,len,min,&index);
        checkStatus(status);        
        return index;
    }
    void MinMax(const Ipp32f * src, int len, Ipp32f * min, Ipp32f * max) {    
        IppStatus status = ippsMinMax_32f(src,len,min,max);
        checkStatus(status);            
    }
    void MinMaxIndex(const Ipp32f * src, int len, Ipp32f *min, int * min_index, Ipp32f *max, int * max_index) {
        IppStatus status = ippsMinMaxIndx_32f(src,len,min,min_index,max,max_index);
        checkStatus(status);            
    }
    void ReplaceNAN(Ipp32f * p, Ipp32f v, int len) {
        IppStatus status = ippsReplaceNAN_32f_I(p,len,v);
        checkStatus(status);            
    }
    Ipp32f Mean(const Ipp32f * src, int len)
    {
        Ipp32f mean;
        IppStatus status = ippsMean_32f(src,len,&mean,ippAlgHintNone);
        checkStatus(status);            
        return mean;
    }
    Ipp32f StdDev(const Ipp32f * src, int len)
    {
        Ipp32f mean;
        IppStatus status = ippsStdDev_32f(src,len,&mean,ippAlgHintNone);
        checkStatus(status);            
        return mean;
    }
    void MeanStdDev(const Ipp32f * src, int len, Ipp32f * mean, Ipp32f * dev)
    {        
        IppStatus status = ippsMeanStdDev_32f(src,len,mean,dev,ippAlgHintNone);
        checkStatus(status);                    
    }
    Ipp32f NormInf(const Ipp32f * src, int len) {
        Ipp32f norm;
        IppStatus status = ippsNorm_Inf_32f(src,len,&norm);
        checkStatus(status);                    
        return norm;
    }
    Ipp32f NormL1(const Ipp32f * src, int len) {
        Ipp32f norm;
        IppStatus status = ippsNorm_L1_32f(src,len,&norm);
        checkStatus(status);                    
        return norm;
    }
    Ipp32f NormL2(const Ipp32f * src, int len) {
        Ipp32f norm;
        IppStatus status = ippsNorm_L2_32f(src,len,&norm);
        checkStatus(status);                    
        return norm;
    }
    Ipp32f NormDiffInf(const Ipp32f * src1,const Ipp32f * src2, int len) {
        Ipp32f norm;
        IppStatus status = ippsNormDiff_Inf_32f(src1,src2,len,&norm);
        checkStatus(status);                    
        return norm;
    }
    Ipp32f NormDiffL1(const Ipp32f * src1, const Ipp32f * src2, int len) {
        Ipp32f norm;
        IppStatus status = ippsNormDiff_L1_32f(src1,src2,len,&norm);
        checkStatus(status);                    
        return norm;
    }
    Ipp32f NormDiffL2(const Ipp32f * src1, const Ipp32f * src2, int len) {
        Ipp32f norm;
        IppStatus status = ippsNormDiff_L2_32f(src1,src2,len,&norm);
        checkStatus(status);                    
        return norm;
    }
    Ipp32f DotProduct(const Ipp32f * src1, const Ipp32f * src2, int len) {
        Ipp32f dp;
        IppStatus status = ippsDotProd_32f(src1,src2,len,&dp);
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
    Ipp32f MSE(const Ipp32f * thetaAbs, const Ipp32f * thetaModel, int len) {
        Ipp32f r = 0;
        for(size_t i = 0; i <  len; i++)
            r += std::pow(thetaAbs[i] - thetaModel[i],2.0f);
        return (1.0f / (float)len) * r;
    }
    Ipp32f RMSE(const Ipp32f * thetaAbs, const Ipp32f * thetaModel, int len) {
        Ipp32f r = MSE(thetaAbs,thetaModel,len);
        return std::sqrt(r);
    }
    Ipp32f MeanSquare(const Ipp32f * x, int len) {
        Ipp32f r = 0;
        for(size_t i = 0; i < len; i++) r += std::pow(x[i],2.0);
        return (1.0f / (float)len) * r;
    }
    Ipp32f AmpToDb(Ipp32f a) {
        return std::pow(10.0,a/20.0);
    }
    Ipp32f DbToAmp(Ipp32f a) {
        return 20.0*std::log10(a);
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
        
    // Short
    // Float
    // Double
    // FloatComplex
    // DoubleComplex
    struct FloatArray
    {
        std::shared_ptr<Ipp32f> ptr;
        Ipp32f * array;
        size_t   len;
        int      r,w;

        FloatArray(size_t n) {
            array = ippsMalloc_32f(len = n);
            ptr = std::shared_ptr<Ipp32f>(array,[](Ipp32f* p) { ippsFree(p); });
            assert(array != NULL);
            Zero(array,n);
            r = 0;
            w = 0;
        }
        FloatArray(const FloatArray &x) {
            ptr = x.ptr;
            array = x.array;
            len = x.len;
        }
        ~FloatArray() {
            
        }

        Ipp32f& operator[] (size_t i) { return array[i]; }

        Ipp32f __getitem__(size_t i) { return array[i]; }
        void __setitem__(size_t i, Ipp32f v) { array[i] = v; }

        FloatArray& operator = (const FloatArray & x) {
            ptr.reset();
            ptr = x.ptr;
            array = x.array;
            len = x.len;
            return *this;
        }
        FloatArray& Copy(Ipp32f * p, int length) {
            if(len != length) {
                ptr.reset();
                array = ippsMalloc_32f(length);
                ptr = std::shared_ptr<Ipp32f>(array,[](Ipp32f* p) { ippsFree(p); });
                assert(array != NULL);
                Zero(array,length);
                r = 0;
                w = 0;
            }
            memcpy(array,p,length);
            len = length;
            return *this;
        }
        void resize(size_t n) {
            Ipp32f * p = ippsMalloc_32f(n);
            Move(array,p,n);
            ippsFree(array);
            array  = p;
            len    = n;
        }
        void fill(Ipp32f value) {
            if(array == NULL) return;
            Set(value,array,len);
        }

        void ring_push(Ipp32f value) {
            array[w++] = value;
            w = w % len;
        }
        Ipp32f ring_pop() {
            Ipp32f v = array[r++];
            r = r % len;
            return v;
        }
        Ipp32f ring_linear_pop() {
            Ipp32f v1 = array[r];
            r = (r+1) % len;
            Ipp32f v2 = array[r];            
            Ipp32f frac = v1 - std::floor(v1);
            return v1 + frac*(v2-v1);
        }
        FloatArray operator + (const Ipp32f value) {
            FloatArray r(*this);
            AddC(array,value,r.array,len);
            return r;
        }
        FloatArray operator + (const FloatArray & b) {
            FloatArray r(*this);
            assert(len == b.len);
            Add(array,b.array,r.array,len);
            return r;
        }
        FloatArray operator - (const Ipp32f value) {
            FloatArray r(*this);
            SubC(array,value,r.array,len);
            return r;
        }
        FloatArray operator - (const FloatArray & b) {
            FloatArray r(*this);
            assert(len == b.len);
            Sub(array,b.array,r.array,len);
            return r;
        }        
        FloatArray operator * (const Ipp32f value) {
            FloatArray r(*this);
            MulC(array,value,r.array,len);
            return r;
        }        
        FloatArray operator * (const FloatArray & b) {
            FloatArray r(*this);
            assert(len == b.len);
            Mul(array,b.array,r.array,len);
            return r;
        }
        FloatArray operator / (const Ipp32f value) {
            FloatArray r(*this);
            DivC(array,value,r.array,len);
            return r;
        }        
        FloatArray operator / (const FloatArray & b) {
            FloatArray r(*this);
            assert(len == b.len);
            Div(array,b.array,r.array,len);
            return r;
        }
        
        Ipp32f Sum() {
            Ipp32f r = 0;
            Casino::Sum(array,len,&r);
            return r;
        }
        
        void print() {
            std::cout << "Array[" << len << "]=";
            for(size_t i = 0; i < len-1; i++) std::cout << array[i] << ",";
            std::cout << array[len-1] << std::endl;
        }
        FloatArray& copy(const FloatArray & a) {
            ptr.reset();
            array = ippsMalloc_32f(a.len);
            memcpy(array,a.array,a.len*sizeof(ipp32f));
            return *this; 
       }

        /*
        void RandUniform(Ipp32f low, Ipp32f high, unsigned int seed = -1)
        {
            if(seed == -1) seed = time(NULL);            
            IppsRandUniState_32f state;
            IppStatus status = ippsRandUniformInit_32f(&state,low,high,seed)
            checkStatus(status);
            status = ippsRandUniform_32f(array,len,&state);
        }
        void RandGaussian(Ipp32f mean, Ipp32f stdDev, unsigned int seed = -1)
        {
            if(seed == -1) seed = time(NULL);
            IppsRandUniState_32f state;
            IppStatus status = ippsRandGaussianInit_32f(&state,mean,stdDev,seed)
            checkStatus(status);
            status = ippsRandGaussian_32f(array,len,&state);
        }
        */        
    };

    sample_vector<float> ConvertToSampleVector(const FloatArray &a) {
        sample_vector<float> r(a.len);
        memcpy(r.data(),a.array,a.len*sizeof(float));
        return r;
    }
    /*
    Vector::Vector4f ConvertToVector4f(const FloatArray &a) {
        Vector::Vector4f r(a.len);
        memcpy(r.data(),a.array,a.len*sizeof(float));
        return r;
    }
    Vector::Vector8f ConvertToVector8f(const FloatArray &a) {
        Vector::Vector8f r(a.len);
        memcpy(r.data(),a.array,a.len*sizeof(float));
        return r;
    }
    */
    sample_vector<double> ConvertToSampleVectorDouble(const FloatArray &a) {
        sample_vector<double> r(a.len);
        for(size_t i = 0; i < a.len; i++) r[i] = a.array[i];
        return r;
    }
    /*
    Vector::Vector4f ConvertToVector2d(const FloatArray &a) {
        Vector::Vector4f r(a.len);
        for(size_t i = 0; i < a.len; i++) r[i] = a.array[i];
        return r;
    }
    Vector::Vector8f ConvertToVector4d(const FloatArray &a) {
        Vector::Vector8f r(a.len);
        for(size_t i = 0; i < a.len; i++) r[i] = a.array[i];
        return r;
    }
    */
    struct AutoCorrFloat
    {
        Ipp8u * buffer;
        int bufferLen;
        int      type;
        int srcLen,dstLen;

        AutoCorrFloat(size_t srcLen, size_t dstLen, int algorithm = (int)Algorithm::ALG_FFT) {
            ippsAutoCorrNormGetBufferSize(srcLen,dstLen,ipp32f,algorithm,&bufferLen);
            if(bufferLen > 0) buffer = ippsMalloc_8u(bufferLen);
            assert(buffer != NULL);
            type = algorithm;
            this->srcLen = srcLen;
            this->dstLen = dstLen;
        }
        ~AutoCorrFloat() {
            if(buffer) ippsFree(buffer);
        }
        void Process(Ipp32f * src, Ipp32f * dst) {
            IppStatus status = ippsAutoCorrNorm_32f(src,srcLen,dst,dstLen,type,buffer);
            checkStatus(status);
        }
    };

    struct AutoCorrFloatComplex
    {
        Ipp8u * buffer;
        int bufferLen;
        int      type;
        int srcLen,dstLen;

        AutoCorrFloatComplex(size_t srcLen, size_t dstLen, int algorithm = (int)Algorithm::ALG_FFT) {
            ippsAutoCorrNormGetBufferSize(srcLen,dstLen,ipp32fc,algorithm,&bufferLen);
            if(bufferLen > 0) buffer = ippsMalloc_8u(bufferLen);
            assert(buffer != NULL);
            type = algorithm;
            this->srcLen = srcLen;
            this->dstLen = dstLen;
        }
        ~AutoCorrFloatComplex() {
            if(buffer) ippsFree(buffer);
        }
        void Process(Ipp32fc * src, Ipp32fc * dst) {
            IppStatus status = ippsAutoCorrNorm_32fc(src,srcLen,dst,dstLen,type,buffer);
            checkStatus(status);
        }
    };

    
    struct AutoCorrDouble
    {
        Ipp8u * buffer;
        int bufferLen;
        int      type;
        int srcLen,dstLen;

        AutoCorrDouble(size_t srcLen, size_t dstLen, int algorithm = (int)Algorithm::ALG_FFT) {
            ippsAutoCorrNormGetBufferSize(srcLen,dstLen,ipp64f,algorithm,&bufferLen);
            if(bufferLen > 0) buffer = ippsMalloc_8u(bufferLen);
            assert(buffer != NULL);
            type = algorithm;
            this->srcLen = srcLen;
            this->dstLen = dstLen;
        }
        ~AutoCorrDouble() {
            if(buffer) ippsFree(buffer);
        }
        void Process(Ipp64f * src, Ipp64f * dst) {
            IppStatus status = ippsAutoCorrNorm_64f(src,srcLen,dst,dstLen,type,buffer);
            checkStatus(status);
        }
    };

    struct AutoCorrDoubleComplex
    {
        Ipp8u * buffer;
        int bufferLen;
        int      type;
        int srcLen,dstLen;

        AutoCorrDoubleComplex(size_t srcLen, size_t dstLen, int algorithm = (int)Algorithm::ALG_FFT) {
            ippsAutoCorrNormGetBufferSize(srcLen,dstLen,ipp64fc,algorithm,&bufferLen);
            if(bufferLen > 0) buffer = ippsMalloc_8u(bufferLen);
            assert(buffer != NULL);
            type = algorithm;
            this->srcLen = srcLen;
            this->dstLen = dstLen;
        }
        ~AutoCorrDoubleComplex() {
            if(buffer) ippsFree(buffer);
        }
        void Process(Ipp64fc * src, Ipp64fc * dst) {
            IppStatus status = ippsAutoCorrNorm_64fc(src,srcLen,dst,dstLen,type,buffer);
            checkStatus(status);
        }
    };

    void acorr(size_t srcLen, Ipp32f * src, size_t dstLen, Ipp32f * dst, int algorithm = (int)Algorithm::ALG_FFT)
    {
        AutoCorrFloat a(srcLen,dstLen,algorithm);
        a.Process(src,dst);
    }
    
    void acorr(size_t srcLen, Ipp32fc * src, size_t dstLen, Ipp32fc * dst, int algorithm = (int)Algorithm::ALG_FFT)
    {
        AutoCorrFloatComplex a(srcLen,dstLen,algorithm);
        a.Process(src,dst);
    }
    
    void acorr(size_t srcLen, Ipp64f * src, size_t dstLen, Ipp64f * dst, int algorithm = (int)Algorithm::ALG_FFT)
    {
        AutoCorrDouble a(srcLen,dstLen,algorithm);
        a.Process(src,dst);
    }
    
    void acorr(size_t srcLen, Ipp64fc * src, size_t dstLen, Ipp64fc * dst, int algorithm = (int)Algorithm::ALG_FFT)
    {
        AutoCorrDoubleComplex a(srcLen,dstLen,algorithm);
        a.Process(src,dst);
    }
    
    
    struct CrossCorrFloat
    {
        Ipp8u * buffer;
        int bufferLen;
        int      type;
        int src1Len,src2Len,dstLen,lowLag;

        CrossCorrFloat(size_t src1Len, size_t src2Len, size_t dstLen, int lowLag, int algorithm = (int)Algorithm::ALG_FFT) {
            ippsCrossCorrNormGetBufferSize(src1Len,src2Len,dstLen,lowLag,ipp32f,algorithm,&bufferLen);
            if(bufferLen > 0) buffer = ippsMalloc_8u(bufferLen);
            assert(buffer != NULL);
            type = algorithm;
            this->lowLag = lowLag;
            this->src1Len = src1Len;
            this->src2Len = src2Len;
            this->dstLen = dstLen;
        }
        ~CrossCorrFloat() {
            if(buffer) ippsFree(buffer);
        }
        void Process(Ipp32f * src1, Ipp32f * src2, Ipp32f * dst) {
            IppStatus status = ippsCrossCorrNorm_32f(src1,src1Len,src2,src2Len,dst,dstLen,lowLag,type,buffer);
            checkStatus(status);
        }
    };

    

    struct CrossCorrFloatComplex
    {
        Ipp8u * buffer;
        int bufferLen;
        int      type;
        int src1Len,src2Len,dstLen,lowLag;

        CrossCorrFloatComplex(size_t src1Len, size_t src2Len, size_t dstLen, int lowLag, int algorithm = (int)Algorithm::ALG_FFT) {
            ippsCrossCorrNormGetBufferSize(src1Len,src2Len,dstLen,lowLag,ipp32fc,algorithm,&bufferLen);
            if(bufferLen > 0) buffer = ippsMalloc_8u(bufferLen);
            assert(buffer != NULL);
            type = algorithm;
            this->src1Len = src1Len;
            this->src2Len = src2Len;
            this->dstLen = dstLen;
            this->lowLag = lowLag;
        }
        ~CrossCorrFloatComplex() {
            if(buffer) ippsFree(buffer);
        }
        void Process(Ipp32fc * src1, Ipp32fc * src2, Ipp32fc * dst) {
            IppStatus status = ippsCrossCorrNorm_32fc(src1,src1Len,src2,src2Len,dst,dstLen,lowLag,type,buffer);
            checkStatus(status);
        }
    };

    void xcorr(size_t srcLen, Ipp32f * src1, size_t srcLen2, Ipp32f* src2, size_t dstLen, Ipp32f * dst, int lowLag,int algorithm = (int)Algorithm::ALG_FFT)
    {
        CrossCorrFloat c(srcLen,srcLen2,dstLen,lowLag,algorithm);
        c.Process(src1,src2,dst);
    }        
    void xcorr(size_t srcLen, Ipp32fc * src1, size_t srcLen2, Ipp32fc* src2, size_t dstLen, Ipp32fc * dst, int lowLag,int algorithm = (int)Algorithm::ALG_FFT)
    {
        CrossCorrFloatComplex c(srcLen,srcLen2,dstLen,lowLag,algorithm);
        c.Process(src1,src2,dst);
    }        
    
    struct CrossCorrDouble
    {
        Ipp8u * buffer;
        int bufferLen;
        int      type;
        int src1Len,src2Len,dstLen,lowLag;

        CrossCorrDouble(size_t src1Len,  size_t src2Len, size_t dstLen, int lowLag, int algorithm = (int)Algorithm::ALG_FFT) {
            ippsCrossCorrNormGetBufferSize(src1Len,src2Len,dstLen,lowLag,ipp64f,algorithm,&bufferLen);
            if(bufferLen > 0) buffer = ippsMalloc_8u(bufferLen);
            assert(buffer != NULL);
            type = algorithm;
            this->src1Len = src1Len;
            this->src2Len = src2Len;
            this->dstLen = dstLen;
            this->lowLag = lowLag;
        }
        ~CrossCorrDouble() {
            if(buffer) ippsFree(buffer);
        }
        void Process(Ipp64f * src1, Ipp64f * src2, Ipp64f * dst) {
            IppStatus status = ippsCrossCorrNorm_64f(src1,src1Len,src2,src2Len,dst,dstLen,lowLag,type,buffer);
            checkStatus(status);
        }
    };

    struct CrossCorrDoubleComplex
    {
        Ipp8u * buffer;
        int bufferLen;
        int      type;
        int src1Len,src2Len,dstLen,lowLag;

        CrossCorrDoubleComplex(size_t src1Len, size_t src2Len, size_t dstLen, int lowLag, int algorithm = (int)Algorithm::ALG_FFT) {
            ippsCrossCorrNormGetBufferSize(src1Len,src2Len,dstLen,lowLag,ipp64fc,algorithm,&bufferLen);
            if(bufferLen > 0) buffer = ippsMalloc_8u(bufferLen);
            assert(buffer != NULL);
            type = algorithm;
            this->src1Len = src1Len;
            this->src2Len = src2Len;
            this->dstLen = dstLen;
            this->lowLag = lowLag;
        }
        ~CrossCorrDoubleComplex() {
            if(buffer) ippsFree(buffer);
        }
        void Process(Ipp64fc * src1, Ipp64fc * src2, Ipp64fc * dst) {
            IppStatus status = ippsCrossCorrNorm_64fc(src1,src1Len,src2,src2Len,dst,dstLen,lowLag,type,buffer);
            checkStatus(status);
        }
    };

    void xcorr(size_t srcLen, Ipp64f * src1, size_t srcLen2, Ipp64f* src2, size_t dstLen, Ipp64f * dst, int lowLag,int algorithm = (int)Algorithm::ALG_FFT)
    {
        CrossCorrDouble c(srcLen,srcLen2,dstLen,lowLag,algorithm);
        c.Process(src1,src2,dst);
    }        
    void xcorr(size_t srcLen, Ipp64fc * src1, size_t srcLen2, Ipp64fc* src2, size_t dstLen, Ipp64fc * dst, int lowLag,int algorithm = (int)Algorithm::ALG_FFT)
    {
        CrossCorrDoubleComplex c(srcLen,srcLen2,dstLen,lowLag,algorithm);
        c.Process(src1,src2,dst);
    }        

    struct ConvolveFloat
    {
        Ipp8u * buffer;
        int bufferLen;
        int      type;
        int src1Len,dstLen,src2Len;

        ConvolveFloat(size_t src1Len,size_t src2Len, int algorithm = (int)Algorithm::ALG_AUTO) {
            ippsConvolveGetBufferSize(src1Len,src2Len,ipp32f,algorithm,&bufferLen);
            if(bufferLen > 0) buffer = ippsMalloc_8u(bufferLen);
            assert(buffer != NULL);
            type = algorithm;
            this->src1Len = src1Len;
            this->src2Len = src2Len;        
        }
        ~ConvolveFloat() {
            if(buffer) ippsFree(buffer);
        }
        void Process(Ipp32f * src1, Ipp32f * src2, Ipp32f * dst) {
            IppStatus status = ippsConvolve_32f(src1,src1Len,src2,src2Len,dst,type,buffer);
            checkStatus(status);
        }
    };
    
    struct ConvolveDouble
    {
        Ipp8u * buffer;
        int bufferLen;
        int      type;
        int src1Len,dstLen,src2Len;

        ConvolveDouble(size_t src1Len,  size_t src2Len, int algorithm = (int)Algorithm::ALG_FFT) {
            ippsConvolveGetBufferSize(src1Len,src2Len,ipp64f,algorithm,&bufferLen);
            if(bufferLen > 0) buffer = ippsMalloc_8u(bufferLen);
            assert(buffer != NULL);
            type = algorithm;
            this->src1Len = src1Len;
            this->src2Len = src2Len;            
        }
        ~ConvolveDouble() {
            if(buffer) ippsFree(buffer);
        }
        void Process(Ipp64f * src1, Ipp64f * src2, Ipp64f * dst) {
            IppStatus status = ippsConvolve_64f(src1,src1Len,src2,src2Len,dst,type,buffer);
            checkStatus(status);
        }
    };
    void conv(size_t src1Len, Ipp32f * src1, size_t src2Len, Ipp32f * src2, Ipp32f * dst, int algorithm = (int)Algorithm::ALG_AUTO)
    {
        ConvolveFloat c(src1Len,src2Len,algorithm);
        c.Process(src1,src2,dst);
    }
    void conv(size_t src1Len, Ipp64f * src1, size_t src2Len, Ipp64f * src2, Ipp64f * dst, int algorithm = (int)Algorithm::ALG_AUTO)
    {
        ConvolveDouble c(src1Len,src2Len,algorithm);
        c.Process(src1,src2,dst);
    }

    struct ConvolutionFilterFloat
    {
        ConvolveFloat *filter;
        std::shared_ptr<Ipp32f> h;
        size_t            len,block;
        std::vector<Ipp32f> ola,temp;

        ConvolutionFilterFloat(size_t n, Ipp32f * impulse, size_t block_size)
        {
            filter = new ConvolveFloat(block_size,block_size);
            block = block_size;
            Ipp32f *hi = ippsMalloc_32f(len);
            h = std::shared_ptr<Ipp32f>(hi,[](Ipp32f * p) { ippsFree(p); }); 
            len = n;
            assert(filter != NULL);
            assert(hi != nullptr);
            ola.resize(block_size);
            temp.resize(block_size + n -1);
        }
        ~ConvolutionFilterFloat() {

        }
        void ProcessBlock(Ipp32f * signal, Ipp32f * dest) {
            filter->Process(h.get(),signal,dest);
            for(size_t i = 0; i < block; i++) {
                dest[i] = temp[i] + ola[i];
                ola[i]  = temp[i+block];
            }
        }
    };

    struct ConvolutionFilterDouble
    {
        ConvolveDouble *filter;
        std::shared_ptr<Ipp64f> h;
        size_t            len,block;
        std::vector<Ipp64f> ola,temp;

        ConvolutionFilterDouble(size_t n, Ipp64f * impulse, size_t block_size)
        {
            filter = new ConvolveDouble(block_size,block_size);
            block = block_size;
            Ipp64f *hi = ippsMalloc_64f(len);
            h = std::shared_ptr<Ipp64f>(hi,[](Ipp64f * p) { ippsFree(p); }); 
            len = n;
            assert(filter != NULL);
            assert(hi != nullptr);
            ola.resize(block_size);
            temp.resize(block_size + n -1);
        }
        ~ConvolutionFilterDouble() {

        }
        void ProcessBlock(Ipp64f * signal, Ipp64f * dest) {
            filter->Process(h.get(),signal,dest);
            for(size_t i = 0; i < block; i++) {
                dest[i] = temp[i] + ola[i];
                ola[i]  = temp[i+block];
            }
        }
    };
    
    struct CFFTFC
    {
        Ipp8u * pBuffer;
        Ipp8u * pSpec;
        Ipp8u * pSpecBuffer;
        IppsFFTSpec_C_32fc * fft;
                
        CFFTFC(size_t n) {
            int size,specbuffer,spec;            
            int flag = IPP_FFT_DIV_FWD_BY_N;
            ippsFFTGetSize_C_32fc(n,flag,ippAlgHintNone,&spec,&specbuffer,&size);            
            pSpec = ippsMalloc_8u(spec);
            pSpecBuffer = specbuffer > 0? ippsMalloc_8u(specbuffer) : NULL;
            pBuffer = ippsMalloc_8u(size);
            int order = std::log2(n);
            IppStatus status = ippsFFTInit_C_32fc(&fft,order, flag,ippAlgHintNone,pSpec,pSpecBuffer);            
            checkStatus(status);
        }
        ~CFFTFC() {
            if(pSpec) ippsFree(pSpec);
            if(pSpecBuffer) ippsFree(pSpecBuffer);
            if(pBuffer) ippsFree(pBuffer);
        }
        void Forward(const Ipp32fc* pSrc, Ipp32fc * pDst)
        {                               
            IppStatus status = ippsFFTFwd_CToC_32fc(pSrc,pDst, fft, pBuffer);                        
            checkStatus(status);
        }
        void Inverse(const Ipp32fc* pSrc, Ipp32fc * pDst)
        {                               
            IppStatus status = ippsFFTInv_CToC_32fc(pSrc,pDst, fft, pBuffer);
            checkStatus(status);
        }
    };
    
    void fft(size_t n, const Ipp32fc * pSrc, Ipp32fc * pDst)
    {
        CFFTFC f(n);
        f.Forward(pSrc,pDst);
    }
    void ifft(size_t n, const Ipp32fc * pSrc, Ipp32fc * pDst)
    {
        CFFTFC f(n);
        f.Inverse(pSrc,pDst);
    }

    struct CFFTF
    {
        Ipp8u * pBuffer;
        Ipp8u * pSpec;
        Ipp8u * pSpecBuffer;
        IppsFFTSpec_C_32f * fft;
                        
        CFFTF(size_t n) {
            int size,specbuffer,spec;            
            int flag = IPP_FFT_DIV_FWD_BY_N;
            ippsFFTGetSize_C_32f(n,flag,ippAlgHintNone,&spec,&specbuffer,&size);            
            pSpec = ippsMalloc_8u(spec);
            pSpecBuffer = specbuffer > 0? ippsMalloc_8u(specbuffer) : NULL;
            pBuffer = ippsMalloc_8u(size);
            int order = std::log2(n);
            IppStatus status = ippsFFTInit_C_32f(&fft,order, flag,ippAlgHintNone,pSpec,pSpecBuffer);            
            checkStatus(status);
        }
        ~CFFTF() {
            if(pSpec) ippsFree(pSpec);
            if(pSpecBuffer) ippsFree(pSpecBuffer);
            if(pBuffer) ippsFree(pBuffer);
        }
        void Forward(const Ipp32f* pSrcRe, const Ipp32f* pSrcIm, Ipp32f * pDstRe, Ipp32f * pDstIm)
        {                   
            ippsFFTFwd_CToC_32f(pSrcRe,pSrcIm,pDstRe,pDstIm,fft,pBuffer);
        }
        void Inverse(const Ipp32f* pSrcRe, const Ipp32f* pSrcIm, Ipp32f * pDstRe, Ipp32f * pDstIm)
        {                   
            ippsFFTInv_CToC_32f(pSrcRe,pSrcIm,pDstRe,pDstIm,fft,pBuffer);
        }
    };

    void fft(size_t n, const Ipp32f * pSrcRe, const Ipp32f* pSrcIm, Ipp32f * pDstRe, Ipp32f * pDstIm)
    {
        CFFTF f(n);
        f.Forward(pSrcRe,pSrcIm,pDstRe,pDstIm);
    }
    void ifft(size_t n, const Ipp32f * pSrcRe, const Ipp32f * pSrcIm, Ipp32f * pDstRe, Ipp32f * pDstIm)    
    {
        CFFTF f(n);
        f.Inverse(pSrcRe, pSrcIm, pDstRe, pDstIm);
    }

    struct RFFTF
    {
        Ipp8u * pBuffer;
        Ipp8u * pSpec;
        Ipp8u * pSpecBuffer;
        IppsFFTSpec_R_32f * fft;
        
                
        RFFTF(size_t n) {
            int size,specbuffer,spec;            
            int flag = IPP_FFT_DIV_FWD_BY_N;
            ippsFFTGetSize_R_32f(n,flag,ippAlgHintNone,&spec,&specbuffer,&size);            
            pSpec = ippsMalloc_8u(spec);
            pSpecBuffer = specbuffer > 0? ippsMalloc_8u(specbuffer) : NULL;
            pBuffer = ippsMalloc_8u(size);
            int order = std::log2(n);
            IppStatus status = ippsFFTInit_R_32f(&fft,order, flag,ippAlgHintNone,pSpec,pSpecBuffer);            
            checkStatus(status); 
        }
        ~RFFTF() {
            if(pSpec) ippsFree(pSpec);
            if(pSpecBuffer) ippsFree(pSpecBuffer);
            if(pBuffer) ippsFree(pBuffer);            
        }
        void Forward(const Ipp32f* pSrc, Ipp32f * pDst)
        {                   
            IppStatus status;
            status = ippsFFTFwd_RToPack_32f(pSrc,pDst,fft,pBuffer);
            checkStatus(status);            
        }
        void Inverse(const Ipp32f* pSrc, Ipp32f * pDst)
        {                   
            IppStatus status = ippsFFTInv_PackToR_32f(pSrc,pDst,fft,pBuffer);
            checkStatus(status);            
        }
        void Unpack(int len, Ipp32f * pSrc, Ipp32fc * pDst) {
            ippsConjPack_32fc(pSrc,pDst,len);
        }
    };

    void fft(size_t n, const Ipp32f * pSrc, Ipp32f * pPacked, Ipp32fc * pDst)
    {
        RFFTF f(n);
        f.Forward(pSrc,pPacked);
        f.Unpack(n,pPacked,pDst);
    }
    void ifft(size_t n, const Ipp32f * pSrcPacked, Ipp32f * pDst)
    {
        RFFTF f(n);
        f.Inverse(pSrcPacked,pDst);
    }
    
    
    
    struct CDFTF
    {
        Ipp8u * pBuffer;        
        Ipp8u * pSpecBuffer;
        IppsDFTSpec_C_32fc * fft;
                        
        CDFTF(size_t n) {
            int size,specbuffer,spec;            
            int flag = IPP_FFT_DIV_FWD_BY_N;
            ippsDFTGetSize_C_32f(n,flag,ippAlgHintNone,&spec,&specbuffer,&size);            
            fft = (IppsDFTSpec_C_32fc*)ippsMalloc_8u(spec);
            pSpecBuffer = specbuffer > 0? ippsMalloc_8u(specbuffer) : NULL;
            pBuffer = ippsMalloc_8u(size);
            int order = n;
            IppStatus status = ippsDFTInit_C_32fc(order,flag,ippAlgHintNone,fft,pSpecBuffer);            
            checkStatus(status);
        }
        ~CDFTF() {
            if(fft) ippsFree(fft);
            if(pSpecBuffer) ippsFree(pSpecBuffer);
            if(pBuffer) ippsFree(pBuffer);
        }
        void Forward(const Ipp32fc* pSrc, Ipp32fc * pDst)
        {                   
            ippsDFTFwd_CToC_32fc(pSrc,pDst,fft,pBuffer);
        }
        void Inverse(const Ipp32fc* pSrc, Ipp32fc * pDst)
        {                   
            ippsDFTInv_CToC_32fc(pSrc,pDst,fft,pBuffer);
        }
    };

    void dft(size_t n, const Ipp32fc * pSrc, Ipp32fc * pDst)
    {
        CDFTF d(n);
        d.Forward(pSrc,pDst);
    }
    void idft(size_t n, const Ipp32fc * pSrc, Ipp32fc * pDst)
    {
        CDFTF d(n);
        d.Inverse(pSrc,pDst);
    }

    struct RDFTF
    {
        Ipp8u * pBuffer;        
        Ipp8u * pSpecBuffer;
        IppsDFTSpec_R_32f * fft;
        
                
        RDFTF(size_t n) {
            int size,specbuffer,spec;            
            int flag = IPP_FFT_DIV_FWD_BY_N;
            ippsDFTGetSize_R_32f(n,flag,ippAlgHintNone,&spec,&specbuffer,&size);            
            fft = (IppsDFTSpec_R_32f*)ippsMalloc_8u(spec);
            pSpecBuffer = specbuffer > 0? ippsMalloc_8u(specbuffer) : NULL;
            pBuffer = ippsMalloc_8u(size);
            int order = n;
            IppStatus status = ippsDFTInit_R_32f(order,flag,ippAlgHintNone,fft,pSpecBuffer);            
            checkStatus(status); 
        }
        ~RDFTF() {
            if(fft) ippsFree(fft);
            if(pSpecBuffer) ippsFree(pSpecBuffer);
            if(pBuffer) ippsFree(pBuffer);            
        }
        void Forward(const Ipp32f* pSrc, Ipp32f * pDst)
        {                   
            IppStatus status;
            status = ippsDFTFwd_RToPack_32f(pSrc,pDst,fft,pBuffer);
            checkStatus(status);            
        }
        void Inverse(const Ipp32f* pSrc, Ipp32f * pDst)
        {                   
            IppStatus status = ippsDFTInv_PackToR_32f(pSrc,pDst,fft,pBuffer);
            checkStatus(status);            
        }
        void Unpack(int len, Ipp32f * pSrc, Ipp32fc * pDst) {
            ippsConjPack_32fc(pSrc,pDst,len);
        }
    };

    void dft(size_t n, const Ipp32f * pSrc, Ipp32f * pDst, Ipp32fc * unpacked)
    {
        RDFTF d(n);
        d.Forward(pSrc,pDst);
        d.Unpack(n,pDst,unpacked);
    }
    void idft(size_t n, const Ipp32f * pSrc, Ipp32f * pDst)
    {
        RDFTF d(n);
        d.Inverse(pSrc,pDst);
    }
    
    struct DCTF
    {
        Ipp8u * pBuffer[2];        
        Ipp8u * pSpec[2];
        Ipp8u * pSpecBuffer[2];
        IppsDCTFwdSpec_32f * forward;
        IppsDCTInvSpec_32f *inverse;
                        
        DCTF(size_t n) {
            int size,specbuffer,spec;            
            int flag = IPP_FFT_DIV_FWD_BY_N;

            ippsDCTFwdGetSize_32f(n,ippAlgHintNone,&spec,&specbuffer,&size);            
            pSpec[0]       = ippsMalloc_8u(spec);
            pSpecBuffer[0] = specbuffer > 0? ippsMalloc_8u(specbuffer) : NULL;
            pBuffer[0]     = ippsMalloc_8u(size);

            int order = n;
            IppStatus status = ippsDCTFwdInit_32f(&forward,n,ippAlgHintNone,pSpec[0],pSpecBuffer[0]);            
            checkStatus(status);

            ippsDCTInvGetSize_32f(n,ippAlgHintNone,&spec,&specbuffer,&size);            
            pSpec[1]       = ippsMalloc_8u(spec);
            pSpecBuffer[1] = specbuffer > 0? ippsMalloc_8u(specbuffer) : NULL;
            pBuffer[1]     = ippsMalloc_8u(size);

            status = ippsDCTInvInit_32f(&inverse,n,ippAlgHintNone,pSpec[1],pSpecBuffer[1]);            
            checkStatus(status);
        }
        ~DCTF() {
            if(forward) ippsFree(forward);
            if(inverse) ippsFree(inverse);
            if(pSpecBuffer[0]) ippsFree(pSpecBuffer[0]);
            if(pBuffer[0]) ippsFree(pBuffer[0]);
            if(pSpecBuffer[1]) ippsFree(pSpecBuffer[1]);
            if(pBuffer[1]) ippsFree(pBuffer[1]);
        }
        void Forward(const Ipp32f* pSrc, Ipp32f * pDst)
        {                   
            IppStatus status = ippsDCTFwd_32f(pSrc,pDst,forward,pSpecBuffer[0]);
            checkStatus(status);
        }
        void Inverse(const Ipp32f* pSrc, Ipp32f * pDst)
        {                   
            IppStatus status = ippsDCTInv_32f(pSrc,pDst,inverse,pSpecBuffer[1]);
            checkStatus(status);
        }
    };


    void dct(size_t n, const Ipp32f * pSrc, Ipp32f * pDst)
    {
        DCTF d(n);
        d.Forward(pSrc,pDst);
    }
    void idct(size_t n, const Ipp32f * pSrc, Ipp32f * pDst)
    {
        DCTF d(n);
        d.Inverse(pSrc,pDst);
    }

    void Goertzal(int len, const Ipp32fc * src, Ipp32fc * dst, Ipp32f freq) {
        IppStatus status = ippsGoertz_32fc(src,len,dst,freq);
    }
    void Goertzal(int len, const Ipp32f * src, Ipp32fc * dst, Ipp32f freq) {
        IppStatus status = ippsGoertz_32f(src,len,dst,freq);
    }
    void Goertzal(int len, const Ipp64fc * src, Ipp64fc * dst, Ipp64f freq) {
        IppStatus status = ippsGoertz_64fc(src,len,dst,freq);
    }
    void Goertzal(int len, const Ipp64f * src, Ipp64fc * dst, Ipp64f freq) {
        IppStatus status = ippsGoertz_64f(src,len,dst,freq);
    }

    struct HilbertTransformF
    {
        Ipp8u * pSpec;
        Ipp8u * pBuffer;
        size_t  size;
        HilbertTransformF(size_t n) {
            int specSize,bufferSize;
            size = n;
            IppStatus status = ippsHilbertGetSize_32f32fc(n,ippAlgHintNone,&specSize,&bufferSize);
            checkStatus(status);
            pSpec = ippsMalloc_8u(specSize);
            pBuffer= ippsMalloc_8u(bufferSize);
            status = ippsHilbertInit_32f32fc(n,ippAlgHintNone,(IppsHilbertSpec*)pSpec,pBuffer);
            checkStatus(status);
        }
        ~HilbertTransformF() {
            if(pSpec) ippsFree(pSpec);
            if(pBuffer) ippsFree(pBuffer);
        }
        void Execute(const Ipp32f * pSrc, Ipp32fc * pDst) {
            IppStatus status = ippsHilbert_32f32fc(pSrc,pDst,(IppsHilbertSpec*)pSpec,pBuffer);
            checkStatus(status);            
        }
    };
    

    void UpSample(const Ipp16s * pSrc, int srcLen, Ipp16s* pDst, int *dstLen, int fact, int * pPhase) {
        IppStatus status = ippsSampleUp_16s(pSrc,srcLen,pDst,dstLen,fact,pPhase);
        checkStatus(status);
    }
    void UpSample(const Ipp32f * pSrc, int srcLen, Ipp32f* pDst, int *dstLen, int fact, int * pPhase) {
        IppStatus status = ippsSampleUp_32f(pSrc,srcLen,pDst,dstLen,fact,pPhase);
        checkStatus(status);
    }
    void UpSample(const Ipp64f * pSrc, int srcLen, Ipp64f* pDst, int *dstLen, int fact, int * pPhase) {
        IppStatus status = ippsSampleUp_64f(pSrc,srcLen,pDst,dstLen,fact,pPhase);
        checkStatus(status);
    }
    void UpSample(const Ipp32fc * pSrc, int srcLen, Ipp32fc* pDst, int *dstLen, int fact, int * pPhase) {
        IppStatus status = ippsSampleUp_32fc(pSrc,srcLen,pDst,dstLen,fact,pPhase);
        checkStatus(status);
    }
    void UpSample(const Ipp64fc * pSrc, int srcLen, Ipp64fc* pDst, int *dstLen, int fact, int * pPhase) {
        IppStatus status = ippsSampleUp_64fc(pSrc,srcLen,pDst,dstLen,fact,pPhase);
        checkStatus(status);
    }
    void DownSample(const Ipp16s * pSrc, int srcLen, Ipp16s* pDst, int *dstLen, int fact, int * pPhase) {
        IppStatus status = ippsSampleDown_16s(pSrc,srcLen,pDst,dstLen,fact,pPhase);
        checkStatus(status);
    }
    void DownSample(const Ipp32f * pSrc, int srcLen, Ipp32f* pDst, int *dstLen, int fact, int * pPhase) {
        IppStatus status = ippsSampleDown_32f(pSrc,srcLen,pDst,dstLen,fact,pPhase);
        checkStatus(status);
    }
    void DownSample(const Ipp64f * pSrc, int srcLen, Ipp64f* pDst, int *dstLen, int fact, int * pPhase) {
        IppStatus status = ippsSampleDown_64f(pSrc,srcLen,pDst,dstLen,fact,pPhase);
        checkStatus(status);
    }
    void DownSample(const Ipp32fc * pSrc, int srcLen, Ipp32fc* pDst, int *dstLen, int fact, int * pPhase) {
        IppStatus status = ippsSampleDown_32fc(pSrc,srcLen,pDst,dstLen,fact,pPhase);
        checkStatus(status);
    }
    void DownSample(const Ipp64fc * pSrc, int srcLen, Ipp64fc* pDst, int *dstLen, int fact, int * pPhase) {
        IppStatus status = ippsSampleDown_64fc(pSrc,srcLen,pDst,dstLen,fact,pPhase);
        checkStatus(status);
    }
    
    void DoubleToFloat(size_t len, double * src, float * dst) {
        for(size_t i = 0; i < len; i++) dst[i] = (float)src[i];
    }
    void FloatToDouble(size_t len, float * src, double * dst) {
        for(size_t i = 0; i < len; i++) dst[i] = (double)src[i];
    }

    // warning - do not know if this works yet
    struct FIRMR32F
    {
        Ipp8u * pSpec;
        Ipp8u * pBuffer;
        size_t  size;
        Ipp32f *pTaps;

        FIRMR32F(size_t n, int up, int down, int upPhase=0,int downPhase=0) {
            int specSize,bufferSize;
            IppStatus status = ippsFIRMRGetSize(n,up,down,ipp32f,&specSize,&bufferSize);
            checkStatus(status);
            pSpec = ippsMalloc_8u(specSize);
            pBuffer= ippsMalloc_8u(bufferSize);
            //pDlySrc= ippsMalloc_32f(((n + up -1) / up);
            //pDlyDst= ippsMalloc_32f(((n + down -1) / down);
            pTaps   = ippsMalloc_32f(n);
            status = ippsFIRMRInit_32f(pTaps,n,up,upPhase,down,downPhase,(IppsFIRSpec_32f*)pSpec);
            checkStatus(status);
        }
        ~FIRMR32F() {
            if(pSpec) ippsFree(pSpec);
            if(pBuffer) ippsFree(pBuffer);
            //if(pDlySrc) ippFree(pDlySrc);
            //if(pDlyDst) ippFree(pDlyDst);
            if(pTaps) ippFree(pTaps);
        }
        void Execute(const Ipp32f * pSrc, Ipp32f* pDst, int numIters) {
            IppStatus status = ippsFIRMR_32f(pSrc,pDst,numIters,(IppsFIRSpec_32f*)pSpec,NULL,NULL,pBuffer);
            checkStatus(status);
        }
    };    


    struct FIRSR32F
    {
        Ipp8u * pSpec;
        Ipp8u * pBuffer;
        size_t  size,nTaps;
        Ipp32f *pTaps;

        enum {
            LP,
            HP,
            BP,
            BS
        };
        int type = LP;

        FIRSR32F(size_t n, Ipp32f * taps, IppAlgType algType = ippAlgFFT) {
            int specSize,bufferSize;
            IppStatus status = ippsFIRSRGetSize(n,ipp32f,&specSize,&bufferSize);
            checkStatus(status);
            pSpec = ippsMalloc_8u(specSize);
            pBuffer= ippsMalloc_8u(bufferSize);
            //pDlySrc= ippsMalloc_32f(((n + up -1) / up);
            //pDlyDst= ippsMalloc_32f(((n + down -1) / down);
            nTaps   = n;
            pTaps   = ippsMalloc_32f(n);
            memcpy(pTaps,taps,n*sizeof(Ipp32f));
            status = ippsFIRSRInit_32f(pTaps,n,algType,(IppsFIRSpec_32f*)pSpec);
            checkStatus(status);
        }
        FIRSR32F(int type, size_t n, Ipp32f fc1, Ipp32f fc2=0, IppAlgType algType = ippAlgFFT) {
            int specSize,bufferSize;
            IppStatus status = ippsFIRSRGetSize(n,ipp32f,&specSize,&bufferSize);
            checkStatus(status);
            pSpec = ippsMalloc_8u(specSize);
            pBuffer= ippsMalloc_8u(bufferSize);
            //pDlySrc= ippsMalloc_32f(((n + up -1) / up);
            //pDlyDst= ippsMalloc_32f(((n + down -1) / down);
            pTaps   = ippsMalloc_32f(n);
            makeFilter(type,fc1,fc2);
            status = ippsFIRSRInit_32f(pTaps,n,algType,(IppsFIRSpec_32f*)pSpec);
            checkStatus(status);
        }
        ~FIRSR32F() {
            if(pSpec) ippsFree(pSpec);
            if(pBuffer) ippsFree(pBuffer);
            //if(pDlySrc) ippFree(pDlySrc);
            //if(pDlyDst) ippFree(pDlyDst);
            if(pTaps) ippFree(pTaps);
        }
        // bandpass/bandstop - fc1 = lower, fc2 = upper
        void makeFilter(int filterType, Ipp32f fc1, Ipp32f fc2=0) {
            type = filterType;
            switch(type) {
                case LP: makeLowPass(fc1);  break;
                case HP: makeHighPass(fc1);  break;
                case BP: makeBandPass(fc1,fc2);  break;
                case BS: makeBandStop(fc1,fc2);  break;
            }
        }
        void makeLowPass(float fc) {            
            for(size_t i = 0; i < nTaps; i++)
            {   
                double x= fc*(double)(i-nTaps/2);            
                pTaps[i] = fc * std::sin(M_PI*x)/(M_PI*x);
            }
            WinHamming(pTaps,pTaps,nTaps);
        }
        void makeHighPass(float fc) {            
            for(size_t i = 0; i < nTaps; i++)
            {   
                double x= fc*(double)(i-nTaps/2);            
                pTaps[i] = (std::sin(M_PI*x)/(M_PI*x)) - fc * std::sin(M_PI*x)/(M_PI*x);
            }            
            WinHamming(pTaps,pTaps,nTaps);
        }
        void makeBandPass(float fc1, float fc2) {
            for(size_t i = 0; i < nTaps; i++)
            {   
                double xl= fc1*(double)(i-nTaps/2);            
                double xu= fc2*(double)(i-nTaps/2);            
                pTaps[i] = (xu * std::sin(M_PI*xu)/(M_PI*xu)) - (xl * std::sin(M_PI*xl)/(M_PI*xl)) ;
            }
            WinHamming(pTaps,pTaps,nTaps);
        }
        void makeBandStop(float fc1, float fc2) {
            for(size_t i = 0; i < nTaps; i++)
            {   
                double xl= fc1*(double)(i-nTaps/2);            
                double xu= fc2*(double)(i-nTaps/2);            
                pTaps[i] = (xu * std::sin(M_PI*xu)/(M_PI*xu)) - (xl * std::sin(M_PI*xl)/(M_PI*xl)) ;
            }
            WinHamming(pTaps,pTaps,nTaps);
        }        
        void Execute(const Ipp32f * pSrc, Ipp32f* pDst, int numIters) {
            IppStatus status = ippsFIRSR_32f(pSrc,pDst,numIters,(IppsFIRSpec_32f*)pSpec,NULL,NULL,pBuffer);
            checkStatus(status);
        }
    };
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
        void Execute(const Ipp32f * pSrc, Ipp32f * pRef, Ipp32f* pDst, int len, float mu) {
            IppStatus status = ippsFIRLMS_32f(pSrc,pRef,pDst,len,mu,pState);
            checkStatus(status);
        }
    };
    
    double kaiserBeta(double As)
    {
        if(As > 50.0)
            return 0.1102 * (As - 8.7);
        else if(As >= 21.)
            return 0.5842 * std::pow(As - 21.0, 0.4) + 0.07886 * (As - 21.0);
        else
            return 0.;
    }

    int kaiserTapsEstimate(double delta, double As)
    {
        return int((As - 8.0) / (2.285 * delta) + 0.5);
    }
    struct Resample32F
    {
        IppsResamplingPolyphaseFixed_32f * pSpec;
        size_t N;
        Ipp64f Time;
        int Outlen;
        Resample32F(int inRate, int outRate, Ipp32f rollf=0.9f, Ipp32f as = 80.0f) {
            int size,height,length;
            Ipp32f alpha = kaiserBeta(as);
            Ipp32f delta = (1.0f - rollf) * M_PI;
            int n = kaiserTapsEstimate(delta,as);
            IppStatus status = ippsResamplePolyphaseFixedGetSize_32f(inRate,outRate,n,&size,&length,&height,ippAlgHintFast);
            checkStatus(status);
            pSpec = (IppsResamplingPolyphaseFixed_32f*)ippsMalloc_8u(size);
            status = ippsResamplePolyphaseFixedInit_32f(inRate,outRate,n,rollf,alpha,(IppsResamplingPolyphaseFixed_32f*)pSpec,ippAlgHintFast);
            checkStatus(status);            
        }
        ~Resample32F() {
            if(pSpec) ippsFree(pSpec);
        }
        void Execute(const Ipp32f * pSrc, int len, Ipp32f* pDst, Ipp64f norm=0.98) {
            IppStatus status = ippsResamplePolyphaseFixed_32f(pSrc,len,pDst,norm,&Time,&Outlen,(IppsResamplingPolyphaseFixed_32f*)pSpec);
            checkStatus(status);
        }
    };

    struct IIR32F
    {
        Ipp8u * pBuffer;
        IppsIIRState_32f * pState;
        
        // B0,B1..Border,A0,A1..Aorder = 2*(order+1)
        IIR32F(size_t n, int order, const Ipp32f * taps) {
            int bufferSize;
            IppStatus status = ippsIIRGetStateSize_32f(n,&bufferSize);
            checkStatus(status);
            pBuffer = ippsMalloc_8u(bufferSize);            
            status = ippsIIRInit_32f(&pState,taps,order,NULL,pBuffer);
            checkStatus(status);
        }
        ~IIR32F() {
            if(pBuffer) ippsFree(pBuffer);        
        }
        void setCoefficients(int n, int order, const Ipp32f* taps) {
            if(pBuffer) ippsFree(pBuffer);
            int bufferSize;
            IppStatus status = ippsIIRGetStateSize_32f(n,&bufferSize);
            checkStatus(status);
            pBuffer = ippsMalloc_8u(bufferSize);            
            status = ippsIIRInit_32f(&pState,taps,order,NULL,pBuffer);
            checkStatus(status);
        }
    };
    struct IIRBiquad32F
    {
        Ipp8u * pBuffer;
        IppsIIRState_32f * pState;
        
        // B0,B1..Border,A0,A1..Aorder = 2*(order+1)
        IIRBiquad32F(size_t n, int numBiquads, const Ipp32f * taps) {
            int bufferSize;
            IppStatus status = ippsIIRGetStateSize_32f(n,&bufferSize);
            checkStatus(status);
            pBuffer = ippsMalloc_8u(bufferSize);            
            status = ippsIIRInit_BiQuad_32f(&pState,taps,numBiquads,NULL,pBuffer);
            checkStatus(status);
        }
        ~IIRBiquad32F() {
            if(pBuffer) ippsFree(pBuffer);        
        }
        void setCoefficients(int order, int numBiquads, const Ipp32f * taps) {
            if(pBuffer) ippsFree(pBuffer);
            int bufferSize;
            IppStatus status = ippsIIRGetStateSize_32f(numBiquads,&bufferSize);
            checkStatus(status);
            pBuffer = ippsMalloc_8u(bufferSize);            
            status = ippsIIRInit_BiQuad_32f(&pState,taps,numBiquads,NULL,pBuffer);
            checkStatus(status);
        }
    };

    /*
    struct STFTF32
    {
        RFFTF fft;    
        size_t blockSize,windowSize,hopSize;

        STFTF32(size_t window_size, size_t hop_size) : fft(window_size) {            
            windowSize= window_size;
            hopSize = hop_size;
        }
        std::vector<std::vector<std::complex<Ipp32f>>> 
        stft(std::vector<float> data) {
            std::vector<std::vector<std::complex<float>>> output;
            std::vector<float> _data(windowSize);
            int data_size = data.size();    
            size_t result_size = (data_size / hopSize);                       
            size_t idx = 0;
            for (size_t pos = 0; pos < data_size; pos += hopSize) {
                for (size_t i = 0; i < windowSize; ++i) {
                    if (pos + i < data_size)
                        _data[i] = data[pos + i];
                    else
                        _data[i] = 0;
                }
                WinHann(_data.data(),_data.data(),windowSize);
                fft.set_input(_data.data());
                fft.Execute();
                fft.Normalize();
                temp = fft.get_output();
                output.push_back(temp);
                idx++;
            }            
            return output;
        }

        std::vector<std::vector<float>> 
        istft(complex_vector<float> data, size_t window_size, size_t hop_size) {
            int data_size = data.size();
            size_t result_size = data_size * hop_size + (window_size - hop_size);
            std::vector<std::vector<float>> output;
            std::vector<float> temp;
            std::vector<float> frame(window_size);
            std::vector<float> result(result_size);
            std::vector<std::complex<float>> slice(window_size);
                                    
            for (size_t i = 0; i < data_size; ++i) {                
                memcpy(slice.data(), data.data() + i * window_size, sizeof(fftw_complex) * window_size);
                fft.set_input(slice);
                fft.Execute();
                fft.Normalize();
                result = fft.get_output();
                output.push_back(result);
            }                
            return output;
        }
    };
    */
}