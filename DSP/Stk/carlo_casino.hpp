// going to specialize everything
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


// this is just std::vector
#include "carlo_samples.hpp"
#include "carlo_samplesdsp.hpp"

// this is MKL stuff
//#include "mkl_array.hpp"
//#include "mkl_vector.hpp"
//#include "mkl_matrix.hpp"

// we can do SIMD kernels with these
//#include "vector2d.hpp"
//#include "vector4d.hpp"
//#include "vector4f.hpp"
//#include "vector8f.hpp"

// Eigen
//#include "TinyEigen.hpp"
//#include "SampleVector.h"

namespace Casino
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
}

#include "carlo_casinodsp.hpp"