#pragma once
#include <cassert>
#include <cstdlib>
#include <memory>
#include <iostream>
#include <complex>
#include <algorithm>
#include <mkl.h>

namespace Sample {
    
    template<typename T, int align=256>
    struct Array
    {
    private:
        std::shared_ptr<T> array;
        
        void allocate(size_t n) {
            array.reset();
            N = n;
            _array = (T*)mkl_malloc(n*sizeof(T),256);            
            assert(_array != nullptr);
            array = std::shared_ptr<T>(_array, [](T * ptr) { mkl_free(ptr); });
        }
    public:
        
        T* _array;
        size_t N;

        Array(size_t n, T val = T()) {
            allocate(n);
            fill(val);
        }
        Array(T * ptr, size_t n) {
            _array = ptr;
            array  = std::shared_ptr<T>(_array, [](T * ptr) { mkl_free(ptr); });
        }
        Array(const Array & a) {        
            allocate(a.N);
            copy(a);
        }

        
        T& operator[](size_t i) { return _array[i-1]; }

        Array& operator = (const Array& a) {
            array.reset();
            array = a.array;
            _array = a._array;
            N     = a.N;
            return *this;
        }


        Array operator + (const Array& a) {
            Array<T> r(N);
            for(size_t i = 0; i < N; i++) r._array[i] = _array[i] + a._array[i];
            return r;
        }
        Array operator - (const Array& a) {
            Array<T> r(N);
            for(size_t i = 0; i < N; i++) r._array[i] = _array[i] - a._array[i];
            return r;
        }
        Array operator * (const Array& a) {
            Array<T> r(N);
            for(size_t i = 0; i < N; i++) r._array[i] = _array[i] * a._array[i];
            return r;
        }
        Array operator / (const Array& a) {
            Array<T> r(N);
            for(size_t i = 0; i < N; i++) r._array[i] = _array[i] / a._array[i];
            return r;
        }
        
        
        void copy(const Array & a) {            
            memcpy(_array,a._array,N*sizeof(T));
        }
        void resize(size_t n) {
            array.reset();
            allocate(n);
        }    
        void fill(const T& val) {
            for(size_t i = 0; i < N; i++) _array[i] = val;
        }
        size_t size() const { return N; }
        
        void print() 
        {
            std::cout << "Array[" << N << "]=";
            for(size_t i = 0; i < N-1; i++)
                std::cout << _array[i] << ",";
            std::cout << _array[N-1] << std::endl;
        }
    };

    template<>
    Array<float> Array<float>::operator + (const Array<float> &b) {
        Array<float> r(N);        
        vsAdd(N,_array,b._array,r._array);
        return r;
    }
    template<>
    Array<double> Array<double>::operator + (const Array<double> &b) {
        Array<double> r(N);
        vdAdd(N,_array,b._array,r._array);
        return r;
    }
    template<> 
    Array<std::complex<float>> Array<std::complex<float>>::operator + (const Array<std::complex<float>> &b) {
        Array<std::complex<float>> r(N);
        vcAdd(N,(MKL_Complex8*)_array,(MKL_Complex8*)b._array,(MKL_Complex8*)r._array);
        return r;
    }
    template<>
    Array<std::complex<double>> Array<std::complex<double>>::operator + (const Array<std::complex<double>> &b) {
        Array<std::complex<double>> r(N);        
        vzAdd(N,(MKL_Complex16*)_array,(MKL_Complex16*)b._array,(MKL_Complex16*)r._array);
        return r;
    }  

    template<>
    Array<float> Array<float>::operator - (const Array<float> &b) {
        Array<float> r(N);        
        vsSub(N,_array,b._array,r._array);
        return r;
    }
    template<>
    Array<double> Array<double>::operator - (const Array<double> &b) {
        Array<double> r(N);
        vdSub(N,_array,b._array,r._array);
        return r;
    }
    template<> 
    Array<std::complex<float>> Array<std::complex<float>>::operator - (const Array<std::complex<float>> &b) {
        Array<std::complex<float>> r(N);
        vcSub(N,(MKL_Complex8*)_array,(MKL_Complex8*)b._array,(MKL_Complex8*)r._array);
        return r;
    }
    template<>
    Array<std::complex<double>> Array<std::complex<double>>::operator - (const Array<std::complex<double>> &b) {
        Array<std::complex<double>> r(N);        
        vzSub(N,(MKL_Complex16*)_array,(MKL_Complex16*)b._array,(MKL_Complex16*)r._array);
        return r;
    }  

    template<>
    Array<float> Array<float>::operator * (const Array<float> &b) {
        Array<float> r(N);        
        vsMul(N,_array,b._array,r._array);
        return r;
    }
    template<>
    Array<double> Array<double>::operator * (const Array<double> &b) {
        Array<double> r(N);
        vdMul(N,_array,b._array,r._array);
        return r;
    }
    template<> 
    Array<std::complex<float>> Array<std::complex<float>>::operator * (const Array<std::complex<float>> &b) {
        Array<std::complex<float>> r(N);
        vcMul(N,(MKL_Complex8*)_array,(MKL_Complex8*)b._array,(MKL_Complex8*)r._array);
        return r;
    }
    template<>
    Array<std::complex<double>> Array<std::complex<double>>::operator * (const Array<std::complex<double>> &b) {
        Array<std::complex<double>> r(N);        
        vzMul(N,(MKL_Complex16*)_array,(MKL_Complex16*)b._array,(MKL_Complex16*)r._array);
        return r;
    }  

    template<>
    Array<float> Array<float>::operator / (const Array<float> &b) {
        Array<float> r(N);        
        vsDiv(N,_array,b._array,r._array);
        return r;
    }
    template<>
    Array<double> Array<double>::operator / (const Array<double> &b) {
        Array<double> r(N);
        vdDiv(N,_array,b._array,r._array);
        return r;
    }
    template<> 
    Array<std::complex<float>> Array<std::complex<float>>::operator / (const Array<std::complex<float>> &b) {
        Array<std::complex<float>> r(N);
        vcDiv(N,(MKL_Complex8*)_array,(MKL_Complex8*)b._array,(MKL_Complex8*)r._array);
        return r;
    }
    template<>
    Array<std::complex<double>> Array<std::complex<double>>::operator / (const Array<std::complex<double>> &b) {
        Array<std::complex<double>> r(N);        
        vzDiv(N,(MKL_Complex16*)_array,(MKL_Complex16*)b._array,(MKL_Complex16*)r._array);
        return r;
    }  


    inline float sum(const Array<float> & a) {        
        return cblas_sasum(a.N,a._array,1);        
    }

    inline double sum(const Array<double> & a) {        
        return cblas_dasum(a.N,a._array,1);        
    }

    inline float sum(const Array<std::complex<float>> &a) {        
        return (cblas_scasum(a.N,(MKL_Complex8*)a._array,1));        
    }

    inline double sum(const Array<std::complex<double>> & a) {        
        return (cblas_dzasum(a.N,(MKL_Complex16*)a._array,1));    
    }  
    
    inline Array<float> sqr(const Array<float> & a) {
        Array<float> r(a.N);        
        vsSqr(a.N,a._array,r._array);
        return r;
    }
    
    inline Array<double> sqr(const Array<double> & a) {
        Array<double> r(a.N);
        vdSqr(a.N,a._array,r._array);
        return r;
    }
    
    inline Array<float> abs(const Array<float> & a) {
        Array<float> r(a.N);        
        vsAbs(a.N,a._array,r._array);
        return r;
    }
    
    inline Array<double> abs(const Array<double> & a) {
        Array<double> r(a.N);
        vdAbs(a.N,a._array,r._array);
        return r;
    }
    inline Array<float> abs(const Array<std::complex<float>> & a) {
        Array<float> r(a.N);        
        vcAbs(a.N,(MKL_Complex8*)a._array,r._array);
        return r;
    }
    
    inline Array<double> abs(const Array<std::complex<double>> & a) {
        Array<double> r(a.N);
        vzAbs(a.N,(MKL_Complex16*)a._array,r._array);
        return r;
    }
    
    
    inline  Array<std::complex<float>> mulByConj(const Array<std::complex<float>> & a, const Array<std::complex<float>> & b) {
        Array<std::complex<float>> r(a.N);
        vcMulByConj(a.N,(MKL_Complex8*)a._array,(MKL_Complex8*)b._array,(MKL_Complex8*)r._array);
        return r;
    }
    inline  Array<std::complex<double>> mulByConj(const Array<std::complex<double>> & a, const Array<std::complex<double>> & b) {
        Array<std::complex<double>> r(a.N);
        vzMulByConj(a.N,(MKL_Complex16*)a._array,(MKL_Complex16*)b._array,(MKL_Complex16*)r._array);
        return r;
    }

    inline  Array<std::complex<float>> conj(const Array<std::complex<float>> & a) {
        Array<std::complex<float>> r(a.N);
        vcConj(a.N,(MKL_Complex8*)a._array,(MKL_Complex8*)r._array);
        return r;
    }
    inline  Array<std::complex<double>> conj(const Array<std::complex<double>> & a) {
        Array<std::complex<double>> r(a.N);
        vzConj(a.N,(MKL_Complex16*)a._array,(MKL_Complex16*)r._array);
        return r;
    }

    inline  Array<std::complex<float>> arg(const Array<std::complex<float>> & a) {
        Array<std::complex<float>> r(a.N);
        vcConj(a.N,(MKL_Complex8*)a._array,(MKL_Complex8*)r._array);
        return r;
    }
    inline  Array<std::complex<double>> arg(const Array<std::complex<double>> & a) {
        Array<std::complex<double>> r(a.N);
        vzConj(a.N,(MKL_Complex16*)a._array,(MKL_Complex16*)r._array);
        return r;
    }

    inline Array<float> fmod(const Array<float> & a, const Array<float> & b) {
        Array<float> r(a.N);        
        vsFmod(a.N,a._array,b._array,r._array);
        return r;
    }

    inline Array<double> fmod(const Array<double> & a, const Array<double> & b) {
        Array<double> r(a.N);        
        vdFmod(a.N,a._array,b._array,r._array);
        return r;
    }

    inline Array<float> remainder(const Array<float> & a, const Array<float> & b) {
        Array<float> r(a.N);        
        vsRemainder(a.N,a._array,b._array,r._array);
        return r;
    }

    inline Array<double> remainder(const Array<double> & a, const Array<double> & b) {
        Array<double> r(a.N);        
        vdRemainder(a.N,a._array,b._array,r._array);
        return r;
    }

    
    inline Array<float> recip(const Array<float> & a) {
        Array<float> r(a.N);        
        vsInv(a.N,a._array,r._array);
        return r;
    }
    
    inline Array<double> recip(const Array<double> & a) {
        Array<double> r(a.N);
        vdInv(a.N,a._array,r._array);
        return r;
    }

    inline Array<float> sqrt(const Array<float> & a) {
        Array<float> r(a.N);        
        vsSqrt(a.N,a._array,r._array);
        return r;
    }
    
    inline Array<double> sqrt(const Array<double> & a) {
        Array<double> r(a.N);
        vdSqrt(a.N,a._array,r._array);
        return r;
    }
    inline Array<std::complex<float>> sqrt(const Array<std::complex<float>> & a) {
        Array<std::complex<float>> r(a.N);        
        vcSqrt(a.N,(MKL_Complex8*)a._array,(MKL_Complex8*)r._array);
        return r;
    }
    
    inline Array<std::complex<double>> sqrt(const Array<std::complex<double>> & a) {
        Array<std::complex<double>> r(a.N);
        vzSqrt(a.N,(MKL_Complex16*)a._array,(MKL_Complex16*)r._array);
        return r;
    }

    inline Array<float> rsqrt(const Array<float> & a) {
        Array<float> r(a.N);        
        vsInvSqrt(a.N,a._array,r._array);
        return r;
    }
    
    inline Array<double> rsqrt(const Array<double> & a) {
        Array<double> r(a.N);
        vdInvSqrt(a.N,a._array,r._array);
        return r;
    }
    
    inline Array<float> cbrt(const Array<float> & a) {
        Array<float> r(a.N);        
        vsCbrt(a.N,a._array,r._array);
        return r;
    }
    
    inline Array<double> cbrt(const Array<double> & a) {
        Array<double> r(a.N);
        vdCbrt(a.N,a._array,r._array);
        return r;
    }
    inline Array<float> rcbrt(const Array<float> & a) {
        Array<float> r(a.N);        
        vsInvCbrt(a.N,a._array,r._array);
        return r;
    }
    
    inline Array<double> rcbrt(const Array<double> & a) {
        Array<double> r(a.N);
        vdInvCbrt(a.N,a._array,r._array);
        return r;
    }

    inline Array<float> pow(const Array<float> & a, const Array<float> & b) {
        Array<float> r(a.N);        
        vsPow(a.N,a._array,b._array,r._array);
        return r;
    }
    inline Array<double> pow(const Array<double> & a, const Array<double> & b) {
        Array<double> r(a.N);        
        vdPow(a.N,a._array,b._array,r._array);
        return r;
    }
    inline Array<std::complex<float>> pow(const Array<std::complex<float>> & a, const Array<std::complex<float>> & b) {
        Array<std::complex<float>> r(a.N);        
        vcPow(a.N,(MKL_Complex8*)a._array,(MKL_Complex8*)b._array,(MKL_Complex8*)r._array);
        return r;
    }
    inline Array<std::complex<double>> pow(const Array<std::complex<double>> & a, const Array<std::complex<double>> & b) {
        Array<std::complex<double>> r(a.N);        
        vzPow(a.N,(MKL_Complex16*)a._array,(MKL_Complex16*)b._array,(MKL_Complex16*)r._array);
        return r;
    }

    inline Array<float> powx(const Array<float> & a, const float b) {
        Array<float> r(a.N);        
        vsPowx(a.N,a._array,b,r._array);
        return r;
    }
    inline Array<double> pow(const Array<double> & a, const double b) {
        Array<double> r(a.N);        
        vdPowx(a.N,a._array,b,r._array);
        return r;
    }
    /*
    inline Array<std::complex<float>> powx(const Array<std::complex<float>> & a, const std::complex<float> & b) {
        Array<std::complex<float>> r(a.N);        
        vcPowx(a.N,(MKL_Complex8*)a._array,static_cast<MKL_Complex8>(b),(MKL_Complex8*)r._array);
        return r;
    }
    inline Array<std::complex<double>> powx(const Array<std::complex<double>> & a, const std::complex<double> & b) {
        Array<std::complex<double>> r(a.N);        
        vzPowx(a.N,(MKL_Complex16*)a._array,(MKL_Complex16)b._array,(MKL_Complex16*)r._array);
        return r;
    }
    */

    
    inline Array<float> exp(const Array<float> & a) {
        Array<float> r(a.N);        
        vsExp(a.N,a._array,r._array);
        return r;
    }    
    inline Array<double> exp(const Array<double> & a) {
        Array<double> r(a.N);
        vdExp(a.N,a._array,r._array);
        return r;
    }
    inline Array<float> exp(const Array<std::complex<float>> & a) {
        Array<float> r(a.N);        
        vcExp(a.N,(MKL_Complex8*)a._array,(MKL_Complex8*)r._array);
        return r;
    }
    
    inline Array<double> exp(const Array<std::complex<double>> & a) {
        Array<double> r(a.N);
        vzExp(a.N,(MKL_Complex16*)a._array,(MKL_Complex16*)r._array);
        return r;
    }

    inline Array<float> exp2(const Array<float> & a) {
        Array<float> r(a.N);        
        vsExp2(a.N,a._array,r._array);
        return r;
    }
    inline Array<double> exp2(const Array<double> & a) {
        Array<double> r(a.N);
        vdExp2(a.N,a._array,r._array);
        return r;
    }
    inline Array<float> exp10(const Array<float> & a) {
        Array<float> r(a.N);        
        vsExp10(a.N,a._array,r._array);
        return r;
    }    
    inline Array<double> exp10(const Array<double> & a) {
        Array<double> r(a.N);
        vdExp10(a.N,a._array,r._array);
        return r;
    }
    inline Array<float> expm1(const Array<float> & a) {
        Array<float> r(a.N);        
        vsExpm1(a.N,a._array,r._array);
        return r;
    }    
    inline Array<double> expm1(const Array<double> & a) {
        Array<double> r(a.N);
        vdExpm1(a.N,a._array,r._array);
        return r;
    }

    inline Array<float> ln(const Array<float> & a) {
        Array<float> r(a.N);        
        vsLn(a.N,a._array,r._array);
        return r;
    }    
    inline Array<double> ln(const Array<double> & a) {
        Array<double> r(a.N);
        vdLn(a.N,a._array,r._array);
        return r;
    }
    inline Array<float> ln(const Array<std::complex<float>> & a) {
        Array<float> r(a.N);        
        vcLn(a.N,(MKL_Complex8*)a._array,(MKL_Complex8*)r._array);
        return r;
    }    
    inline Array<double> ln(const Array<std::complex<double>> & a) {
        Array<double> r(a.N);
        vzLn(a.N,(MKL_Complex16*)a._array,(MKL_Complex16*)r._array);
        return r;
    }
    inline Array<float> log2(const Array<float> & a) {
        Array<float> r(a.N);        
        vsLog2(a.N,a._array,r._array);
        return r;
    }    
    inline Array<double> log2(const Array<double> & a) {
        Array<double> r(a.N);
        vdLog2(a.N,a._array,r._array);
        return r;
    }
    inline Array<float> log10(const Array<float> & a) {
        Array<float> r(a.N);        
        vsLog10(a.N,a._array,r._array);
        return r;
    }    
    inline Array<double> log10(const Array<double> & a) {
        Array<double> r(a.N);
        vdLog10(a.N,a._array,r._array);
        return r;
    }
    inline Array<float> log1p(const Array<float> & a) {
        Array<float> r(a.N);        
        vsLog1p(a.N,a._array,r._array);
        return r;
    }    
    inline Array<double> log1p(const Array<double> & a) {
        Array<double> r(a.N);
        vdLog1p(a.N,a._array,r._array);
        return r;
    }
    inline Array<float> logb(const Array<float> & a) {
        Array<float> r(a.N);        
        vsLogb(a.N,a._array,r._array);
        return r;
    }    
    inline Array<double> logb(const Array<double> & a) {
        Array<double> r(a.N);
        vdLogb(a.N,a._array,r._array);
        return r;
    }

    inline Array<float> sin(const Array<float> & a) {
        Array<float> r(a.N);        
        vsSin(a.N,a._array,r._array);
        return r;
    }    
    inline Array<double> sin(const Array<double> & a) {
        Array<double> r(a.N);
        vdSin(a.N,a._array,r._array);
        return r;
    }
    inline Array<float> sin(const Array<std::complex<float>> & a) {
        Array<float> r(a.N);        
        vcSin(a.N,(MKL_Complex8*)a._array,(MKL_Complex8*)r._array);
        return r;
    }    
    inline Array<double> sin(const Array<std::complex<double>> & a) {
        Array<double> r(a.N);
        vzSin(a.N,(MKL_Complex16*)a._array,(MKL_Complex16*)r._array);
        return r;
    }

    inline Array<float> cos(const Array<float> & a) {
        Array<float> r(a.N);        
        vsCos(a.N,a._array,r._array);
        return r;
    }    
    inline Array<double> cos(const Array<double> & a) {
        Array<double> r(a.N);
        vdCos(a.N,a._array,r._array);
        return r;
    }
    inline Array<float> cos(const Array<std::complex<float>> & a) {
        Array<float> r(a.N);        
        vcCos(a.N,(MKL_Complex8*)a._array,(MKL_Complex8*)r._array);
        return r;
    }    
    inline Array<double> cos(const Array<std::complex<double>> & a) {
        Array<double> r(a.N);
        vzCos(a.N,(MKL_Complex16*)a._array,(MKL_Complex16*)r._array);
        return r;
    }
    inline Array<float> tan(const Array<float> & a) {
        Array<float> r(a.N);        
        vsTan(a.N,a._array,r._array);
        return r;
    }    
    inline Array<double> tan(const Array<double> & a) {
        Array<double> r(a.N);
        vdTan(a.N,a._array,r._array);
        return r;
    }
    inline Array<float> tan(const Array<std::complex<float>> & a) {
        Array<float> r(a.N);        
        vcTan(a.N,(MKL_Complex8*)a._array,(MKL_Complex8*)r._array);
        return r;
    }    
    inline Array<double> tan(const Array<std::complex<double>> & a) {
        Array<double> r(a.N);
        vzTan(a.N,(MKL_Complex16*)a._array,(MKL_Complex16*)r._array);
        return r;
    }
    inline void sincos(const Array<float> & a, Array<float> & sins, Array<float> & coss) {
        sins.resize(a.size());
        coss.resize(a.size());
        vsSinCos(a.N,a._array,sins._array,coss._array);    
    }    
    inline void sincos(const Array<double> & a, Array<double> & sins, Array<double> & coss) {
        sins.resize(a.size());
        coss.resize(a.size());
        vdSinCos(a.N,a._array,sins._array,coss._array);    
    }    
    

    inline Array<float> asin(const Array<float> & a) {
        Array<float> r(a.N);        
        vsAsin(a.N,a._array,r._array);
        return r;
    }    
    inline Array<double> asin(const Array<double> & a) {
        Array<double> r(a.N);
        vdAsin(a.N,a._array,r._array);
        return r;
    }
    inline Array<float> asin(const Array<std::complex<float>> & a) {
        Array<float> r(a.N);        
        vcAsin(a.N,(MKL_Complex8*)a._array,(MKL_Complex8*)r._array);
        return r;
    }    
    inline Array<double> asin(const Array<std::complex<double>> & a) {
        Array<double> r(a.N);
        vzAsin(a.N,(MKL_Complex16*)a._array,(MKL_Complex16*)r._array);
        return r;
    }

    inline Array<float> acos(const Array<float> & a) {
        Array<float> r(a.N);        
        vsAcos(a.N,a._array,r._array);
        return r;
    }    
    inline Array<double> acos(const Array<double> & a) {
        Array<double> r(a.N);
        vdAcos(a.N,a._array,r._array);
        return r;
    }
    inline Array<float> acos(const Array<std::complex<float>> & a) {
        Array<float> r(a.N);        
        vcAcos(a.N,(MKL_Complex8*)a._array,(MKL_Complex8*)r._array);
        return r;
    }    
    inline Array<double> acos(const Array<std::complex<double>> & a) {
        Array<double> r(a.N);
        vzAcos(a.N,(MKL_Complex16*)a._array,(MKL_Complex16*)r._array);
        return r;
    }
    inline Array<float> atan(const Array<float> & a) {
        Array<float> r(a.N);        
        vsAtan(a.N,a._array,r._array);
        return r;
    }    
    inline Array<double> atan(const Array<double> & a) {
        Array<double> r(a.N);
        vdAtan(a.N,a._array,r._array);
        return r;
    }
    inline Array<float> atan(const Array<std::complex<float>> & a) {
        Array<float> r(a.N);        
        vcAtan(a.N,(MKL_Complex8*)a._array,(MKL_Complex8*)r._array);
        return r;
    }    
    inline Array<double> atan(const Array<std::complex<double>> & a) {
        Array<double> r(a.N);
        vzAtan(a.N,(MKL_Complex16*)a._array,(MKL_Complex16*)r._array);
        return r;
    }
    inline Array<float> atan2(const Array<float> & a, const Array<float> & b) {
        Array<float> r(a.N);        
        vsAtan2(a.N,a._array,b._array, r._array);
        return r;
    }    
    inline Array<double> atan2(const Array<double> & a, const Array<double> & b) {
        Array<double> r(a.N);
        vdAtan2(a.N,a._array,b._array,r._array);
        return r;
    }

    inline Array<float> sinh(const Array<float> & a) {
        Array<float> r(a.N);        
        vsSinh(a.N,a._array,r._array);
        return r;
    }    
    inline Array<double> sinh(const Array<double> & a) {
        Array<double> r(a.N);
        vdSinh(a.N,a._array,r._array);
        return r;
    }
    inline Array<float> sinh(const Array<std::complex<float>> & a) {
        Array<float> r(a.N);        
        vcSinh(a.N,(MKL_Complex8*)a._array,(MKL_Complex8*)r._array);
        return r;
    }    
    inline Array<double> sinh(const Array<std::complex<double>> & a) {
        Array<double> r(a.N);
        vzSinh(a.N,(MKL_Complex16*)a._array,(MKL_Complex16*)r._array);
        return r;
    }

    inline Array<float> cosh(const Array<float> & a) {
        Array<float> r(a.N);        
        vsCosh(a.N,a._array,r._array);
        return r;
    }    
    inline Array<double> cosh(const Array<double> & a) {
        Array<double> r(a.N);
        vdCosh(a.N,a._array,r._array);
        return r;
    }
    inline Array<float> cosh(const Array<std::complex<float>> & a) {
        Array<float> r(a.N);        
        vcCosh(a.N,(MKL_Complex8*)a._array,(MKL_Complex8*)r._array);
        return r;
    }    
    inline Array<double> cosh(const Array<std::complex<double>> & a) {
        Array<double> r(a.N);
        vzCosh(a.N,(MKL_Complex16*)a._array,(MKL_Complex16*)r._array);
        return r;
    }
    inline Array<float> tanh(const Array<float> & a) {
        Array<float> r(a.N);        
        vsTanh(a.N,a._array,r._array);
        return r;
    }    
    inline Array<double> tanh(const Array<double> & a) {
        Array<double> r(a.N);
        vdTanh(a.N,a._array,r._array);
        return r;
    }
    inline Array<float> tahn(const Array<std::complex<float>> & a) {
        Array<float> r(a.N);        
        vcTanh(a.N,(MKL_Complex8*)a._array,(MKL_Complex8*)r._array);
        return r;
    }    
    inline Array<double> tanh(const Array<std::complex<double>> & a) {
        Array<double> r(a.N);
        vzTanh(a.N,(MKL_Complex16*)a._array,(MKL_Complex16*)r._array);
        return r;
    }

    inline Array<float> asinh(const Array<float> & a) {
        Array<float> r(a.N);        
        vsAsinh(a.N,a._array,r._array);
        return r;
    }    
    inline Array<double> asinh(const Array<double> & a) {
        Array<double> r(a.N);
        vdAsinh(a.N,a._array,r._array);
        return r;
    }
    inline Array<float> asinh(const Array<std::complex<float>> & a) {
        Array<float> r(a.N);        
        vcAsinh(a.N,(MKL_Complex8*)a._array,(MKL_Complex8*)r._array);
        return r;
    }    
    inline Array<double> asinh(const Array<std::complex<double>> & a) {
        Array<double> r(a.N);
        vzAsinh(a.N,(MKL_Complex16*)a._array,(MKL_Complex16*)r._array);
        return r;
    }

    inline Array<float> acosh(const Array<float> & a) {
        Array<float> r(a.N);        
        vsAcosh(a.N,a._array,r._array);
        return r;
    }    
    inline Array<double> acosh(const Array<double> & a) {
        Array<double> r(a.N);
        vdAcosh(a.N,a._array,r._array);
        return r;
    }
    inline Array<float> acosh(const Array<std::complex<float>> & a) {
        Array<float> r(a.N);        
        vcAcosh(a.N,(MKL_Complex8*)a._array,(MKL_Complex8*)r._array);
        return r;
    }    
    inline Array<double> acosh(const Array<std::complex<double>> & a) {
        Array<double> r(a.N);
        vzAcosh(a.N,(MKL_Complex16*)a._array,(MKL_Complex16*)r._array);
        return r;
    }
    inline Array<float> atanh(const Array<float> & a) {
        Array<float> r(a.N);        
        vsAtanh(a.N,a._array,r._array);
        return r;
    }    
    inline Array<double> atanh(const Array<double> & a) {
        Array<double> r(a.N);
        vdAtanh(a.N,a._array,r._array);
        return r;
    }
    inline Array<float> atanh(const Array<std::complex<float>> & a) {
        Array<float> r(a.N);        
        vcAtanh(a.N,(MKL_Complex8*)a._array,(MKL_Complex8*)r._array);
        return r;
    }    
    inline Array<double> atanh(const Array<std::complex<double>> & a) {
        Array<double> r(a.N);
        vzAtanh(a.N,(MKL_Complex16*)a._array,(MKL_Complex16*)r._array);
        return r;
    }
    
    inline Array<float> floor(const Array<float> & a) {
        Array<float> r(a.N);        
        vsFloor(a.N,a._array,r._array);
        return r;
    }    
    inline Array<double> floor(const Array<double> & a) {
        Array<double> r(a.N);
        vdFloor(a.N,a._array,r._array);
        return r;
    }
    inline Array<float> ceil(const Array<float> & a) {
        Array<float> r(a.N);        
        vsCeil(a.N,a._array,r._array);
        return r;
    }    
    inline Array<double> ceil(const Array<double> & a) {
        Array<double> r(a.N);
        vdCeil(a.N,a._array,r._array);
        return r;
    }
    inline Array<float> trunc(const Array<float> & a) {
        Array<float> r(a.N);        
        vsTrunc(a.N,a._array,r._array);
        return r;
    }    
    inline Array<double> trunc(const Array<double> & a) {
        Array<double> r(a.N);
        vdTrunc(a.N,a._array,r._array);
        return r;
    }
    inline Array<float> round(const Array<float> & a) {
        Array<float> r(a.N);        
        vsRound(a.N,a._array,r._array);
        return r;
    }    
    inline Array<double> round(const Array<double> & a) {
        Array<double> r(a.N);
        vdRound(a.N,a._array,r._array);
        return r;
    }

    inline Array<float> frac(const Array<float> & a) {
        Array<float> r(a.N);        
        vsFrac(a.N,a._array,r._array);
        return r;
    }    
    inline Array<double> frac(const Array<double> & a) {
        Array<double> r(a.N);
        vdFrac(a.N,a._array,r._array);
        return r;
    }


    template<typename T>
    T sum(const Array<T> & a) {
        T x = T(0);
        for(size_t i = 0; i < a.N; i++) x += a._array[i];
        return x;
    }
    template<typename T>
    T min(const Array<T> & a) {
        T min = a._array[0];
        for(size_t i = 1; i < a.N; i++) 
            if(a._array[i] < min)
                min = a._array[i];
        return min;            
    }
    template<typename T>
    T max(const Array<T> & a) {
        T max = a._array[0];
        for(size_t i = 1; i < a.N; i++) 
            if(a._array[i] > max)
                max = a._array[i];
        return max;            
    }
    template<typename T>
    size_t min_element(const Array<T> & a) {
        T min = a._array[0];
        size_t mini = 0;
        for(size_t i = 1; i < a.N; i++)
            if(a._array[i] < min) 
            {
                min = a._array[i];
                mini= i;
            }
        return mini+1;
    }
    template<typename T>
    size_t max_element(const Array<T> & a) {
        T max = a._array[0];
        size_t maxi = 0;
        for(size_t i = 1; i < a.N; i++)
            if(a._array[i] > max) 
            {
                max = a._array[i];
                maxi= i;
            }
        return maxi+1;
    }

} // MKL