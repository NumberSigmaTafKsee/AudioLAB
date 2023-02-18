#pragma once
#include <cassert>
#include <cstdlib>
#include <memory>
#include <iostream>
#include <complex>
#include <algorithm>
#include <mkl.h>

namespace MKL {
    
    template<typename T, int align=64>
    struct Array2D
    {
    private:
        std::shared_ptr<T> array;
        
        void allocate(size_t m, size_t n) {
            array.reset();
            M = m;
            N = n;
            _array = (T*)mkl_malloc(M*N*sizeof(T),align);            
            assert(_array != nullptr);
            array = std::shared_ptr<T>(_array, [](T * ptr) { mkl_free(ptr); });
        }
    public:
        
        T* _array;
        size_t M,N;

        Array2D(size_t n, T val = T()) {
            allocate(n);
            fill(val);
        }
        Array2D(T * ptr, size_t n) {
            _array = ptr;
            array  = std::shared_ptr<T>(_array, [](T * ptr) { mkl_free(ptr); });
        }
        Array2D(const Array2D & a) {        
            allocate(a.N);
            copy(a);
        }

        
        T& operator[](size_t i) { return _array[i-1]; }
        T& operator()(size_t i, size_t j) { return _array[i*N + j]; }

        Array2D& operator = (const Array2D& a) {
            array.reset();
            array = a.array;
            _array = a._array;
            M     = a.M;
            N     = a.N;
            return *this;
        }


        Array2D operator + (const Array2D& a) {
            Array2D<T> r(M,N);
            for(size_t i = 0; i < M*N; i++) r._array[i] = _array[i] + a._array[i];
            return r;
        }
        Array2D operator - (const Array2D& a) {
            Array2D<T> r(M*N);
            for(size_t i = 0; i < M*N; i++) r._array[i] = _array[i] - a._array[i];
            return r;
        }
        Array2D operator * (const Array2D& a) {
            Array2D<T> r(M*N);
            for(size_t i = 0; i < M*N; i++) r._array[i] = _array[i] * a._array[i];
            return r;
        }
        Array2D operator / (const Array2D& a) {
            Array2D<T> r(M*N);
            for(size_t i = 0; i < M*N; i++) r._array[i] = _array[i] / a._array[i];
            return r;
        }
        
        void copy(const Array2D & a) {            
            memcpy(_array,a._array,N*sizeof(T));
        }
        void resize(size_t m, size_t n) {
            array.reset();
            allocate(m,n);
        }    
        void fill(const T& val) {
            for(size_t i = 0; i < M*N; i++) _array[i] = val;
        }
        size_t size() const { return M*N; }
        
        void print() 
        {
            std::cout << "Array2D[" << M << "," << N << "]=";
            for(size_t j = 0; j < M; j++)
                for(size_t i = 0; i < N; i++)
                    std::cout << _array[j*N+i] << ",";
                std::cout << std::endl;
        }
            std::cout << _array[N-1] << std::endl;
        }
    };

    template<>
    Array2D<float> Array2D<float>::operator + (const Array2D<float> &b) {
        Array2D<float> r(M,N);        
        vsAdd(M*N,_array,b._array,r._array);
        return r;
    }
    template<>
    Array2D<double> Array2D<double>::operator + (const Array2D<double> &b) {
        Array2D<double> r(M,N);
        vdAdd(M*N,_array,b._array,r._array);
        return r;
    }
    template<> 
    Array2D<std::complex<float>> Array2D<std::complex<float>>::operator + (const Array2D<std::complex<float>> &b) {
        Array2D<std::complex<float>> r(M,N);
        vcAdd(M*N,(MKL_Complex8*)_array,(MKL_Complex8*)b._array,(MKL_Complex8*)r._array);
        return r;
    }
    template<>
    Array2D<std::complex<double>> Array2D<std::complex<double>>::operator + (const Array2D<std::complex<double>> &b) {
        Array2D<std::complex<double>> r(M,N);        
        vzAdd(M*N,(MKL_Complex16*)_array,(MKL_Complex16*)b._array,(MKL_Complex16*)r._array);
        return r;
    }  

    template<>
    Array2D<float> Array2D<float>::operator - (const Array2D<float> &b) {
        Array2D<float> r(M,N);        
        vsSub(M*N,_array,b._array,r._array);
        return r;
    }
    template<>
    Array2D<double> Array2D<double>::operator - (const Array2D<double> &b) {
        Array2D<double> r(M,N);
        vdSub(M,N,_array,b._array,r._array);
        return r;
    }
    template<> 
    Array2D<std::complex<float>> Array2D<std::complex<float>>::operator - (const Array2D<std::complex<float>> &b) {
        Array2D<std::complex<float>> r(M,N);
        vcSub(M,N,(MKL_Complex8*)_array,(MKL_Complex8*)b._array,(MKL_Complex8*)r._array);
        return r;
    }
    template<>
    Array2D<std::complex<double>> Array2D<std::complex<double>>::operator - (const Array2D<std::complex<double>> &b) {
        Array2D<std::complex<double>> r(M,N);        
        vzSub(M*N,(MKL_Complex16*)_array,(MKL_Complex16*)b._array,(MKL_Complex16*)r._array);
        return r;
    }  

    template<>
    Array2D<float> Array2D<float>::operator * (const Array2D<float> &b) {
        Array2D<float> r(M,N);        
        vsMul(M*N,_array,b._array,r._array);
        return r;
    }
    template<>
    Array2D<double> Array2D<double>::operator * (const Array2D<double> &b) {
        Array2D<double> r(M,N);
        vdMul(M*N,_array,b._array,r._array);
        return r;
    }
    template<> 
    Array2D<std::complex<float>> Array2D<std::complex<float>>::operator * (const Array2D<std::complex<float>> &b) {
        Array2D<std::complex<float>> r(M,N);
        vcMul(M*N,(MKL_Complex8*)_array,(MKL_Complex8*)b._array,(MKL_Complex8*)r._array);
        return r;
    }
    template<>
    Array2D<std::complex<double>> Array2D<std::complex<double>>::operator * (const Array2D<std::complex<double>> &b) {
        Array2D<std::complex<double>> r(M,N);        
        vzMul(M*N,(MKL_Complex16*)_array,(MKL_Complex16*)b._array,(MKL_Complex16*)r._array);
        return r;
    }  

    template<>
    Array2D<float> Array2D<float>::operator / (const Array2D<float> &b) {
        Array2D<float> r(M,N);        
        vsDiv(M*N,_array,b._array,r._array);
        return r;
    }
    template<>
    Array2D<double> Array2D<double>::operator / (const Array2D<double> &b) {
        Array2D<double> r(M,N);
        vdDiv(M*N,_array,b._array,r._array);
        return r;
    }
    template<> 
    Array2D<std::complex<float>> Array2D<std::complex<float>>::operator / (const Array2D<std::complex<float>> &b) {
        Array2D<std::complex<float>> r(M,N);
        vcDiv(M*N,(MKL_Complex8*)_array,(MKL_Complex8*)b._array,(MKL_Complex8*)r._array);
        return r;
    }
    template<>
    Array2D<std::complex<double>> Array2D<std::complex<double>>::operator / (const Array2D<std::complex<double>> &b) {
        Array2D<std::complex<double>> r(M,N);        
        vzDiv(M*N,(MKL_Complex16*)_array,(MKL_Complex16*)b._array,(MKL_Complex16*)r._array);
        return r;
    }  


    inline float sum(const Array2D<float> & a) {        
        return cblas_sasum(a.M*a.N,a._array,1);        
    }

    inline double sum(const Array2D<double> & a) {        
        return cblas_dasum(a.M*a.N,a._array,1);        
    }

    inline float sum(const Array2D<std::complex<float>> &a) {        
        return (cblas_scasum(a.M*a.N,(MKL_Complex8*)a._array,1));        
    }

    inline double sum(const Array2D<std::complex<double>> & a) {        
        return (cblas_dzasum(a.M*a.N,(MKL_Complex16*)a._array,1));    
    }  
    
    inline Array2D<float> sqr(const Array2D<float> & a) {
        Array2D<float> r(a.M,a.N);        
        vsSqr(a.M*a.N,a._array,r._array);
        return r;
    }
    
    inline Array2D<double> sqr(const Array2D<double> & a) {
        Array2D<double> r(a.M,a.N);
        vdSqr(a.M*a.N,a._array,r._array);
        return r;
    }
    
    inline Array2D<float> abs(const Array2D<float> & a) {
        Array2D<float> r(a.M,a.N);        
        vsAbs(a.M*a.N,a._array,r._array);
        return r;
    }
    
    inline Array2D<double> abs(const Array2D<double> & a) {
        Array2D<double> r(a.M,a.N);
        vdAbs(a.M*a.N,a._array,r._array);
        return r;
    }
    inline Array2D<float> abs(const Array2D<std::complex<float>> & a) {
        Array2D<float> r(a.M,a.N);        
        vcAbs(a.M*a.N,(MKL_Complex8*)a._array,r._array);
        return r;
    }
    
    inline Array2D<double> abs(const Array2D<std::complex<double>> & a) {
        Array2D<double> r(a.M,a.N);
        vzAbs(a.M*a.N,(MKL_Complex16*)a._array,r._array);
        return r;
    }
    
    
    inline  Array2D<std::complex<float>> mulByConj(const Array2D<std::complex<float>> & a, const Array2D<std::complex<float>> & b) {
        Array2D<std::complex<float>> r(a.M,a.N);
        vcMulByConj(a.M*a.N,(MKL_Complex8*)a._array,(MKL_Complex8*)b._array,(MKL_Complex8*)r._array);
        return r;
    }
    inline  Array2D<std::complex<double>> mulByConj(const Array2D<std::complex<double>> & a, const Array2D<std::complex<double>> & b) {
        Array2D<std::complex<double>> r(a.M,a.N);
        vzMulByConj(a.M*a.N,(MKL_Complex16*)a._array,(MKL_Complex16*)b._array,(MKL_Complex16*)r._array);
        return r;
    }

    inline  Array2D<std::complex<float>> conj(const Array2D<std::complex<float>> & a) {
        Array2D<std::complex<float>> r(a.M,a.N);
        vcConj(a.M*a.N,(MKL_Complex8*)a._array,(MKL_Complex8*)r._array);
        return r;
    }
    inline  Array2D<std::complex<double>> conj(const Array2D<std::complex<double>> & a) {
        Array2D<std::complex<double>> r(a.M,a.N);
        vzConj(a.M*a.N,(MKL_Complex16*)a._array,(MKL_Complex16*)r._array);
        return r;
    }

    inline  Array2D<std::complex<float>> arg(const Array2D<std::complex<float>> & a) {
        Array2D<std::complex<float>> r(a.M,a.N);
        vcConj(a.M*a.N,(MKL_Complex8*)a._array,(MKL_Complex8*)r._array);
        return r;
    }
    inline  Array2D<std::complex<double>> arg(const Array2D<std::complex<double>> & a) {
        Array2D<std::complex<double>> r(a.M,a.N);
        vzConj(a.M*a.N,(MKL_Complex16*)a._array,(MKL_Complex16*)r._array);
        return r;
    }

    inline Array2D<float> fmod(const Array2D<float> & a, const Array2D<float> & b) {
        Array2D<float> r(a.M,a.N);        
        vsFmod(a.M*a.N,a._array,b._array,r._array);
        return r;
    }

    inline Array2D<double> fmod(const Array2D<double> & a, const Array2D<double> & b) {
        Array2D<double> r(a.M,a.N);        
        vdFmod(a.M*a.N,a._array,b._array,r._array);
        return r;
    }

    inline Array2D<float> remainder(const Array2D<float> & a, const Array2D<float> & b) {
        Array2D<float> r(a.M,a.N);        
        vsRemainder(a.M*a.N,a._array,b._array,r._array);
        return r;
    }

    inline Array2D<double> remainder(const Array2D<double> & a, const Array2D<double> & b) {
        Array2D<double> r(a.N);        
        vdRemainder(a.N,a._array,b._array,r._array);
        return r;
    }

    
    inline Array2D<float> recip(const Array2D<float> & a) {
        Array2D<float> r(a.N);        
        vsInv(a.N,a._array,r._array);
        return r;
    }
    
    inline Array2D<double> recip(const Array2D<double> & a) {
        Array2D<double> r(a.N);
        vdInv(a.N,a._array,r._array);
        return r;
    }

    inline Array2D<float> sqrt(const Array2D<float> & a) {
        Array2D<float> r(a.N);        
        vsSqrt(a.N,a._array,r._array);
        return r;
    }
    
    inline Array2D<double> sqrt(const Array2D<double> & a) {
        Array2D<double> r(a.N);
        vdSqrt(a.N,a._array,r._array);
        return r;
    }
    inline Array2D<std::complex<float>> sqrt(const Array2D<std::complex<float>> & a) {
        Array2D<std::complex<float>> r(a.N);        
        vcSqrt(a.N,(MKL_Complex8*)a._array,(MKL_Complex8*)r._array);
        return r;
    }
    
    inline Array2D<std::complex<double>> sqrt(const Array2D<std::complex<double>> & a) {
        Array2D<std::complex<double>> r(a.N);
        vzSqrt(a.N,(MKL_Complex16*)a._array,(MKL_Complex16*)r._array);
        return r;
    }

    inline Array2D<float> rsqrt(const Array2D<float> & a) {
        Array2D<float> r(a.N);        
        vsInvSqrt(a.N,a._array,r._array);
        return r;
    }
    
    inline Array2D<double> rsqrt(const Array2D<double> & a) {
        Array2D<double> r(a.N);
        vdInvSqrt(a.N,a._array,r._array);
        return r;
    }
    
    inline Array2D<float> cbrt(const Array2D<float> & a) {
        Array2D<float> r(a.N);        
        vsCbrt(a.N,a._array,r._array);
        return r;
    }
    
    inline Array2D<double> cbrt(const Array2D<double> & a) {
        Array2D<double> r(a.N);
        vdCbrt(a.N,a._array,r._array);
        return r;
    }
    inline Array2D<float> rcbrt(const Array2D<float> & a) {
        Array2D<float> r(a.N);        
        vsInvCbrt(a.N,a._array,r._array);
        return r;
    }
    
    inline Array2D<double> rcbrt(const Array2D<double> & a) {
        Array2D<double> r(a.N);
        vdInvCbrt(a.N,a._array,r._array);
        return r;
    }

    inline Array2D<float> pow(const Array2D<float> & a, const Array2D<float> & b) {
        Array2D<float> r(a.N);        
        vsPow(a.N,a._array,b._array,r._array);
        return r;
    }
    inline Array2D<double> pow(const Array2D<double> & a, const Array2D<double> & b) {
        Array2D<double> r(a.N);        
        vdPow(a.N,a._array,b._array,r._array);
        return r;
    }
    inline Array2D<std::complex<float>> pow(const Array2D<std::complex<float>> & a, const Array2D<std::complex<float>> & b) {
        Array2D<std::complex<float>> r(a.N);        
        vcPow(a.N,(MKL_Complex8*)a._array,(MKL_Complex8*)b._array,(MKL_Complex8*)r._array);
        return r;
    }
    inline Array2D<std::complex<double>> pow(const Array2D<std::complex<double>> & a, const Array2D<std::complex<double>> & b) {
        Array2D<std::complex<double>> r(a.N);        
        vzPow(a.N,(MKL_Complex16*)a._array,(MKL_Complex16*)b._array,(MKL_Complex16*)r._array);
        return r;
    }

    inline Array2D<float> powx(const Array2D<float> & a, const float b) {
        Array2D<float> r(a.N);        
        vsPowx(a.N,a._array,b,r._array);
        return r;
    }
    inline Array2D<double> pow(const Array2D<double> & a, const double b) {
        Array2D<double> r(a.N);        
        vdPowx(a.N,a._array,b,r._array);
        return r;
    }
    /*
    inline Array2D<std::complex<float>> powx(const Array2D<std::complex<float>> & a, const std::complex<float> & b) {
        Array2D<std::complex<float>> r(a.N);        
        vcPowx(a.N,(MKL_Complex8*)a._array,static_cast<MKL_Complex8>(b),(MKL_Complex8*)r._array);
        return r;
    }
    inline Array2D<std::complex<double>> powx(const Array2D<std::complex<double>> & a, const std::complex<double> & b) {
        Array2D<std::complex<double>> r(a.N);        
        vzPowx(a.N,(MKL_Complex16*)a._array,(MKL_Complex16)b._array,(MKL_Complex16*)r._array);
        return r;
    }
    */

    
    inline Array2D<float> exp(const Array2D<float> & a) {
        Array2D<float> r(a.N);        
        vsExp(a.N,a._array,r._array);
        return r;
    }    
    inline Array2D<double> exp(const Array2D<double> & a) {
        Array2D<double> r(a.N);
        vdExp(a.N,a._array,r._array);
        return r;
    }
    inline Array2D<float> exp(const Array2D<std::complex<float>> & a) {
        Array2D<float> r(a.N);        
        vcExp(a.N,(MKL_Complex8*)a._array,(MKL_Complex8*)r._array);
        return r;
    }
    
    inline Array2D<double> exp(const Array2D<std::complex<double>> & a) {
        Array2D<double> r(a.N);
        vzExp(a.N,(MKL_Complex16*)a._array,(MKL_Complex16*)r._array);
        return r;
    }

    inline Array2D<float> exp2(const Array2D<float> & a) {
        Array2D<float> r(a.N);        
        vsExp2(a.N,a._array,r._array);
        return r;
    }
    inline Array2D<double> exp2(const Array2D<double> & a) {
        Array2D<double> r(a.N);
        vdExp2(a.N,a._array,r._array);
        return r;
    }
    inline Array2D<float> exp10(const Array2D<float> & a) {
        Array2D<float> r(a.N);        
        vsExp10(a.N,a._array,r._array);
        return r;
    }    
    inline Array2D<double> exp10(const Array2D<double> & a) {
        Array2D<double> r(a.N);
        vdExp10(a.N,a._array,r._array);
        return r;
    }
    inline Array2D<float> expm1(const Array2D<float> & a) {
        Array2D<float> r(a.N);        
        vsExpm1(a.N,a._array,r._array);
        return r;
    }    
    inline Array2D<double> expm1(const Array2D<double> & a) {
        Array2D<double> r(a.N);
        vdExpm1(a.N,a._array,r._array);
        return r;
    }

    inline Array2D<float> ln(const Array2D<float> & a) {
        Array2D<float> r(a.N);        
        vsLn(a.N,a._array,r._array);
        return r;
    }    
    inline Array2D<double> ln(const Array2D<double> & a) {
        Array2D<double> r(a.N);
        vdLn(a.N,a._array,r._array);
        return r;
    }
    inline Array2D<float> ln(const Array2D<std::complex<float>> & a) {
        Array2D<float> r(a.N);        
        vcLn(a.N,(MKL_Complex8*)a._array,(MKL_Complex8*)r._array);
        return r;
    }    
    inline Array2D<double> ln(const Array2D<std::complex<double>> & a) {
        Array2D<double> r(a.N);
        vzLn(a.N,(MKL_Complex16*)a._array,(MKL_Complex16*)r._array);
        return r;
    }
    inline Array2D<float> log2(const Array2D<float> & a) {
        Array2D<float> r(a.N);        
        vsLog2(a.N,a._array,r._array);
        return r;
    }    
    inline Array2D<double> log2(const Array2D<double> & a) {
        Array2D<double> r(a.N);
        vdLog2(a.N,a._array,r._array);
        return r;
    }
    inline Array2D<float> log10(const Array2D<float> & a) {
        Array2D<float> r(a.N);        
        vsLog10(a.N,a._array,r._array);
        return r;
    }    
    inline Array2D<double> log10(const Array2D<double> & a) {
        Array2D<double> r(a.N);
        vdLog10(a.N,a._array,r._array);
        return r;
    }
    inline Array2D<float> log1p(const Array2D<float> & a) {
        Array2D<float> r(a.N);        
        vsLog1p(a.N,a._array,r._array);
        return r;
    }    
    inline Array2D<double> log1p(const Array2D<double> & a) {
        Array2D<double> r(a.N);
        vdLog1p(a.N,a._array,r._array);
        return r;
    }
    inline Array2D<float> logb(const Array2D<float> & a) {
        Array2D<float> r(a.N);        
        vsLogb(a.N,a._array,r._array);
        return r;
    }    
    inline Array2D<double> logb(const Array2D<double> & a) {
        Array2D<double> r(a.N);
        vdLogb(a.N,a._array,r._array);
        return r;
    }

    inline Array2D<float> sin(const Array2D<float> & a) {
        Array2D<float> r(a.N);        
        vsSin(a.N,a._array,r._array);
        return r;
    }    
    inline Array2D<double> sin(const Array2D<double> & a) {
        Array2D<double> r(a.N);
        vdSin(a.N,a._array,r._array);
        return r;
    }
    inline Array2D<float> sin(const Array2D<std::complex<float>> & a) {
        Array2D<float> r(a.N);        
        vcSin(a.N,(MKL_Complex8*)a._array,(MKL_Complex8*)r._array);
        return r;
    }    
    inline Array2D<double> sin(const Array2D<std::complex<double>> & a) {
        Array2D<double> r(a.N);
        vzSin(a.N,(MKL_Complex16*)a._array,(MKL_Complex16*)r._array);
        return r;
    }

    inline Array2D<float> cos(const Array2D<float> & a) {
        Array2D<float> r(a.N);        
        vsCos(a.N,a._array,r._array);
        return r;
    }    
    inline Array2D<double> cos(const Array2D<double> & a) {
        Array2D<double> r(a.N);
        vdCos(a.N,a._array,r._array);
        return r;
    }
    inline Array2D<float> cos(const Array2D<std::complex<float>> & a) {
        Array2D<float> r(a.N);        
        vcCos(a.N,(MKL_Complex8*)a._array,(MKL_Complex8*)r._array);
        return r;
    }    
    inline Array2D<double> cos(const Array2D<std::complex<double>> & a) {
        Array2D<double> r(a.N);
        vzCos(a.N,(MKL_Complex16*)a._array,(MKL_Complex16*)r._array);
        return r;
    }
    inline Array2D<float> tan(const Array2D<float> & a) {
        Array2D<float> r(a.N);        
        vsTan(a.N,a._array,r._array);
        return r;
    }    
    inline Array2D<double> tan(const Array2D<double> & a) {
        Array2D<double> r(a.N);
        vdTan(a.N,a._array,r._array);
        return r;
    }
    inline Array2D<float> tan(const Array2D<std::complex<float>> & a) {
        Array2D<float> r(a.N);        
        vcTan(a.N,(MKL_Complex8*)a._array,(MKL_Complex8*)r._array);
        return r;
    }    
    inline Array2D<double> tan(const Array2D<std::complex<double>> & a) {
        Array2D<double> r(a.N);
        vzTan(a.N,(MKL_Complex16*)a._array,(MKL_Complex16*)r._array);
        return r;
    }
    inline void sincos(const Array2D<float> & a, Array2D<float> & sins, Array2D<float> & coss) {
        sins.resize(a.size());
        coss.resize(a.size());
        vsSinCos(a.N,a._array,sins._array,coss._array);    
    }    
    inline void sincos(const Array2D<double> & a, Array2D<double> & sins, Array2D<double> & coss) {
        sins.resize(a.size());
        coss.resize(a.size());
        vdSinCos(a.N,a._array,sins._array,coss._array);    
    }    
    

    inline Array2D<float> asin(const Array2D<float> & a) {
        Array2D<float> r(a.N);        
        vsAsin(a.N,a._array,r._array);
        return r;
    }    
    inline Array2D<double> asin(const Array2D<double> & a) {
        Array2D<double> r(a.N);
        vdAsin(a.N,a._array,r._array);
        return r;
    }
    inline Array2D<float> asin(const Array2D<std::complex<float>> & a) {
        Array2D<float> r(a.N);        
        vcAsin(a.N,(MKL_Complex8*)a._array,(MKL_Complex8*)r._array);
        return r;
    }    
    inline Array2D<double> asin(const Array2D<std::complex<double>> & a) {
        Array2D<double> r(a.N);
        vzAsin(a.N,(MKL_Complex16*)a._array,(MKL_Complex16*)r._array);
        return r;
    }

    inline Array2D<float> acos(const Array2D<float> & a) {
        Array2D<float> r(a.N);        
        vsAcos(a.N,a._array,r._array);
        return r;
    }    
    inline Array2D<double> acos(const Array2D<double> & a) {
        Array2D<double> r(a.N);
        vdAcos(a.N,a._array,r._array);
        return r;
    }
    inline Array2D<float> acos(const Array2D<std::complex<float>> & a) {
        Array2D<float> r(a.N);        
        vcAcos(a.N,(MKL_Complex8*)a._array,(MKL_Complex8*)r._array);
        return r;
    }    
    inline Array2D<double> acos(const Array2D<std::complex<double>> & a) {
        Array2D<double> r(a.N);
        vzAcos(a.N,(MKL_Complex16*)a._array,(MKL_Complex16*)r._array);
        return r;
    }
    inline Array2D<float> atan(const Array2D<float> & a) {
        Array2D<float> r(a.N);        
        vsAtan(a.N,a._array,r._array);
        return r;
    }    
    inline Array2D<double> atan(const Array2D<double> & a) {
        Array2D<double> r(a.N);
        vdAtan(a.N,a._array,r._array);
        return r;
    }
    inline Array2D<float> atan(const Array2D<std::complex<float>> & a) {
        Array2D<float> r(a.N);        
        vcAtan(a.N,(MKL_Complex8*)a._array,(MKL_Complex8*)r._array);
        return r;
    }    
    inline Array2D<double> atan(const Array2D<std::complex<double>> & a) {
        Array2D<double> r(a.N);
        vzAtan(a.N,(MKL_Complex16*)a._array,(MKL_Complex16*)r._array);
        return r;
    }
    inline Array2D<float> atan2(const Array2D<float> & a, const Array2D<float> & b) {
        Array2D<float> r(a.N);        
        vsAtan2(a.N,a._array,b._array, r._array);
        return r;
    }    
    inline Array2D<double> atan2(const Array2D<double> & a, const Array2D<double> & b) {
        Array2D<double> r(a.N);
        vdAtan2(a.N,a._array,b._array,r._array);
        return r;
    }

    inline Array2D<float> sinh(const Array2D<float> & a) {
        Array2D<float> r(a.N);        
        vsSinh(a.N,a._array,r._array);
        return r;
    }    
    inline Array2D<double> sinh(const Array2D<double> & a) {
        Array2D<double> r(a.N);
        vdSinh(a.N,a._array,r._array);
        return r;
    }
    inline Array2D<float> sinh(const Array2D<std::complex<float>> & a) {
        Array2D<float> r(a.N);        
        vcSinh(a.N,(MKL_Complex8*)a._array,(MKL_Complex8*)r._array);
        return r;
    }    
    inline Array2D<double> sinh(const Array2D<std::complex<double>> & a) {
        Array2D<double> r(a.N);
        vzSinh(a.N,(MKL_Complex16*)a._array,(MKL_Complex16*)r._array);
        return r;
    }

    inline Array2D<float> cosh(const Array2D<float> & a) {
        Array2D<float> r(a.N);        
        vsCosh(a.N,a._array,r._array);
        return r;
    }    
    inline Array2D<double> cosh(const Array2D<double> & a) {
        Array2D<double> r(a.N);
        vdCosh(a.N,a._array,r._array);
        return r;
    }
    inline Array2D<float> cosh(const Array2D<std::complex<float>> & a) {
        Array2D<float> r(a.N);        
        vcCosh(a.N,(MKL_Complex8*)a._array,(MKL_Complex8*)r._array);
        return r;
    }    
    inline Array2D<double> cosh(const Array2D<std::complex<double>> & a) {
        Array2D<double> r(a.N);
        vzCosh(a.N,(MKL_Complex16*)a._array,(MKL_Complex16*)r._array);
        return r;
    }
    inline Array2D<float> tanh(const Array2D<float> & a) {
        Array2D<float> r(a.N);        
        vsTanh(a.N,a._array,r._array);
        return r;
    }    
    inline Array2D<double> tanh(const Array2D<double> & a) {
        Array2D<double> r(a.N);
        vdTanh(a.N,a._array,r._array);
        return r;
    }
    inline Array2D<float> tahn(const Array2D<std::complex<float>> & a) {
        Array2D<float> r(a.N);        
        vcTanh(a.N,(MKL_Complex8*)a._array,(MKL_Complex8*)r._array);
        return r;
    }    
    inline Array2D<double> tanh(const Array2D<std::complex<double>> & a) {
        Array2D<double> r(a.N);
        vzTanh(a.N,(MKL_Complex16*)a._array,(MKL_Complex16*)r._array);
        return r;
    }

    inline Array2D<float> asinh(const Array2D<float> & a) {
        Array2D<float> r(a.N);        
        vsAsinh(a.N,a._array,r._array);
        return r;
    }    
    inline Array2D<double> asinh(const Array2D<double> & a) {
        Array2D<double> r(a.N);
        vdAsinh(a.N,a._array,r._array);
        return r;
    }
    inline Array2D<float> asinh(const Array2D<std::complex<float>> & a) {
        Array2D<float> r(a.N);        
        vcAsinh(a.N,(MKL_Complex8*)a._array,(MKL_Complex8*)r._array);
        return r;
    }    
    inline Array2D<double> asinh(const Array2D<std::complex<double>> & a) {
        Array2D<double> r(a.N);
        vzAsinh(a.N,(MKL_Complex16*)a._array,(MKL_Complex16*)r._array);
        return r;
    }

    inline Array2D<float> acosh(const Array2D<float> & a) {
        Array2D<float> r(a.N);        
        vsAcosh(a.N,a._array,r._array);
        return r;
    }    
    inline Array2D<double> acosh(const Array2D<double> & a) {
        Array2D<double> r(a.N);
        vdAcosh(a.N,a._array,r._array);
        return r;
    }
    inline Array2D<float> acosh(const Array2D<std::complex<float>> & a) {
        Array2D<float> r(a.N);        
        vcAcosh(a.N,(MKL_Complex8*)a._array,(MKL_Complex8*)r._array);
        return r;
    }    
    inline Array2D<double> acosh(const Array2D<std::complex<double>> & a) {
        Array2D<double> r(a.N);
        vzAcosh(a.N,(MKL_Complex16*)a._array,(MKL_Complex16*)r._array);
        return r;
    }
    inline Array2D<float> atanh(const Array2D<float> & a) {
        Array2D<float> r(a.N);        
        vsAtanh(a.N,a._array,r._array);
        return r;
    }    
    inline Array2D<double> atanh(const Array2D<double> & a) {
        Array2D<double> r(a.N);
        vdAtanh(a.N,a._array,r._array);
        return r;
    }
    inline Array2D<float> atanh(const Array2D<std::complex<float>> & a) {
        Array2D<float> r(a.N);        
        vcAtanh(a.N,(MKL_Complex8*)a._array,(MKL_Complex8*)r._array);
        return r;
    }    
    inline Array2D<double> atanh(const Array2D<std::complex<double>> & a) {
        Array2D<double> r(a.N);
        vzAtanh(a.N,(MKL_Complex16*)a._array,(MKL_Complex16*)r._array);
        return r;
    }
    
    inline Array2D<float> floor(const Array2D<float> & a) {
        Array2D<float> r(a.N);        
        vsFloor(a.N,a._array,r._array);
        return r;
    }    
    inline Array2D<double> floor(const Array2D<double> & a) {
        Array2D<double> r(a.N);
        vdFloor(a.N,a._array,r._array);
        return r;
    }
    inline Array2D<float> ceil(const Array2D<float> & a) {
        Array2D<float> r(a.N);        
        vsCeil(a.N,a._array,r._array);
        return r;
    }    
    inline Array2D<double> ceil(const Array2D<double> & a) {
        Array2D<double> r(a.N);
        vdCeil(a.N,a._array,r._array);
        return r;
    }
    inline Array2D<float> trunc(const Array2D<float> & a) {
        Array2D<float> r(a.N);        
        vsTrunc(a.N,a._array,r._array);
        return r;
    }    
    inline Array2D<double> trunc(const Array2D<double> & a) {
        Array2D<double> r(a.N);
        vdTrunc(a.N,a._array,r._array);
        return r;
    }
    inline Array2D<float> round(const Array2D<float> & a) {
        Array2D<float> r(a.N);        
        vsRound(a.N,a._array,r._array);
        return r;
    }    
    inline Array2D<double> round(const Array2D<double> & a) {
        Array2D<double> r(a.N);
        vdRound(a.N,a._array,r._array);
        return r;
    }

    inline Array2D<float> frac(const Array2D<float> & a) {
        Array2D<float> r(a.N);        
        vsFrac(a.N,a._array,r._array);
        return r;
    }    
    inline Array2D<double> frac(const Array2D<double> & a) {
        Array2D<double> r(a.N);
        vdFrac(a.N,a._array,r._array);
        return r;
    }


    template<typename T>
    T sum(const Array2D<T> & a) {
        T x = T(0);
        for(size_t i = 0; i < a.N; i++) x += a._array[i];
        return x;
    }
    template<typename T>
    T min(const Array2D<T> & a) {
        T min = a._array[0];
        for(size_t i = 1; i < a.N; i++) 
            if(a._array[i] < min)
                min = a._array[i];
        return min;            
    }
    template<typename T>
    T max(const Array2D<T> & a) {
        T max = a._array[0];
        for(size_t i = 1; i < a.N; i++) 
            if(a._array[i] > max)
                max = a._array[i];
        return max;            
    }
    template<typename T>
    size_t min_element(const Array2D<T> & a) {
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
    size_t max_element(const Array2D<T> & a) {
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