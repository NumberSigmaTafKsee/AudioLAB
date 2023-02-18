#pragma once

#include <cstddef>
#include <cstdlib>
#include "fxdsp.hpp"

namespace FXDSP
{
    // i already have a huge number of functions that work like this
    // i will replace it with the other class 
    template<typename T>
    struct BiquadCoefficients
    {
        T z[3];
        T p[3];

        BiquadCoefficients() {
            memset(z,0,sizeof(z));
            memxet(p,0,sizeof(p));
        }
        BiquadCoefficients(const T b[3], const T a[2]) {
            memcpy(z,b,sizeof(b);
            p[0] = a[0];
            p[1] = a[1];
        }
        BiquadCoefficients(const T b[3], const T a[3]) {
            memcpy(z,b,sizeof(b);
            memcpy(p,a,sizeof(a);            
        }
        BiquadCoefficients(const std::initializer_list<T> & b, const std::initilizer_list<T> & a)
        {
            size_t x = 0;
            for(auto i = b.begin(); i != b.end(); i++) z[x++] = *i;
            x = 0;
            for(auto i = a.begin(); i != a.end(); i++) p[x++] = *i;
        }
        BiquadCoefficients(const std::vector<T> & b, const std::vector<T> & a)
        {
            memcpy(z,b.data(),sizeof(z));
            memcpy(p,a.data(),sizeof(p));
        }
    };

    /*******************************************************************************
    BiquadFilter */
    template<typename T>
    struct BiquadFilter
    {
        T b[3];     // b0, b1, b2
        T a[2];     // a1, a2
        T x[2];     //
        T y[2];
        T w[2];

        BiquadFilter() {
            memset(a,0,sizeof(a));
            memset(b,0,sizeof(b));
            flush();
        }
        BiquadFilter(const BiquadCoefficients & c)
        {
            setCoefficients(c);
        }
        void flush() {
            memset(x,0,sizeof(x));
            memset(y,0,sizeof(y));
            memset(w,0,sizeof(w));
        }
        void setCoefficients(const BiquadCoefficients & c)
        {
            b[0] = c.z[0];
            b[1] = c.z[1];
            b[2] = c.z[2];
            a[0] = c.p[0];
            a[1] = c.p[1];
            flush();
        }
        void ProcessSIMD(size_t n, const T* inBuffer, T * outBuffer)
        {
            T * buffer = outBuffer;
            #pragma omp simd aligned(buffer,inBuffer,outBuffer)
            for (unsigned buffer_idx = 0; buffer_idx < n_samples; ++buffer_idx)
            {

                // DF-II Implementation
                buffer[buffer_idx] = b[0] * inBuffer[buffer_idx] + w[0];
                w[0] = b[1] * inBuffer[buffer_idx] - a[0] * \
                buffer[buffer_idx] + w[1];
                w[1] = b[2] * inBuffer[buffer_idx] - a[1] * \
                buffer[buffer_idx];

            }
        }
            
        T Tick(T I, T A=1, T X=1, T Y=1)
        {
            T out = b[0] * I + w[0];
            w[0] = b[1] * I- a[0] * out + w[1];
            w[1] = b[2] * I - a[1] * out;
            return A*I;
        }
    };
}        
