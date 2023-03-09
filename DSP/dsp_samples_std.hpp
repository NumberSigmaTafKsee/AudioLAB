#pragma once

#include <vector>
#include <algorithm>
#include <cmath>
#include <algorithm>
#include <functional>
#include <iostream>
#include <complex>
#include <vector>
#include <new>
#include <chrono>
#include <random>
#include <cassert>

#include "Core/core_allocator.hpp"
#define OMPSIMD

namespace Casino
{
    /* swig will not generate these :(
    template<typename T> using sample_vector = std::vector<T,Allocator::aligned_allocator<T,64>>;
    template<typename T> using sample_matrix = std::vector<sample_vector<T>>;
    
    template<typename T>
    using complex_vector = sample_vector<std::complex<T>>;

    template<typename T>
    using complex_matrix = sample_matrix<std::complex<T>>;    
    */
    template<typename T>
    struct sample_vector : public std::vector<T,Allocator::aligned_allocator<T,64>>
    {
        using base = std::vector<T,Allocator::aligned_allocator<T,64>>;
        sample_vector() = default;
        sample_vector(size_t i) : base(i) {}

        using base::operator [];        
        using base::size;
        using base::resize;
        using base::max_size;
        using base::capacity;
        using base::empty;
        using base::reserve;
        using base::shrink_to_fit;
        using base::at;
        using base::front;
        using base::back;
        using base::data;
        using base::assign;
        using base::push_back;
        using base::pop_back;
        using base::insert;
        using base::erase;
        using base::swap;
        using base::clear;
        using base::emplace;
        using base::emplace_back;

        void fill(T v) {
            for(size_t i = 0; i < size(); i++) (*this)[i] = v;
        }
        void print() {
            std::cout << "VECTOR[" << size() << "]";
            for(size_t i = 0; i < size(); i++)
                std::cout << (*this)[i] << ",";
            std::cout << std::endl;
        }
    };
    template<typename T>
    struct sample_matrix : public std::vector<std::vector<T,Allocator::aligned_allocator<T,64>>>
    {
        using base = std::vector<std::vector<T,Allocator::aligned_allocator<T,64>>>;
        sample_matrix() = default;
        sample_matrix(size_t i, size_t j) : base(i) {
            for(size_t x = 0; x < i; x++)
                (*this)[i]->resize(j);
        }
    };
    template<typename T>
    struct complex_vector : public std::vector<std::complex<T>,Allocator::aligned_allocator<std::complex<T>,64>>
    {
        using base = std::vector<std::complex<T>,Allocator::aligned_allocator<std::complex<T>,64>>;
        complex_vector() = default;
        complex_vector(size_t i) : base(i) {}
    };
    template<typename T>
    struct complex_matrix : public std::vector<std::vector<std::complex<T>,Allocator::aligned_allocator<std::complex<T>,64>>>
    {
        using base = std::vector<std::vector<std::complex<T>,Allocator::aligned_allocator<std::complex<T>,64>>>;
        complex_matrix() = default;
        complex_matrix(size_t i, size_t j) : base(i) {
            for(size_t x = 0; x < i; x++)
                (*this)[i]->resize(j);
        }
    };

    template<typename T>
    struct StereoVector
    {
        sample_vector<T> samples[2];

        StereoVector(size_t n) {
            for(size_t i = 0; i < n; i++) samples[i].resize(n);
        }

        sample_vector<T>& operator[](size_t channel) { return samples[channel]; }

        sample_vector<T> get_left_channel() { return samples[0]; }
        sample_vector<T> get_right_channel() { return samples[1]; }

        void set_left_channel(sample_vector<T> & v) { samples[0] = v; }
        void set_right_channel(sample_vector<T> & v) { samples[1] = v; }
    };

    template<typename T>
    std::ostream& operator << (std::ostream & o, const sample_matrix<T> & m )
    {
        for(size_t i = 0; i < m.rows(); i++)
        {
            for(size_t j = 0; j < m.cols(); j++)
                o << m(i,j) << ",";
            o << std::endl;
        }
        return o;
    }

    template<typename T>
    T get_stride(size_t ch, size_t num_channels, size_t pos, sample_vector<T> & samples)
    {
        return samples[pos*num_channels + ch];
    }
        
    template<typename T>
    void set_stride(size_t ch, size_t num_channels, size_t pos, sample_vector<T> & samples, T sample)
    {
        samples[pos*num_channels + ch] = sample;
    }

    template<typename T>
    sample_vector<T> get_left_channel(const sample_vector<T> & in) {
        sample_vector<T> r(in.size()/2);
        size_t x = 0;
        #pragma omp simd
        for(size_t i = 0; i < in.size(); i+=2) r[x++] = in[i];
        return r;
    }
    template<typename T>
    sample_vector<T> get_right_channel(const sample_vector<T> & in) {
        sample_vector<T> r(in.size()/2);
        size_t x = 0;
        #pragma omp simd
        for(size_t i = 1; i < in.size(); i+=2) r[x++] = in[i];
        return r;
    }
    template<typename T>
    sample_vector<T> get_channel(size_t ch, const sample_vector<T> & in) {
        sample_vector<T> r(in.size()/2);
        size_t x = 0;
        #pragma omp simd
        for(size_t i = ch; i < in.size(); i+=2) r[x++] = in[i];
        return r;
    }
    template<typename T>
    void set_left_channel(const sample_vector<T> & left, sample_vector<T> & out) {
        size_t x = 0;
        #pragma omp simd
        for(size_t i = 0; i < out.size(); i+=2) out[i] = left[x++];
    }
    template<typename T>
    void set_right_channel(const sample_vector<T> & right, sample_vector<T> & out) {
        size_t x = 0;
        #pragma omp simd
        for(size_t i = 1; i < out.size(); i+=2) out[i] = right[x++];
    }
    template<typename T>
    void set_channel(size_t ch, const sample_vector<T> & in, sample_vector<T> & out) {
        size_t x = 0;
        #pragma omp simd
        for(size_t i = ch; i < out.size(); i+=2) out[i] = in[x++];
    }
    template<typename T>
    sample_vector<T> interleave(size_t n, size_t channels, const sample_vector<sample_vector<T>> & in) {
        sample_vector<T> r(n*channels);
        #pragma omp simd
        for(size_t i = 0; i < channels; i++)
            for(size_t j = 0; j < n; j++)
                r[j*channels + i] = in[i][j];
        return r;
    }
    template<typename T>
    sample_vector<T> interleave(size_t n, size_t channels, const sample_vector<T*> & in) {
        sample_vector<T> r(n*channels);
        #pragma omp simd
        for(size_t i = 0; i < channels; i++)
            for(size_t j = 0; j < n; j++)
                r[j*channels + i] = in[i][j];
        return r;
    }
    template<typename T>
    sample_vector<sample_vector<T>> deinterleave(size_t n, size_t channels, const sample_vector<T> & in) {
        sample_vector<sample_vector<T>> r(n);
        #pragma omp simd
        for(size_t i = 0; i < channels; i++)
        {
            r[i].resize(n);            
            for(size_t j = 0; j < n; j++)
                r[i][j] = in[j*channels + i];
        }
        return r;
    }

    template<typename T>
    T get_stride(size_t ch, size_t num_channels, size_t pos, T * samples)
    {
        return samples[pos*num_channels + ch];
    }
    template<typename T>
    void set_stride(size_t ch, size_t num_channels, size_t pos, T * samples, T sample)
    {
        samples[pos*num_channels + ch] = sample;
    }

    template<typename T>
    sample_vector<T> get_left_channel(size_t n, const T* in) {
        sample_vector<T> r(n/2);
        size_t x = 0;
        #pragma omp simd
        for(size_t i = 0; i < n; i+=2) r[x++] = in[i];
        return r;
    }
    template<typename T>
    sample_vector<T> get_right_channel(size_t n, const T* & in) {
        sample_vector<T> r(n/2);
        size_t x = 0;
        #pragma omp simd
        for(size_t i = 1; i < n; i+=2) r[x++] = in[i];
        return r;
    }
    template<typename T>
    sample_vector<T> get_channel(size_t ch, size_t n, T* in) {
        sample_vector<T> r(n/2);
        size_t x = 0;
        #pragma omp simd
        for(size_t i = ch; i < n; i+=2) r[x++] = in[i];
        return r;
    }

    template<typename T>
    void set_left_channel(size_t n, const T* left, T* out) {
        size_t x = 0;
        #pragma omp simd
        for(size_t i = 0; i < n; i+=2) out[i] = left[x++];
    }
    template<typename T>
    void set_right_channel(size_t n, const T* right, T* out) {
        size_t x = 0;
        #pragma omp simd
        for(size_t i = 1; i < n; i+=2) out[i] = right[x++];
    }
    template<typename T>
    void set_channel(size_t ch, size_t n, const T* in, T* out) {
        size_t x = 0;
        #pragma omp simd
        for(size_t i = ch; i < n; i+=2) out[i] = in[x++];
    }

    template<typename T>
    void fill_left_channel(size_t n, const T v, T* out) {
        size_t x = 0;
        #pragma omp simd
        for(size_t i = 0; i < n; i+=2) out[i] = v;
    }
    template<typename T>
    void fill_right_channel(size_t n, const T v, T* out) {
        size_t x = 0;
        #pragma omp simd
        for(size_t i = 1; i < n; i+=2) out[i] = v;
    }
    template<typename T>
    void fill_channel(size_t ch, size_t n, const T v, T* out) {
        size_t x = 0;
        #pragma omp simd
        for(size_t i = ch; i < n; i+=2) out[i] = v;
    }

    template<typename T>
    void fill_left_channel(sample_vector<T> & out, const T v) {
        size_t x = 0;
        #pragma omp simd
        for(size_t i = 0; i < out.size(); i+=2) out[i] = v;
    }
    template<typename T>
    void fill_right_channel(sample_vector<T> & out, const T v) {
        size_t x = 0;
        #pragma omp simd
        for(size_t i = 1; i < out.size(); i+=2) out[i] = v;
    }
    template<typename T>
    void fill_channel(size_t ch, sample_vector<T> & out, const T v) {
        size_t x = 0;
        #pragma omp simd
        for(size_t i = ch; i < out.size(); i+=2) out[i] = v;
    }


    template<typename T>
    sample_vector<T> interleave(size_t n, size_t channels, const T** & in) {
        sample_vector<T> r(n*channels);
        #pragma omp simd
        for(size_t i = 0; i < channels; i++)
            for(size_t j = 0; j < n; j++)
                r[j*channels + i] = in[i][j];
        return r;
    }
    template<typename T>
    sample_vector<sample_vector<T>> deinterleave(size_t n, size_t channels, const T* & in) {
        sample_vector<sample_vector<T>> r(n);
        #pragma omp simd
        for(size_t i = 0; i < channels; i++)
        {
            r[i].resize(n);
            for(size_t j = 0; j < n; j++)
                r[i][j] = in[j*channels + i];
        }
        return r;
    }

    template<typename T>
    bool equal_vector (sample_vector<T> & a, sample_vector<T> & b) {        
        return std::equal(a.begin(),a.end(),b.begin());
    }

    template<typename T>
    void copy_vector(sample_vector<T> & dst, sample_vector<T> & src) {        
        std::copy(src.begin(),src.end(),dst.begin());
    }
    template<typename T>
    void copy_vector(sample_vector<T> & dst, size_t n, T * src) {        
        std::copy(&src[0],&src[n-1],dst.begin());
    }
    template<typename T>
    sample_vector<T> slice_vector(size_t start, size_t end, sample_vector<T> & src) {
        sample_vector<T> r(end-start);        
        std::copy(src.begin()+start,src.begin()+end,r.begin());
        return r;
    }

    template<typename T>
    void copy_buffer(size_t n, T * dst, T * src) {                        
        memcpy(dst,src,n*sizeof(T));
    }

    template<typename T>
    sample_vector<T> slice_buffer(size_t start, size_t end, T * ptr) {
        sample_vector<T> r(end-start);
        std::copy(&ptr[start],&ptr[end],r.begin());
        return r;
    }

    template<typename T>
    void split_stereo(size_t n, const T* input, T * left, T * right)
    {
        size_t x=0;
        for(size_t i = 0; i < n; i+=2)
        {
            left[x] = input[i];
            right[x++] = input[i+1];
        }
    }

    template<typename T>
    void split_stereo(const sample_vector<T> & input, sample_vector<T> & left, sample_vector<T> & right) {
        size_t x = input.size();
        left.resize(x/2);
        right.resize(x/2);
        split_stereo(x,input.data(),left.data(),right.data());
    }

    template<typename T>
    void swap(sample_vector<T> & left, sample_vector<T> & right) {
        std::swap(left,right);
    }

    template<typename T>
    bool contains(const sample_vector<T> & v, const T val) {
        return std::find(v.begin(),v.end(),val) != v.end();
    }

    template<typename T>
    sample_vector<T> mix(const sample_vector<T> & a, const sample_vector<T> & b, T m = 0.5)
    {
        assert(a.size() == b.size());
        sample_vector<T> r(a.size());
        T max = -99999;
        #pragma omp simd
        for(size_t i = 0; i < r.size(); i++) 
        {            
            r[i] = m*(a[i]+b[i]);
        }                
        return r;
    }
    template<typename T>
    sample_vector<T> normalize(const sample_vector<T> & a) {
        sample_vector<T> r(a);        
        auto max = std::max_element(r.begin(),r.end());
        if(*max > 0) 
            #pragma omp simd
            for(size_t i = 0; i < r.size(); i++) r[i] /= *max;
        return r;
    }
    template<class A, class B>
    sample_vector<B> convert(const sample_vector<A> & v) {
        sample_vector<B> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = B(v[i]);
        return r;
    }
    template<class T>
    sample_vector<T> kernel(const sample_vector<T> & v, T (*f)(T value)) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = f(v[i]);
        return r;
    }
    template<class T>
    sample_vector<T> kernel(const sample_vector<T> & v, std::function<T (T)> func) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = func(v[i]);
        return r;
    }
    template<class T>
    void inplace_add(const sample_vector<T> & a, sample_vector<T> & r, std::function<T (T)> func) {        
        #pragma omp simd
        for(size_t i = 0; i < a.size(); i++) r[i] += func(a[i]);        
    }
    template<class T>
    void inplace_sub(const sample_vector<T> & a, sample_vector<T> & r, std::function<T (T)> func) {        
        #pragma omp simd
        for(size_t i = 0; i < a.size(); i++) r[i] -= func(a[i]);        
    }
    template<class T>
    void inplace_mul(const sample_vector<T> & a, sample_vector<T> & r, std::function<T (T)> func) {        
        #pragma omp simd
        for(size_t i = 0; i < a.size(); i++) r[i] *= func(a[i]);        
    }
    template<class T>
    void inplace_div(const sample_vector<T> & a, sample_vector<T> & r, std::function<T (T)> func) {        
        #pragma omp simd
        for(size_t i = 0; i < a.size(); i++) r[i] /= func(a[i]);
    }


    template<class T>
    void fill(sample_vector<T> & in, T x)
    {
        std::fill(in.begin(),in.end(),x);        
    }
    template<class T>
    void zeros(sample_vector<T> & in)
    {
        fill(in,T(0));
    }
    template<class T>
    void ones(sample_vector<T> & in)
    {
        fill(in,T(1));
    }
    

    template<class T>
    sample_vector<T> cos(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::cos(v[i]);
        return r;
    }
    template<class T>
    sample_vector<T> sin(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::sin(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> tan(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::tan(v[i]);
        return r;
    }

    template<class T>
    sample_vector<T> acos(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::acos(v[i]);
        return r;
    }
    template<class T>
    sample_vector<T> asin(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::asin(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> atan(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::atan(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> atan2(const sample_vector<T> & v, const T value) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::atan2(v[i], value);
        return r;
    }    
    template<class T>
    sample_vector<T> atan2(const sample_vector<T> & v, const sample_vector<T> value) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::atan2(v[i], value[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> cosh(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::cosh(v[i]);
        return r;
    }
    template<class T>
    sample_vector<T> sinh(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::sinh(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> tanh(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::tanh(v[i]);
        return r;
    }

    template<class T>
    sample_vector<T> acosh(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::acosh(v[i]);
        return r;
    }
    template<class T>
    sample_vector<T> asinh(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::asinh(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> atanh(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::atanh(v[i]);
        return r;
    }    

    template<class T>
    sample_vector<T> exp(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::exp(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> log(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::log(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> log10(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::log10(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> exp2(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::exp2(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> expm1(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::expm1(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> ilogb(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::ilogb(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> log2(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::log2(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> log1p(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::log1p(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> logb(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::logb(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> scalbn(const sample_vector<T> & v, const sample_vector<int> & x) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::scalbn(v[i],x[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> scalbn(const sample_vector<T> & v, const int x) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::scalbn(v[i],x);
        return r;
    }    
    template<class T>
    sample_vector<T> scalbln(const sample_vector<T> & v, const sample_vector<long int> & x) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::scalbln(v[i],x[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> scalbln(const sample_vector<T> & v, const long int x) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::scalbln(v[i],x);
        return r;
    }    
    template<class T>
    sample_vector<T> pow(const sample_vector<T> & v, const sample_vector<T> & x) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::pow(v[i],x[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> pow(const sample_vector<T> & v, const T x) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::pow(v[i],x);
        return r;
    }    
    template<class T>
    sample_vector<T> pow(const T x, const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::pow(x,v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> sqrt(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::sqrt(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> cbrt(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::cbrt(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> hypot(const sample_vector<T> & v, const sample_vector<T> & x) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::hypot(v[i],x[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> hypot(const sample_vector<T> & v, const T x) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::hypot(v[i],x);
        return r;
    }    
    template<class T>
    sample_vector<T> hypot(const T x, const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::hypot(x,v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> erf(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::erf(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> erfc(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::erfc(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> tgamma(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::tgamma(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> lgamma(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::lgamma(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> ceil(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::ceil(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> floor(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::floor(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> fmod(const sample_vector<T> & v, const sample_vector<T> & x) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::fmod(v[i],x[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> fmod(const sample_vector<T> & v, const T x) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::fmod(v[i],x);
        return r;
    }    
    template<class T>
    sample_vector<T> fmod(const T x, const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::fmod(x,v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> trunc(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::trunc(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> round(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::round(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<long int> lround(const sample_vector<T> & v) {
        sample_vector<long int> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::lround(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<long long int> llround(const sample_vector<T> & v) {
        sample_vector<long long int> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::llround(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> nearbyint(const sample_vector<T> & v) {
        sample_vector<long long int> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::nearbyint(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> remainder(const sample_vector<T> & a, const sample_vector<T> & b) {
        sample_vector<long long int> r(a.size());
        #pragma omp simd
        for(size_t i = 0; i < a.size(); i++) r[i] = std::remainder(a[i],b[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> copysign(const sample_vector<T> & a, const sample_vector<T> & b) {
        sample_vector<long long int> r(a.size());
        #pragma omp simd
        for(size_t i = 0; i < a.size(); i++) r[i] = std::copysign(a[i],b[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> fdim(const sample_vector<T> & a, const sample_vector<T> & b) {
        sample_vector<long long int> r(a.size());
        #pragma omp simd
        for(size_t i = 0; i < a.size(); i++) r[i] = std::fdim(a[i],b[i]);
        return r;
    }    
    #undef fmax
    template<class T>
    sample_vector<T> fmax(const sample_vector<T> & a, const sample_vector<T> & b) {
        sample_vector<long long int> r(a.size());
        #pragma omp simd
        for(size_t i = 0; i < a.size(); i++) r[i] = std::fmax(a[i],b[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> fmin(const sample_vector<T> & a, const sample_vector<T> & b) {
        sample_vector<long long int> r(a.size());
        #pragma omp simd
        for(size_t i = 0; i < a.size(); i++) r[i] = std::fmin(a[i],b[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> fma(const sample_vector<T> & a, const sample_vector<T> & b, const sample_vector<T> & c) {
        sample_vector<long long int> r(a.size());
        #pragma omp simd
        for(size_t i = 0; i < a.size(); i++) r[i] = std::fma(a[i],b[i],c[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> fabs(const sample_vector<T> & a) {
        sample_vector<T> r(a.size());
        #pragma omp simd
        for(size_t i = 0; i < a.size(); i++) r[i] = std::fabs(a[i]);
        return r;
    }    

    template<class T>
    sample_vector<T> operator +(const sample_vector<T> & a, const sample_vector<T> & b) {
        sample_vector<T> r(a.size());
        #pragma omp simd
        for(size_t i = 0; i < a.size(); i++) r[i] = a[i] + b[i];
        return r;
    }   
    template<class T>
    sample_vector<T> operator -(const sample_vector<T> & a, const sample_vector<T> & b) {
        sample_vector<T> r(a.size());
        #pragma omp simd
        for(size_t i = 0; i < a.size(); i++) r[i] = a[i] - b[i];
        return r;
    }   
    template<class T>
    sample_vector<T> operator *(const sample_vector<T> & a, const sample_vector<T> & b) {
        sample_vector<T> r(a.size());
        #pragma omp simd
        for(size_t i = 0; i < a.size(); i++) r[i] = a[i] * b[i];
        return r;
    }   
    template<class T>
    sample_vector<T> operator /(const sample_vector<T> & a, const sample_vector<T> & b) {
        sample_vector<T> r(a.size());
        #pragma omp simd
        for(size_t i = 0; i < a.size(); i++) r[i] = a[i] / b[i];
        return r;
    }   
    template<class T>
    sample_vector<T> operator %(const sample_vector<T> & a, const sample_vector<T> & b) {
        sample_vector<T> r(a.size());
        #pragma omp simd
        for(size_t i = 0; i < a.size(); i++) r[i] = std::fmod(a[i],b[i]);
        return r;
    }   
    template<class T>
    sample_vector<T> operator +(const sample_vector<T> & a, const T& b) {
        sample_vector<T> r(a.size());
        #pragma omp simd
        for(size_t i = 0; i < a.size(); i++) r[i] = a[i] + b;
        return r;
    }   
    template<class T>
    sample_vector<T> operator -(const sample_vector<T> & a, const T& b) {
        sample_vector<T> r(a.size());
        #pragma omp simd
        for(size_t i = 0; i < a.size(); i++) r[i] = a[i] - b;
        return r;
    }   
    template<class T>
    sample_vector<T> operator *(const sample_vector<T> & a, const T& b) {
        sample_vector<T> r(a.size());
        #pragma omp simd
        for(size_t i = 0; i < a.size(); i++) r[i] = a[i] * b;
        return r;
    }   
    template<class T>
    sample_vector<T> operator / (const sample_vector<T> & a, const T& b) {
        sample_vector<T> r(a.size());
        #pragma omp simd
        for(size_t i = 0; i < a.size(); i++) r[i] = a[i] / b;
        return r;
    }   
    template<class T>
    sample_vector<T> operator %(const sample_vector<T> & a, const T& b) {
        sample_vector<T> r(a.size());
        #pragma omp simd
        for(size_t i = 0; i < a.size(); i++) r[i] = std::fmod(a[i],b);
        return r;
    }   
    template<class T>
    sample_vector<T> operator +(const T & a, const sample_vector<T>& b) {
        sample_vector<T> r(a.size());
        #pragma omp simd
        for(size_t i = 0; i < a.size(); i++) r[i] = a + b[i];
        return r;
    }   
    template<class T>
    sample_vector<T> operator -(const T & a, const sample_vector<T>& b) {
        sample_vector<T> r(a.size());
        #pragma omp simd
        for(size_t i = 0; i < a.size(); i++) r[i] = a - b[i];
        return r;
    }   
    template<class T>
    sample_vector<T> operator *(const T & a, const sample_vector<T>& b) {
        sample_vector<T> r(a.size());
        #pragma omp simd
        for(size_t i = 0; i < a.size(); i++) r[i] = a * b[i];
        return r;
    }   
    template<class T>
    sample_vector<T> operator /(const T & a, const sample_vector<T>& b) {
        sample_vector<T> r(a.size());
        #pragma omp simd
        for(size_t i = 0; i < a.size(); i++) r[i] = a / b[i];
        return r;
    }   

    template<typename T>
    sample_vector<T> linspace(int start, int end)
    {        
        T delta = start >= end? -1.0 : 1.0;
        int m = std::floor((end-start)/delta);
        if(start > end) m++;
        sample_vector<T> r(m);
        #pragma omp simd
        for(size_t i = 0; i < m; i++) r[i] = start + i*delta;
        return r;
    }

    template<typename T>
    sample_vector<T> regspace(int start, int end, T delta=(T)0.0)
    {
        if(delta==0.0) delta = start <= end? 1.0 : -1.0;                
        int m = std::floor((end-start)/delta);
        if(start > end) m++;
        sample_vector<T> r(m);
        #pragma omp simd
        for(size_t i = 0; i < m; i++) r[i] = start + i*delta;
        return r;
    }

    // convolution
    // will be enormously slow         
    // fft is much faster
    // useful for very small vector only

    template<typename T>
    sample_vector<T> linear_convolution(T * h, T * x, int M, int N)
    {
        sample_vector<T> y(M+N-1);        
            
        for(int i=0;i<M+N-1;i++)
        {  
            y[i]=0;                 
            for(int j=0;j<=i;j++)
                y[i]+=x[j]*h[i-j];
        }
        return y;
    }
    

    // linear or acyclic
    template<typename T>
    sample_vector<T> linear_convolution(sample_vector<T> & h, sample_vector<T> & x)
    {
        int M=h.size();
        int N=x.size();
        sample_vector<T> y(M+N-1);        
            
                          
        for(int i=0;i<M+N-1;i++)
        {  
            y[i]=0;                    
            for(int j=0;j<=i;j++)
                y[i]+=x[j]*h[i-j];
        }
        return y;
    }

        
    template<typename T>
    sample_vector<T> circular_convolution(T * h, T * x, int M, int N)
    {        
        sample_vector<T> y(M);                
        for(int n=0;n < M;n++)
        { 
            y[n]=0;                                  
            for(int k=0;k < N;k++)                
            {
                int j = (k>n)? n-k+N : n-k;
                y[n]=y[n]+x[k]*h[j];
            }
        }    
        return y;    
    }   
    template<typename T>
    sample_vector<T> circular_convolution(sample_vector<T>& h, sample_vector<T>& x)
    {
        int M = h.size();
        int N = x.size();
        sample_vector<T> y(M);    
                        
        for(int n=0;n < M;n++)
        { 
            y[n]=0;                                  
            for(int k=0;k < N;k++)                
            {
                int j = (k>n)? n-k+N : n-k;
                y[n]=y[n]+x[k]*h[j];
            }
        }    
        return y;    
    }   

    
    template<typename T>
    void cshiftr(sample_vector<T> & v)
    {
        T t = v[0];
        #pragma omp simd
        for(size_t i = 1; i < v.size(); i++)
            v[i-1] = v[i];
        v[v.size()-1] = t;
    }
    template<typename T>
    void cshiftl(sample_vector<T> & v)
    {
        T t = v[v.size()-1];
        #pragma omp simd
        for(size_t i = v.size()-1; i > 0; i--)
            v[i] = v[i-1];
        v[0] = t;
    }
    template<typename T>
    void cshiftleft(sample_vector<T> & v, size_t steps)
    {        
        for(size_t i = 0; i < steps; i++) cshiftl(v);
    }
    template<typename T>
    void cshiftright(sample_vector<T> & v, size_t steps)
    {
        for(size_t i = 0; i < steps; i++) cshiftr(v);
    }    
    template<typename T>
    void cshift(sample_vector<T> & v, int steps)
    {
        if(steps > 0) cshiftleft(v,steps);
        else if(steps < 0) cshiftright(v,abs(steps));
    }

    #define SQR(x) ((x)*(x))

    /** Window types */
    typedef enum
    {
    win_ones,
    win_rectangle,
    win_hamming,
    win_hanning,
    win_hanningz,
    win_blackman,
    win_blackman_harris,
    win_gaussian,
    win_welch,
    win_parzen,
    win_default = win_hanningz,
    } window_type;

    template<typename T>
    T unwrap2pi (T phase)
    {
        /* mod(phase+pi,-2pi)+pi */
        return phase + 2.0*M_PI * (1. + std::floor(-(phase + M_PI) / 2.0*M_PI));
    }


    template<typename T>
    T mean (const sample_vector<T> & s)
    {
        T tmp = 0.0;
        #if defined(HAVE_INTEL_IPP)
            ippsMean(s.data(), (int)s.size(), &tmp);
            return tmp;
        #elif defined(HAVE_ACCELERATE)
            vDSP_meanv(s.data(), 1, &tmp, s.size());
            return tmp;
        #else
            size_t j;
            #pragma omp simd
            for (j = 0; j < s.size(); j++) {
                tmp += s[j];
            }
            return tmp / (T)(s.size());
        #endif
    }

    template<typename T>
    T sum (const sample_vector<T>& s)
    {
        T tmp = 0.0;
        #if defined(HAVE_INTEL_IPP)
            ippsSum(s.data(), (int)s.size(), &tmp);
        #elif defined(HAVE_ACCELERATE)
            vDSP_sve(s.data(), 1, &tmp, s.size());
        #else
            size_t j;
            #pragma omp simd
            for (j = 0; j < s.size(); j++) {
                tmp += s[j];
            }
        #endif
        return tmp;
    }

    template<typename T>
    T max (const sample_vector<T>& s)
    {
        #if defined(HAVE_INTEL_IPP)
            T tmp = 0.;
            ippsMax( s.data(), (int)s.size(), &tmp);
        #elif defined(HAVE_ACCELERATE)
            T tmp = 0.;
            vDSP_maxv( s.data(), 1, &tmp, s.size() );
        #else
            size_t j;
            T tmp = s.data()[0];
            for (j = 1; j < s.size(); j++) {
                tmp = (tmp > s[j]) ? tmp : s[j];
            }
        #endif
        return tmp;
    }

    template<typename T>
    T min (const sample_vector<T>& s)
    {
        #if defined(HAVE_INTEL_IPP)
            T tmp = 0.;
            ippsMin(s.data(), (int)s.size(), &tmp);
        #elif defined(HAVE_ACCELERATE)
            T tmp = 0.;
            vDSP_minv(s.data(), 1, &tmp, s.size());
        #else
            size_t j;
            T tmp = s[0];
            #pragma omp simd
            for (j = 1; j < s.size(); j++) {
                tmp = (tmp < s[j]) ? tmp : s[j];
            }
        #endif
        return tmp;
    }

    /* use std::min_element
    template<typename T>
    size_t min_elem (const sample_vector<T>& s)
    {
        #ifndef HAVE_ACCELERATE
            size_t j, pos = 0.;
            T tmp = s[0];
            #pragma omp simd
            for (j = 0; j < s.size(); j++) {
                pos = (tmp < s[j]) ? pos : j;
                tmp = (tmp < s[j]) ? tmp : s[j];
            }
        #else
            T tmp = 0.;
            vDSP_Length pos = 0;
            vDSP_minvi(s, 1, &tmp, &pos, s.size());
        #endif
        return (size_t)pos;
    }

    template<typename T>
    size_t max_elem (const sample_vector<T>& s)
    {
    #ifndef HAVE_ACCELERATE
    size_t j, pos = 0;
    T tmp = 0.0;
    #pragma omp simd
    for (j = 0; j < s.size(); j++) {
        pos = (tmp > s[j]) ? pos : j;
        tmp = (tmp > s[j]) ? tmp : s[j];
    }
    #else
    T tmp = 0.;
    vDSP_Length pos = 0;
    vDSP_maxvi(s.data(), 1, &tmp, &pos, s.size());
    #endif
    return (size_t)pos;
    }
    */
    template<typename T>
    void swap(T & a, T & b) {
        T x = a;
        a = b;
        b = x;
    }
    template<typename T>
    void shift (sample_vector<T>& s)
    {
        size_t half = s.size() / 2, start = half, j;
        // if length is odd, middle element is moved to the end
        if (2 * half < s.size()) start ++;
        #ifndef HAVE_BLAS
            for (j = 0; j < half; j++) {
                swap(s[j], s[j + start]);
            }
        #else
            cblas_swap(half, s.data(), 1, s.data() + start, 1);
        #endif
        if (start != half) {
            for (j = 0; j < half; j++) {
            swap(s[j + start - 1], s[j + start]);
            }
        }
    }

    template<typename T>
    void ishift (const sample_vector<T> & s)
    {
        size_t half = s.size() / 2, start = half, j;
        // if length is odd, middle element is moved to the beginning
        if (2 * half < s.size()) start ++;
        #ifndef HAVE_BLAS
            #pragma omp simd
            for (j = 0; j < half; j++) {
                swap(s[j], s[j + start]);
            }
        #else
            cblas_swap(half, s->data, 1, s->data + start, 1);
        #endif            
            if (start != half) {
                #pragma omp simd
                for (j = 0; j < half; j++) {
                    swap(s[half], s[j]);
                }
            }
        }

    
    template<typename T>
    void clamp(const sample_vector<T>& in, T absmax) {
        size_t i;
        absmax = fabs(absmax);
        #pragma omp simd  
        for (i = 0; i < in->length; i++) in[i] = std::clamp(in[i],-absmax,absmax);  
    }


    template<typename T>
    T level_lin (const sample_vector<T>& f)
    {
        T energy = 0.;
        #ifndef HAVE_BLAS
            size_t j;
            #pragma omp simd
            for (j = 0; j < f.size(); j++) {
                energy += SQR (f[j]);
            }
        #else
            energy = cblas_dot(f.size(), f.data(), 1, f.data(), 1);
        #endif
        return energy / f.size();
    }

    template<typename T>
    T local_hfc (const sample_vector<T>& v)
    {
        T hfc = 0.;
        size_t j;
        #pragma omp simd
        for (j = 0; j < v->length; j++) {
            hfc += (j + 1) * v[j];
        }
        return hfc;
        }

    template<typename T>
    void min_removal (const sample_vector<T>& v)
    {
        T v_min = min(v);
        v += -v_min;  
    }

    template<typename T>
    T alpha_norm (const sample_vector<T>& o, T alpha)
    {
        size_t j;
        T tmp = 0.;
        #pragma omp simd
        for (j = 0; j < o.size(); j++) {
            tmp += std::pow(std::fabs(o[j]), alpha);
        }
        return std::pow(tmp / o.size(), 1. / alpha);
    }

    template<typename T>
    void alpha_normalise (const sample_vector<T>& o, T alpha)
    {
        size_t j;
        T norm = alpha_norm (o, alpha);
        o /= norm;  
    }


    template<typename T>
    T median (const sample_vector<T>& input) {
        size_t n = input->length;
        T * arr = input.data();
        size_t low, high ;
        size_t median;
        size_t middle, ll, hh;

        low = 0 ; high = n-1 ; median = (low + high) / 2;
        
        for (;;) {
            if (high <= low) /* One element only */
            return arr[median] ;

            if (high == low + 1) {  /* Two elements only */
            if (arr[low] > arr[high])
                std::swap(arr[low], arr[high]) ;
            return arr[median] ;
            }

            /* Find median of low, middle and high items; swap into position low */
            middle = (low + high) / 2;
            if (arr[middle] > arr[high])    std::swap(arr[middle], arr[high]);
            if (arr[low]    > arr[high])    std::swap(arr[low],    arr[high]);
            if (arr[middle] > arr[low])     std::swap(arr[middle], arr[low]) ;

            /* Swap low item (now in position middle) into position (low+1) */
            std::swap(arr[middle], arr[low+1]) ;

            /* Nibble from each end towards middle, swapping items when stuck */
            ll = low + 1;
            hh = high;
            for (;;) {
            do ll++; while (arr[low] > arr[ll]) ;
            do hh--; while (arr[hh]  > arr[low]) ;

            if (hh < ll)
                break;

            std::swap(arr[ll], arr[hh]) ;
            }

            /* Swap middle item (in position low) back into correct position */
            std::swap(arr[low], arr[hh]) ;

            /* Re-set active partition */
            if (hh <= median)
            low = ll;
            if (hh >= median)
            high = hh - 1;
        }
    }


    template<typename T>
    T moving_thres (const sample_vector<T>& vec, const sample_vector<T>& tmpvec,
        size_t post, size_t pre, size_t pos)
    {
        size_t k;
        T *medar = (T *) tmpvec;
        size_t win_length = post + pre + 1;
        size_t length = vec.size();
        /* post part of the buffer does not exist */
        if (pos < post + 1) {
            for (k = 0; k < post + 1 - pos; k++)
            medar[k] = 0.;            /* 0-padding at the beginning */
            #pragma omp simd
            for (k = post + 1 - pos; k < win_length; k++)
            medar[k] = vec[k + pos - post];
            /* the buffer is fully defined */
        } else if (pos + pre < length) {
            #pragma omp simd
            for (k = 0; k < win_length; k++)
            medar[k] = vec[k + pos - post];
            /* pre part of the buffer does not exist */
        } else {
            #pragma omp simd
            for (k = 0; k < length - pos + post; k++)
            medar[k] = vec[k + pos - post];
            #pragma omp simd
            for (k = length - pos + post; k < win_length; k++)
            medar[k] = 0.;            /* 0-padding at the end */
        }
        return median (tmpvec);
    }

    template<typename T>
    void sample_vector_adapt_thres(const sample_vector<T>& vec, const sample_vector<T>& tmp,
        size_t post, size_t pre) {
        size_t length = vec.size(), j;  
        for (j=0;j<length;j++) {
            vec[j] -= moving_thres(vec, tmp, post, pre, j);
        }
    }

    template<typename T>
    T quadratic_peak_pos (const sample_vector<T>& x, size_t pos) {
        T s0, s1, s2; 
        size_t x0, x2;
        T half = .5, two = 2.;
        if (pos == 0 || pos == x.size() - 1) return pos;
        x0 = (pos < 1) ? pos : pos - 1;
        x2 = (pos + 1 < x.size()) ? pos + 1 : pos;
        if (x0 == pos) return (x[pos] <= x[x2]) ? pos : x2;
        if (x2 == pos) return (x[pos] <= x[x0]) ? pos : x0;
        s0 = x[x0];
        s1 = x[pos];
        s2 = x[x2];
        return pos + half * (s0 - s2 ) / (s0 - two * s1 + s2);
    }


    template<typename T>
    T quadratic_peak_mag (const sample_vector<T>& x, T pos) {
        T x0, x1, x2;
        size_t index = (size_t)(pos - .5) + 1;
        if (pos >= x.size() || pos < 0.) return 0.;
        if ((T)index == pos) return x[index];
        x0 = x[index - 1];
        x1 = x[index];
        x2 = x[index + 1];
        return x1 - .25 * (x0 - x2) * (pos - index);
    }

    template<typename T>
    size_t peakpick(const sample_vector<T>& onset, size_t pos) {
        size_t tmp=0;
        tmp = (onset[pos] > onset[pos-1]
            &&  onset[pos] > onset[pos+1]
            &&  onset[pos] > 0.);
        return tmp;
    }

    template<typename T>
    T quadfrac (T s0, T s1, T s2, T pf)
    {
        T tmp =
            s0 + (pf / 2.) * (pf * (s0 - 2. * s1 + s2) - 3. * s0 + 4. * s1 - s2);
        return tmp;
    }

    template<typename T>
    T freqtomidi (T freq)
    {
        T midi;
        if (freq < 2. || freq > 100000.) return 0.; // avoid nans and infs
        /* log(freq/A-2)/log(2) */
        midi = freq / 6.875;
        midi = std::log (midi) / 0.6931471805599453;
        midi *= 12;
        midi -= 3;
        return midi;
    }

    template<typename T>
    T miditofreq (T midi)
    {
        T freq;
        if (midi > 140.) return 0.; // avoid infs
        freq = (midi + 3.) / 12.;
        freq = EXP (freq * 0.6931471805599453);
        freq *= 6.875;
        return freq;
    }

    template<typename T>
    T bintofreq (T bin, T samplerate, T fftsize)
    {
        T freq = samplerate / fftsize;
        return freq * MAX(bin, 0);
    }

    template<typename T>
    T bintomidi (T bin, T samplerate, T fftsize)
    {
        T midi = bintofreq (bin, samplerate, fftsize);
        return freqtomidi (midi);
    }

    template<typename T>
    T freqtobin (T freq, T samplerate, T fftsize)
    {
        T bin = fftsize / samplerate;
        return MAX(freq, 0) * bin;
    }

    template<typename T>
    T miditobin (T midi, T samplerate, T fftsize)
    {
        T freq = miditofreq (midi);
        return freqtobin (freq, samplerate, fftsize);
    }


    size_t is_power_of_two (size_t a)
    {
        if ((a & (a - 1)) == 0) {
            return 1;
        } else {
            return 0;
        }
    }


    size_t next_power_of_two (size_t a)
    {
        size_t i = 1;
        while (i < a) i <<= 1;
            return i;
        }

    size_t power_of_two_order (size_t a)
    {
        int order = 0;
        int temp = next_power_of_two(a);
        while (temp >>= 1) {
            ++order;
        }
        return order;
    }

    template<typename T>
    T db_spl (const sample_vector<T>& o)
    {
        return 10. * LOG10 (level_lin (o));
    }

    template<typename T>
    size_t silence_detection (const sample_vector<T>& o, T threshold)
    {
        return (db_spl (o) < threshold);
    }

    template<typename T>
    T level_detection (const sample_vector<T>& o, T threshold)
    {
        T db_spl = db_spl (o);
        if (db_spl < threshold) {
            return 1.;
        } else {
            return db_spl;
        }
    }

    template<typename T>
    T zero_crossing_rate (sample_vector<T>& input)
    {
        size_t j;
        size_t zcr = 0;
        for (j = 1; j < input->length; j++) {
            // previous was strictly negative
            if (input->data[j - 1] < 0.) {
            // current is positive or null
            if (input->data[j] >= 0.) {
                zcr += 1;
            }
            // previous was positive or null
            } else {
            // current is strictly negative
            if (input->data[j] < 0.) {
                zcr += 1;
            }
            }
        }
        return zcr / (T) input->length;
    }

    template<typename T>
    void autocorr (const sample_vector<T>& input, sample_vector<T>& output)
    {
        size_t i, j, length = input->length;
        T *data, *acf;
        T tmp = 0;
        data = input->data;
        acf = output->data;
        #pragma omp simd
        for (i = 0; i < length; i++) {
            tmp = 0.;
            for (j = i; j < length; j++) {
                tmp += data[j - i] * data[j];
            }
            acf[i] = tmp / (T) (length - i);
        }
    }

    template<typename T>
    sample_vector<T> create_window (size_t n, window_type wintype) {
        sample_vector<T> win(n);
        T * w = win.data();
        size_t i, size =n;
    
        switch(wintype) {
        case win_ones:
            win.ones();      
            break;
        case win_rectangle:
            win.fill(.5);
            break;
        case win_hamming:
            #pragma omp simd
            for (i=0;i<size;i++)
                w[i] = 0.54 - 0.46 * std::cos(2.0*M_PI * i / (size));
            break;
        case win_hanning:
            #pragma omp simd
            for (i=0;i<size;i++)
                w[i] = 0.5 - (0.5 * std::cos(2.0*M_PI * i / (size)));
            break;
        case win_hanningz:
            #pragma omp simd
            for (i=0;i<size;i++)
                w[i] = 0.5 * (1.0 - std::cos(2.0*M_PI * i / (size)));
            break;
        case win_blackman:
            #pragma omp simd
            for (i=0;i<size;i++)
                w[i] = 0.42
                - 0.50 * std::cos(    2.0*M_PI*i/(size-1.0))
                + 0.08 * std::cos(2.0*2.0*M_PI*i/(size-1.0));
            break;
        case win_blackman_harris:
        #pragma omp simd
            for (i=0;i<size;i++)
                w[i] = 0.35875
                - 0.48829 * std::cos(    2.0*M_PI*i/(size-1.0))
                + 0.14128 * std::cos(2.0*2.0*M_PI*i/(size-1.0))
                - 0.01168 * std::cos(3.0*2.0*M_PI*i/(size-1.0));
            break;
        case win_gaussian:
            {
            T a, b, c = 0.5;
            size_t n;
            #pragma omp simd
            for (n = 0; n < size; n++)
            {
                a = (n-c*(size-1))/(SQR(c)*(size-1));
                b = -c*SQR(a);
                w[n] = std::exp(b);
            }
            }
            break;
        case win_welch:
            #pragma omp simd
            for (i=0;i<size;i++)
                w[i] = 1.0 - SQR((2.*i-size)/(size+1.0));
            break;
        case win_parzen:
            #pragma omp simd
            for (i=0;i<size;i++)
                w[i] = 1.0 - std::fabs((2.f*i-size)/(size+1.0f));
            break;
        default:
            break;
        }  
    }

    template<typename T>
    T hztomel (T freq)
    {
        const T lin_space = 3./200.;
        const T split_hz = 1000.;
        const T split_mel = split_hz * lin_space;
        const T log_space = 27./std::log(6400/1000.);
        if (freq < 0) {
            std::cerr << ("hztomel: input frequency should be >= 0\n");
            return 0;
        }
        if (freq < split_hz)
        {
            return freq * lin_space;
        } else {
            return split_mel + log_space * std::log (freq / split_hz);
        }

    }

    template<typename T>
    T meltohz (T mel)
    {
        const T lin_space = 200./3.;
        const T split_hz = 1000.;
        const T split_mel = split_hz / lin_space;
        const T logSpacing = std::pow(6400/1000., 1/27.);
        if (mel < 0) {
            std::cerr << "meltohz: input mel should be >= 0\n";
            return 0;
        }
        if (mel < split_mel) {
            return lin_space * mel;
        } else {
            return split_hz * std::pow(logSpacing, mel - split_mel);
        }
    }

    template<typename T>
    T hztomel_htk (T freq)
    {
        const T split_hz = 700.;
        const T log_space = 1127.;
        if (freq < 0) {
            std::cerr << "hztomel_htk: input frequency should be >= 0\n";
            return 0;
        }
        return log_space * std::log(1 + freq / split_hz);
    }

    template<typename T>
    T meltohz_htk (T mel)
    {
        const T split_hz = 700.;
        const T log_space = 1./1127.;
        if (mel < 0) {
            std::cerr << "meltohz_htk: input frequency should be >= 0\n";
            return 0;
        }
        return split_hz * (std::exp( mel * log_space) - 1.);
    }

}