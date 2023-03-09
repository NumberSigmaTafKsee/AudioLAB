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

//#include "samples/Allocator.hpp"
#include "carlo_mkl.hpp"
#include "SndFile.hpp"

namespace Casino
{
    
    struct wav_data {
        int64_t frames;
        size_t  size;
        int     samplerate;
        int     channels;
        int     format;
        int     sections;
    };
    MKL::Vector<float> load_wave_float(const std::string& file, wav_data & info)
    {
        SndFileReaderFloat r(file.c_str());
        MKL::Vector<float> wav(r.size());
        r.read(r.size(),wav.data());
        info.frames = r.frames();
        info.size   = r.size();
        info.samplerate = r.samplerate();
        info.channels = r.channels();
        info.format   = r.format();
        info.sections = r.sections();
        return wav;
    }
    void save_wave_float(const std::string& file, MKL::Vector<float> & samples, wav_data& info)
    {
        SndFileWriterFloat r(file.c_str(),info.format,info.channels,info.samplerate);
        r.write(samples.size(),samples.data());
    }
    MKL::Vector<double> load_wave_double(const std::string& file, wav_data & info)
    {
        SndFileReaderDouble r(file.c_str());
        MKL::Vector<double> wav(r.size());
        r.read(r.size(),wav.data());
        info.frames = r.frames();
        info.size   = r.size();
        info.samplerate = r.samplerate();
        info.channels = r.channels();
        info.format   = r.format();
        info.sections = r.sections();
        return wav;
    }
    void save_wave_double(const std::string& file, MKL::Vector<double> & samples, wav_data& info)
    {
        SndFileWriterDouble r(file.c_str(),info.format,info.channels,info.samplerate);
        r.write(samples.size(),samples.data());
    }


    template<typename T>
    T get_stride(size_t ch, size_t num_channels, size_t pos, MKL::Vector<T> & samples)
    {
        return samples[pos*num_channels + ch];
    }
        
    template<typename T>
    void set_stride(size_t ch, size_t num_channels, size_t pos, MKL::Vector<T> & samples, T sample)
    {
        samples[pos*num_channels + ch] = sample;
    }

    template<typename T>
    MKL::Vector<T> get_left_channel(const MKL::Vector<T> & in) {
        MKL::Vector<T> r(in.size()/2);
        size_t x = 0;
        #pragma omp simd
        for(size_t i = 0; i < in.size(); i+=2) r[x++] = in[i];
        return r;
    }
    template<typename T>
    MKL::Vector<T> get_right_channel(const MKL::Vector<T> & in) {
        MKL::Vector<T> r(in.size()/2);
        size_t x = 0;
        #pragma omp simd
        for(size_t i = 1; i < in.size(); i+=2) r[x++] = in[i];
        return r;
    }
    template<typename T>
    MKL::Vector<T> get_channel(size_t ch, const MKL::Vector<T> & in) {
        MKL::Vector<T> r(in.size()/2);
        size_t x = 0;
        #pragma omp simd
        for(size_t i = ch; i < in.size(); i+=2) r[x++] = in[i];
        return r;
    }
    template<typename T>
    void set_left_channel(const MKL::Vector<T> & left, MKL::Vector<T> & out) {
        size_t x = 0;
        #pragma omp simd
        for(size_t i = 0; i < out.size(); i+=2) out[i] = left[x++];
    }
    template<typename T>
    void set_right_channel(const MKL::Vector<T> & right, MKL::Vector<T> & out) {
        size_t x = 0;
        #pragma omp simd
        for(size_t i = 1; i < out.size(); i+=2) out[i] = right[x++];
    }
    template<typename T>
    void set_channel(size_t ch, const MKL::Vector<T> & in, MKL::Vector<T> & out) {
        size_t x = 0;
        #pragma omp simd
        for(size_t i = ch; i < out.size(); i+=2) out[i] = in[x++];
    }
    template<typename T>
    MKL::Vector<T> interleave(size_t n, size_t channels, const MKL::Vector<MKL::Vector<T>> & in) {
        MKL::Vector<T> r(n*channels);
        #pragma omp parallel
        {       
        for(size_t i = 0; i < channels; i++)
            #pragma omp simd
            for(size_t j = 0; j < n; j++)
                r[j*channels + i] = in[i][j];
        }
        return r;
    }
    template<typename T>
    MKL::Vector<T> interleave(size_t n, size_t channels, const MKL::Vector<T*> & in) {
        MKL::Vector<T> r(n*channels);
        for(size_t i = 0; i < channels; i++)
            for(size_t j = 0; j < n; j++)
                r[j*channels + i] = in[i][j];
        return r;
    }
    template<typename T>
    MKL::Vector<MKL::Vector<T>> deinterleave(size_t n, size_t channels, const MKL::Vector<T> & in) {
        MKL::Vector<MKL::Vector<T>> r(n);
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
    MKL::Vector<T> get_left_channel(size_t n, const T* in) {
        MKL::Vector<T> r(n/2);
        size_t x = 0;
        #pragma omp simd
        for(size_t i = 0; i < n; i+=2) r[x++] = in[i];
        return r;
    }
    template<typename T>
    MKL::Vector<T> get_right_channel(size_t n, const T* & in) {
        MKL::Vector<T> r(n/2);
        size_t x = 0;
        #pragma omp simd
        for(size_t i = 1; i < n; i+=2) r[x++] = in[i];
        return r;
    }
    template<typename T>
    MKL::Vector<T> get_channel(size_t ch, size_t n, T* in) {
        MKL::Vector<T> r(n/2);
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
    MKL::Vector<T> interleave(size_t n, size_t channels, const T** & in) {
        MKL::Vector<T> r(n*channels);
        #pragma omp simd
        for(size_t i = 0; i < channels; i++)
            for(size_t j = 0; j < n; j++)
                r[j*channels + i] = in[i][j];
        return r;
    }
    template<typename T>
    MKL::Vector<MKL::Vector<T>> deinterleave(size_t n, size_t channels, const T* & in) {
        MKL::Vector<MKL::Vector<T>> r(n);
        #pragma omp parallel
        {
            for(size_t i = 0; i < channels; i++)
            {
                r[i].resize(n);
                #pragma omp simd
                for(size_t j = 0; j < n; j++)
                    r[i][j] = in[j*channels + i];
            }
        }
        return r;
    }

    template<typename T>
    bool equal_vector (MKL::Vector<T> & a, MKL::Vector<T> & b) {
        return std::equal(a.begin(),a.end(),b.end());
    }

    template<typename T>
    void copy_vector(MKL::Vector<T> & dst, MKL::Vector<T> & src) {
        std::copy(src.begin(),src.end(),dst.begin());
    }
    template<typename T>
    void copy_vector(MKL::Vector<T> & dst, size_t n, T * src) {
        std::copy(&src[0],&src[n-1],dst.begin());
    }
    template<typename T>
    MKL::Vector<T> slice_vector(size_t start, size_t end, MKL::Vector<T> & src) {
        MKL::Vector<T> r(end-start);
        std::copy(src.begin()+start,src.begin()+end,r.begin());
        return r;
    }

    template<typename T>
    void copy_buffer(size_t n, T * dst, T * src) {
        memcpy(dst,src,n*sizeof(T));
    }

    template<typename T>
    MKL::Vector<T> slice_buffer(size_t start, size_t end, T * ptr) {
        MKL::Vector<T> r(end-start);    
        std::copy(&ptr[start],&ptr[end],r.begin());
        return r;
    }

    template<typename T>
    void split_stereo(size_t n, const T* input, T * left, T * right)
    {
        size_t x=0;
        #pragma omp simd
        for(size_t i = 0; i < n; i+=2)
        {
            left[x] = input[i];
            right[x++] = input[i+1];
        }
    }

    template<typename T>
    void split_stereo(const MKL::Vector<T> & input, MKL::Vector<T> & left, MKL::Vector<T> & right) {
        size_t x = input.size();
        left.resize(x/2);
        right.resize(x/2);
        split_stereo(x,input.data(),left.data(),right.data());
    }

    template<typename T>
    T insert_front(size_t n, T in, T * buffer) {
        T r = buffer[n-1];
        #pragma omp simd
        for(size_t i=0; i < n-1; i++) buffer[n+1] = buffer[n];
        buffer[0] = in;
        return r;
    }
  
    
    template<typename T>
    void swap(MKL::Vector<T> & left, MKL::Vector<T> & right) {
        std::swap(left,right);
    }

    template<typename T>
    bool is_in(const MKL::Vector<T> & v, const T val) {
        return std::find(v.begin(),v.end(),val) != v.end();
    }

    template<typename T>
    MKL::Vector<T> normalize(const MKL::Vector<T> & a) {
        MKL::Vector<T> r(a);        
        T max = *std::max_element(r.begin(),r.end());        
        if(max > 0) 
            #pragma omp simd
            for(size_t i = 0; i < r.size(); i++) r[i] /= max;
        return r;
    }

    template<typename T>
    MKL::Vector<T> mix(const MKL::Vector<T> & a, const MKL::Vector<T> & b)
    {
        assert(a.size() == b.size());
        MKL::Vector<T> r(a.size());        
        #pragma omp simd
        for(size_t i = 0; i < r.size(); i++) 
        {
            r[i] = a[i]+b[i];        
        }
        return normalize(r);        
    }
    
    template<class A, class B>
    MKL::Vector<B> convert(const MKL::Vector<A> & v) {
        MKL::Vector<B> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = B(v[i]);
    }
    template<class T>
    MKL::Vector<T> kernel(const MKL::Vector<T> & v, T (*f)(T value)) {
        MKL::Vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = f(v[i]);
        return r;
    }
    template<class T>
    MKL::Vector<T> kernel(const MKL::Vector<T> & v, std::function<T (T)> func) {
        MKL::Vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = func(v[i]);
        return r;
    }
    template<class T>
    void inplace_add(const MKL::Vector<T> & a, MKL::Vector<T> & r, std::function<T (T)> func) {        
        #pragma omp simd
        for(size_t i = 0; i < a.size(); i++) r[i] += func(a[i]);        
    }
    template<class T>
    void inplace_sub(const MKL::Vector<T> & a, MKL::Vector<T> & r, std::function<T (T)> func) {        
        #pragma omp simd
        for(size_t i = 0; i < a.size(); i++) r[i] -= func(a[i]);        
    }
    template<class T>
    void inplace_mul(const MKL::Vector<T> & a, MKL::Vector<T> & r, std::function<T (T)> func) {        
        #pragma omp simd
        for(size_t i = 0; i < a.size(); i++) r[i] *= func(a[i]);        
    }
    template<class T>
    void inplace_div(const MKL::Vector<T> & a, MKL::Vector<T> & r, std::function<T (T)> func) {        
        #pragma omp simd
        for(size_t i = 0; i < a.size(); i++) r[i] /= func(a[i]);
    }


    template<class T>
    void fill(MKL::Vector<T> & in, T x)
    {
        #pragma omp simd
        for(size_t i = 0; i < in.size(); i++) in[i] = x;
    }
    template<class T>
    void zeros(MKL::Vector<T> & in)
    {
        fill(in,T(0));
    }
    template<class T>
    void ones(MKL::Vector<T> & in)
    {
        fill(in,T(1));
    }

    
    template<typename A, typename B>
    MKL::Vector<A> vector_cast(MKL::Vector<B> & in) {
        MKL::Vector<A> r(in.size());
        #pragma omp simd
        for(size_t i = 0; i < in.size(); i++)
            r[i] = (A)in[i];
        return r;
    }
    template<typename T>
    MKL::Vector<T> vector_copy(T * ptr, size_t n) {
        MKL::Vector<T> r(n);
        #pragma omp simd
        for(size_t i = 0; i < n; i++)
            r[i] = ptr[i];
        return r;
    }


    template <class T>
    void zeros(std::vector<T> & v) {
        std::fill(v.begin(),v.end(),T(0));
    }

}

/*
#include "carlo_samples.inl"
#include "carlo_vtk.hpp"
#include "carlo_ringbuffer.hpp"
*/