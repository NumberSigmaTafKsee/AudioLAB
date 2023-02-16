#pragma once
#include <algorithm>
#include <functional>
#include <cstdint>
#include <cmath>
#include <iostream>
#include "SndFile.hpp"
#include "SoundObject.hpp"
#include "audio_function_generator.hpp"
#include "audio_functions.hpp"

namespace AudioDSP
{    
    struct wav_data {
        int64_t frames;
        size_t  size;
        int     samplerate;
        int     channels;
        int     format;
        int     sections;
    };


    void load_wave(const char* file, wav_data & info, size_t n, float * wav)
    {
        SndFileReaderFloat r(file);        
        r.read(n,wav);
        info.frames = r.frames();
        info.size   = r.size();
        info.samplerate = r.samplerate();
        info.channels = r.channels();
        info.format   = r.format();
        info.sections = r.sections();        
    }    
    
    void save_wave(const char* file, size_t n, float * samples, wav_data& info)
    {
        SndFileWriterFloat r(file,info.format,info.channels,info.samplerate);
        r.write(n,samples);
    }

    void load_wave(const char* file, wav_data & info, size_t n, double * wav)
    {
        SndFileReaderDouble r(file);        
        r.read(n,wav);
        info.frames = r.frames();
        info.size   = r.size();
        info.samplerate = r.samplerate();
        info.channels = r.channels();
        info.format   = r.format();
        info.sections = r.sections();        
    }    
    
    void save_wave(const char* file, size_t n, double * samples, wav_data& info)
    {
        SndFileWriterDouble r(file,info.format,info.channels,info.samplerate);
        r.write(n,samples);
    }
    
    
    template<typename T>
    T get_stride(size_t ch, size_t num_channels, size_t pos, size_t n, T * samples)
    {
        return samples[pos*num_channels + ch];
    }
        
    template<typename T>
    void set_stride(size_t ch, size_t num_channels, size_t pos, size_t n, T * samples, T sample)
    {
        samples[pos*num_channels + ch] = sample;
    }

    template<typename T>
    void get_left_channel(const size_t n, const T * in, T * r) {        
        size_t x = 0;
        #pragma omp simd
        for(size_t i = 0; i < in.size(); i+=2) r[x++] = in[i];        
    }
    template<typename T>
    void get_right_channel(const size_t n, const T * in, T  *r) {        
        size_t x = 0;
        #pragma omp simd
        for(size_t i = 1; i < in.size(); i+=2) r[x++] = in[i];    
    }
    template<typename T>
    void get_channel(size_t ch, const size_t n, T * in, T * r) {        
        size_t x = 0;
        #pragma omp simd
        for(size_t i = ch; i < in.size(); i+=2) r[x++] = in[i];        
    }
    template<typename T>
    void set_left_channel(const size_t n, T * left, T * out) {
        size_t x = 0;
        #pragma omp simd
        for(size_t i = 0; i < out.size(); i+=2) out[i] = left[x++];
    }
    template<typename T>
    void set_right_channel(const size_t n, T * right, T * out) {
        size_t x = 0;
        #pragma omp simd
        for(size_t i = 1; i < out.size(); i+=2) out[i] = right[x++];
    }
    template<typename T>
    void set_channel(size_t ch, size_t stride, const size_t n, T * in,  T * out) {
        size_t x = 0;
        #pragma omp simd
        for(size_t i = ch; i < out.size(); i+=stride) out[i] = in[x++];
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
    void fill_left_channel(size_t n, T * out, const T v) {
        size_t x = 0;
        #pragma omp simd
        for(size_t i = 0; i < out.size(); i+=2) out[i] = v;
    }
    template<typename T>
    void fill_right_channel(size_t n, T * out, const T v) {
        size_t x = 0;
        #pragma omp simd
        for(size_t i = 1; i < out.size(); i+=2) out[i] = v;
    }
    template<typename T>
    void fill_channel(size_t ch, size_t n, T * out, const T v) {
        size_t x = 0;
        #pragma omp simd
        for(size_t i = ch; i < out.size(); i+=2) out[i] = v;
    }


    
    template<typename T>
    void copy_vector(size_t n, T * dst, const T * src) {        
        memcpy(dst,src, n*sizeof(T));
    }
        
    template<typename T>
    void copy_buffer(size_t n, T * dst, const T * src) {                        
        memcpy(dst,src,n*sizeof(T));
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
    void merge_stereo(size_t n, const T * left, const T * right, T * out)
    {
        size_t x=0;
        for(size_t i = 0; i < n; i+=2)
        {
            out[i] = left[x];
            out[i+1] = right[x++];            
        }
    }
    
    template<typename T>
    void swap_channels(size_t n, T * left, T * right) {
        for(size_t i = 0; i < n; i++)
        {
            T x = left[i];
            left[i] = right[i];
            right[i] = x;
        }
    }


    template<typename T>
    void mix(const size_t n, const T * a, const T * b, T * r, T m = 0.5)
    {        
        #pragma omp simd
        for(size_t i = 0; i < r.size(); i++) 
        {            
            r[i] = m*(a[i]+b[i]);
        }                        
    }
    template<typename T>
    void normalize(const size_t n, const T * a, T * r) {        
        auto max = std::max_element(r,r+n);
        if(*max > 0) 
            #pragma omp simd
            for(size_t i = 0; i < r.size(); i++) r[i] /= *max;     
    }
    
    template<class T>
    void fill(size_t n, T * in, T x)
    {
        std::fill(in,in+n,x);        
    }
    template<class T>
    void zeros(size_t n, T * in)
    {
        fill(in,T(0));
    }
    template<class T>
    void ones(size_t n, T * in)
    {
        fill(in,T(1));
    }
    
    
    
    template<typename T>
    void linear_convolution(const T * h, const T * x, int M, int N, T * y)
    {            
        for(int i=0;i<M+N-1;i++)
        {  
            y[i]=0;                 
            for(int j=0;j<=i;j++)
                y[i]+=x[j]*h[i-j];
        }    
    }
    
        
    template<typename T>
    void circular_convolution(const T * h, const T * x, int M, int N, T * y)
    {                
        for(int n=0;n < M;n++)
        { 
            y[n]=0;                                  
            for(int k=0;k < N;k++)                
            {
                int j = (k>n)? n-k+N : n-k;
                y[n]=y[n]+x[k]*h[j];
            }
        }            
    }   
    
    template<typename T>
    void print_vector(size_t n, T * ptr) {
        for(size_t i = 0; i < n; i++) std::cout << ptr[i] << ",";
        std::cout << std::endl;
    }

    template<typename T>
    void generate_noise(size_t n, T * r)
    {        
        std::generate(r,r+n,[]() { return Oscillators::Generators::function_noise(); });        
    }
    template<typename T>
    void generate_sin(DspFloatType freq, DspFloatType sampleRate, size_t n, T * r)
    {        
        Oscillators::Generators::SineGenerator sine(freq,sampleRate);
        std::generate(r,r+n,[&sine]() { return sine(); });        
    }
    template<typename T>
    void generate_cos(DspFloatType freq, DspFloatType sampleRate, size_t n, T * r)
    {        
        Oscillators::Generators::CosGenerator cose(freq,sampleRate);
        std::generate(r,r+n,[&cose]() { return cose(); });     
    }
    template<typename T>
    void generate_tan(DspFloatType freq, DspFloatType sampleRate, size_t n, T * r)
    {        
        Oscillators::Generators::TanGenerator tane(freq,sampleRate);
        std::generate(r,r+n,[&tane]() { return tane(); });        
    }
    template<typename T>
    void generate_phasor(DspFloatType freq, DspFloatType sampleRate, size_t n, T * r)
    {        
        Oscillators::Generators::PhasorGenerator phasore(freq,sampleRate);
        std::generate(r,r+n,[&phasore]() { return phasore(); });        
    }
    template<typename T>
    void generate_square(DspFloatType freq, DspFloatType sampleRate, size_t n, T * r)
    {     
        Oscillators::Generators::SineGenerator squaree(freq,sampleRate);
        std::generate(r,r+n,[&squaree]() { return squaree(); });     
    }
    template<typename T>
    void generate_saw(DspFloatType freq, DspFloatType sampleRate, size_t n, T * r)
    {        
        Oscillators::Generators::SawGenerator sawe(freq,sampleRate);
        std::generate(r,r+n,[&sawe]() { return sawe(); });     
    }
    template<typename T>
    void generate_triangle(DspFloatType freq, DspFloatType sampleRate, size_t n, T * r)
    {        
        Oscillators::Generators::TriangleGenerator trie(freq,sampleRate);
        std::generate(r,r+n,[&trie]() { return trie(); });     
    }
    /*
    template<typename T>
    void generate_function(Oscillators::Generators::FunctionGenerator::Type type, DspFloatType freq, DspFloatType sampleRate, size_t n)
    {
        void r(n);
        Oscillators::Generators::FunctionGenerator func(type,freq,sampleRate);
        std::generate(r,r+n,[&func]() { return func(); });
    }
    */
    template<typename T>
    void oscillator(OscillatorProcessor & osc, size_t n, T * r)
    {        
        std::generate(r,r+n,[&osc]() { return osc.Tick(); });     
    }
    template<typename T>
    void generator(GeneratorProcessor & osc, size_t n, T * r)
    {        
        std::generate(r,r+n,[&osc]() { return osc.Tick(); });     
    }
    template<typename T>
    void filter(FilterProcessor & filt, size_t n, const T * in, T *r)
    {        
        memcpy(r.data(),in,n*sizeof(T));        
        std::for_each(r,r+n,[&filt](T & x) { x = filt.Tick(x); });     
    }
    template<typename T>
    void function(FunctionProcessor & func, size_t n, const T * in, T *r)
    {        
        memcpy(r.data(),in,n*sizeof(T));        
        std::for_each(r,r+n,[&func](T & x) { x = func.Tick(x); });     
    }


    template<typename T>
    T minIndex(size_t n, const T * in)
    {
        return *std::min_element(in,in+n);
    }
    template<typename T>
    T maxIndex(size_t n, const T * in)
    {
        return *std::max_element(in,in+n);
    }
    template<typename T>
    int min(size_t n, const T * in)
    {
        return std::min(in,in+n);
    }
    template<typename T>
    int max(size_t n, const T * in)
    {
        return std::max(in,in+n);
    }
}
