#pragma once
#include "FIRSincFilter.hpp"

struct Upsampler
{
    
    size_t upfactor;
    Upsampler(size_t x) {
        upfactor = x;
    }

    void ProcessBlock(size_t n, DspFloatType * in, DspFloatType * out) {
        #pragma omp simd
        for(size_t i = 0; i < n; i++) {
            size_t x = 0;
            out[x++] = in[i];
            for(; x < upfactor; x++)
                out[x] = 0;
        }        
    }
};

struct Downsampler
{
    FIRSincFilter filter;    
    size_t downfactor;
    Downsampler(size_t x) : filter(65,20000.0),0,"lp","hamming") {
        downfactor = x;        
    }
    void ProcessBlock(size_t n, DspFloatType * in, DspFloatType * out) {
        #pragma omp simd
        size_t q = 0;
        for(size_t i = 0; i < n; i++) {
            DspFloatType x = filter.filter(in[q += downfactor]);
            out[i] = x;
        }        
    }
};
