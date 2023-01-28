#pragma once 
#include "include/samples/kfr_sample_dsp.hpp"
#include "SimpleResampler.hpp"
#include "Decimators.hpp"

template<typename T>
sample_vector<T> Upsample(int factor, sample_vector<T> in) {
    sample_vector<T> r(factor * in.size());
    memset(r.data(),0,factor*in.size()*sizeof(T));
    for(size_t i = 0; i < in.size(); i++)
        r[i*factor] = in[i];
    return r;
}

template<typename T>
sample_vector<T> Downsample(int factor, sample_vector<T> in) {
    sample_vector<T> r(in.size()/factor);    
    memset(r.data(),0,in.size()/factor*sizeof(T));        
    for(size_t i = 0; i < r.size(); i++)
        r[i] = in[i*factor];
    return r;
}


