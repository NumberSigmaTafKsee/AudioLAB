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

/* doesn't work like this
sample_vector<float> Upsampler(float sample_rate, unsigned int factor, size_t n, float * p)
{
    SimpleResampler resample;
    resample.setup(sample_rate,factor);
    sample_vector<float> out(factor * n);
    resample.up(n,p,out.data());
    return out;
}

sample_vector<float> Downsampler(float sample_rate, unsigned int factor, size_t n, float * p)
{
    SimpleResampler resample;
    resample.setup(sample_rate,factor);
    // need to filter it first
    sample_vector<float> out(n/factor);
    resample.down(n,p,out.data());
    return out;
}
*/

/* doesn't work this way the resampler needs to remain around
template<typename T>
kfr::univector<T> resample(kfr::sample_rate_conversion_quality quality, const kfr::univector<T> & input, size_t output_sr, size_t input_sr) {    
    auto r = kfr::resampler<T>(quality,output_sr, input_sr);    
    kfr::univector<T> output(input.size() * output_sr/input_sr + r.get_delay());
    r.process(output,input);
    return output;    
}

template<typename T>
kfr::univector<T> kfr_upsampler(kfr::resample_quality quality, int factor,  sample_vector<T> in, T sr=44100.0f)
{    
    sample_vector<T> out(factor*in.size());
    auto r = kfr::resampler<T>(quality,factor*sr,sr);
    r.process(out,in);
    return out;
}
template<typename T>
kfr::univector<T> kfr_downsampler(kfr::resample_quality quality, int factor, sample_vector<T> in, T sr = 44100.0f)
{    
    sample_vector<T> out(in.size()/factor);
    auto r = kfr::resampler<T>(quality,sr,factor*sr);
    r.process(out,in);    
    return out;
}
template<typename T>
kfr::univector<T> downsampler_draft(int factor,  sample_vector<T> in)
{
    return kfr_downsampler(kfr::resample_quality::draft,factor,in);
}
template<typename T>
kfr::univector<T> downsampler_low(int factor,  sample_vector<T> in)
{
    return kfr_downsampler(kfr::resample_quality::low,factor,in);
}
template<typename T>
kfr::univector<T> downsampler_normal(int factor, sample_vector<T> in)
{
    return kfr_downsampler(kfr::resample_quality::normal,factor,in);
}
template<typename T>
kfr::univector<T> downsampler_high(int factor, sample_vector<T> in)
{
    return kfr_downsampler(kfr::resample_quality::high,factor,in);
}
template<typename T>
kfr::univector<T> downsampler_perfect(int factor, sample_vector<T> in)
{
    return kfr_downsampler(kfr::resample_quality::perfect,factor,in);
}

template<typename T>
kfr::univector<T> upsampler_draft(int factor, sample_vector<T> in)
{
    return kfr_upsampler(kfr::resample_quality::draft,factor,in);
}
template<typename T>
kfr::univector<T> upsampler_low(int factor, sample_vector<T> in)
{
    return kfr_upsampler(kfr::resample_quality::low,factor,in);
}
template<typename T>
kfr::univector<T> upsampler_normal(int factor, sample_vector<T> in)
{
    return kfr_upsampler(kfr::resample_quality::normal,factor,in);
}
template<typename T>
kfr::univector<T> upsampler_high(int factor, sample_vector<T> in)
{
    return kfr_upsampler(kfr::resample_quality::high,factor,in);
}
template<typename T>
kfr::univector<T> upsampler_perfect(int factor, sample_vector<T> in)
{
    return kfr_upsampler(kfr::resample_quality::perfect,factor,in);
}
*/