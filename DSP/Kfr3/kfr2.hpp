#pragma once

#include <vector>
#include <cmath>
#include <complex>

// will be huge problem with c++
#define KFR_NO_C_COMPLEX_TYPES
#include <kfr/capi.h>
#include <kfr/dft.hpp>
#include <kfr/dsp.hpp>

#include "kfr3_samples.hpp"


namespace kfr3
{
    template<typename T>
    std::vector<T> 
    cheby1_lowpass(T cutoff, T sample_rate) {
        kfr::zpk<T> filt = kfr::iir_lowpass(kfr::chebyshev1<T>(order,w),cutoff,sample_rate);        
        std::vector<kfr::biquad_params<T>> bqs = kfr::to_sos<T>(filt);
        std::vector<T> r(bqs.size()*6);
        size_t x = 0;
        for(size_t i = 0; i < bqs.size(); i)
        {            
            r[x++] = bqs[i].b0;
            r[x++] = bqs[i].b1;
            r[x++] = bqs[i].b2;
            r[x++] = bqs[i].a0;
            r[x++] = bqs[i].a1;
            r[x++] = bqs[i].a2;
        }            
        return r;
    }

    template<typename T>
    std::vector<T> 
    cheby1_highpass(T cutoff, T sample_rate) {
        kfr::zpk<T> filt = kfr::iir_highpass(kfr::chebyshev1<T>(order,w),cutoff,sample_rate);        
        std::vector<T> r(bqs.size()*6);
        size_t x = 0;
        for(size_t i = 0; i < bqs.size(); i)
        {            
            r[x++] = bqs[i].b0;
            r[x++] = bqs[i].b1;
            r[x++] = bqs[i].b2;
            r[x++] = bqs[i].a0;
            r[x++] = bqs[i].a1;
            r[x++] = bqs[i].a2;
        }            
        return r;
    }

    template<typename T>
    std::vector<T> 
    cheby1_bandstop(T locutoff, T hicutoff, T sample_rate) {
        kfr::zpk<T> filt = kfr::iir_bandstop(kfr::chebyshev1<T>(order,w),locutoff,hicutoff,sample_rate);        
        std::vector<T> r(bqs.size()*6);
        size_t x = 0;
        for(size_t i = 0; i < bqs.size(); i)
        {            
            r[x++] = bqs[i].b0;
            r[x++] = bqs[i].b1;
            r[x++] = bqs[i].b2;
            r[x++] = bqs[i].a0;
            r[x++] = bqs[i].a1;
            r[x++] = bqs[i].a2;
        }            
        return r;
    }

    template<typename T>
    std::vector<T> 
    cheby1_bandpass(T locutoff, T hicutoff, T sample_rate) {
        kfr::zpk<T> filt = kfr::iir_bandpass(kfr::chebyshev1<T>(order,w),locutoff,hicutoff,sample_rate);        
        std::vector<T> r(bqs.size()*6);
        size_t x = 0;
        for(size_t i = 0; i < bqs.size(); i)
        {            
            r[x++] = bqs[i].b0;
            r[x++] = bqs[i].b1;
            r[x++] = bqs[i].b2;
            r[x++] = bqs[i].a0;
            r[x++] = bqs[i].a1;
            r[x++] = bqs[i].a2;
        }            
        return r;
    }

    template<typename T>
    sample_vector<T> convolve(size_t n, const T *a, const T * b) {
        kfr::univector<T> v1(n);
        kfr::univector<T> v2(n);
        memcpy(v1.data(),a,sizeof(T)*n);
        memcpy(v2.data(),b,sizeof(T)*n);
        kfr::univector<T> r = (kfr::convolve(v1,v2));
        sample_vector<T> output(r.size());
        memcpy(output.data(),r.data(),n*sizeof(T));
        return sample_vector<T>(output);
    }

    template<typename T>
    sample_vector<T> convolve(const sample_vector<T> &a, const sample_vector<T> &b) {        
        assert(a.size() == b.size());
        kfr::univector<T> r = (kfr::convolve(a,b));        
        return sample_vector<T>(r);
    }
    
    // there isn't lag in kfr
    template<typename T>
    sample_vector<T> xcorr(size_t n, const T *src1, const T * src2) {
        kfr::univector<T> v1(n);
        kfr::univector<T> v2(n);
        memcpy(v1.data(),src1.data(),sizeof(T)*n);
        memcpy(v2.data(),src2.data(),sizeof(T)*n);
        kfr::univector<T> r = (kfr::correlate(v1,v2));
        sample_vector<T> output(r.size());
        memcpy(output.data(),r.data(),n*sizeof(T));
        return sample_vector<T>(output);
    }
    
    template<typename T>
    sample_vector<T> xcorr(const sample_vector<T> &a, const sample_vector<T> &b) {        
        assert(a.size() == b.size());
        kfr::univector<T> r = (kfr::correlate(a,b));        
        return sample_vector<T>(r);
    }
    
    // there isn't lag in kfr
    template<typename T>
    sample_vector<T> acorr(size_t n, const T *src) {
        kfr::univector<T> v1(n);        
        memcpy(v1.data(),src.data(),sizeof(T)*n);        
        kfr::univector<T> r = (kfr::autocorrelate(v1));
        sample_vector<T> output(r.size());
        memcpy(output.data(),r.data(),n*sizeof(T));
        return sample_vector<T>(output);
    }
    
    template<typename T>
    sample_vector<T> acorr(const sample_vector<T> &a) {        
        kfr::univector<T> r = (kfr::autocorrelate(a));        
        return sample_vector<T>(r);
    }

    
    template<typename T>
    sample_vector<T> make_univec(kfr::univector<T> & r) {        
        sample_vector<T> x(s);
        memcpy(s.data(),r.data(),s*sizeof(T));
        return x;
    }

    template<typename T>
    sample_vector<T> make_window_hann(size_t s) {
        return make_univec(kfr::univector<T>(kfr::window_hann<T>(s)));
    }
    template<typename T>
    sample_vector<T> make_window_hamming(size_t s) {
        return make_univec(kfr::univector<T>(kfr::window_hamming<T>(s)));
    }
    template<typename T>
    sample_vector<T> make_window_blackman(size_t s, const T alpha=T(0.16), window_symmetry symmetry = kfr::window_symmetry::symmetric) {
        return make_univec(kfr::univector<T>(kfr::window_blackman<T>(s,alpha,symmetry)));
    }
    template<typename T>
    sample_vector<T> make_window_blackman_harris(size_t s, window_symmetry symmetry = kfr::window_symmetry::symmetric) {
        return make_univec(kfr::univector<T>(kfr::window_blackman_harris<T>(s,symmetry)));
    }
    template<typename T>
    sample_vector<T> make_window_gaussian(size_t s, const T alpha=T(0.25)) {
        return make_univec(kfr::univector<T>(kfr::window_gaussian<T>(s,alpha)));
    }
    template<typename T>
    sample_vector<T> make_window_triangular(size_t s) {
        return make_univec(kfr::univector<T>(kfr::window_triangular<T>(s)));
    }
    template<typename T>
    sample_vector<T> make_window_bartlett(size_t s) {
        return make_univec(kfr::univector<T>(kfr::window_bartlett<T>(s)));
    }
    template<typename T>
    sample_vector<T> make_window_cosine(size_t s) {
        return make_univec(kfr::univector<T>(kfr::window_cosine<T>(s)));
    }
    template<typename T>
    sample_vector<T> make_window_bartlett_hann(size_t s) {
        return make_univec(kfr::univector<T>(kfr::window_bartlett_hann<T>(s)));
    }
    template<typename T>
    sample_vector<T> make_window_bohman(size_t s) {
        return make_univec(kfr::univector<T>(kfr::window_bohman<T>(s)));
    }
    template<typename T>
    sample_vector<T> make_window_lanczos(size_t s) {
        return make_univec(kfr::univector<T>(kfr::window_lanczos<T>(s)));
    }
    template<typename T>
    sample_vector<T> make_window_flattop(size_t s) {
        return make_univec(kfr::univector<T>(kfr::window_flattop<T>(s)));
    }
    template<typename T>
    sample_vector<T> make_window_rectangular(size_t s) {
        return make_univec(kfr::univector<T>(kfr::window_rectangular<T>(s)));
    }
    template<typename T>
    sample_vector<T> make_window_kaiser(size_t s, const T beta = T(0.5)) {
        return make_univec(kfr::univector<T>(kfr::window_kaiser<T>(s,beta)));
    }

    template<typename T>
    T energy_to_loudness(T energy) {
        return kfr::energy_to_loudness(energy);
    }
    template<typename T>
    T loudness_to_energy(T loudness) {
        return kfr::loudness_to_energy(loudness);
    }    
    template<typename T> T normalize_frequency(T f, T sample_rate) {
        return f/sample_rate;
    }
    
    template<typename T>
    T amp_to_dB(const T & in) {
        return kfr::amp_to_dB(in);
    }        

        
    template<typename T>
    sample_vector<T> dcremove(const sample_vector<T> & input, T cutoff) {
         kfr::univector<T> out(input.size());
        auto x = kfr::dcremove(input,cutoff);
        x.apply(input,out);
        return make_univec(out);
    }
    template<typename T>
    sample_vector<T> sinewave(size_t n, T freq, T sample_rate, T phase=(T)0) {
        kfr::univector<T> r(n);                
        for(size_t i = 0; i < n; i++)
        {
            r[i] = kfr::sine(2*M_PI*phase);
            phase += freq/sample_rate;
            if(phase > (T)1.0) phase-=(T)1.0;
        }
        return make_univev(r);
    }
    template<typename T>
    sample_vector<T> squarewave(size_t n, T freq, T sample_rate, T phase=(T)0) {
        kfr::univector<T> r(n);
        for(size_t i = 0; i < n; i++)
        {
            r[i] = kfr::square(2*M_PI*phase);
            phase += freq/sample_rate;
            if(phase > (T)1.0) phase-=(T)1.0;
        }
        return make_univec(r);
    }
    template<typename T>
    kfr::univector<T> trianglewave(size_t n, T freq, T sample_rate, T phase=(T)0) {
        kfr::univector<T> r(n);
        for(size_t i = 0; i < n; i++)
        {
            r[i] = kfr::triangle(2*M_PI*phase);
            phase += freq/sample_rate;
            if(phase > (T)1.0) phase-=(T)1.0;
        }
        return r;
    }
    template<typename T>
    sample_vector<T> sawtoothwave(size_t n, T freq, T sample_rate, T phase=(T)0) {
        kfr::univector<T> r(n);
        for(size_t i = 0; i < n; i++)
        {
            r[i] = kfr::sawtooth(2*M_PI*phase);
            phase += freq/sample_rate;
            if(phase > (T)1.0) phase-=(T)1.0;
        }
        return make_univec(r);
    }
    template<typename T>
    sample_vector<T> generate_sin(size_t n, T start, T step) {
        kfr::univector<T> r(n);
        r = kfr::gen_sin(start,step);
        return make_univec(r);
    }

    template<typename T>
    sample_vector<T> generate_linear(size_t n, T start, T step) {
        kfr::univector<T> r(n);
        r = kfr::gen_linear(start,step);
        return make_univec(r);
    }

    template<typename T>
    sample_vector<T> generate_exp(size_t n, T start, T step) {
        kfr::univector<T> r(n);
        r = kfr::gen_exp(start,step);
        return make_univec(r);
    }

    template<typename T>
    sample_vector<T> generate_exp2(size_t n, T start, T step) {
        kfr::univector<T> r(n);
        r = kfr::gen_exp2(start,step);
        return make_univec(r);        
    }

    template<typename T>
    sample_vector<T> generate_cossin(size_t n, T start, T step) {
        kfr::univector<T> r(n);
        r = kfr::gen_cossin(start,step);
        return make_univec(r);
    }

    template<typename T>
    void deinterleave(sample_matrix<T> & out, const T * ptr, size_t size, size_t channels) {                        
        kfr::deinterleave(p.v.data(),ptr,channels,size);
    }
    
    template<typename T>
    sample_vector<T> interleave(const sample_matrix<T> & input) {        
        return make_univec(kfr::interleave(input.v));
    }

    template<typename T>
    sample_vector<T> generate_sequence(int n, int start, int step) {
        sample_vector<T> r(n);
        for(size_t i = 0; i < n; i++)
            r[i] = start + step*n;
        return r;
    }

    std::string dB_to_string(const double value, double min=-140.0f) { return kfr::dB_to_string(value,min); }
    std::string dB_to_utf8string(const double value, double min=-140.0f) { return kfr::dB_to_utf8string(value,min); }
}