#pragma once

#include <vector>
#include <cmath>
#include <complex>

#include <kfr/base.hpp>
#include <kfr/dft.hpp>
#include <kfr/dsp.hpp>
#include <kfr/io.hpp>

// will be huge problem with c++
#define KFR_NO_C_COMPLEX_TYPES
#include <kfr/capi.h>

#include "kfr3_samples.hpp"
#include "kfr3_capi.hpp"

namespace kfr3
{
    template<typename T>
    sample_vector<T> dcremove(const sample_vector<T> & input, T cutoff) {
         kfr::univector<T> out(input.size());
        auto x = kfr::dcremove(input,cutoff);
        x.apply(input,out);
        return make_univec(out);
    }

    template<typename T>
    struct filter_coefficients
    {
        T z[3];
        T p[3];        
    };

    // bessel wont preserve group delay
    // it is ok for audio
    template<typename T>
    std::vector<filter_coefficients<T>>
    bessel_lowpass(int order, T cutoff, T sample_rate) {
        kfr::zpk<T> filt = kfr::iir_lowpass(kfr::bessel<T>(order),cutoff,sample_rate);        
        sample_vector<kfr::biquad_params<T>> bqs = kfr::to_sos<T>(filt);
        std::vector<filter_coefficients<T>> r(bqs.size());        
        for(size_t i = 0; i < bqs.size(); i)
        {            
            r[i].z[0] = bqs[i].b0;
            r[i].z[1] = bqs[i].b1;
            r[i].z[2] = bqs[i].b2;
            r[i].p[0] = bqs[i].a0;
            r[i].p[1] = bqs[i].a1;
            r[i].p[2] = bqs[i].a2;
        }            
        return r;
    }

    template<typename T>
    std::vector<filter_coefficients<T>>
    bessel_highpass(int order,T cutoff, T sample_rate) {
        kfr::zpk<T> filt = kfr::iir_highpass(kfr::bessel<T>(order),cutoff,sample_rate);        
        sample_vector<kfr::biquad_params<T>> bqs = kfr::to_sos<T>(filt);
        std::vector<filter_coefficients<T>> r(bqs.size());
        for(size_t i = 0; i < bqs.size(); i)
        {            
            r[i].z[0] = bqs[i].b0;
            r[i].z[1] = bqs[i].b1;
            r[i].z[2] = bqs[i].b2;
            r[i].p[0] = bqs[i].a0;
            r[i].p[1] = bqs[i].a1;
            r[i].p[2] = bqs[i].a2;
        }            
        return r;
    }

    template<typename T>
    std::vector<filter_coefficients<T>>
    bessel_bandstop(int order, T locutoff, T hicutoff, T sample_rate) {
        kfr::zpk<T> filt = kfr::iir_bandstop(kfr::bessel<T>(order),locutoff,hicutoff,sample_rate);        
        sample_vector<kfr::biquad_params<T>> bqs = kfr::to_sos<T>(filt);
        std::vector<filter_coefficients<T>> r(bqs.size());
        for(size_t i = 0; i < bqs.size(); i)
        {            
            r[i].z[0] = bqs[i].b0;
            r[i].z[1] = bqs[i].b1;
            r[i].z[2] = bqs[i].b2;
            r[i].p[0] = bqs[i].a0;
            r[i].p[1] = bqs[i].a1;
            r[i].p[2] = bqs[i].a2;
        }            
        return r;
    }

    template<typename T>
    std::vector<filter_coefficients<T>>
    bessel_bandpass(int order, T locutoff, T hicutoff, T sample_rate) {
        kfr::zpk<T> filt = kfr::iir_bandpass(kfr::bessel<T>(order),locutoff,hicutoff,sample_rate);        
        sample_vector<kfr::biquad_params<T>> bqs = kfr::to_sos<T>(filt);
        std::vector<filter_coefficients<T>> r(bqs.size());
        for(size_t i = 0; i < bqs.size(); i)
        {            
            r[i].z[0] = bqs[i].b0;
            r[i].z[1] = bqs[i].b1;
            r[i].z[2] = bqs[i].b2;
            r[i].p[0] = bqs[i].a0;
            r[i].p[1] = bqs[i].a1;
            r[i].p[2] = bqs[i].a2;
        }            
        return r;
    }

    template<typename T>
    std::vector<filter_coefficients<T>>
    butterworth_lowpass(int order, T cutoff, T sample_rate) {
        kfr::zpk<T> filt = kfr::iir_lowpass(kfr::butterworth<T>(order),cutoff,sample_rate);        
        sample_vector<kfr::biquad_params<T>> bqs = kfr::to_sos<T>(filt);
        std::vector<filter_coefficients<T>> r(bqs.size());        
        for(size_t i = 0; i < bqs.size(); i)
        {            
            r[i].z[0] = bqs[i].b0;
            r[i].z[1] = bqs[i].b1;
            r[i].z[2] = bqs[i].b2;
            r[i].p[0] = bqs[i].a0;
            r[i].p[1] = bqs[i].a1;
            r[i].p[2] = bqs[i].a2;
        }            
        return r;
    }

    template<typename T>
    std::vector<filter_coefficients<T>>
    butterworth_highpass(int order,T cutoff, T sample_rate) {
        kfr::zpk<T> filt = kfr::iir_highpass(kfr::butterworth<T>(order),cutoff,sample_rate);        
        sample_vector<kfr::biquad_params<T>> bqs = kfr::to_sos<T>(filt);
        std::vector<filter_coefficients<T>> r(bqs.size());
        for(size_t i = 0; i < bqs.size(); i)
        {            
            r[i].z[0] = bqs[i].b0;
            r[i].z[1] = bqs[i].b1;
            r[i].z[2] = bqs[i].b2;
            r[i].p[0] = bqs[i].a0;
            r[i].p[1] = bqs[i].a1;
            r[i].p[2] = bqs[i].a2;
        }            
        return r;
    }

    template<typename T>
    std::vector<filter_coefficients<T>>
    butterworth_bandstop(int order, T locutoff, T hicutoff, T sample_rate) {
        kfr::zpk<T> filt = kfr::iir_bandstop(kfr::butterworth<T>(order),locutoff,hicutoff,sample_rate);        
        sample_vector<kfr::biquad_params<T>> bqs = kfr::to_sos<T>(filt);
        std::vector<filter_coefficients<T>> r(bqs.size());
        for(size_t i = 0; i < bqs.size(); i)
        {            
            r[i].z[0] = bqs[i].b0;
            r[i].z[1] = bqs[i].b1;
            r[i].z[2] = bqs[i].b2;
            r[i].p[0] = bqs[i].a0;
            r[i].p[1] = bqs[i].a1;
            r[i].p[2] = bqs[i].a2;
        }            
        return r;
    }

    template<typename T>
    std::vector<filter_coefficients<T>>
    butterworth_bandpass(int order, T locutoff, T hicutoff, T sample_rate) {
        kfr::zpk<T> filt = kfr::iir_bandpass(kfr::butterworth<T>(order),locutoff,hicutoff,sample_rate);        
        sample_vector<kfr::biquad_params<T>> bqs = kfr::to_sos<T>(filt);
        std::vector<filter_coefficients<T>> r(bqs.size());
        for(size_t i = 0; i < bqs.size(); i)
        {            
            r[i].z[0] = bqs[i].b0;
            r[i].z[1] = bqs[i].b1;
            r[i].z[2] = bqs[i].b2;
            r[i].p[0] = bqs[i].a0;
            r[i].p[1] = bqs[i].a1;
            r[i].p[2] = bqs[i].a2;
        }            
        return r;
    }

    template<typename T>
    std::vector<filter_coefficients<T>>
    cheby1_lowpass(int order, T w, T cutoff, T sample_rate) {
        kfr::zpk<T> filt = kfr::iir_lowpass(kfr::chebyshev1<T>(order,w),cutoff,sample_rate);        
        sample_vector<kfr::biquad_params<T>> bqs = kfr::to_sos<T>(filt);
        std::vector<filter_coefficients<T>> r(bqs.size());        
        for(size_t i = 0; i < bqs.size(); i)
        {            
            r[i].z[0] = bqs[i].b0;
            r[i].z[1] = bqs[i].b1;
            r[i].z[2] = bqs[i].b2;
            r[i].p[0] = bqs[i].a0;
            r[i].p[1] = bqs[i].a1;
            r[i].p[2] = bqs[i].a2;
        }            
        return r;
    }

    template<typename T>
    std::vector<filter_coefficients<T>>
    cheby1_highpass(int order, T w, T cutoff, T sample_rate) {
        kfr::zpk<T> filt = kfr::iir_highpass(kfr::chebyshev1<T>(order,w),cutoff,sample_rate);        
        sample_vector<kfr::biquad_params<T>> bqs = kfr::to_sos<T>(filt);
        std::vector<filter_coefficients<T>> r(bqs.size());
        for(size_t i = 0; i < bqs.size(); i)
        {            
            r[i].z[0] = bqs[i].b0;
            r[i].z[1] = bqs[i].b1;
            r[i].z[2] = bqs[i].b2;
            r[i].p[0] = bqs[i].a0;
            r[i].p[1] = bqs[i].a1;
            r[i].p[2] = bqs[i].a2;
        }            
        return r;
    }

    template<typename T>
    std::vector<filter_coefficients<T>>
    cheby1_bandstop(int order, T w, T locutoff, T hicutoff, T sample_rate) {
        kfr::zpk<T> filt = kfr::iir_bandstop(kfr::chebyshev1<T>(order,w),locutoff,hicutoff,sample_rate);        
        sample_vector<kfr::biquad_params<T>> bqs = kfr::to_sos<T>(filt);
        std::vector<filter_coefficients<T>> r(bqs.size());
        for(size_t i = 0; i < bqs.size(); i)
        {            
            r[i].z[0] = bqs[i].b0;
            r[i].z[1] = bqs[i].b1;
            r[i].z[2] = bqs[i].b2;
            r[i].p[0] = bqs[i].a0;
            r[i].p[1] = bqs[i].a1;
            r[i].p[2] = bqs[i].a2;
        }            
        return r;
    }

    template<typename T>
    std::vector<filter_coefficients<T>>
    cheby1_bandpass(int order, T w, T locutoff, T hicutoff, T sample_rate) {
        kfr::zpk<T> filt = kfr::iir_bandpass(kfr::chebyshev1<T>(order,w),locutoff,hicutoff,sample_rate);        
        sample_vector<kfr::biquad_params<T>> bqs = kfr::to_sos<T>(filt);
        std::vector<filter_coefficients<T>> r(bqs.size());
        for(size_t i = 0; i < bqs.size(); i)
        {            
            r[i].z[0] = bqs[i].b0;
            r[i].z[1] = bqs[i].b1;
            r[i].z[2] = bqs[i].b2;
            r[i].p[0] = bqs[i].a0;
            r[i].p[1] = bqs[i].a1;
            r[i].p[2] = bqs[i].a2;
        }            
        return r;
    }

    template<typename T>
    std::vector<filter_coefficients<T>>
    cheby2_lowpass(int order, T w, T cutoff, T sample_rate) {
        kfr::zpk<T> filt = kfr::iir_lowpass(kfr::chebyshev2<T>(order,w),cutoff,sample_rate);        
        sample_vector<kfr::biquad_params<T>> bqs = kfr::to_sos<T>(filt);
        std::vector<filter_coefficients<T>> r(bqs.size());        
        for(size_t i = 0; i < bqs.size(); i)
        {            
            r[i].z[0] = bqs[i].b0;
            r[i].z[1] = bqs[i].b1;
            r[i].z[2] = bqs[i].b2;
            r[i].p[0] = bqs[i].a0;
            r[i].p[1] = bqs[i].a1;
            r[i].p[2] = bqs[i].a2;
        }            
        return r;
    }

    template<typename T>
    std::vector<filter_coefficients<T>>
    cheby2_highpass(int order, T w, T cutoff, T sample_rate) {
        kfr::zpk<T> filt = kfr::iir_highpass(kfr::chebyshev2<T>(order,w),cutoff,sample_rate);        
        sample_vector<kfr::biquad_params<T>> bqs = kfr::to_sos<T>(filt);
        std::vector<filter_coefficients<T>> r(bqs.size());
        for(size_t i = 0; i < bqs.size(); i)
        {            
            r[i].z[0] = bqs[i].b0;
            r[i].z[1] = bqs[i].b1;
            r[i].z[2] = bqs[i].b2;
            r[i].p[0] = bqs[i].a0;
            r[i].p[1] = bqs[i].a1;
            r[i].p[2] = bqs[i].a2;
        }            
        return r;
    }

    template<typename T>
    std::vector<filter_coefficients<T>>
    cheby2_bandstop(int order, T w, T locutoff, T hicutoff, T sample_rate) {
        kfr::zpk<T> filt = kfr::iir_bandstop(kfr::chebyshev2<T>(order,w),locutoff,hicutoff,sample_rate);        
        sample_vector<kfr::biquad_params<T>> bqs = kfr::to_sos<T>(filt);
        std::vector<filter_coefficients<T>> r(bqs.size());
        for(size_t i = 0; i < bqs.size(); i)
        {            
            r[i].z[0] = bqs[i].b0;
            r[i].z[1] = bqs[i].b1;
            r[i].z[2] = bqs[i].b2;
            r[i].p[0] = bqs[i].a0;
            r[i].p[1] = bqs[i].a1;
            r[i].p[2] = bqs[i].a2;
        }            
        return r;
    }

    template<typename T>
    std::vector<filter_coefficients<T>>
    cheby2_bandpass(int order, T w, T locutoff, T hicutoff, T sample_rate) {
        kfr::zpk<T> filt = kfr::iir_bandpass(kfr::chebyshev2<T>(order,w),locutoff,hicutoff,sample_rate);        
        sample_vector<kfr::biquad_params<T>> bqs = kfr::to_sos<T>(filt);
        std::vector<filter_coefficients<T>> r(bqs.size());
        for(size_t i = 0; i < bqs.size(); i)
        {            
            r[i].z[0] = bqs[i].b0;
            r[i].z[1] = bqs[i].b1;
            r[i].z[2] = bqs[i].b2;
            r[i].p[0] = bqs[i].a0;
            r[i].p[1] = bqs[i].a1;
            r[i].p[2] = bqs[i].a2;
        }            
        return r;
    }
    
    template<typename T>
    sample_vector<T> fir_lowpass_window_taps(size_t taps, T fc, T alpha=3.0)
    {
        sample_vector<T> r(taps);             
        auto window = kfr::to_handle(kfr::window_kaiser(taps,alpha));
        fir_lowpass(r,fc,window,true);
        return r;
    }
    template<typename T>
    sample_vector<T> fir_lowpass_hann_taps(size_t taps, T fc)
    {
        sample_vector<T> r(taps);
        auto window = kfr::to_handle(kfr::window_hann(taps));
        fir_lowpass(r,fc,window,true);
        return r;
    }
    template<typename T>
    sample_vector<T> fir_lowpass_hamming_taps(size_t taps, T fc)
    {
        sample_vector<T> r(taps);
        auto window = kfr::to_handle(kfr::window_hamming(taps));
        fir_lowpass(r,fc,window,true);
        return r;
    }
    template<typename T>
    sample_vector<T> fir_lowpass_blackman_taps(size_t taps, T fc)
    {
        sample_vector<T> r(taps);
        auto window = kfr::to_handle(kfr::window_blackman(taps));
        fir_lowpass(r,fc,window,true);
        return r;
    }
    template<typename T>
    sample_vector<T> fir_lowpass_rectangular_taps(size_t taps, T fc)
    {
        sample_vector<T> r(taps);
        auto window = kfr::to_handle(kfr::window_rectangular(taps));
        fir_lowpass(r,fc,window,true);
        return r;
    }
    template<typename T>
    sample_vector<T> fir_lowpass_triangular_taps(size_t taps, T fc)
    {
        sample_vector<T> r(taps);
        auto window = kfr::to_handle(kfr::window_triangular(taps));
        fir_lowpass(r,fc,window,true);
        return r;
    }

    template<typename T>
    sample_vector<T> fir_highpass_window_taps(size_t taps, T fc, T alpha=3.0)
    {
        sample_vector<T> r(taps);
        auto window = kfr::to_handle(kfr::window_kaiser(taps,alpha));
        fir_highpass(r,fc,window,true);
        return r;
    }
    template<typename T>
    sample_vector<T> fir_highpass_hann_taps(size_t taps, T fc)
    {
        sample_vector<T> r(taps);
        auto window = kfr::to_handle(kfr::window_hann(taps));
        fir_highpass(r,fc,window,true);
        return r;
    }
    template<typename T>
    sample_vector<T> fir_highpass_hamming_taps(size_t taps, T fc)
    {
        sample_vector<T> r(taps);
        auto window = kfr::to_handle(kfr::window_hamming(taps));
        fir_highpass(r,fc,window,true);
        return r;
    }
    template<typename T>
    sample_vector<T> fir_highpass_blackman_taps(size_t taps, T fc)
    {
        sample_vector<T> r(taps);
        auto window = kfr::to_handle(kfr::window_blackman(taps));
        fir_highpass(r,fc,window,true);
        return r;
    }
    template<typename T>
    sample_vector<T> fir_highpass_rectangular_taps(size_t taps, T fc)
    {
        sample_vector<T> r(taps);
        auto window = kfr::to_handle(kfr::window_rectangular(taps));
        fir_highpass(r,fc,window,true);
        return r;
    }
    template<typename T>
    sample_vector<T> fir_highpass_triangular_taps(size_t taps, T fc)
    {
        sample_vector<T> r(taps);
        auto window = kfr::to_handle(kfr::window_triangular(taps));
        fir_highpass(r,fc,window,true);
        return r;
    }

    template<typename T>
    sample_vector<T> fir_bandpass_window_taps(size_t taps, T fc1, T fc2, T alpha=3.0)
    {
        sample_vector<T> r(taps);
        auto window = kfr::to_handle(kfr::window_kaiser(taps,alpha));
        fir_bandpass(r,fc1,fc2,window,true);
        return r;
    }
    template<typename T>
    sample_vector<T> fir_bandpass_hann_taps(size_t taps, T fc1, T fc2)
    {
        sample_vector<T> r(taps);
        auto window = kfr::to_handle(kfr::window_hann(taps));
        fir_bandpass(r,fc1,fc2,window,true);
        return r;
    }
    template<typename T>
    sample_vector<T> fir_bandpass_hamming_taps(size_t taps, T fc1, T fc2)
    {
        sample_vector<T> r(taps);
        auto window = kfr::to_handle(kfr::window_hamming(taps));
        fir_bandpass(r,fc1,fc2,window,true);
        return r;
    }
    template<typename T>
    sample_vector<T> fir_bandpass_blackman_taps(size_t taps, T fc1, T fc2)
    {
        sample_vector<T> r(taps);
        auto window = kfr::to_handle(kfr::window_blackman(taps));
        fir_bandpass(r,fc1,fc2,window,true);
        return r;
    }
    template<typename T>
    sample_vector<T> fir_bandpass_rectangular_taps(size_t taps, T fc1, T fc2)
    {
        sample_vector<T> r(taps);
        auto window = kfr::to_handle(kfr::window_rectangular(taps));
        fir_bandpass(r,fc1,fc2,window,true);
        return r;
    }
    template<typename T>
    sample_vector<T> fir_bandpass_triangular_taps(size_t taps, T fc1, T fc2)
    {
        sample_vector<T> r(taps);
        auto window = kfr::to_handle(kfr::window_triangular(taps));
        fir_bandpass(r,fc1,fc2,window,true);
        return r;
    }

    template<typename T>
    sample_vector<T> fir_bandstop_window_taps(size_t taps, T fc1, T fc2, T alpha=3.0)
    {
        sample_vector<T> r(taps);
        auto window = kfr::to_handle(kfr::window_kaiser(taps,alpha));
        fir_bandstop(r,fc1,fc2,window,true);
        return r;
    }
    template<typename T>
    sample_vector<T> fir_bandstop_hann_taps(size_t taps, T fc1, T fc2)
    {
        sample_vector<T> r(taps);
        auto window = kfr::to_handle(kfr::window_hann(taps));
        fir_bandstop(r,fc1,fc2,window,true);
        return r;
    }
    template<typename T>
    sample_vector<T> fir_bandstop_hamming_taps(size_t taps, T fc1, T fc2)
    {
        sample_vector<T> r(taps);
        auto window = kfr::to_handle(kfr::window_hamming(taps));
        fir_bandstop(r,fc1,fc2,window,true);
        return r;
    }
    template<typename T>
    sample_vector<T> fir_bandstop_blackman_taps(size_t taps, T fc1, T fc2)
    {
        sample_vector<T> r(taps);
        auto window = kfr::to_handle(kfr::window_blackman(taps));
        fir_bandstop(r,fc1,fc2,window,true);
        return r;
    }
    template<typename T>
    sample_vector<T> fir_bandstop_rectangular_taps(size_t taps, T fc1, T fc2)
    {
        sample_vector<T> r(taps);
        auto window = kfr::to_handle(kfr::window_rectangular(taps));
        fir_bandstop(r,fc1,fc2,window,true);
        return r;
    }
    template<typename T>
    sample_vector<T> fir_bandstop_triangular_taps(size_t taps, T fc1, T fc2)
    {
        sample_vector<T> r(taps);
        auto window = kfr::to_handle(kfr::window_triangular(taps));
        fir_bandstop(r,fc1,fc2,window,true);
        return r;
    }

    template<typename T>
    struct Resampler 
    {
        kfr::samplerate_converter<T> * sc;

        using quality = kfr::sample_rate_conversion_quality;

        Resampler(quality q, int64_t interp_factor, int64_t decimation_factor, T scale=1.0, T cutoff=0.5)
        {
            sc = new kfr::samplerate_converter<T>(q,interp_factor,decimation_factor,scale,cutoff);
        }
        ~Resampler()
        {
            if(sc) delete sc;
        }

        size_t ProcessBlock(size_t n, T * in, T * out)
        {
            kfr::univector<T> I(n);
            kfr::univector<T> O(n);
            memcpy(I.data(),in,n*sizeof(T));
            memcpy(O.data(),out,n*sizeof(T));
            size_t r = sc->process(O,I);
            memcpy(out,O.data(),n*sizeof(T));
            return r;
        }
        size_t process(sample_vector<T> & in, sample_vector<T> & out)
        {
            return sc->process(in,out);
        }
    };

    template<typename T>
    size_t resample(Resampler<T>& r, sample_vector<T> & in, sample_vector<T> & out)
    {
        return r.process(in,out);
    }      
    template<typename T>
    size_t resample(Resampler<T> & r, size_t n, T * in, T * out)
    {
        return r.ProcessBlock(n,in,out);
    }  

}