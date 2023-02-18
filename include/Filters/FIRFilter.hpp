#pragma once

#include <stdlib.h>
#include <cmath>
#include <string.h>
#include <string.h>
#include <vector>
#include <StdSamples/stdsamples.hpp>

namespace Filters::FIR
{
    class FIR_filter
    {
    public:
        FIR_filter( int taps=0, DspFloatType f1=0, DspFloatType f2=0, const char* type="", 
                    const char* window="");
        ~FIR_filter();

        AudioDSP::sample_vector<DspFloatType> getCoefficients();

        DspFloatType filter(DspFloatType new_sample);

		DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {
			return filter(I);
		}
		void ProcessSIMD(size_t n, DspFloatType * in, DspFloatType * out);
		void ProcessBlock(size_t n, DspFloatType * in, DspFloatType * out) { ProcessSIMD(n,in,out); }
		void ProcessInplace(size_t n, DspFloatType * in) { ProcessSIMD(n,in,in); }
		
    private:
        AudioDSP::sample_vector<DspFloatType> lowPass_coefficient( int taps, DspFloatType f);
        AudioDSP::sample_vector<DspFloatType> highPass_coefficient(int taps, DspFloatType f);
        AudioDSP::sample_vector<DspFloatType> bandPass_coefficient(int taps, DspFloatType f1, DspFloatType f2);
        AudioDSP::sample_vector<DspFloatType> bandStop_coefficient(int taps, DspFloatType f1, DspFloatType f2);

        AudioDSP::sample_vector<DspFloatType> window_hamming(int taps);
        AudioDSP::sample_vector<DspFloatType> window_triangle(int taps);
        AudioDSP::sample_vector<DspFloatType> window_hanning(int taps);
        AudioDSP::sample_vector<DspFloatType> window_blackman(int taps);

        AudioDSP::sample_vector<DspFloatType> h;       // FIR coefficients
        AudioDSP::sample_vector<DspFloatType> samples; // FIR delay

        int idx;        // Round robin index
        int taps;       // Number of taps of the filter
    };

    static DspFloatType sinc(const DspFloatType x)
    {
        if (x == 0)
            return 1;

        return sin(M_PI * x) / (M_PI * x);
    }

    FIR_filter::FIR_filter( int taps, DspFloatType f1, DspFloatType f2, const char* type,
                            const char* window): h(taps), samples(taps)
    {
        this->idx     = 0;
        this->taps    = taps;

        AudioDSP::sample_vector<DspFloatType> h;  // Buffer FIR coefficients
        AudioDSP::sample_vector<DspFloatType> w;  // Buffer window coefficients

        // Calculate the coefficient corresponding to the filter type
        if (!strcmp(type, "lp")) {
            h = lowPass_coefficient(taps, f1);
        }
        else if (!strcmp(type, "hp")) {
            h = highPass_coefficient(taps, f1);
        }
        else if (!strcmp(type, "bp")) {
            h = bandPass_coefficient(taps, f1, f2);
        }
        else if (!strcmp(type, "sb")) {
            h = bandStop_coefficient(taps, f1, f2);
        }

        // Calculate the window to improve the FIR filter
        if (!strcmp(window, "hamming")) {
            w = window_hamming(taps);
        }
        else if (!strcmp(window, "triangle")) {
            w = window_triangle(taps);
        }
        else if (!strcmp(window, "hanning")) {
            w = window_hanning(taps);
        }
        else if (!strcmp(window, "blackman")) {
            w = window_blackman(taps);
        }

        if (!strcmp(window, "")) {
            this->h = h;
        }
        else
        {
			DspFloatType * wx = w.data();
			DspFloatType * hx = h.data();
			#pragma omp simd aligned(hx,wx)
            for(int n = 0; n < taps; n++) {
                this->h[n] = hx[n] * wx[n];
            }
        }
    }

    FIR_filter::~FIR_filter()
    {

    }

    AudioDSP::sample_vector<DspFloatType> FIR_filter::getCoefficients()
    {
        return this->h;
    }

    AudioDSP::sample_vector<DspFloatType> FIR_filter::lowPass_coefficient(int taps, DspFloatType f)
    {
        AudioDSP::sample_vector<int>    nv(taps);
        AudioDSP::sample_vector<DspFloatType> hv(taps);
        int * n = nv.data();
        DspFloatType * h = hv.data();

        #pragma omp simd aligned(n)
        for(int i = 0; i < taps; i++) {
            n[i] = i - int(taps/2);
        }

        #pragma omp simd aligned(h,n)
        for(int i = 0; i < taps; i++) {
            h[i] = 2.0*f*sinc(2.0*f*n[i]);
        }

        return hv;
    }

    AudioDSP::sample_vector<DspFloatType> FIR_filter::highPass_coefficient(int taps, DspFloatType f)
    {
        AudioDSP::sample_vector<int>    nv(taps);
        AudioDSP::sample_vector<DspFloatType> hv(taps);
		int * n = nv.data();
        DspFloatType * h = hv.data();
        
        #pragma omp simd aligned(n)
        for(int i = 0; i < taps; i++) {
            n[i] = i - int(taps/2);
        }

        #pragma omp simd aligned(h,n)
        for(int i = 0; i < taps; i++) {
            h[i] = sinc(n[i]) - 2.0*f*sinc(2.0*f*n[i]);
        }

        return hv;
    }

    AudioDSP::sample_vector<DspFloatType> FIR_filter::bandPass_coefficient(int taps, DspFloatType f1, DspFloatType f2)
    {
        AudioDSP::sample_vector<int>    nv(taps);
        AudioDSP::sample_vector<DspFloatType> hv(taps);
		int * n = nv.data();
        DspFloatType * h = hv.data();
        
        #pragma omp simd aligned(n,h)
        for(int i = 0; i < taps; i++) {
            n[i] = i - int(taps/2);
        }

        #pragma omp simd aligned(n,h)
        for(int i = 0; i < taps; i++) {
            h[i] = 2.0*f1*sinc(2.0*f1*n[i]) - 2.0*f2*sinc(2.0*f2*n[i]);
        }

        return hv;
    }

    AudioDSP::sample_vector<DspFloatType> FIR_filter::bandStop_coefficient(int taps, DspFloatType f1, DspFloatType f2)
    {
        AudioDSP::sample_vector<int>    nv(taps);
        AudioDSP::sample_vector<DspFloatType> hv(taps);
        int * n = nv.data();
        DspFloatType * h = hv.data();

        #pragma omp simd aligned(n,h)
        for(int i = 0; i < taps; i++) {
            n[i] = i - int(taps/2);
        }

        #pragma omp simd aligned(n,h)
        for(int i = 0; i < taps; i++) {
            h[i] = 2.0*f1*sinc(2.0*f1*n[i]) - 2.0*f2*sinc(2.0*f2*n[i]) + sinc(n[i]);
        }

        return hv;
    }

    AudioDSP::sample_vector<DspFloatType> FIR_filter::window_hamming(int taps)
    {        
        AudioDSP::sample_vector<DspFloatType> wv(taps);        
        DspFloatType * w = wv.data();

        DspFloatType alpha   = 0.54;
        DspFloatType beta    = 0.46;

        #pragma omp simd aligned(w)
        for(int i = 0; i < taps; i++) {
            w[i] = alpha - beta * std::cos(2.0 * M_PI * i / (taps - 1));
        }

        return wv;
    }

    AudioDSP::sample_vector<DspFloatType> FIR_filter::window_hanning(int taps)
    {
        AudioDSP::sample_vector<DspFloatType> wv(taps);		
        DspFloatType * w = wv.data();
        #pragma omp simd aligned(w)
        for(int i = 0; i < taps; i++) {
            w[i] =  std::sin(((DspFloatType) M_PI * i) / (taps - 1)) *
                    std::sin(((DspFloatType) M_PI * i) / (taps - 1));
        }

        return wv;
    }

    AudioDSP::sample_vector<DspFloatType> FIR_filter::window_triangle(int taps)
    {
        AudioDSP::sample_vector<DspFloatType> wv(taps);
		DspFloatType * w = wv.data();
        DspFloatType l = taps;

        #pragma omp simd aligned(w)
        for(int i = 0; i < taps; i++) {
            w[i] = 1 - std::abs((i - (((DspFloatType)(taps-1)) / 2.0)) / (((DspFloatType)l) / 2.0));
        }

        return wv;
    }

    AudioDSP::sample_vector<DspFloatType> FIR_filter::window_blackman(int taps)
    {
        AudioDSP::sample_vector<DspFloatType> wv(taps);
		DspFloatType * w = wv.data();
        DspFloatType alpha0 = 0.42;
        DspFloatType alpha1 = 0.5;
        DspFloatType alpha2 = 0.08;

		
        #pragma omp simd aligned(w)
        for(int i = 0; i < taps; i++) {
            w[i] = alpha0 - alpha1 * std::cos(2.0 * M_PI * i / (taps - 1))
                        - alpha2 * std::cos(4.0 * M_PI * i / (taps - 1));
        }

        return wv;
    }

    DspFloatType FIR_filter::filter(DspFloatType new_sample)
    {
        DspFloatType result = 0;

        // Save the new sample
        this->samples[this->idx] = new_sample;

        // Calculate the output
        DspFloatType * s = this->samples.data();
        DspFloatType * h = this->h.data();
        #pragma omp simd aligned(s,h)
        for(int n = 0; n < this->taps; n++)
            result += this->samples[(this->idx + n) % this->taps] * this->h[n];

        // Increase the round robin index
        this->idx = (this->idx + 1) % this->taps;

        return result;
    }
    void FIR_filter::ProcessSIMD(size_t n, DspFloatType * in, DspFloatType * out)
    {
		// Calculate the output
		DspFloatType * s = this->samples.data();
		DspFloatType * h = this->h.data();
		
		#pragma omp simd aligned(in,out,s,h)
		for(size_t i = 0; i < n; i++)
		{
			DspFloatType result = 0;
			const DspFloatType new_sample = in[i];
			// Save the new sample
			this->samples[this->idx] = new_sample;

			
			for(int n = 0; n < this->taps; n++)
				result += this->samples[(this->idx + n) % this->taps] * this->h[n];

			// Increase the round robin index
			this->idx = (this->idx + 1) % this->taps;

			out[i] = result;
		}
	}
}
