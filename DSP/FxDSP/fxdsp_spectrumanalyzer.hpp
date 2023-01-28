#pragma once

#include <cstddef>
#include <cstdlib>
#include "fxdsp.hpp"

namespace FXDSP
{
    template<typename T>
    static inline void
    calculate_bin_frequencies(T* dest, unsigned fft_length, T sample_rate)
    {
        T freq_step = sample_rate / fft_length;
        #pragma omp simd
        for(unsigned bin = 0; bin < fft_length / 2; ++bin)
        {
            dest[bin] = bin * freq_step;
        }
    }

    template<typename T, class FFT>
    struct SpectrumAnalyzer
    {
        unsigned        fft_length;
        unsigned        bins;
        T           sample_rate;
        T           mag_sum;
        T*          frequencies;
        T*          real;
        T*          imag;
        T*          mag;
        T*          phase;
        T*          root_moment;
        FFT*         fft;
        Window_t            window_type;
        WindowFunction<T>*  window;

        SpectrumAnalyzer(unsigned fft_length, T sample_rate)
        {
            fft = new FFT(fft_length);
            window = new WindowFunction<T>(fft_length, BLACKMAN);

            frequencies = new T[fft_length/2];
            real = new T[fft_length/2];
            imag = new T[fft_length/2];
            mag = new T[fft_length/2];
            phase = new T[fft_length/2];
            root_moment = new T[fft_length/2];

            if ((NULL != window) && (NULL != fft) && (NULL != frequencies) && (NULL != real) \
                && (NULL != imag) && (NULL != mag) && (NULL != phase) && (NULL != root_moment))
            {
                this->fft_length = fft_length;
                this->bins = fft_length / 2;
                this->sample_rate = sample_rate;
                this->mag_sum = 0.0;                                
                this->window_type = BLACKMAN;
                calculate_bin_frequencies(frequencies, fft_length, sample_rate);
                *(root_moment) = 0.0;
            }
            else
            {
                throw std::runtime_error("Out of memory");
            }
        }

        void analyze(T * signal)
        {
            T scratch[fft_length];
            window->Process(scratch, signal, fft_length);
            fft->R2C(scratch, real, imag);
            VectorRectToPolar(mag, phase, real, imag, bins);
            mag_sum = VectorSum(mag, bins);
            root_moment[0] = 0.0;
        }
        T centroid() {
            T num[bins];
            VectorVectorMultiply(num, mag, frequencies, bins);
            return VectorSum(num, bins) / mag_sum;
        }
        T spread() {
            T mu = centroid(analyzer);
            T num[bins];
            if (root_moment[0] == 0.0)
            {
                VectorScalarAdd(root_moment, frequencies, -mu, bins);
            }
            VectorPower(num, root_moment, 2, bins);
            return VectorSum(num, bins) / mag_sum;
        }
        T skewness() {
            T mu = centroid(analyzer);
            T num[bins];
            if (root_moment[0] == 0.0)
            {
                VectorScalarAdd(root_moment, frequencies, -mu, bins);
            }
            VectorPower(num, root_moment, 3, bins);
            return VectorSum(num, bins) / mag_sum;
        }
        T kurtosis() {
            T mu = centroid(analyzer);
            T num[bins];
            if (root_moment[0] == 0.0)
            {
                VectorScalarAdd(root_moment, frequencies, -mu, bins);
            }
            VectorPower(num, root_moment, 4, bins);
            return VectorSum(num, bins) / mag_sum;
        }
    };

    using SpectrumAnalyzerFloat = SpectrumAnalyzer<float,FFTFloat>;
    using SpectrumAnalyzerDouble= SpectrumAnalyzer<double,FFTDouble>;
}