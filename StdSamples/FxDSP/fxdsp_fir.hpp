#pragma once

#include <cstddef>
#include <cstdlib>

#include "fxdsp.hpp"
#include "fxdsp_fft.hpp"


namespace FXDSP
{

/** The kernel length at which to use FFT convolution vs direct */
/* So this is a pretty good value for now */
#define USE_FFT_CONVOLUTION_LENGTH (128)


/** Convolution Algorithm to use */
typedef enum _ConvolutionMode
{
    /** Choose the best algorithm based on filter size */
    BEST    = 0,

    /** Use direct convolution */
    DIRECT  = 1,

    /** Use FFT Convolution (Better for longer filter kernels */
    FFT     = 2

} ConvolutionMode_t;


template<typename T, class FFT=FFTFloat, class Split=FFTSplitComplex>
struct FIRFilter
{
    T*                  kernel;
    const T*            kernel_end;
    T*                  overlap;
    unsigned            kernel_length;
    unsigned            overlap_length;
    ConvolutionMode_t   conv_mode;
    FFT*                fft_config;
    Spli                fft_kernel;
    unsigned            fft_length;

    FIRFilter(const T* filter_kernel, unsigned length, ConvolutionMode_t convolution_mode)
    {
            // Array lengths and sizes
        unsigned kernel_length = length;                    // IN SAMPLES!
        unsigned overlap_length = kernel_length - 1;        // IN SAMPLES!
        
        kernel  = new T[kernel_length];
        overlap = new T[overlap_length]; 
        
        // Initialize Buffers
        CopyBuffer(kernel, filter_kernel, kernel_length);
        ClearBuffer(overlap, overlap_length);

        // Set up the struct
        kernel = kernel;
        kernel_end = filter_kernel + (kernel_length - 1);
        overlap = overlap;
        kernel_length = kernel_length;
        overlap_length = overlap_length;
        fft_config = NULL;
        fft_kernel.realp = NULL;
        fft_kernel.imagp = NULL;

        if (((convolution_mode == BEST) &&
            (kernel_length < USE_FFT_CONVOLUTION_LENGTH)) ||
            (convolution_mode == DIRECT))
        {
            conv_mode = DIRECT;
        }

        else
        {
            conv_mode = FFT;
        }
    }
    ~FIRFilter() {
        if(kernel) delete [] kernel;
        if(overlap) delete [] overlap;
    }

    void flush() { 
            ClearBuffer(overlap, overlap_length);
    }
    void setKernel(T * filter_kernel)
    {
        CopyBuffer(kernel,filter_kernel,kernel_length);
    }
    void ProcessBlock(size_t n, T * inBuffer, T * outBuffer)
    {
        // Do direct convolution
        if (conv_mode == DIRECT)
        {
            unsigned resultLength = n_samples + (kernel_length - 1);
            // Temporary buffer to store full result of filtering..
            T buffer[resultLength];

            Convolve((T*)inBuffer, n_samples,
                            kernel, kernel_length, buffer);

            // Add in the overlap from the last block
            VectorVectorAdd(buffer, overlap, buffer, overlap_length);
            CopyBuffer(overlap, buffer + n_samples, overlap_length);
            CopyBuffer(outBuffer, buffer, n_samples);
        }
        // Otherwise do FFT Convolution
        else
        {
            // Configure the FFT on the first run, that way we can figure out how
            // long the input blocks are going to be. This makes the filter more
            // complicated internally in order to make the convolution transparent.
            // Calculate length of FFT
            if(fft_config == 0)
            {
                // Calculate FFT Length
                fft_length = next_pow2(n_samples + kernel_length - 1);
                fft_config = new FFT(fft_length);

                // fft kernel buffers
                T padded_kernel[fft_length];

                // Allocate memory for filter kernel
                fft_kernel.realp = (T*) malloc(fft_length * sizeof(T));
                fft_kernel.imagp = fft_kernel.realp +(fft_length / 2);

                // Write zero padded kernel to buffer
                CopyBuffer(padded_kernel, kernel, kernel_length);
                ClearBuffer((padded_kernel + kernel_length), (fft_length - kernel_length));

                // Calculate FFT of filter kernel
                fft_config->R2C(padded_kernel, fft_kernel);
            }

            // Buffer for transformed input
            T buffer[fft_length];

            // Convolve
            fft_config->convolve((T*)inBuffer, n_samples, fft_kernel, buffer);

            // Add in the overlap from the last block
            VectorVectorAdd(buffer, overlap, buffer, overlap_length);
            CopyBuffer(overlap, buffer + n_samples, overlap_length);
            CopyBuffer(outBuffer, buffer, n_samples);
        }            
    }
};

using FIRFilterFloat    = FIRFilter<>;
using FIRFilterDouble   = FIRFilter<double,FFTDouble,FFTSplitComplexD>;
}