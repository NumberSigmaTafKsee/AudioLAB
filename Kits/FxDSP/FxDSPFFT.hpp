#pragma once

#include "FxDSP.hpp"

namespace DSPFX
{

    struct FxFFT 
    {
        #ifdef DSPFLOATDOUBLE
        FFTConfigD * config;
        #else
        FFTConfig * config;
        #endif

        FxFFT(size_t n) {
            #ifdef DSPFLOATDOUBLE
            config = FFTInitD(n);
            #else
            config = FFTInit(n);
            #endif
            assert(config != NULL);
        }
        ~FxFFT() {
            #ifdef DSPFLOATDOUBLE
            if(config) FFTFreeD(config);
            #else
            if(config) FFTFree(config);
            #endif
        }
        void R2C(const DspFloatType * in, DspFloatType * real, DspFloatType * imag) {
            #ifdef DSPFLOATDOUBLE
            Error_t err = FFT_R2CD(config, in, real, imag);
            #else
            Error_t err = FFT_R2C(config, in, real, imag);
            #endif
        }
        void C2R(const DspFloatType * real, const DspFloatType * imag,DspFloatType * out) {
            #ifdef DSPFLOATDOUBLE
            Error_t err = IFFT_C2RD(config, real, imag, out);
            #else
            Error_t err = IFFT_C2R(config, real, imag, out);
            #endif
        }
        void convolve(DspFloatType * in1, size_t len1, DspFloatType *in2, size_t len2, DspFloatType * out) {
            #ifdef DSPFLOATDOUBLE
            FFTConvolveD(config,in1,len1,in2,len2,out);
            #else
            FFTConvolve(config,in1,len1,in2,len2,out);
            #endif
        }
        void filterConvolve(DspFloatType * in, size_t len, DspFloatType * real, DspFloatType * imag, DspFloatType * out)
        {
            #ifdef DSPFLOATDOUBLE
            FFTSplitComplexD split;
            #else
            FFTSplitComplex split;
            #endif

            split.realp = real;
            split.imagp = imag;

            #ifdef DSPFLOATDOUBLE
            FFTFilterConvolveD(config,in,len,split,out);
            #else
            FFTFilterConvolve(config,in,len,split,out);
            #endif
        }                     
    };

}
