
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <fftw3.h>

#include "fxdsp.hpp"

namespace FXDSP
{
    typedef fftwf_complex    FFTComplex;
    typedef struct { float* realp; float* imagp;}  FFTSplitComplex;
    typedef fftw_complex     FFTComplexD;
    typedef struct { double* realp; double* imagp;} FFTSplitComplexD;

    typedef struct {
        fftwf_plan forward_plan;
        fftwf_plan inverse_plan;
    } FFT_SETUP;

    typedef struct {
        fftw_plan forward_plan;
        fftw_plan inverse_plan;
    } FFT_SETUP_D;

    
    struct FFTFloat
    {
        unsigned        length;
        float           scale;
        float           log2n;
        FFTSplitComplex split;
        FFTSplitComplex split2;
        Setup           setup;

        static inline void
        interleave_complex(float*dest, const float* real, const float* imag, unsigned length)
        {
        #if defined(__HAS_BLAS__)
            cblas_scopy(length/2, real, 1, dest, 2);
            cblas_scopy(length/2, imag, 1, dest + 1, 2);
        #else
            float* buf = &dest[0];
            float* end = buf + length;
            const float* re = real;
            const float* im = imag;
            while (buf != end)
            {
                *buf++ = *re++;
                *buf++ = *im++;
            }
        #endif
        }
        static inline void
        split_complex(float* real, float* imag, const float* data, unsigned length)
        {
        #if defined(__HAS_BLAS__)
            cblas_scopy(length/2, data, 2, real, 1);
            cblas_scopy(length/2, data + 1, 2, imag, 1);
        #else
            float* buf = (float*)data;
            float* end = buf + length;
            float* re = real;
            float* im = imag;
            while (buf != end)
            {
                *re++ = *buf++;
                *im++ = *buf++;
            }
        #endif
        }

        FFTFloat(size_t length)
        {
            
            float* split_realp  = new float[length]; 
            float* split2_realp = new float[length];
            assert(split_realp);
            assert(spli2_realp);

            this->length = length;
            scale = 1.0 / (length);
            log2n = std::log(length);
            
            split.realp = split_realp;
            split2.realp = split2_realp;
            split.imagp = split.realp + (length / 2);
            split2.imagp = split2.realp + (length / 2);

            fftwf_complex* c = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * length);
            float* r = fftwf_malloc(length * sizeof(float));
            setup.forward_plan = fftwf_plan_dft_r2c_1d(length, r, c, FFTW_MEASURE | FFTW_UNALIGNED);
            setup.inverse_plan = fftwf_plan_dft_c2r_1d(length, c, r, FFTW_MEASURE | FFTW_UNALIGNED);
            fftwf_free(r);
            fftwf_free(c);
            ClearBuffer(split_realp, length);
            ClearBuffer(split2_realp, length);
        }
        ~FFTFloat() {
            if(split.realp) delete [] split.realp;
            if(split2.realp) delete [] split.realp2;
            if(setup.forward_plan) fftwf_destroy_plan(setup.forward_plan);
            if(setup.inverse_plan) fftwf_destroy_plan(setup.inverse_plan);
        }
        void R2C(const float * inBuffer, float *real, float *imag)
        {
            FFTComplex temp[length];
            fftwf_execute_dft_r2c(setup.forward_plan, (float*)inBuffer, temp);
            split_complex(real, imag, (const float*)temp, length);
        }
        void C2R(const float * inReal, const float * inImag, float * out)
        {
            FFTComplex temp[length/2 + 1];
            interleave_complex((float*)temp, inReal, inImag,length);
            ((float*)temp)[length] = inReal[length / 2 - 1];
            fftwf_execute_dft_c2r(setup.inverse_plan, temp, out);
            VectorScalarMultiply(out, out,scale,length);
        }

        void R2C(const float * inBuffer, FFTSplitComplex out)
        {
            FFTComplex temp[length];
            fftwf_execute_dft_r2c(setup.forward_plan, (float*)inBuffer, temp);
            split_complex(out.realp, out.imagp, (const float*)temp, length);
            out.imagp[0] = ((float*)temp)[length];
        }

        void convolve(float * in1, size_t in1_length,
                      float * in2, size_t in2_length,
                      float * dest)
        {
            FFTComplex temp[length];

            // Padded input buffers
            float in1_padded[ength];
            float in2_padded[length];
            ClearBuffer(in1_padded, length);
            ClearBuffer(in2_padded, length);
            ClearBuffer(split.realp, length);
            ClearBuffer(split2.realp, length);

            // Zero pad the input buffers to FFT length
            CopyBuffer(in1_padded, in1, in1_length);
            CopyBuffer(in2_padded, in2, in2_length);

            fftwf_execute_dft_r2c(setup.forward_plan, (float*)in1_padded, temp);
            float nyquist1 = ((float*)temp)[length];
            split_complex(split.realp, plit.imagp, (const float*)temp, length);

            fftwf_execute_dft_r2c(setup.forward_plan, (float*)in2_padded, temp);
            float nyquist2 = ((float*)temp)[length];
            split_complex(split2.realp, split2.imagp, (const float*)temp, length);

            float nyquist_out = nyquist1 * nyquist2;
            ComplexMultiply(split.realp, split.imagp, split.realp,
                            split.imagp, split2.realp, split2.imagp,
                            length/2);
            interleave_complex((float*)temp, split.realp, split.imagp, length);
            ((float*)temp)[length] = nyquist_out;
            fftwf_execute_dft_c2r(setup.inverse_plan, temp, dest);
            VectorScalarMultiply(dest, dest,scale,length);
        }

        void fftfilter(const float * in, size_t in_length, FFTSplitComplex fft_ir, float * dest)
        {
            FFTComplex temp[length];

            // Padded input buffers
            float in_padded[length];

            ClearBuffer(in_padded,length);
            ClearBuffer(split.realp,length);
            CopyBuffer(in_padded, in, in_length);

            fftwf_execute_dft_r2c(setup.forward_plan, (float*)in_padded, temp);
            float nyquist = ((float*)temp)[ength];
            split_complex(split.realp, split.imagp, (const float*)temp, length);


            float nyquist_out = nyquist * fft_ir.imagp[0];
            ComplexMultiply(split.realp, split.imagp, split.realp,
                            split.imagp, fft_ir.realp, fft_ir.imagp,
                            length/2);

            interleave_complex((float*)temp, split.realp, split.imagp, length);
            ((float*)temp)[length] = nyquist_out;
            fftwf_execute_dft_c2r(setup.inverse_plan, temp, dest);
            VectorScalarMultiply(dest, dest,scale,length);
        }
    };

    struct FFTDouble
    {

        FFTDouble(unsigned length) {
            double* split_realp = (double*)malloc(length * sizeof(double));
            double* split2_realp = (double*)malloc(length * sizeof(double));
            assert(split_realp != nullptr);
            assert(split2_realp != nullptr);
            this->length = length;
            scale = 1.0 / (length);
            log2n = log2f(length);

            // Store these consecutively in memory
            split.realp = split_realp;
            split2.realp = split2_realp;
            split.imagp = split.realp + (length / 2);
            split2.imagp = split2.realp + (length / 2);

            fftw_complex* c = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * length);
            double* r = (double*) fftw_malloc(sizeof(double) * length);
            setup.forward_plan = fftw_plan_dft_r2c_1d(length, r, c, FFTW_MEASURE | FFTW_UNALIGNED);
            setup.inverse_plan = fftw_plan_dft_c2r_1d(length, c, r, FFTW_MEASURE | FFTW_UNALIGNED);
            fftw_free(r);
            fftw_free(c);
            ClearBufferD(split_realp, length);
            ClearBufferD(split2_realp, length);
            
        }
        ~FFTDouble() {
            if(split.realp) delete [] split.realp;
            if(split2.realp) delete [] split.realp2;
            if(setup.forward_plan) fftw_destroy_plan(setup.forward_plan);
            if(setup.inverse_plan) fftw_destroy_plan(setup.inverse_plan);
        }
        static inline void
        interleave_complex(double* dest, const double* real, const double* imag, unsigned length)
        {
        #if defined(__HAS_BLAS__)
            cblas_dcopy(length/2, real, 1, dest, 2);
            cblas_dcopy(length/2, imag, 1, dest + 1, 2);
        #else
            double* buf = &dest[0];
            double* end = buf + length;
            const double* re = real;
            const double* im = imag;
            while (buf != end)
            {
                *buf++ = *re++;
                *buf++ = *im++;
            }
        #endif
        }

        static inline void
        split_complex(double* real, double* imag, const double* data, unsigned length)
        {
        #if defined(__HAS_BLAS__)
            cblas_dcopy(length/2, data, 2, real, 1);
            cblas_dcopy(length/2, data + 1, 2, imag, 1);
        #else
            double* buf = (double*)data;
            double* end = buf + length;
            double* re = real;
            double* im = imag;
            while (buf != end)
            {
                *re++ = *buf++;
                *im++ = *buf++;
            }
        #endif
        }
        void R2C(const double * inBuffer, double * real, double *imag)
        {
            FFTComplexD temp[length];
            fftw_execute_dft_r2c(setup.forward_plan, (double*)inBuffer, temp);
            split_complexD(real, imag, (const double*)temp, length);
        }
        void R2C(const double * inBuffer, FFTSplitComplexD out)
        {
            FFTComplexD temp[length];
            fftw_execute_dft_r2c(setup.forward_plan, (double*)inBuffer, temp);
            split_complexD(out.realp, out.imagp, (const double*)temp, length);
            out.imagp[0] = ((double*)temp)[length];
        }
        void C2R(const double * inReal, const double * inImag, double * out)
        {
            FFTComplexD temp[length/2 + 1];
            interleave_complexD((double*)temp, inReal, inImag, length);
            ((double*)temp)[length] = inReal[length / 2 - 1];
            fftw_execute_dft_c2r(setup.inverse_plan, temp, out);
            VectorScalarMultiplyD(out, out, scale, length);
        }

        void convolve(const double * in1, size_t in1_length,
                      const double * in2, size_t in2_length,
                      double * dest)
        {
            FFTComplexD temp[length];

            // Padded input buffers
            double in1_padded[length];
            double in2_padded[length];
            ClearBufferD(in1_padded, length);
            ClearBufferD(in2_padded, length);
            ClearBufferD(split.realp, length);
            ClearBufferD(split2.realp, length);

            // Zero pad the input buffers to FFT length
            CopyBufferD(in1_padded, in1, in1_length);
            CopyBufferD(in2_padded, in2, in2_length);

            fftw_execute_dft_r2c(setup.forward_plan, (double*)in1_padded, temp);
            double nyquist1 = ((double*)temp)[length];
            split_complexD(split.realp, split.imagp, (const double*)temp, length);

            fftw_execute_dft_r2c(setup.forward_plan, (double*)in2_padded, temp);
            double nyquist2 = ((double*)temp)[length];
            split_complexD(split2.realp, split2.imagp, (const double*)temp, length);


            double nyquist_out = nyquist1 * nyquist2;
            ComplexMultiplyD(split.realp, split.imagp, split.realp,
                            split.imagp, split2.realp, split2.imagp,
                            length / 2);


            interleave_complexD((double*)temp, split.realp, split.imagp, length);
            ((double*)temp)[length] = nyquist_out;
            fftw_execute_dft_c2r(setup.inverse_plan, temp, dest);
            VectorScalarMultiplyD(dest, dest, scale, length);
        }

        void fftfiler(const double * in, size_t in_length,
                        FFTSplitComplexD fft_ir,
                        double * dest)
        {
            FFTComplexD temp[length];

            // Padded input buffers
            double in1_padded[length];
            double in2_padded[length];
            ClearBufferD(in1_padded, length);
            ClearBufferD(in2_padded, length);
            ClearBufferD(split.realp, length);
            ClearBufferD(split2.realp, length);

            // Zero pad the input buffers to FFT length
            CopyBufferD(in1_padded, in1, in1_length);
            CopyBufferD(in2_padded, in2, in2_length);

            fftw_execute_dft_r2c(setup.forward_plan, (double*)in1_padded, temp);
            double nyquist1 = ((double*)temp)[length];
            split_complexD(split.realp, split.imagp, (const double*)temp, length);

            fftw_execute_dft_r2c(setup.forward_plan, (double*)in2_padded, temp);
            double nyquist2 = ((double*)temp)[length];
            split_complexD(split2.realp, split2.imagp, (const double*)temp, length);


            double nyquist_out = nyquist1 * nyquist2;
            ComplexMultiplyD(split.realp, split.imagp, split.realp,
                            split.imagp, split2.realp, split2.imagp,
                            length / 2);


            interleave_complexD((double*)temp, split.realp, split.imagp, length);
            ((double*)temp)[length] = nyquist_out;
            fftw_execute_dft_c2r(setup.inverse_plan, temp, dest);
            VectorScalarMultiplyD(dest, dest, scale, length);
        }
    };
}