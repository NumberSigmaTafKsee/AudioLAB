#pragma once

#include <cmath>
#include <cstring>
#include <cfloat>

#if defined(USE_BLAS)
#include <cblas.h>
#endif

/* Scalar for converting int to float samples (1/32767.0) */
#define INT16_TO_FLOAT_SCALAR (0.00003051850947599719f)

/* 1/ln(2) */
#define INV_LN2 (1.442695040888963f)

/* 2*pi */
#define TWO_PI (6.283185307179586f)

/* pi/2 */
#define PI_OVER_TWO (1.5707963267948966f)

/* 1/(TWO_PI) */
#define INVERSE_TWO_PI (0.159154943091895f)

/* ln(10.0)/20.0 */
#define LOG_TEN_OVER_TWENTY (0.11512925464970228420089957273422)

/* 20.0/ln(10.0) */
#define TWENTY_OVER_LOG_TEN (8.6858896380650365530225783783321)

#define SQRT_TWO_OVER_TWO (0.70710678118654757273731092936941422522068023681641)

/* Utility Function Macros ****************************************************/

/* Limit value value to the range (l, u) */
#define LIMIT(value,lower,upper) ((value) < (lower) ? (lower) : \
                                 ((value) > (upper) ? (upper) : (value)))


/* Linearly interpolate between y0 and y1
  y(x) = (1 - x) * y(0) + x * y(1) | x in (0, 1)
 Function-style signature:
    float LIN_INTERP(float x, float y0, float y1)
    double LIN_INTERP(double x, double y0, double y1)
 Parameters:
    x:  x-value at which to calculate y.
    y0: y-value at x = 0.
    y1: y-value at x = 1.
 Returns:
    Interpolated y-value at specified x.
 */
#define LIN_INTERP(x, y0, y1) ((y0) + (x) * ((y1) - (y0)))


/* Convert frequency from Hz to Radians per second
 Function-style signature:
    float HZ_TO_RAD(float frequency)
    double HZ_TO_RAD(double frequency)
 Parameters:
    frequency:  Frequency in Hz.
 Returns:
    Frequency in Radians per second.
 */
#define HZ_TO_RAD(f) (TWO_PI * (f))


/* Convert frequency from Radians per second to Hz
 Function-style signature:
    float RAD_TO_HZ(float frequency)
    double RAD_TO_HZ(double frequency)
 Parameters:
 frequency:  Frequency in Radians per second.
 Returns:
 Frequency in Hz.
 */
#define RAD_TO_HZ(omega) (INVERSE_TWO_PI * (omega))


/* Fast exponentiation function
 y = e^x
 Function-style signature:
    float F_EXP(float x)
    double F_EXP(double x)
 Parameters:
    x:  Value to exponentiate.
 Returns:
    e^x.
 */
#define F_EXP(x) ((362880 + (x) * (362880 + (x) * (181440 + (x) * \
                  (60480 + (x) * (15120 + (x) * (3024 + (x) * (504 + (x) * \
                  (72 + (x) * (9 + (x) ))))))))) * 2.75573192e-6)


/* Decibel to Amplitude Conversion */
#define DB_TO_AMP(x) ((x) > -150.0 ? expf((x) * LOG_TEN_OVER_TWENTY) : 0.0f)
#define DB_TO_AMPD(x) ((x) > -150.0 ? exp((x) * LOG_TEN_OVER_TWENTY) : 0.0)

/* Amplitude to Decibel Conversion */
#define AMP_TO_DB(x) (((x) < 0.0000000298023223876953125) ? -150.0 : \
                      (logf(x) * TWENTY_OVER_LOG_TEN))
#define AMP_TO_DBD(x) (((x) < 0.0000000298023223876953125) ? -150.0 : \
                       (log(x) * TWENTY_OVER_LOG_TEN))

/* Smoothed Absolute Value */
#define SMOOTH_ABS(x) (sqrt(((x) * (x)) + 0.025) - sqrt(0.025))


namespace FXDSP
{
    /* 32 bit "pointer cast" union */
    typedef union
    {
        float f;
        int i;
    } f_pcast32;



    int next_pow2(int x)
    {
        if (x < 0)
            return 0;
        --x;
        x |= x >> 1;
        x |= x >> 2;
        x |= x >> 4;
        x |= x >> 8;
        x |= x >> 16;
        return x+1;
    }

    inline float f_abs(float f)
    {
        // flip the sign bit...
        int i = ((*(int*)&f) & 0x7fffffff);
        return (*(float*)&i);
    }





    /* f_max **********************************************************************/
    inline float f_max(float x, float a)
    {
        x -= a;
        x += fabs(x);
        x *= 0.5;
        x += a;
        return x;
    }



    /* f_min **********************************************************************/
    inline float f_min(float x, float b)
    {
        x = b - x;
        x += fabs(x);
        x *= 0.5;
        x = b - x;
        return x;
    }



    /* f_clamp ********************************************************************/
    inline float f_clamp(float x, float a, float b)
    {
        const float x1 = fabs(x - a);
        const float x2 = fabs(x - b);
        x = x1 + a + b;
        x -= x2;
        x *= 0.5;
        return x;
    }



    /* f_pow2 *********************************************************************/
    inline float f_pow2(float x)
    {
        f_pcast32 *px;
        f_pcast32 tx;
        f_pcast32 lx;
        float dx;

        px = (f_pcast32 *)&x;          // store address of float as long pointer
        tx.f = (x-0.5f) + (3<<22);      // temporary value for truncation
        lx.i = tx.i - 0x4b400000;       // integer power of 2
        dx = x - (float)lx.i;           // float remainder of power of 2

        x = 1.0f + dx * (0.6960656421638072f +  // cubic apporoximation of 2^x
            dx * (0.224494337302845f +       // for x in the range [0, 1]
            dx * (0.07944023841053369f)));
        (*px).i += (lx.i << 23);                // add int power of 2 to exponent

        return (*px).f;
    }

    /* f_tanh *********************************************************************/
    inline float f_tanh(float x)
    {
        double xa = f_abs(x);
        double x2 = xa * xa;
        double x3 = xa * x2;
        double x4 = x2 * x2;
        double x7 = x3 * x4;
        double res = (1.0 - 1.0 / (1.0 + xa + x2 + 0.58576695 * x3 + 0.55442112 * x4 + 0.057481508 * x7));
        return (x > 0 ? res : -res);
    }


    /* int16ToFloat ***************************************************************/
    inline float int16ToFloat(signed short sample)
    {
        return (float)(sample * INT16_TO_FLOAT_SCALAR);
    }

    /* floatToInt16 ***************************************************************/
    inline signed short floatToInt16(float sample)
    {
        return (signed short)(sample * 32767.0);
    }


    /* ratioToDb ******************************************************************/    
    inline float AmpToDb(float amplitude)
    {
        return 20.0*log10f(amplitude);
    }

    inline double AmpToDb(double amplitude)
    {
        return 20.0*log10(amplitude);
    }

    /* dbToRatio ******************************************************************/
    inline float DbToAmp(float dB)
    {
        return (dB > -150.0f ? expf(dB * LOG_TEN_OVER_TWENTY): 0.0f);
    }

    inline double DbToAmp(double dB)
    {
        return (dB > -150.0 ? exp(dB * LOG_TEN_OVER_TWENTY): 0.0);
    }

    void RectToPolar(float real, float imag, float* outMag, float* outPhase)
    {
        *outMag = sqrtf(powf(real, 2.0) + powf(imag, 2.0));
        double phase = atanf(imag/real);
        if (phase < 0)
            phase += M_PI;
        *outPhase = phase;
    }

    void RectToPolar(double real, double imag, double* outMag, double* outPhase)
    {
        *outMag = sqrt(pow(real, 2.0) + pow(imag, 2.0));
        double phase = atan(imag/real);
        if (phase < 0)
            phase += M_PI;
        *outPhase = phase;
    }

    void PolarToRect(float mag, float phase, float* outReal, float* outImag)
    {
        *outReal = mag * cosf(phase);
        *outImag = mag * sinf(phase);
    }

    void PolarToRect(double mag, double phase, double* outReal, double* outImag)
    {
        *outReal = mag * cos(phase);
        *outImag = mag * sin(phase);
    }

    /*******************************************************************************
    FloatBufferToInt16 */
    Error_t
    FloatBufferToInt16(signed short* dest, const float* src, unsigned length)
    {
    #ifdef __APPLE__
        // Use the Accelerate framework if we have it
        float scale = (float)INT16_MAX;
        float temp[length];
        vDSP_vsmul(src, 1, &scale, temp, 1, length);
        vDSP_vfix16(temp,1,dest,1,length);

    #else
        // Otherwise do it manually
        unsigned i;        
        const unsigned end = 4 * (length / 4);
        #pragma omp simd aligned(dest,src)
        for (i = 0; i < end; i+=4)
        {
            dest[i] = floatToInt16(*src++);
            dest[i + 1] = floatToInt16(*src++);
            dest[i + 2] = floatToInt16(*src++);
            dest[i + 3] = floatToInt16(*src++);
        }
        #pragma omp simd aligned(dest,src)
        for (i = end; i < length; ++i)
        {
            dest[i] = floatToInt16(*src++);
        }
    #endif
        return NOERR;
    }

    /*******************************************************************************
    DoubleBufferToInt16 */
    Error_t
    DoubleBufferToInt16(signed short* dest, const double* src, unsigned length)
    {
    #ifdef __APPLE__
        // Use the Accelerate framework if we have it
        double scale = (float)INT16_MAX;
        double temp[length];
        vDSP_vsmulD(src, 1, &scale, temp, 1, length);
        vDSP_vfix16D(temp,1,dest,1,length);

    #else
        // Otherwise do it manually
        unsigned i;
        const unsigned end = 4 * (length / 4);
        #pragma omp simd aligned(dest,src)
        for (i = 0; i < end; i+=4)
        {
            dest[i] = floatToInt16(*src++);
            dest[i + 1] = floatToInt16(*src++);
            dest[i + 2] = floatToInt16(*src++);
            dest[i + 3] = floatToInt16(*src++);
        }
        #pragma omp simd aligned(dest,src)
        for (i = end; i < length; ++i)
        {
            dest[i] = floatToInt16(*src++);
        }
    #endif
        return NOERR;
    }

    /*******************************************************************************
    Int16BufferToFloat */
    Error_t
    Int16BufferToFloat(float* dest, const signed short* src, unsigned length)
    {
    #ifdef __APPLE__
        // Use the Accelerate framework if we have it
        float temp[length];
        float scale = 1.0 / (float)INT16_MAX;
        vDSP_vflt16(src,1,temp,1,length);
        vDSP_vsmul(temp, 1, &scale, dest, 1, length);

    #else
        // Otherwise do it manually
        unsigned i;
        #pragma omp simd aligned(dest,src)
        const unsigned end = 4 * (length / 4);
        for (i = 0; i < end; i+=4)
        {
            dest[i] = int16ToFloat(*src++);
            dest[i + 1] = int16ToFloat(*src++);
            dest[i + 2] = int16ToFloat(*src++);
            dest[i + 3] = int16ToFloat(*src++);
        }
        #pragma omp simd aligned(dest,src)
        for (i = end; i < length; ++i)
        {
            dest[i] = int16ToFloat(*src++);
        }
    #endif
        return NOERR;
    }

    /*******************************************************************************
    Int16BufferToDouble */
    Error_t
    Int16BufferToDouble(double* dest, const signed short* src, unsigned length)
    {
    #ifdef __APPLE__
        // Use the Accelerate framework if we have it
        double temp[length];
        double scale = 1.0 / (double)INT16_MAX;
        vDSP_vflt16D(src,1,temp,1,length);
        vDSP_vsmulD(temp, 1, &scale, dest, 1, length);

    #else
        // Otherwise do it manually
        unsigned i;
        const unsigned end = 4 * (length / 4);
        #pragma omp simd aligned(dest,src)
        for (i = 0; i < end; i+=4)
        {
            dest[i] = int16ToFloat(*src++);
            dest[i + 1] = int16ToFloat(*src++);
            dest[i + 2] = int16ToFloat(*src++);
            dest[i + 3] = int16ToFloat(*src++);
        }
        #pragma omp simd aligned(dest,src)
        for (i = end; i < length; ++i)
        {
            dest[i] = int16ToFloat(*src++);
        }
    #endif
        return NOERR;
    }


    /*******************************************************************************
    DoubleToFloat */
    Error_t
    DoubleToFloat(float* dest, const double* src, unsigned length)
    {
    #ifdef __APPLE__
        vDSP_vdpsp(src, 1, dest, 1, length);
    #else
        // Otherwise do it manually
        unsigned i;
        const unsigned end = 4 * (length / 4);
        #pragma omp simd aligned(dest,src)
        for (i = 0; i < end; i+=4)
        {
            dest[i] = (float)src[i];
            dest[i + 1] = (float)src[i + 1];
            dest[i + 2] = (float)src[i + 2];
            dest[i + 3] = (float)src[i + 3];
        }
        #pragma omp simd aligned(dest,src)
        for (i = end; i < length; ++i)
        {
            dest[i] = (float)src[i];
        }
    #endif
        return NOERR;
    }


    /*******************************************************************************
    Float To Double */
    Error_t
    FloatToDouble(double* dest, const float* src, unsigned length)
    {
    #ifdef __APPLE__
        vDSP_vspdp(src, 1, dest, 1, length);
    #else
        // Otherwise do it manually
        unsigned i;
        const unsigned end = 4 * (length / 4);
        for (i = 0; i < end; i+=4)
        {
            dest[i] = (double)(src[i]);
            dest[i + 1] = (double)(src[i + 1]);
            dest[i + 2] = (double)(src[i + 2]);
            dest[i + 3] = (double)(src[i + 3]);
        }
        for (i = end; i < length; ++i)
        {
            dest[i] = (double)(src[i]);
        }
    #endif
        return NOERR;
    }


    /*******************************************************************************
    FillBuffer */
    template<typename T>
    void FillBuffer(T *dest, unsigned length, T value)
    {        
        unsigned i = 0;        
        #pragma omp simd aligned(dest)
        for (i = 0; i < length; ++i)
        {
            dest[i] = value;     
        }    
    }
    /*******************************************************************************
    ClearBuffer */
    template<typename T>
    void  ClearBuffer(T *dest, unsigned length)
    {
        memset(dest, 0, length * sizeof(T));
    }

    /*******************************************************************************
    CopyBuffer */
    template<typename T>
    void CopyBuffer(T* dest, const T* src, unsigned length)
    {
        if (src != dest)
        {
            memmove(dest, src, length * sizeof(T));
        }    
    }

    /*******************************************************************************
    CopyBufferStride */
    template<typename T>
    void  CopyBufferStride(T*         dest,
                    unsigned       dest_stride,
                    const T*   src,
                    unsigned       src_stride,
                    unsigned       length)
    {
    #if defined(USE_BLAS)        
        cblas_scopy(length, src, src_stride, dest, dest_stride);
    #else
        #pragma omp simd aligned(dest,src)
        for (unsigned i = 0; i < length; ++i)
        {
            dest[i * dest_stride] = src[i * src_stride];
        }
    #endif        
    }

    /*******************************************************************************
    SplitToInterleaved */
    template<typaname T>
    void  SplitToInterleaved(T* dest, const T* real, const T* imag, unsigned length)
    {
    #if defined(USE_BLAS)
        cblas_scopy(length, real, 1, dest, 2);
        cblas_scopy(length, imag, 1, dest + 1, 2);
    #else
        unsigned i;
        unsigned i2;      
        #pragma omp simd aligned(dest,real,imag)
        for (i = 0; i < length; ++i)
        {
            i2 = i * 2;
            dest[i2]     = real[i];
            dest[i2 + 1] = imag[i];
        }
    #endif        
    }

    template<typaname T>
    void InterleavedToSplit(T*       real,
                    T*       imag,
                    const T* input,
                    unsigned     length)
    {
    #if defined(USE_BLAS)
        cblas_scopy(length, input, 2, real, 1);
        cblas_scopy(length, input + 1, 2, imag, 1);
    #else
        unsigned i;
        unsigned i2;
        #pragma omp simd aligned(real,imag,input)
        for (i = 0; i < length; ++i)
        {
            i2 = i * 2;
            real[i] = input[i2];
            imag[i] = input[i2 + 1];
        }

    #endif        
    }

    /*******************************************************************************
     VectorAbs */
     template<typaname T>
    void VectorAbs(T *dest, const T *in, unsigned length)
    {
        // Otherwise do it manually
        unsigned i = 0;
        #pragma omp simd
        for (unsigned i = 0; i < length; ++i)
        {
            dest[i] = fabsf(in[i]);
        }
    }

    template<typename T>
    void VectorNegate(T          *dest,
                const T    *in,
                unsigned       length)
    {
        unsigned i;
        #pragma omp simd aligned(dest,in)
        for (i = 0; i < length; ++i)
        {
            dest[i] = -in[i];
        }            
    }

    /*******************************************************************************
    VectorSum */
    template<typename T>
    T VectorSum(const T* src, unsigned length)
    {
        float res = 0.0;
        #pragma omp simd aligned(src)
        for (unsigned i = 0; i < length; ++i)
        {
            res += src[i];
        }
        return res;
    }

    /*******************************************************************************
     VectorMax */
    template<typename T>
    T VectorMax(const T* src, unsigned length)
    {
        float res = FLT_MIN;
        #pragma omp simd aligned(src)
        for (unsigned i = 0; i < length; ++i)
        {
            if (src[i] > res)
            {
                res = src[i];
            }
        }
        return res;
    }

    /*******************************************************************************
    VectorMaxVI */
    template<typename T>
    void VectorMaxVI(T* value, unsigned* index, const T* src, unsigned length)
    {
        float res = FLT_MIN;
        #pragma omp simd aligned(src,value,index)
        for (unsigned i = 0; i < length; ++i)
        {
            if (src[i] > res)
            {
                *value = res = src[i];
                *index = i;
            }
        }
    }

    /*******************************************************************************
    VectorMin */
    template<typename T>
    T VectorMin(const T* src, unsigned length)
    {
        float res = FLT_MAX;
        #pragma omp simd aligned(src)
        for (unsigned i = 0; i < length; ++i)
        {
            if (src[i] < res)
            {
                res = src[i];
            }
        }    
        return res;
    }

    /*******************************************************************************
    VectorMinVI */
    template<typename T>
    void VectorMinVI(T* value, unsigned* index, const T* src, unsigned length)
    {

        float res = src[0];
        #pragma omp simd aligned(src,value,index)
        for (unsigned i = 0; i < length; ++i)
        {
            if (src[i] < res)
            {
                *value = res = src[i];
                *index = i;
            }
        }
    }

    /*******************************************************************************
    VectorVectorAdd */
    template<typename T>
    void VectorVectorAdd(T         *dest,
                    const T   *in1,
                    const T   *in2,
                    unsigned      length)
    {     
        unsigned i;
        #pragma omp simd aligned(dest,in1,in2)
        for (i = 0; i < length; ++i)
        {
            dest[i] = in1[i] + in2[i];
        }
    }

    /*******************************************************************************
    VectorVectorSub */
    template<typename T>
    void VectorVectorSub(T         *dest,
                    const T   *in1,
                    const T   *in2,
                    unsigned      length)
    {
        unsigned i;
        #pragma omp simd aligned(dest,in1,in2)
        for (i = 0; i < length; ++i)
        {
            dest[i] = in2[i] - in1[i];
        }
    }

    /*******************************************************************************
    VectorScalarAdd */
    template<typename T>
    void VectorScalarAdd(T           *dest,
                    const T     *in1,
                    const T     scalar,
                    unsigned        length)
    {
        unsigned i;
        #pragma omp simd aligned(dest,in1)
        for (i = 0; i < length; ++i)
        {
            dest[i] = in1[i] + scalar;
        }
    }

    /*******************************************************************************
    VectorVectorMultiply */
    template<typename T>
    void VectorVectorMultiply(T    *dest,
                        const T    *in1,
                        const T    *in2,
                        unsigned       length)
    {
        unsigned i;
        #pragma omp simd aligned(dest,in1,in2)
        for (i = 0; i < length; ++i)
        {
            dest[i] = in1[i] * in2[i];
        }
    }

    /*******************************************************************************
    VectorScalarMultiply */
    template<typename T>
    void VectorScalarMultiply(T          *dest,
                        const T    *in1,
                        const T    scalar,
                        unsigned       length)
    {    
        unsigned i;
        #pragma omp simd aligned(dest,in1)
        for (i = 0; i < length; ++i)
        {
            dest[i] = in1[i] * scalar;
        }
    }

    /*******************************************************************************
    VectorVectorMix */

    template<typename T>
    void VectorVectorMix( T  *dest,
                    const T  *in1,
                    const T  *scalar1,
                    const T  *in2,
                    const T  *scalar2,
                    unsigned     length)
    {    
        unsigned i;
        #pragma omp simd aligned(dest,in1,in2)
        for (i = 0; i < length; ++i)
        {
            dest[i] = (in1[i] * (*scalar1)) + (in2[i] * (*scalar2));
        }        
    }

    template<typename T>
    void VectorVectorSumScale(T        *dest,
                        const T  *in1,
                        const T  *in2,
                        const T  *scalar,
                        unsigned     length)
    {    
        unsigned i;
        #pragma omp simd aligned(dest,in1,in2)
        for (i = 0; i < length; ++i)
        {
            dest[i] = (in1[i] + in2[i]) * (*scalar);
        }    
    }

    /*******************************************************************************
    VectorPower */
    template<typename T>
    void VectorPower(T* dest, const T* in, T power, unsigned length)
    {    
        unsigned i;
        #pragma omp simd aligned(dest,in)
        for (i = 0; i < length; ++i)
        {
            dest[i] = powf(in[i], power);
        }
    }

    /*******************************************************************************
    Convolve */
    template<typename T>
    void Convolve(T     *in1,
            unsigned    in1_length,
            T           *in2,
            unsigned    in2_length,
            T           *dest)
    {

        unsigned resultLength = in1_length + (in2_length - 1);            
        unsigned i;
        #pragma omp simd aligned(dest,in1,in2)
        for (i = 0; i <resultLength; ++i)
        {
            unsigned kmin, kmax, k;
            dest[i] = 0;

            kmin = (i >= (in2_length - 1)) ? i - (in2_length - 1) : 0;
            kmax = (i < in1_length - 1) ? i : in1_length - 1;            
            for (k = kmin; k <= kmax; k++)
            {
                dest[i] += in1[k] * in2[i - k];
            }
        }
    }

    /*******************************************************************************
    VectorDbConvert */
    template<typename T>
    void VectorDbConvert(T* dest,
                    const T* in,
                    unsigned length)
    {
        #pragma omp simd aligned(dest,in)
        for (unsigned i = 0; i < length; ++i)
        {
            dest[i] = AMP_TO_DB(in[i]);
        }
    }

    /*******************************************************************************
    ComplexMultiply */
    template<typename T>
    void ComplexMultiply(T*          re,
                    T*          im,
                    const T*    re1,
                    const T*    im1,
                    const T*    re2,
                    const T*    im2,
                    unsigned        length)
    {
        #pragma omp simd aligned(re,im,re1,im1,re2,im2)
        for (unsigned i = 0; i < length; ++i)
        {
            float ire1 = re1[i];
            float iim1 = im1[i];
            float ire2 = re2[i];
            float iim2 = im2[i];
            re[i] = (ire1 * ire2 - iim1 * iim2);
            im[i] = (iim1 * ire2 + iim2 * ire1);
        }
    }

    /*******************************************************************************
    VectorRectToPolar */
    template<typename T>
    void VectorRectToPolar(T*        magnitude,
                    T*        phase,
                    const T*  real,
                    const T*  imaginary,
                    unsigned      length)
    {
        unsigned i;
        // might not do anything
        #pragma omp simd         aligned(magnitude,phase,real,imaginary)
        for (i = 0; i < length; ++i)
        {
            RectToPolar(real[i], imaginary[i], &magnitude[i], &phase[i]);
        }
    }

    /*******************************************************************************
    MeanSquare */
    template<typename T>
    T MeanSquare(const T* data, unsigned length)
    {
        float result = 0.;
        float scratch[length];
        VectorPower(scratch, data, 2, length);
        result = VectorSum(scratch, length) / length;
        return result;
    }

}
