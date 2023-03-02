
/*
 * Dsp.c
 * Hamilton Kibbe
 * Copyright 2013 Hamilton Kibbe
 */

#pragma once

#include <cmath>
#include <cstring>
#include <cstddef>
#include <cstdlib>
#include <cfloat>

#include <fftw3.h>



/* Function Declarations ******************************************************/

/* Define log2 and log2f for MSVC */
#if !defined(log2) || !defined(log2f)
#define _USE_FXDSP_LOG

double log2(double n);

float log2f(float n);

#endif


/**  Find the nearest power of two
 * @param x     number to process
 * @return      Absolute value of f.
 */
int next_pow2(int x);


/**  Fast absolute value
 * @details     Fast fabs() implementation
 * @param f     Value to process
 * @return      Absolute value of f.
 */
float f_abs(float f);


/**  Max of two floats
 * @details branchless max() implementation
 * @param   x first value to compare,
 * @param   a second value to compare.
 * @return  the maximum value of the two arguments.
 */
float f_max(float x, float a);


/**  Min of two floats
 * @details branchless min() implementation
 * @param   x first value to compare,
 * @param   a second value to compare.
 * @return  the minimum value of the two arguments.
 */
float f_min(float x, float b);


/**  Clamp values to range
 * @details branchless LIMIT() implementation
 * @param x value to clamp
 * @param a lower bound
 * @param b upper bound
 * @return  val clamped to range (a, b)
 */
float f_clamp(float x, float a, float b);


/** Calculate pow(2, x)
 * @details fast, branchless pow(2,x) approximation
 * @param x     power of 2 to calculate.
 * @return      2^x
 */
float f_pow2(float x);


/** Calculate tanh_x
* @details fast tanh approximation
* @param x     input
* @return      ~tanh(x)
*/
float f_tanh(float x);


/** Convert signed sample to float
 *
 * @details convert a signed 16 bit sample to a 32 bit float sample in the range
 * [-1.0, 1,0]
 *
 * @param sample    The sample to convert.
 * @return          The sample as a float.
 */
float int16ToFloat(signed short sample);


/** Convert a float sample to signed
 *
 * @details convert a 32 bit float sample in the range [-1.0, 1,0] to a signed
 * 16 bit sample.
 *
 * @param sample    The sample to convert.
 * @return          The sample as a 16-bit signed int.
 */
signed short floatToInt16(float sample);


/** Convert an amplitude to dB
 * @details     Convert a voltage amplitude to dB.
 * @param amp   The amplitude to convert.
 * @return      Amplitude value in dB.
 */
float AmpToDb(float ratio);

double AmpToDbD(double ratio);


/** Convert a value in dB to an amplitude
 * @details convert a dBFS value to a voltage amplitude
 * @param dB        The value in dB to convert.
 * @return          dB value as a voltage amplitude.
 */
float DbToAmp(float dB);

double DbToAmpD(double dB);


/** Convert complex value to magnitude/phase
 * @param real      Real part of input.
 * @param imag      Imaginary part of input.
 * @param outMag    Magnitude output.
 * @param outPhase  Phase output.
 */
void RectToPolar(float real, float imag, float* outMag, float* outPhase);

void RectToPolarD(double real, double imag, double* outMag, double* outPhase);


/** Convert magnitude/phase to complex
 * @param mag       Magnitude input.
 * @param phase     Phase input.
 * @param outReal   Real part output.
 * @param outImag   Imaginary part output.
 */
void PolarToRect(float mag, float phase, float* outReal, float* outImag);

void PolarToRectD(double mag, double phase, double* outReal, double* outImag);

/** Error codes */
typedef enum Error
{
    /** No Error (0) */
    NOERR,

    /** Generic Error (1) */
    ERROR,

    /** Malloc failure... */
    NULL_PTR_ERROR,

    /** invalid value... */
    VALUE_ERROR,

    /** Number of defined error codes */
    N_ERRORS
}Error_t;


typedef enum _bias_t
{
    /** Pass positive signals, clamp netagive signals to 0 */
    FORWARD_BIAS,

    /** Pass negative signals, clamp positive signals to 0 */
    REVERSE_BIAS,

    /** Full-wave rectification. */
    FULL_WAVE
}bias_t;



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


struct FFTConfig
{
    unsigned        length;
    float           scale;
    float           log2n;
    FFTSplitComplex split;
    FFTSplitComplex split2;
    FFT_SETUP        setup;
};

struct FFTConfigD
{
    unsigned                length;
    double                  scale;
    double                  log2n;
    FFTSplitComplexD        split;
    FFTSplitComplexD        split2;
    FFT_SETUP_D             setup;
};

/** Filter types */
typedef enum Filter_t
{
    /** Lowpass */
    LOWPASS,

    /** Highpass */
    HIGHPASS,

    /** Bandpass */
    BANDPASS,

    /** Allpass */
    ALLPASS,

    /** Notch */
    NOTCH,

    /** Peaking */
    PEAK,

    /** Low Shelf */
    LOW_SHELF,

    /** High Shelf */
    HIGH_SHELF,

    /** Number of Filter types */
    N_FILTER_TYPES
}Filter_t;


/** The kernel length at which to use FFT convolution vs direct */
/* So this is a pretty good value for now */



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


/** Boltzman's constant */
static const float BOLTZMANS_CONSTANT = 1.38e-23;

/** Magnitude of electron charge */
static const float Q = 1.609e-19;


typedef enum
{
    FULL_SCALE,
    K_12,
    K_14,
    K_20
} MeterScale;



    /** Optocoupler types */
typedef enum _Opto_t
{
    /** Light-Dependent-Resistor output. Based
     on Vactrol VTL series. datasheet:
     http://pdf.datasheetcatalog.com/datasheet/perkinelmer/VT500.pdf

     Midpoint Delay values:
       Turn-on delay:   ~10ms
       Turn-off delay:  ~133ms
     */
    OPTO_LDR,

    /** TODO: Add Phototransistor/voltage output opto model*/
    OPTO_PHOTOTRANSISTOR
} Opto_t;


/** Resampling Factor constants */
typedef enum factor
{
    /** 2x resampling */
    X2 = 0,

    /** 4x resampling */
    X4,

    /** 8x resampling */
    X8,

    /** 16x resampling */
    /*X16,*/

    /** number of resampling factors */
    N_FACTORS
} ResampleFactor_t;


typedef enum  TapeSpeed
{
    TS_3_75IPS,
    TS_7_5IPS,
    TS_15IPS,
    TS_30IPS
}TapeSpeed;


/** Window function type */
typedef enum _Window_t
{
    /** Rectangular window */
    BOXCAR,

    /** Hann window */
    HANN,

    /** Hamming window */
    HAMMING,

    /** Blackman window */
    BLACKMAN,

    /** Tukey window */
    TUKEY,

    /** Cosine window */
    COSINE,

    /** Lanczos window */
    LANCZOS,

    /** Bartlett window */
    BARTLETT,

    /** Gauss window */
    GAUSSIAN,

    /** Bartlett-Hann window */
    BARTLETT_HANN,

    /** Kaiser window */
    KAISER,

    /** Nuttall window */
    NUTTALL,

    /** Blackaman-Harris window */
    BLACKMAN_HARRIS,

    /** Blackman-Nuttall window */
    BLACKMAN_NUTTALL,

    /** Flat top window */
    FLATTOP,

    /** Poisson window */
    POISSON,

    /** The number of window types */
    N_WINDOWTYPES
} Window_t;



/*******************************************************************************
 BiquadFilter */
struct BiquadFilter
{
    float b[3];     // b0, b1, b2
    float a[2];     // a1, a2
    float x[2];     //
    float y[2];
    float w[2];
};

struct BiquadFilterD
{
    double b[3];     // b0, b1, b2
    double a[2];     // a1, a2
    double  x[2];     //
    double  y[2];
    double  w[2];
};


/* FIRFilter ***********************************************************/
struct FIRFilter
{
    float*              kernel;
    const float*        kernel_end;
    float*              overlap;
    unsigned            kernel_length;
    unsigned            overlap_length;
    ConvolutionMode_t   conv_mode;
    FFTConfig*          fft_config;
    FFTSplitComplex     fft_kernel;
    unsigned            fft_length;
};

struct FIRFilterD
{
    double*             kernel;
    const double*       kernel_end;
    double*             overlap;
    unsigned            kernel_length;
    unsigned            overlap_length;
    ConvolutionMode_t   conv_mode;
    FFTConfigD*         fft_config;
    FFTSplitComplexD    fft_kernel;
    unsigned            fft_length;
};

/* LadderFilter ********************************************************/
struct LadderFilter
{
    float y[4];
    float w[4];
    float Vt;           // transistor treshold voltage [V]
    float sample_rate;
    float cutoff;
    float resonance;
};


/* RBJFilter ***********************************************************/
struct RBJFilter
{
    BiquadFilter* biquad;
    Filter_t type;
    float omega;
    float Q;
    float cosOmega;
    float sinOmega;
    float alpha;
    float A;
    float dbGain;
    float b[3];
    float a[3];
    float sampleRate;
};

struct RBJFilterD
{
    BiquadFilterD* biquad;
    Filter_t type;
    double omega;
    double Q;
    double cosOmega;
    double sinOmega;
    double alpha;
    double A;
    double dbGain;
    double b[3];
    double a[3];
    double sampleRate;
};


/* LRFilter ***************************************************************/
struct LRFilter
{
    RBJFilter*  filterA;
    RBJFilter*  filterB;
    Filter_t    type;
    float       cutoff;
    float       Q;
    float       sampleRate;
};

struct LRFilterD
{
    RBJFilterD* filterA;
    RBJFilterD* filterB;
    Filter_t    type;
    double      cutoff;
    double      Q;
    double      sampleRate;
};


/*******************************************************************************
 CircularBuffer */
struct CircularBuffer
{
    unsigned    length;
    unsigned    wrap;
    float*      buffer;
    unsigned    read_index;
    unsigned    write_index;
    unsigned    count;
};


/*******************************************************************************
 CircularBufferD */
struct CircularBufferD
{
    unsigned    length;
    unsigned    wrap;
    double*     buffer;
    unsigned    read_index;
    unsigned    write_index;
    unsigned    count;
};

/* Upsampler **********************************************************/
struct Decimator
{
    unsigned factor;
    FIRFilter** polyphase;
};

struct DecimatorD
{
    unsigned factor;
    FIRFilterD** polyphase;
};



/*******************************************************************************
 DiodeRectifier */
struct DiodeRectifier
{
    bias_t  bias;
    float   threshold;
    float   vt;
    float   scale;
    float   abs_coeff;
    float*  scratch;
};

struct DiodeRectifierD
{
    bias_t  bias;
    double  threshold;
    double  vt;
    double  scale;
    double abs_coeff;
    double* scratch;
};




/*******************************************************************************
 Diode */
struct DiodeSaturator
{
    bias_t  bias;
    float   amount;
};

struct DiodeSaturatorD
{
    bias_t  bias;
    double  amount;
};



/*******************************************************************************
 MultibandFilter */
struct MultibandFilter
{
    LRFilter*   LPA;
    LRFilter*   HPA;
    LRFilter*   LPB;
    LRFilter*   HPB;
    RBJFilter*  APF;
    float       lowCutoff;
    float       highCutoff;
    float       sampleRate;
};

struct MultibandFilterD
{
    LRFilterD*  LPA;
    LRFilterD*  HPA;
    LRFilterD*  LPB;
    LRFilterD*  HPB;
    RBJFilterD* APF;
    double      lowCutoff;
    double      highCutoff;
    double      sampleRate;
};


/* OnePoleFilter ********************************************************/
struct OnePole
{
    float a0;
    float b1;
    float y1;
    float cutoff;
    float sampleRate;
    Filter_t type;

};

struct OnePoleD
{
    double a0;
    double b1;
    double y1;
    double cutoff;
    double sampleRate;
    Filter_t type;
};



struct PolySaturator
{
    float a;
    float b;
    float n;
};


struct PolySaturatorD
{
    double a;
    double b;
    double n;
};



/*******************************************************************************
 RMSEstimator */
struct RMSEstimator
{
    float   avgTime;
    float   sampleRate;
    float   avgCoeff;
    float   RMS;
};

struct RMSEstimatorD
{
    double  avgTime;
    double  sampleRate;
    double  avgCoeff;
    double  RMS;
};

/* Upsampler **********************************************************/
struct Upsampler
{
    unsigned factor;
    FIRFilter** polyphase;
};

struct UpsamplerD
{
    unsigned factor;
    FIRFilterD** polyphase;
};


/* Static Function Prototypes */
static float
modZeroBessel(float x);

static double
modZeroBesselD(double x);

static float
chebyshev_poly(int n, float x);

static double
chebyshev_polyD(int n, double x);


/* Implementations */

struct WindowFunction
{
    float*      window;
    unsigned    length;
    Window_t    type;
};

struct WindowFunctionD
{
    double*     window;
    unsigned    length;
    Window_t    type;
};






/*******************************************************************************
 TapeSaturator */
struct Tape
{
    PolySaturator*  polysat;
    TapeSpeed       speed;
    float           sample_rate;
    float           saturation;
    float           hysteresis;
    float           flutter;
    float           pos_peak;
    float           neg_peak;
    float*          flutter_mod;
    unsigned        flutter_mod_length;
};




/*******************************************************************************
 Static Function Prototypes */

static void
calculate_bin_frequencies(float* dest, unsigned fft_length, float sample_rate);

static void
calculate_bin_frequenciesD(double* dest, unsigned fft_length, double sample_rate);


/*******************************************************************************
 Structs */

struct SpectrumAnalyzer
{
    unsigned        fft_length;
    unsigned        bins;
    float           sample_rate;
    float           mag_sum;
    float*          frequencies;
    float*          real;
    float*          imag;
    float*          mag;
    float*          phase;
    float*          root_moment;
    FFTConfig*      fft;
    Window_t        window_type;
    WindowFunction* window;
};

struct SpectrumAnalyzerD
{
    unsigned            fft_length;
    unsigned            bins;
    double              sample_rate;
    double              mag_sum;
    double*             frequencies;
    double*             real;
    double*             imag;
    double*             mag;
    double*             phase;
    double*             root_moment;
    FFTConfigD*         fft;
    Window_t            window_type;
    WindowFunctionD*    window;
};


/* Channel id numbers */
enum
{
    LEFT = 0,
    RIGHT,
    CENTER,
    LEFT_SURROUND,
    RIGHT_SURROUND,
    N_CHANNELS
};
double CHANNEL_GAIN[N_CHANNELS] =
{
    1.0,    /* LEFT */
    1.0,    /* RIGHT */
    1.0,    /* CENTER */
    1.41,   /* LEFT_SURROUND */
    1.41    /* RIGHT_SURROUND */
};


struct KWeightingFilter
{
    BiquadFilter*   pre_filter;
    BiquadFilter*   rlb_filter;
};


struct KWeightingFilterD
{
    BiquadFilterD*  pre_filter;
    BiquadFilterD*  rlb_filter;
};


struct BS1770Meter
{
    KWeightingFilter**  filters;
    Upsampler**         upsamplers;
    CircularBuffer**    buffers;
    unsigned            n_channels;
    unsigned            sample_count;
    unsigned            gate_len;
    unsigned            overlap_len;
};



struct BS1770MeterD
{
    KWeightingFilterD** filters;
    UpsamplerD**        upsamplers;
    CircularBufferD**   buffers;
    unsigned            n_channels;
    unsigned            sample_count;
    unsigned            gate_len;
    unsigned            overlap_len;
};




/** Create a new BiquadFilter
 *
 * @details Allocates memory and returns an initialized BiquadFilter.
 *          Play nice and call BiquadFilterFree on the filter when
 *          you're done with it.
 *
 * @param bCoeff    Numerator coefficients [b0, b1, b2]
 * @param aCoeff    Denominator coefficients [a1, a2]
 * @return          An initialized BiquadFilter
 */
BiquadFilter* BiquadFilterInit(const float* bCoeff, const float* aCoeff);

BiquadFilterD* BiquadFilterInitD(const double *bCoeff, const double *aCoeff);


/** Free memory associated with a BiquadFilter
 *
 * @details release all memory allocated by BiquadFilterInit for the
 *          supplied filter.
 * @param filter    BiquadFilter to free.
 * @return          Error code, 0 on success
 */
Error_t BiquadFilterFree(BiquadFilter* filter);

Error_t BiquadFilterFreeD(BiquadFilterD* filter);


/** Flush filter state buffers
 *
 * @param filter    BiquadFilter to flush.
 * @return          Error code, 0 on success
 */
Error_t BiquadFilterFlush(BiquadFilter* filter);

Error_t BiquadFilterFlushD(BiquadFilterD* filter);


/** Filter a buffer of samples
 * @details Uses a DF-II biquad implementation to filter input samples
 *
 * @param filter    The BiquadFilter to use.
 * @param outBuffer The buffer to write the output to.
 * @param inBuffer  The buffer to filter.
 * @param n_samples The number of samples to filter.
 * @return          Error code, 0 on success
 */
Error_t BiquadFilterProcess(BiquadFilter*   filter,
                    float*          outBuffer,
                    const float*    inBuffer,
                    unsigned        n_samples);

Error_t BiquadFilterProcessD(BiquadFilterD  *filter,
                     double         *outBuffer,
                     const double   *inBuffer,
                     unsigned       n_samples);

/** Filter a single samples
 * @details Uses a DF-II biquad implementation to filter input sample
 *
 * @param filter    The BiquadFilter to use.
 * @param in_sample The sample to process.
 * @return          Filtered sample.
 */
float BiquadFilterTick(BiquadFilter* filter, float in_sample);

double BiquadFilterTickD(BiquadFilterD* filter, double in_sample);


/** Update the filter kernel for a given filter
 *
 * @param filter    The filter to update
 * @param bCoeff    Numerator coefficients [b0, b1, b2]
 * @param aCoeff    Denominator coefficients [a1, a2]
 */
Error_t BiquadFilterUpdateKernel(BiquadFilter*  filter,
                         const float*   bCoeff,
                         const float*   aCoeff);

Error_t BiquadFilterUpdateKernelD(BiquadFilterD *filter,
                          const double  *bCoeff,
                          const double  *aCoeff);


KWeightingFilter* KWeightingFilterInit(float sample_rate);

KWeightingFilterD* KWeightingFilterInitD(double sample_rate);

Error_t KWeightingFilterProcess(KWeightingFilter*   filter,
                        float*              dest,
                        const float*        src,
                        unsigned            length);

Error_t KWeightingFilterProcessD(KWeightingFilterD* filter,
                         double*            dest,
                         const double*      src,
                         unsigned           length);


Error_t KWeightingFilterFlush(KWeightingFilter* filter);

Error_t KWeightingFilterFlushD(KWeightingFilterD* filter);


Error_t KWeightingFilterFree(KWeightingFilter* filter);

Error_t KWeightingFilterFreeD(KWeightingFilterD* filter);


BS1770Meter* BS1770MeterInit(unsigned n_channels, float sample_rate);

BS1770MeterD* BS1770MeterInitD(unsigned n_channels, double sample_rate);

Error_t BS1770MeterProcess(BS1770Meter*     meter,
                   float*           loudness,
                   float**          peaks,
                   const float**    samples,
                   unsigned         n_samples);

Error_t BS1770MeterProcessD(BS1770MeterD*   meter,
                    double*         loudness,
                    double**        peaks,
                    const double**  samples,
                    unsigned        n_samples);

Error_t BS1770MeterFree(BS1770Meter* meter);

Error_t BS1770MeterFreeD(BS1770MeterD* meter);

/** Create a new FIRFilter
 *
 * @details Allocates memory and returns an initialized FIRFilter.
 *			Play nice and call FIRFilterFree on the filter when you're
 *          done with it.
 *
 * @param filter_kernel     The filter coefficients. These are copied to the
 *                          filter so there is no need to keep them around.
 * @param length            The number of coefficients in filter_kernel.
 * @param convolution_mode  Convolution algorithm. Either BEST, FFT, or DIRECT.
 * @return                  An initialized FIRFilter
 */
FIRFilter*
FIRFilterInit(const float*      filter_kernel,
              unsigned          length,
              ConvolutionMode_t convolution_mode);
FIRFilterD*
FIRFilterInitD(const double*        filter_kernel,
               unsigned             length,
               ConvolutionMode_t    convolution_mode);


/** Free memory associated with a FIRFilter
 *
 * @details release all memory allocated by FIRFilterInit for the
 *			supplied filter.
 *
 * @param filter	FIRFilter to free
 * @return			Error code, 0 on success
 */
Error_t
FIRFilterFree(FIRFilter* filter);

Error_t
FIRFilterFreeD(FIRFilterD* filter);


/** Flush filter state buffer
 *
 * @param filter	FIRFilter to flush
 * @return			Error code, 0 on success
 */
Error_t
FIRFilterFlush(FIRFilter* filter);

Error_t
FIRFilterFlushD(FIRFilterD* filter);


/** Filter a buffer of samples
 *
 * @details Uses either FFT or direct-form convolution to filter the samples.
 *
 * @param filter	The FIRFilter to use
 * @param outBuffer	The buffer to write the output to
 * @param inBuffer	The buffer to filter
 * @param n_samples The number of samples to filter
 * @return			Error code, 0 on success
 */
Error_t
FIRFilterProcess(FIRFilter*     filter,
                 float*         outBuffer,
                 const float*   inBuffer,
				 unsigned       n_samples);


Error_t
FIRFilterProcessD(FIRFilterD*   filter,
                  double*       outBuffer,
                  const double* inBuffer,
                  unsigned      n_samples);


/** Update the filter kernel for a given filter
 *
 * @details New kernel must be the same length as the old one!
 *
 * @param filter		The FIRFilter to use
 * @param filter_kernel	The new filter kernel to use
 * @return			Error code, 0 on success
 */
Error_t
FIRFilterUpdateKernel(FIRFilter*    filter,
					  const float*  filter_kernel);

Error_t
FIRFilterUpdateKernelD(FIRFilterD*    filter,
                       const double*  filter_kernel);



/** Create a new DiodeSaturator
 *
 * @details Allocates memory and returns an initialized DiodeSaturator.
 *          call DiodeSaturatorFree to release allocated memory
 *
 * @param amount    Clipping amount
 * @param bias      Diode bias, FORWARD_BIAS or REVERSE_BIAS
 *                  Forward-bias will clip positive signals and leave negative
 *                  signals untouched.
 * @return          An initialized DiodeSaturator
 */
DiodeSaturator*
DiodeSaturatorInit(bias_t bias, float amount);

DiodeSaturatorD*
DiodeSaturatorInitD(bias_t bias, double amount);


/** Free memory associated with a DiodeSaturator
 *
 * @details release all memory allocated by DiodeSaturatorInit for the
 *          given diode.
 *
 * @param diode     DiodeSaturator to free
 * @return          Error code, 0 on success
 */
Error_t
DiodeSaturatorFree(DiodeSaturator* saturator);

Error_t
DiodeSaturatorFreeD(DiodeSaturatorD* saturator);


/** Update DiodeSaturator clipping amount
 *
 * @details Update the diode model's clipping amount
 *
 * @param diode     DiodeSaturator to update
 * @param amount    New diode clipping amount [0 1]
 * @return          Error code, 0 on success
 */
Error_t
DiodeSaturatorSetAmount(DiodeSaturator* saturator, float amount);

Error_t
DiodeSaturatorSetAmountD(DiodeSaturatorD* saturator, double amount);


/** Process a buffer of samples
 * @details Uses a diode saturator model to process input samples
 *
 * @param diode     The DiodeSaturator to use.
 * @param outBuffer The buffer to write the output to.
 * @param inBuffer  The buffer to filter.
 * @param n_samples The number of samples to filter.
 * @return          Error code, 0 on success
 */
Error_t
DiodeSaturatorProcess(DiodeSaturator*   saturator,
                      float*            out_buffer,
                      const float*      in_buffer,
                      unsigned          n_samples);

Error_t
DiodeSaturatorProcessD(DiodeSaturatorD* saturator,
                       double*          out_buffer,
                       const double*    in_buffer,
                       unsigned         n_samples);


/** Process a single sample
 * @details Uses a diode saturator model to process an input sample
 *
 * @param diode     The DiodeSaturator to use.
 * @param in_sample The sample to process.
 * @return          A processed sample.
 */
float
DiodeSaturatorTick(DiodeSaturator* saturator, float in_sample);

double
DiodeSaturatorTickD(DiodeSaturatorD* saturator, double in_sample);


/** Update DiodeRectifier threshold voltage
 *
 * @details Update the diode model's threshold voltage
 *
 * @param diode     DiodeRectifier instance to update
 * @param threshold	New diode threshold [0 1]
 * @return			Error code, 0 on success
 */
Error_t
DiodeRectifierSetThreshold(DiodeRectifier* diode, float threshold);

Error_t
DiodeRectifierSetThresholdD(DiodeRectifierD* diode, double threshold);


/** Process a buffer of samples
 * @details Uses a diode rectifier to process input samples.
 *
 * @param diode     The DiodeRectifier instance to use.
 * @param outBuffer	The buffer to write the output to.
 * @param inBuffer	The buffer to filter.
 * @param n_samples The number of samples to filter.
 * @return			Error code, 0 on success
 */
Error_t
DiodeRectifierProcess(DiodeRectifier*   diode,
                      float*            out_buffer,
                      const float*      in_buffer,
                      unsigned          n_samples);

Error_t
DiodeRectifierProcessD(DiodeRectifierD* diode,
                       double*          out_buffer,
                       const double*    in_buffer,
                       unsigned         n_samples);


/** Process a single sample
 * @details Uses a diode rectifier model to process an input sample
 *
 * @param diode     The DiodeRectifier instance to use.
 * @param in_sample	The sample to process.
 * @return			A processed sample.
 */
float
DiodeRectifierTick(DiodeRectifier* diode, float in_sample);

double
DiodeRectifierTickD(DiodeRectifierD* diode, double in_sample);


/** Create a new CircularBuffer
 *
 * @details Allocates memory and returns an initialized CircularBuffer.
 *			Play nice and call CircularBuffer Free on the filter when you're
 *          done with it.
 *
 * @param length		The minimum number of elements in the circular buffer
 */
CircularBuffer*  CircularBufferInit(unsigned length);

CircularBufferD* CircularBufferInitD(unsigned length);


/** Free Heap Memory associated with CircularBuffer
*
* @details Frees memory allocated by CircularBufferInit
*/
Error_t CircularBufferFree(CircularBuffer* cb);

Error_t CircularBufferFreeD(CircularBufferD* cb);


/** Write samples to circular buffer
*/
Error_t CircularBufferWrite(CircularBuffer* cb, const float* src, unsigned n_samples);

Error_t CircularBufferWriteD(CircularBufferD*   cb,
                     const double*      src,
                     unsigned           n_samples);


/** Read samples from circular buffer
 */
Error_t CircularBufferRead(CircularBuffer* cb, float* dest, unsigned n_samples);

Error_t CircularBufferReadD(CircularBufferD* cb, double* dest, unsigned n_samples);


/** Flush circular buffer
 */
Error_t CircularBufferFlush(CircularBuffer* cb);

Error_t CircularBufferFlushD(CircularBufferD* cb);

/** Rewind the read head of the buffer by `n_samples` samples
 */
Error_t CircularBufferRewind(CircularBuffer* cb, unsigned n_samples);

Error_t CircularBufferRewindD(CircularBufferD* cb, unsigned n_samples);


/** Return the number of unread samples in the buffer
 */
unsigned CircularBufferCount(CircularBuffer* cb);

unsigned CircularBufferCountD(CircularBufferD* cb);



/** Create a new Decimator
 *
 * @details Allocates memory and returns an initialized Decimator with
 *          a given decimation factor.
 *
 * @param factor    Decimation factor
 * @return          An initialized Decimator
 */
Decimator* DecimatorInit(ResampleFactor_t factor);

DecimatorD* DecimatorInitD(ResampleFactor_t factor);

/** Free memory associated with a Upsampler
 *
 * @details release all memory allocated by DecimatorInit for the
 *          supplied filter.
 * @param decimator Decimator to free.
 * @return          Error code, 0 on success
 */
Error_t DecimatorFree(Decimator* decimator);

Error_t DecimatorFreeD(DecimatorD* decimator);

/** Flush decimator state buffers
 *
 * @param decimator Upsampler to flush.
 * @return          Error code, 0 on success
 */
Error_t DecimatorFlush(Decimator* decimator);

Error_t DecimatorFlushD(DecimatorD* decimator);

/** Decimate a buffer of samples
 *
 * @details Decimates given buffer using a polyphase decimator
 *
 * @param decimator The Decimator to use
 * @param outBuffer The buffer to write the output to
 * @param inBuffer  The buffer to filter
 * @param n_samples The number of samples to upsample
 * @return          Error code, 0 on success
 */
Error_t DecimatorProcess(Decimator*     decimator,
                 float*         outBuffer,
                 const float*   inBuffer,
                 unsigned       n_samples);

Error_t DecimatorProcessD(DecimatorD*   decimator,
                 double*        outBuffer,
                 const double*  inBuffer,
                 unsigned       n_samples);


/** Create a new DiodeRectifier
 *
 * @details Allocates memory and returns an initialized DiodeRectifier.
 *			Play nice and call DiodeRectifierFree when you're done with it.
 *
 * @param threshold         Normalized voltage threshold
 * @param bias              DiodeRectifier bias, FORWARD_BIAS or REVERSE_BIAS
 *                          Forward-bias will pass positive signals and clamp
 *                          negative signals to 0.
 * @return                  An initialized DiodeRectifier
 */
DiodeRectifier* DiodeRectifierInit(bias_t bias, float threshold);

DiodeRectifierD* DiodeRectifierInitD(bias_t bias, double threshold);


/** Free memory associated with a DiodeRectifier
 *
 * @details release all memory allocated by DiodeRectifierInit for the
 *			given DiodeRectifier.
 *
 * @param DiodeRectifier     DiodeRectifier to free
 * @return			Error code, 0 on success
 */
Error_t DiodeRectifierFree(DiodeRectifier* diode);

Error_t DiodeRectifierFreeD(DiodeRectifierD* diode);


/** Update DiodeRectifier threshold voltage
 *
 * @details Update the diode model's threshold voltage
 *
 * @param diode     DiodeRectifier instance to update
 * @param threshold	New diode threshold [0 1]
 * @return			Error code, 0 on success
 */
Error_t DiodeRectifierSetThreshold(DiodeRectifier* diode, float threshold);

Error_t DiodeRectifierSetThresholdD(DiodeRectifierD* diode, double threshold);


/** Process a buffer of samples
 * @details Uses a diode rectifier to process input samples.
 *
 * @param diode     The DiodeRectifier instance to use.
 * @param outBuffer	The buffer to write the output to.
 * @param inBuffer	The buffer to filter.
 * @param n_samples The number of samples to filter.
 * @return			Error code, 0 on success
 */
Error_t DiodeRectifierProcess(DiodeRectifier*   diode,
                      float*            out_buffer,
                      const float*      in_buffer,
                      unsigned          n_samples);

Error_t DiodeRectifierProcessD(DiodeRectifierD* diode,
                       double*          out_buffer,
                       const double*    in_buffer,
                       unsigned         n_samples);


/** Process a single sample
 * @details Uses a diode rectifier model to process an input sample
 *
 * @param diode     The DiodeRectifier instance to use.
 * @param in_sample	The sample to process.
 * @return			A processed sample.
 */
float DiodeRectifierTick(DiodeRectifier* diode, float in_sample);

double DiodeRectifierTickD(DiodeRectifierD* diode, double in_sample);


static inline void interleave_complex(float*       dest,
                   const float* real,
                   const float* im,
                   unsigned     length);

static inline void interleave_complexD(double*         dest,
                    const double*   real,
                    const double*   im,
                    unsigned        length);


static inline void split_complex(float*        real,
              float*        im,
              const float*  data,
              unsigned      length);


static inline void split_complexD(double*          real,
               double*          im,
               const double*    data,
               unsigned         length);



/** Generate a Boxcar window of length n
 *
 * @details Create an n-point boxcar window in the supplied buffer. Does not
 *          allocate memory, so the user is responsible for ensuring that the
 *          destination buffer can hold the entire window
 *
 * @param n     The length of the window
 * @param dest  Buffer where the window is written. Buffer size must be at
 *              least n * sizeof(float)
 * @return      Error code, 0 on success
 */
Error_t
boxcar(unsigned n, float* dest);
Error_t
boxcarD(unsigned n, double* dest);


/** Generate a Hann window of length n
 *
 * @details Create an n-point Hann window in the supplied buffer. Does not
 *          allocate memory, so the user is responsible for ensuring that the
 *          destination buffer can hold the entire window
 *
 * @param n     The length of the window
 * @param dest  Buffer where the window is written. Buffer size must be at
 *              least n * sizeof(float)
 * @return      Error code, 0 on success
 */
Error_t
hann(unsigned n, float* dest);
Error_t
hannD(unsigned n, double* dest);

/** Generate a Hamming window of length n
 *
 * @details Create an n-point Hamming window in the supplied buffer. Does not
 *          allocate memory, so the user is responsible for ensuring that the
 *          destination buffer can hold the entire window
 *
 * @param n     The length of the window
 * @param dest  Buffer where the window is written. Buffer size must be at
 *              least n * sizeof(float)
 * @return      Error code, 0 on success
 */
Error_t
hamming(unsigned n, float* dest);

Error_t
hammingD(unsigned n, double* dest);


/** Generate a Blackman window of length n for given alpha
 *
 * @details Create an n-point Blackman window in the supplied buffer. Does not
 *          allocate memory, so the user is responsible for ensuring that the
 *          destination buffer can hold the entire window
 *
 * @param n     The length of the window
 * @param a     Alpha value
 * @param dest  Buffer where the window is written. Buffer size must be at
 *              least n * sizeof(float)
 * @return      Error code, 0 on success
 */
Error_t
blackman(unsigned n, float a, float* dest);

Error_t
blackmanD(unsigned n, double a, double* dest);

/** Generate a Tukey window of length n for given alpha
 *
 * @details Create an n-point Tukey window in the supplied buffer. Does not
 *          allocate memory, so the user is responsible for ensuring that the
 *          destination buffer can hold the entire window
 *
 * @param n     The length of the window
 * @param a     Alpha value
 * @param dest  Buffer where the window is written. Buffer size must be at
 *              least n * sizeof(float)
 * @return      Error code, 0 on success
 */
Error_t
tukey(unsigned n, float a, float* dest);

Error_t
tukeyD(unsigned n, double a, double* dest);


/** Generate a cosine window of length n
 *
 * @details Create an n-point Cosine window in the supplied buffer. Does not
 *          allocate memory, so the user is responsible for ensuring that the
 *          destination buffer can hold the entire window
 *
 * @param n     The length of the window
 * @param dest  Buffer where the window is written. Buffer size must be at
 *              least n * sizeof(float)
 * @return      Error code, 0 on success
 */
Error_t
cosine(unsigned n, float* dest);

Error_t
cosineD(unsigned n, double* dest);


/** Generate a Lanczos window of length n
 *
 * @details Creates an n-point Lanczos window in the supplied buffer. Does not
 *          allocate memory, so the user is responsible for ensuring that the
 *          destination buffer can hold the entire window
 *
 * @param n     The length of the window
 * @param dest  Buffer where the window is written. Buffer size must be at
 *              least n * sizeof(float)
 * @return      Error code, 0 on success
 */
Error_t
lanczos(unsigned n, float* dest);

Error_t
lanczosD(unsigned n, double* dest);


/** Generate a Bartlett window of length n
 *
 * @details Creates an n-point Bartlett window in the supplied buffer. Does not
 *          allocate memory, so the user is responsible for ensuring that the
 *          destination buffer can hold the entire window
 *
 * @param n     The length of the window
 * @param dest  Buffer where the window is written. Buffer size must be at
 *              least n * sizeof(float)
 * @return      Error code, 0 on success
 */
Error_t
bartlett(unsigned n, float* dest);

Error_t
bartlettD(unsigned n, double* dest);


/** Generate a Gaussian window of length n for given sigma
 *
 * @details Creates an n-point Gaussian window in the supplied buffer. Does not
 *          allocate memory, so the user is responsible for ensuring that the
 *          destination buffer can hold the entire window
 *
 * @param n     The length of the window
 * @param sigma Sigma value
 * @param dest  Buffer where the window is written. Buffer size must be at
 *              least n * sizeof(float)
 * @return      Error code, 0 on success
 */
Error_t
gaussian(unsigned n, float sigma, float* dest);

Error_t
gaussianD(unsigned n, double sigma, double* dest);


/** Generate a Bartlett-Hann window of length n
 *
 * @details Creates an n-point Bartlett-Hann window in the supplied buffer.
 *          Does notallocate memory, so the user is responsible for ensuring
 *          that the destination buffer can hold the entire window
 *
 * @param n     The length of the window
 * @param dest  Buffer where the window is written. Buffer size must be at
 *              least n * sizeof(float)
 * @return      Error code, 0 on success
 */
Error_t
bartlett_hann(unsigned n, float* dest);

Error_t
bartlett_hannD(unsigned n, double* dest);

/** Generate a Kaiser window of length n for given alpha
 *
 * @details Create an n-point Kaiser window in the supplied buffer. Does not
 *          allocate memory, so the user is responsible for ensuring that the
 *          destination buffer can hold the entire window
 *
 * @param n     The length of the window
 * @param a     Alpha value = (Beta / PI)
 * @param dest  Buffer where the window is written. Buffer size must be at
 *              least n * sizeof(float)
 * @return      Error code, 0 on success
 */
Error_t
kaiser(unsigned n, float a, float* dest);

Error_t
kaiserD(unsigned n, double a, double* dest);


/** Generate a Nuttall window of length n
 *
 * @details Create an n-point Nuttall window in the supplied buffer. Does not
 *          allocate memory, so the user is responsible for ensuring that the
 *          destination buffer can hold the entire window
 *
 * @param n     The length of the window
 * @param dest  Buffer where the window is written. Buffer size must be at
 *              least n * sizeof(float)
 * @return      Error code, 0 on success
 */
Error_t
nuttall(unsigned n, float* dest);

Error_t
nuttallD(unsigned n, double* dest);


/** Generate a Blackman-Harris window of length n
 *
 * @details Create an n-point Blackman-Harris window in the supplied buffer.
 *          Does not allocate memory, so the user is responsible for ensuring
 *          that the destination buffer can hold the entire window
 *
 * @param n     The length of the window
 * @param dest  Buffer where the window is written. Buffer size must be at
 *              least n * sizeof(float)
 * @return      Error code, 0 on success
 */
Error_t
blackman_harris(unsigned n, float* dest);

Error_t
blackman_harrisD(unsigned n, double* dest);


/** Generate a Blackman-Nuttall window of length n
 *
 * @details Create an n-point Blackman-Nuttall window in the supplied buffer.
 *          Does not allocate memory, so the user is responsible for ensuring
 *          that the destination buffer can hold the entire window
 *
 * @param n     The length of the window
 * @param dest  Buffer where the window is written. Buffer size must be at
 *              least n * sizeof(float)
 * @return      Error code, 0 on success
 */
Error_t
blackman_nuttall(unsigned n, float* dest);

Error_t
blackman_nuttallD(unsigned n, double* dest);


/** Generate a flat top window of length n
 *
 * @details Create an n-point flat top window in the supplied buffer. Does not
 *          allocate memory, so the user is responsible for ensuring that the
 *          destination buffer can hold the entire window
 *
 * @param n     The length of the window
 * @param dest  Buffer where the window is written. Buffer size must be at
 *              least n * sizeof(float)
 * @return      Error code, 0 on success
 */
Error_t
flat_top(unsigned n, float* dest);

Error_t
flat_topD(unsigned n, double* dest);


/** Generate a Poisson window of length n and given D
 *
 * @details Create an n-point Poisson    window in the supplied buffer. Does not
 *          allocate memory, so the user is responsible for ensuring that the
 *          destination buffer can hold the entire window
 *
 * @param n     The length of the window
 * @param D     Target decay in dB over 1/2 window length
 * @param dest  Buffer where the window is written. Buffer size must be at
 *              least n * sizeof(float)
 * @return      Error code, 0 on success
 */
Error_t
poisson(unsigned n, float D, float* dest);

Error_t
poissonD(unsigned n, double D, double* dest);




Error_t
chebyshev(unsigned n, float A, float* dest);

Error_t
chebyshevD(unsigned n, double A, double* dest);



/** Create a new WindowFunction
 *
 * @details Allocates memory and returns an initialized WindowFunction.
 *          Play nice and call WindowFunctionFree on the window when
 *          you're done with it.
 *
 * @param n     Number of points in the window.
 * @param type  Type of window function to generate.
 * @return          Error code, 0 on success
 */
WindowFunction*
WindowFunctionInit(unsigned n, Window_t type);

WindowFunctionD*
WindowFunctionInitD(unsigned n, Window_t type);


/** Free memory associated with a WindowFunction
 *
 * @details release all memory allocated by WindowFunctionInit for the
 *          supplied window.
 *
 * @param window    The window to free
 * @return          Error code, 0 on success
 */
Error_t
WindowFunctionFree(WindowFunction* window);

Error_t
WindowFunctionFreeD(WindowFunctionD* window);


/** Window a buffer of samples
 *
 * @details Applies the window to the buffer of samples passed to it
 *
 * @param window    The WindowFunction to use
 * @param outBuffer The buffer to write the output to
 * @param inBuffer  The buffer to filter
 * @param n_samples The number of samples to window
 * @return          Error code, 0 on success
 */
Error_t
WindowFunctionProcess(WindowFunction*   window,
                      float*            outBuffer,
                      const float*      inBuffer,
                      unsigned          n_samples);

Error_t
WindowFunctionProcessD(WindowFunctionD* window,
                       double*          outBuffer,
                       const double*    inBuffer,
                       unsigned         n_samples);


/** Create a new Upsampler
 *
 * @details Allocates memory and returns an initialized Upsampler with
 *          a given upsampling factor. Play nice and call UpsamplerFree
 *          on the filter whenyou're done with it.
 *
 * @param factor    Upsampling factor
 * @return          An initialized Upsampler
 */
Upsampler*
UpsamplerInit(ResampleFactor_t factor);

UpsamplerD*
UpsamplerInitD(ResampleFactor_t factor);


/** Free memory associated with a Upsampler
 *
 * @details release all memory allocated by Upsampler for the
 *          supplied filter.
 * @param upsampler Upsampler to free.
 * @return          Error code, 0 on success
 */
Error_t
UpsamplerFree(Upsampler* upsampler);

Error_t
UpsamplerFreeD(UpsamplerD* upsampler);


/** Flush upsampler state buffers
 *
 * @param upsampler Upsampler to flush.
 * @return          Error code, 0 on success
 */
Error_t
UpsamplerFlush(Upsampler* upsampler);

Error_t
UpsamplerFlushD(UpsamplerD* upsampler);


/** Upsample a buffer of samples
 *
 * @details Upsamples given buffer using sinc interpolation
 *
 * @param upsampler The Upsampler to use
 * @param outBuffer The buffer to write the output to
 * @param inBuffer  The buffer to filter
 * @param n_samples The number of samples to upsample
 * @return          Error code, 0 on success
 */
Error_t
UpsamplerProcess(Upsampler*     upsampler,
                 float*         outBuffer,
                 const float*   inBuffer,
                 unsigned       n_samples);

Error_t
UpsamplerProcessD(UpsamplerD*   upsampler,
                 double*        outBuffer,
                 const double*  inBuffer,
                 unsigned       n_samples);



/** Create a new FFTConfig
 *
 * @details Allocates memory and returns an initialized FFTConfig,
 *      which is used to store the FFT Configuration. Play nice and call
 *          FFTFree on it when you're done.
 *
 * @param length        length of the FFT. should be a power of 2.
 * @return        An initialized FFTConfig.
 */
FFTConfig*
FFTInit(unsigned length);

FFTConfigD*
FFTInitD(unsigned length);

/** Free memory associated with a FFTConfig
 *
 * @details release all memory allocated by FFTInit for the supplied
 *      fft configuration.
 *
 * @param fft       pointer to the FFTConfig to free.
 * @return      Error code, 0 on success
 */
Error_t
FFTFree(FFTConfig* fft);

Error_t
FFTFreeD(FFTConfigD* fft);


/** Calculate Real to Complex Forward FFT
 *
 * @details Calculates the magnitude of the real forward FFT of the data in
 *          inBuffer.
 *
 * @param fft       Pointer to the FFT configuration.
 * @param inBuffer  Input data. should be the same size as the fft.
 * @param real      Allocated buffer where the real part will be written. length
 *                  should be (fft->length/2).
 * @param imag      Allocated buffer where the imaginary part will be written. l
 *                  length should be (fft->length/2).
 * @return          Error code, 0 on success.
 */
Error_t
FFT_R2C(FFTConfig*      fft,
        const float*    inBuffer,
        float*          real,
        float*          imag);

Error_t
FFT_R2CD(FFTConfigD*    fft,
         const double*  inBuffer,
         double*        real,
         double*        imag);

Error_t
FFT_IR_R2C(FFTConfig*       fft,
           const float*     inBuffer,
           FFTSplitComplex  out);

Error_t
FFT_IR_R2CD(FFTConfigD*         fft,
            const double*       inBuffer,
            FFTSplitComplexD    out);


/** Calculate Complex to Real Inverse FFT
 *
 * @details Calculates the inverse FFT of the data in inBuffer.
 *
 * @param fft       Pointer to the FFT configuration.
 * @param inReal    Input real part. Length fft->length/2
 * @param inImag    Input imaginary part. Length fft->length/2
 * @param out       Allocated buffer where the signal will be written. length
 *                  should be fft->length.
 * @return          Error code, 0 on success.
 */
Error_t
IFFT_C2R(FFTConfig*    fft,
         const float*  inReal,
         const float*  inImag,
         float*        out);

Error_t
IFFT_C2RD(FFTConfigD*   fft,
          const double* inreal,
          const double* inImag,
          double*       out);


/** Perform Convolution using FFT*
 * @details convolve in1 with in2 and write results to dest
 * @param in1           First input to convolve.
 * @param in1_length    Length [samples] of in1.
 * @param in2           Second input to convolve.
 * @param in2_length    Length[samples] of second input.
 * @param dest          Output buffer. needs to be of length
 *                      in1_length + in2_length - 1
 * @return              Error code.
 */
Error_t
FFTConvolve(FFTConfig* fft,
            float       *in1,
            unsigned    in1_length,
            float       *in2,
            unsigned    in2_length,
            float       *dest);

Error_t
FFTConvolveD(FFTConfigD*    fft,
             const double*  in1,
             unsigned       in1_length,
             const double*  in2,
             unsigned       in2_length,
             double*        dest);


/** Perform Convolution using FFT*
 * @details Convolve in1 with IFFT(fft_ir) and write results to dest.
 *          This takes an already transformed kernel as the second argument, to
 *          be used in an LTI filter, where the FFT of the kernel can be pre-
 *          calculated.
 * @param in1           First input to convolve.
 * @param in1_length    Length [samples] of in1.
 * @param fft_ir        Second input to convolve (Already FFT'ed).
 * @param dest          Output buffer. needs to be of length
 *                      in1_length + in2_length - 1
 * @return              Error code.
 */
Error_t
FFTFilterConvolve(FFTConfig*        fft,
                  const float*      in,
                  unsigned          in_length,
                  FFTSplitComplex   fft_ir,
                  float*            dest);

Error_t
FFTFilterConvolveD(FFTConfigD*      fft,
                   const double*    in,
                   unsigned         in_length,
                   FFTSplitComplexD fft_ir,
                   double*          dest);

/** Just prints the complex output
 *
 */
Error_t
FFTdemo(FFTConfig* fft, float* buffer);


/** Create a new RBJFilter
 *
 * @details Allocates memory and returns an initialized RBJFilter.
 *			Play nice and call RBJFilterFree on the filter when you're
 *          done with it.
 *
 * @param type			The filter type
 * @param cutoff		The starting cutoff frequency to use
 * @param sampleRate	The sample rate in Samp/s
 * @return 				An initialized RBJFilter
 */
RBJFilter*
RBJFilterInit(Filter_t type, float cutoff, float sampleRate);

RBJFilterD*
RBJFilterInitD(Filter_t type, double cutoff,double sampleRate);


/** Free memory associated with a RBJFilter
 *
 * @details release all memory allocated by RBJFilterInit for the
 *			supplied filter.
 *
 * @param filter	RBJFilter to free
 * @return			Error code, 0 on success
 */
Error_t
RBJFilterFree(RBJFilter* filter);

Error_t
RBJFilterFreeD(RBJFilterD* filter);


/** Update RBJFilter type
 *
 * @details Update the filter type and recalculate filter coefficients.
 *
 * @param filter	RBJFilter to update
 * @param type		New filter type
 * @return			Error code, 0 on success
 */
Error_t
RBJFilterSetType(RBJFilter* filter, Filter_t type);

Error_t
RBJFilterSetTypeD(RBJFilterD* filter, Filter_t type);


/** Update RBJFilter Cutoff
 *
 * @details Update the filter cutoff/center frequency and recalculate filter
 *			coefficients.
 *
 * @param filter	RBJFilter to update
 * @param cutoff	New filter cutoff/center frequency
 * @return			Error code, 0 on success
 */
Error_t
RBJFilterSetCutoff(RBJFilter* filter, float cutoff);

Error_t
RBJFilterSetCutoffD(RBJFilterD* filterD, double cutoff);

/** Update RBJFilter Q
 *
 * @details Update the filter Q and recalculate filter coefficients.
 *
 * @param filter	RBJFilter to update
 * @param Q			New filter Q
 * @return			Error code, 0 on success
 */
Error_t
RBJFilterSetQ(RBJFilter* filter, float Q);

Error_t
RBJFilterSetQD(RBJFilterD* filter, double Q);



/** Update RBJFilter Parameters
 *
 * @details Update the filter Q and recalculate filter coefficients.
 *
 * @param filter	RBJFilter to update
 * @param type		New filter type
 * @param cutoff	New filter cutoff/center frequency
 * @param Q			New filter Q
 * @return			Error code, 0 on success
 */
Error_t
RBJFilterSetParams(RBJFilter*   filter,
                   Filter_t     type,
                   float        cutoff,
                   float        Q);

Error_t
RBJFilterSetParamsD(RBJFilterD* filter,
                   Filter_t     type,
                   double       cutoff,
                   double       Q);



/** Filter a buffer of samples
 * @details Uses an RBJ-style filter to filter input samples
 *
 * @param filter	The RBJFilter to use.
 * @param outBuffer	The buffer to write the output to.
 * @param inBuffer	The buffer to filter.
 * @param n_samples The number of samples to filter.
 * @return			Error code, 0 on success
 */
Error_t
RBJFilterProcess(RBJFilter*     filter,
                float*          outBuffer,
                const float*    inBuffer,
                unsigned        n_samples);

Error_t
RBJFilterProcessD(RBJFilterD*   filter,
                  double*       outBuffer,
                  const double* inBuffer,
                  unsigned      n_samples);


/** Flush filter state buffers
*
* @param filter    RBJFilter to flush.
* @return          Error code, 0 on success
*/
Error_t
RBJFilterFlush(RBJFilter* filter);

Error_t
RBJFilterFlushD(RBJFilterD* filter);


OnePole*
OnePoleInit(float cutoff, float sampleRate, Filter_t type);
    
OnePoleD*
OnePoleInitD(double cutoff, double sampleRate, Filter_t type);
    
OnePole*
OnePoleRawInit(float beta, float alpha);

OnePoleD*
OnePoleRawInitD(double beta, double alpha);

Error_t
OnePoleFree(OnePole *filter);
    
Error_t
OnePoleFreeD(OnePoleD *filter);

Error_t
OnePoleFlush(OnePole *filter);

Error_t
OnePoleFlushD(OnePoleD *filter);

Error_t
OnePoleSetType(OnePole* filter, Filter_t type);

Error_t
OnePoleSetTypeD(OnePoleD* filter, Filter_t type);

    
Error_t
OnePoleSetCutoff(OnePole* filter, float cutoff);
    
Error_t
OnePoleSetCutoffD(OnePoleD* filter, double cutoff);

    
Error_t
OnePoleSetSampleRate(OnePole* filter, float sampleRate);
    
Error_t
OnePoleSetSampleRateD(OnePoleD* filter, double sampleRate);
    
Error_t
OnePoleSetCoefficients(OnePole* filter, float* beta, float* alpha);

Error_t
OnePoleSetCoefficientsD(OnePoleD* filter, double* beta, double* alpha);
  
  Error_t
OnePoleProcess(OnePole*         filter,
                 float*         outBuffer,
                 const float*   inBuffer,
                 unsigned       n_samples);

Error_t
OnePoleProcessD(OnePoleD*   filter,
                 double*        outBuffer,
                 const double*  inBuffer,
                 unsigned       n_samples);
    
    
float
OnePoleTick(OnePole*    filter,
              float         inSample);
    
double
OnePoleTickD(OnePoleD*  filter,
               double       inSample);


    
float
OnePoleAlpha(OnePole* filter);

double
OnePoleAlphaD(OnePoleD* filter);
 
float
OnePoleBeta(OnePole* filter);

double
OnePoleBetaD(OnePoleD* filter);
    


PolySaturator*
PolySaturatorInit(float n);

PolySaturatorD*
PolySaturatorInitD(double n);

Error_t
PolySaturatorFree(PolySaturator* Saturator);

Error_t
PolySaturatorFreeD(PolySaturatorD* Saturator);

Error_t
PolySaturatorSetN(PolySaturator* saturator, float n);

Error_t
PolySaturatorSetND(PolySaturatorD* saturator, double n);

Error_t
PolySaturatorProcess(PolySaturator* saturator,
                     float*         out_buffer,
                     const float*   in_buffer,
                     unsigned       n_samples);

Error_t
PolySaturatorProcessD(PolySaturatorD* saturator,
                      double*         out_buffer,
                      const double*   in_buffer,
                      unsigned        n_samples);


float
PolySaturatorTick(PolySaturator* saturator, float in_sample);

double
PolySaturatorTickD(PolySaturatorD* saturator, double in_sample);


/** Create a new Opto
 *
 * @details Allocates memory and returns an initialized Opto.
 *			Play nice and call OptoFree when you're done with it.
 *
 * @param opto_type         Optocoupler model type.
 * @param delay             Amount of delay in the optocoupler. The halfway
 *                          point is a good approximation of the actual device,
 *                          higher and lower values are based on model
 *                          extrapolation. The effect is most pronounced with
 *                          the LDR(Vactrol) model as the LDR response is slow.
 *                          while this is not a realistic parameter with higher-
 *                          bandwidth models, higher settings of this parameter
 *                          result in an artifically exaggerated effect.
 * @param sample_rate       system sampling rate.
 * @return                  An initialized Opto
 */
/* Opto *******************************************************************/
struct Opto
{
    Opto_t      type;           //model type
    float       sample_rate;
    float       previous;
    float       delay;
    float       on_cutoff;
    float       off_cutoff;
    char        delta_sign;     // sign of signal dv/dt
    OnePole*  lp;
};



struct OptoD
{
    Opto_t      type;           //model type
    double      sample_rate;
    double      previous;
    double      delay;
    double      on_cutoff;
    double      off_cutoff;
    char        delta_sign;     // sign of signal dv/dt
    OnePoleD* lp;
};


Opto*
OptoInit(Opto_t opto_type, float delay, float sample_rate);

OptoD*
OptoInitD(Opto_t opto_type, double delay, double sample_rate);


Error_t
OptoFree(Opto* optocoupler);

Error_t
OptoFreeD(OptoD* optocoupler);



Error_t
OptoSetDelay(Opto* optocoupler, float delay);

Error_t
OptoSetDelayD(OptoD* optocoupler, double delay);



Error_t
OptoProcess(Opto*           optocoupler,
            float*          out_buffer,
            const float*    in_buffer,
            unsigned        n_samples);

Error_t
OptoProcessD(OptoD*         optocoupler,
             double*        out_buffer,
             const double*  in_buffer,
             unsigned       n_samples);

float
OptoTick(Opto* optocoupler, float in_sample);

double
OptoTickD(OptoD* optocoupler, double in_sample);




/** Create a new Tape
 */
Tape*
TapeInit(TapeSpeed speed, float saturation, float hysteresis, float flutter, float sample_rate);



/** Free memory associated with a Tape
 *
 * @details release all memory allocated by TapeInit for the
 *			given saturator.
 *
 * @param saturator Tape to free
 * @return			Error code, 0 on success
 */
Error_t
TapeFree(Tape* tape);

Error_t
TapeSetSpeed(Tape* tape, TapeSpeed speed);

Error_t
TapeSetSaturation(Tape* tape, float saturation);

Error_t
TapeSetHysteresis(Tape* tape, float hysteresis);

Error_t
TapeSetFlutter(Tape* tape, float flutter);

float
TapeGetSaturation(Tape* tape);

float
TapeGetHysteresis(Tape* tape);

Error_t
TapeProcess(Tape*           tape,
            float*          out_buffer,
            const float*    in_buffer,
            unsigned        n_samples);



/** Process a single sample
 * @details Uses a Tape model to process an input sample
 *
 * @param taoe      The Tape to use.
 * @param in_sample	The sample to process.
 * @return			A processed sample.
 */
float
TapeTick(Tape* tape, float in_sample);




