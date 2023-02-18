
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

#define USE_FFTW_FFT
#define E_INV (0.36787944117144233)

/* Macro Constants ************************************************************/

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

#define USE_FFT_CONVOLUTION_LENGTH (128)
// Sqrt(2)/2
#define FILT_Q (0.70710681186548)

/* pow(10, (-12./20.)) */
#define K12_REF (0.25118864315095801309)

/* pow(10, (-12./20.)); */
#define K14_REF (0.19952623149688797355)

/* pow(10, (-20./20.)); */
#define K20_REF (0.1)


#define PREFILTER_FC    (1500.12162162162167078350)
#define PREFILTER_GAIN  (3.99976976976976983380)
#define PREFILTER_Q     (sqrt(2.0)/2.0) // Close enough....
#define RLBFILTER_FC    (40.2802802803)
#define RLBFILTER_Q     (0.92792792793)
#define GATE_LENGTH_S   (0.4)
#define GATE_OVERLAP    (0.75)

// The higher this value is, the more the saturation amount is affected by tape
// speed
#define SPEED_SATURATION_COEFF (0.0)

#define N_FLUTTER_COMPONENTS (11)

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

/* LadderFilter ********************************************************/
struct LadderFilterD
{
    double y[4];
    double w[4];
    double Vt;           // transistor treshold voltage [V]
    double sample_rate;
    double cutoff;
    double resonance;
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





// REF: http://www.manquen.net/audio/docs/Flutter%20database%2002-8-28.htm
static const float flutterRateBase[N_FLUTTER_COMPONENTS] =
{
    1.01,   // Outer Tension Sense
    2.52,   // Inner Tension Sense
    0.80,   // Pre-stabilizer
    1.01,   // Left lifter
    3.11,   // Record scrape idler
    3.33,   // Capstan
    0.81,   // Pinch roller
    1.01,   // Right lifter
    0.80,   // Motion Sensor
    2.52,   // Inner tension sense
    1.01    // Outer tension sense
};


static inline float
calculate_n(float saturation, TapeSpeed speed)
{
    // Tape speed dependent saturation.
    float n = ((50 * (1-SPEED_SATURATION_COEFF)) + \
            ((unsigned)speed * 50 * SPEED_SATURATION_COEFF)) \
    * powf((1.0075 - saturation), 2.);
    printf("N: %1.5f\n", n);
    return n;
}

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


static const float polyphase2xfilter0[64] =
{
    -0.005650002975f, 0.005204734392f, 0.001785006723f, -0.002015689854f,
    -0.001451580552f, 0.002462821314f, 0.002242911141f, -0.002549842233f,
    -0.003116151551f, 0.00251203333f, 0.004146839026f, -0.002282319823f,
    -0.005326004699f, 0.001805614913f, 0.006656893063f, -0.001012795139f,
    -0.008139784448f, -0.0001834870782f, 0.009799872525f, 0.001914215391f,
    -0.01165371668f, -0.004402380902f, 0.01379979588f, 0.008041893132f,
    -0.01644209214f, -0.01370095555f, 0.02012135088f, 0.02370615676f,
    -0.02660541236f, -0.0471835807f, 0.04586292431f, 0.1908444166f,
    0.2320233136f, 0.1216813177f, -0.01510138717f, -0.04808133841f,
    0.002177739516f, 0.02961835265f, 0.002410362242f, -0.02067817375f,
    -0.004595956765f, 0.01512517128f, 0.005689772312f, -0.01120857149f,
    -0.006160198245f, 0.008245054632f, 0.006234389264f, -0.005928888917f,
    -0.006025022827f, 0.004076622892f, 0.005621950142f, -0.002611333271f,
    -0.005095448345f, 0.001461869455f, 0.004489722662f, -0.0005810860312f,
    -0.00384558877f, -6.323363777e-05f, 0.00320897717f, 0.0005431400496f,
    -0.002466904931f, -0.0003506442299f, 0.003678260138f, 0.007207856979f
};

static const double polyphase2xfilter0D[64] =
{
    -0.005650002975, 0.005204734392, 0.001785006723, -0.002015689854,
    -0.001451580552, 0.002462821314, 0.002242911141, -0.002549842233,
    -0.003116151551, 0.00251203333, 0.004146839026, -0.002282319823,
    -0.005326004699, 0.001805614913, 0.006656893063, -0.001012795139,
    -0.008139784448, -0.0001834870782, 0.009799872525, 0.001914215391,
    -0.01165371668, -0.004402380902, 0.01379979588, 0.008041893132,
    -0.01644209214, -0.01370095555, 0.02012135088, 0.02370615676,
    -0.02660541236, -0.0471835807, 0.04586292431, 0.1908444166,
    0.2320233136, 0.1216813177, -0.01510138717, -0.04808133841,
    0.002177739516, 0.02961835265, 0.002410362242, -0.02067817375,
    -0.004595956765, 0.01512517128, 0.005689772312, -0.01120857149,
    -0.006160198245, 0.008245054632, 0.006234389264, -0.005928888917,
    -0.006025022827, 0.004076622892, 0.005621950142, -0.002611333271,
    -0.005095448345, 0.001461869455, 0.004489722662, -0.0005810860312,
    -0.00384558877, -6.323363777e-05, 0.00320897717, 0.0005431400496,
    -0.002466904931, -0.0003506442299, 0.003678260138, 0.007207856979f
};

static const float polyphase2xfilter1[64] =
{
    0.007207856979f, 0.003678260138f, -0.0003506442299f, -0.002466904931f,
    0.0005431400496f, 0.00320897717f, -6.323363777e-05f, -0.00384558877f,
    -0.0005810860312f, 0.004489722662f, 0.001461869455f, -0.005095448345f,
    -0.002611333271f, 0.005621950142f, 0.004076622892f, -0.006025022827f,
    -0.005928888917f, 0.006234389264f, 0.008245054632f, -0.006160198245f,
    -0.01120857149f, 0.005689772312f, 0.01512517128f, -0.004595956765f,
    -0.02067817375f, 0.002410362242f, 0.02961835265f, 0.002177739516f,
    -0.04808133841f, -0.01510138717f, 0.1216813177f, 0.2320233136f,
    0.1908444166f, 0.04586292431f, -0.0471835807f, -0.02660541236f,
    0.02370615676f, 0.02012135088f, -0.01370095555f, -0.01644209214f,
    0.008041893132f, 0.01379979588f, -0.004402380902f, -0.01165371668f,
    0.001914215391f, 0.009799872525f, -0.0001834870782f, -0.008139784448f,
    -0.001012795139f, 0.006656893063f, 0.001805614913f, -0.005326004699f,
    -0.002282319823f, 0.004146839026f, 0.00251203333f, -0.003116151551f,
    -0.002549842233f, 0.002242911141f, 0.002462821314f, -0.001451580552f,
    -0.002015689854f, 0.001785006723f, 0.005204734392f, -0.005650002975f
};

static const double polyphase2xfilter1D[64] =
{
    0.007207856979, 0.003678260138, -0.0003506442299, -0.002466904931,
    0.0005431400496, 0.00320897717, -6.323363777e-05, -0.00384558877,
    -0.0005810860312, 0.004489722662, 0.001461869455, -0.005095448345,
    -0.002611333271, 0.005621950142, 0.004076622892, -0.006025022827,
    -0.005928888917, 0.006234389264, 0.008245054632, -0.006160198245,
    -0.01120857149, 0.005689772312, 0.01512517128, -0.004595956765,
    -0.02067817375, 0.002410362242, 0.02961835265, 0.002177739516,
    -0.04808133841, -0.01510138717, 0.1216813177, 0.2320233136,
    0.1908444166, 0.04586292431, -0.0471835807, -0.02660541236,
    0.02370615676, 0.02012135088, -0.01370095555, -0.01644209214,
    0.008041893132, 0.01379979588, -0.004402380902, -0.01165371668,
    0.001914215391, 0.009799872525, -0.0001834870782, -0.008139784448,
    -0.001012795139, 0.006656893063, 0.001805614913, -0.005326004699,
    -0.002282319823, 0.004146839026, 0.00251203333, -0.003116151551,
    -0.002549842233, 0.002242911141, 0.002462821314, -0.001451580552,
    -0.002015689854, 0.001785006723, 0.005204734392, -0.005650002975f
};

static const float* polyphase2x[2] =
{
    polyphase2xfilter0,
    polyphase2xfilter1
};

static const double* polyphase2xD[2] =
{
    polyphase2xfilter0D,
    polyphase2xfilter1D
};

/******************************************************************************
 * 4x Polyphase Coefficients
 *****************************************************************************/


/* polyphase coefficients */
static const float polyphase4xfilter0[64] =
{
    -0.00000550193565312270, -0.00023989303973949967, -0.00023999831743559219,
    0.00041857898393486150, -0.00044576332941102142, 0.00037528520356352324,
    -0.00022679507241807901, 0.00001005360016564356, 0.00026568139742544350,
    -0.00058763442732559881, 0.00093780019181280015, -0.00129247065827479230,
    0.00162258415863347420, -0.00189469852631972690, 0.00207245646869515550,
    -0.00211841938609982860, 0.00199613334012931730, -0.00167227201341311420,
    0.00111868364406304500, -0.00031415360657057714, -0.00075432069698170971,
    0.00209096873067570130, -0.00369167469794095570, 0.00554578705976306860,
    -0.00763998795040821650, 0.00996561045112854400, -0.01253265629179869600,
    0.01539915877460620900, -0.01874255508564667800, 0.02307428102189461700,
    -0.03013132631750860600, 0.04986002935449357400, 0.22768231662660360000,
    -0.01055445373187310400,-0.00242922862603297970, 0.00692925189763634630,
    -0.00888569409203613170, 0.00962397381607586060, -0.00963109606102358420,
    0.00915008850822176390, -0.00833349334762865240, 0.00729397293927114360,
    -0.00612312239712655120, 0.00489859007829948190, -0.00368645579334862150,
    0.00254166657955889630, -0.00150780221571094790, 0.00061679932385032938,
    0.00011104639797060849, -0.00066665041677132217, 0.00105123825819202600,
    -0.00127500043910553090, 0.00135538037357331550, -0.00131510562601075520,
    0.00118009349357532760, -0.00097737277931345037, 0.00073317235020506322,
    -0.00047134868059595163, 0.00021240498037212825, 0.00002636709946986627,
    -0.00022811348535030983, 0.00036568444376011244, -0.00035938108509561283,
    -0.00014500093992842017
};

static const double polyphase4xfilter0D[64] =
{
    -0.00000550193565312270, -0.00023989303973949967, -0.00023999831743559219,
    0.00041857898393486150, -0.00044576332941102142, 0.00037528520356352324,
    -0.00022679507241807901, 0.00001005360016564356, 0.00026568139742544350,
    -0.00058763442732559881, 0.00093780019181280015, -0.00129247065827479230,
    0.00162258415863347420, -0.00189469852631972690, 0.00207245646869515550,
    -0.00211841938609982860, 0.00199613334012931730, -0.00167227201341311420,
    0.00111868364406304500, -0.00031415360657057714, -0.00075432069698170971,
    0.00209096873067570130, -0.00369167469794095570, 0.00554578705976306860,
    -0.00763998795040821650, 0.00996561045112854400, -0.01253265629179869600,
    0.01539915877460620900, -0.01874255508564667800, 0.02307428102189461700,
    -0.03013132631750860600, 0.04986002935449357400, 0.22768231662660360000,
    -0.01055445373187310400,-0.00242922862603297970, 0.00692925189763634630,
    -0.00888569409203613170, 0.00962397381607586060, -0.00963109606102358420,
    0.00915008850822176390, -0.00833349334762865240, 0.00729397293927114360,
    -0.00612312239712655120, 0.00489859007829948190, -0.00368645579334862150,
    0.00254166657955889630, -0.00150780221571094790, 0.00061679932385032938,
    0.00011104639797060849, -0.00066665041677132217, 0.00105123825819202600,
    -0.00127500043910553090, 0.00135538037357331550, -0.00131510562601075520,
    0.00118009349357532760, -0.00097737277931345037, 0.00073317235020506322,
    -0.00047134868059595163, 0.00021240498037212825, 0.00002636709946986627,
    -0.00022811348535030983, 0.00036568444376011244, -0.00035938108509561283,
    -0.00014500093992842017
};

static const float polyphase4xfilter1[64] =
{
    -0.00002595797444661836, -0.00033159406371634453, -0.00004003845164712390,
    0.00030958725035557503, -0.00048381356303793505, 0.00057743931777250439,
    -0.00059306657324951866, 0.00052705576764733731, -0.00037460728718835209,
    0.00013324359186715783, 0.00019477642628928376, -0.00060079846081130246,
    0.00106874128368349620, -0.00157464880438700760, 0.00208680496872009270,
    -0.00256637951271162030, 0.00296854274730440780, -0.00324394543205795480,
    0.00334040876665658980, -0.00320460751778328100, 0.00278344742860470220,
    -0.00202471380915492540, 0.00087634989166469565, 0.00071670937456010163,
    -0.00282016652750145290, 0.00552671676047757250, -0.00899126832692335140,
    0.01351254773966352500, -0.01974924724155863900, 0.02942757112618565200,
    -0.04865417259070546600, 0.12300896236825336000, 0.18879981233367804000,
    -0.04448897023269413700, 0.02045323204195124300, -0.01000809975965635500,
    0.00403943112423317280, -0.00022851694087257357, -0.00228829657600720330,
    0.00391514880477506020, -0.00487826635621856080, 0.00532621442446622270,
    -0.00537080289051146410, 0.00510495108219325300, -0.00461045696704973870,
    0.00396098100949913380, -0.00322269618243638470, 0.00245385936712340340,
    -0.00170400440703556350, 0.00101317004140845660, -0.00041140885181395541,
    -0.00008128595117489710, 0.00045457750318969565, -0.00070690655287011350,
    0.00084418158457462541, -0.00087811261976715406, 0.00082428251064974034,
    -0.00070006187684540366, 0.00052250550883928834, -0.00030652573311330681,
    0.00006442703633721771, 0.00018787119054113082, -0.00038361778068667259,
    -0.00007100715719020498
};

static const double polyphase4xfilter1D[64] =
{
    -0.00002595797444661836, -0.00033159406371634453, -0.00004003845164712390,
    0.00030958725035557503, -0.00048381356303793505, 0.00057743931777250439,
    -0.00059306657324951866, 0.00052705576764733731, -0.00037460728718835209,
    0.00013324359186715783, 0.00019477642628928376, -0.00060079846081130246,
    0.00106874128368349620, -0.00157464880438700760, 0.00208680496872009270,
    -0.00256637951271162030, 0.00296854274730440780, -0.00324394543205795480,
    0.00334040876665658980, -0.00320460751778328100, 0.00278344742860470220,
    -0.00202471380915492540, 0.00087634989166469565, 0.00071670937456010163,
    -0.00282016652750145290, 0.00552671676047757250, -0.00899126832692335140,
    0.01351254773966352500, -0.01974924724155863900, 0.02942757112618565200,
    -0.04865417259070546600, 0.12300896236825336000, 0.18879981233367804000,
    -0.04448897023269413700, 0.02045323204195124300, -0.01000809975965635500,
    0.00403943112423317280, -0.00022851694087257357, -0.00228829657600720330,
    0.00391514880477506020, -0.00487826635621856080, 0.00532621442446622270,
    -0.00537080289051146410, 0.00510495108219325300, -0.00461045696704973870,
    0.00396098100949913380, -0.00322269618243638470, 0.00245385936712340340,
    -0.00170400440703556350, 0.00101317004140845660, -0.00041140885181395541,
    -0.00008128595117489710, 0.00045457750318969565, -0.00070690655287011350,
    0.00084418158457462541, -0.00087811261976715406, 0.00082428251064974034,
    -0.00070006187684540366, 0.00052250550883928834, -0.00030652573311330681,
    0.00006442703633721771, 0.00018787119054113082, -0.00038361778068667259,
    -0.00007100715719020498
};


static const float polyphase4xfilter2[64] =
{
    -0.00007100715719020498, -0.00038361778068667259, 0.00018787119054113082, 0.00006442703633721771,
    -0.00030652573311330681, 0.00052250550883928834, -0.00070006187684540366, 0.00082428251064974034,
    -0.00087811261976715406, 0.00084418158457462541, -0.00070690655287011350, 0.00045457750318969565,
    -0.00008128595117489710, -0.00041140885181395541, 0.00101317004140845660, -0.00170400440703556350,
    0.00245385936712340340, -0.00322269618243638470, 0.00396098100949913380, -0.00461045696704973870,
    0.00510495108219325300, -0.00537080289051146410, 0.00532621442446622270, -0.00487826635621856080,
    0.00391514880477506020, -0.00228829657600720330, -0.00022851694087257357, 0.00403943112423317280,
    -0.01000809975965635500, 0.02045323204195124300, -0.04448897023269413700, 0.18879981233367804000,
    0.12300896236825336000, -0.04865417259070546600, 0.02942757112618565200, -0.01974924724155863900,
    0.01351254773966352500, -0.00899126832692335140, 0.00552671676047757250, -0.00282016652750145290,
    0.00071670937456010163, 0.00087634989166469565, -0.00202471380915492540, 0.00278344742860470220,
    -0.00320460751778328100, 0.00334040876665658980, -0.00324394543205795480, 0.00296854274730440780,
    -0.00256637951271162030, 0.00208680496872009270, -0.00157464880438700760, 0.00106874128368349620,
    -0.00060079846081130246, 0.00019477642628928376, 0.00013324359186715783, -0.00037460728718835209,
    0.00052705576764733731, -0.00059306657324951866, 0.00057743931777250439, -0.00048381356303793505,
    0.00030958725035557503, -0.00004003845164712390, -0.00033159406371634453, -0.00002595797444661836
};

static const double polyphase4xfilter2D[64] =
{
    -0.00007100715719020498, -0.00038361778068667259, 0.00018787119054113082, 0.00006442703633721771,
    -0.00030652573311330681, 0.00052250550883928834, -0.00070006187684540366, 0.00082428251064974034,
    -0.00087811261976715406, 0.00084418158457462541, -0.00070690655287011350, 0.00045457750318969565,
    -0.00008128595117489710, -0.00041140885181395541, 0.00101317004140845660, -0.00170400440703556350,
    0.00245385936712340340, -0.00322269618243638470, 0.00396098100949913380, -0.00461045696704973870,
    0.00510495108219325300, -0.00537080289051146410, 0.00532621442446622270, -0.00487826635621856080,
    0.00391514880477506020, -0.00228829657600720330, -0.00022851694087257357, 0.00403943112423317280,
    -0.01000809975965635500, 0.02045323204195124300, -0.04448897023269413700, 0.18879981233367804000,
    0.12300896236825336000, -0.04865417259070546600, 0.02942757112618565200, -0.01974924724155863900,
    0.01351254773966352500, -0.00899126832692335140, 0.00552671676047757250, -0.00282016652750145290,
    0.00071670937456010163, 0.00087634989166469565, -0.00202471380915492540, 0.00278344742860470220,
    -0.00320460751778328100, 0.00334040876665658980, -0.00324394543205795480, 0.00296854274730440780,
    -0.00256637951271162030, 0.00208680496872009270, -0.00157464880438700760, 0.00106874128368349620,
    -0.00060079846081130246, 0.00019477642628928376, 0.00013324359186715783, -0.00037460728718835209,
    0.00052705576764733731, -0.00059306657324951866, 0.00057743931777250439, -0.00048381356303793505,
    0.00030958725035557503, -0.00004003845164712390, -0.00033159406371634453, -0.00002595797444661836
};

static const float polyphase4xfilter3[64] =
{
    -0.00014500093992842017, -0.00035938108509561283, 0.00036568444376011244, -0.00022811348535030983,
    0.00002636709946986627, 0.00021240498037212825, -0.00047134868059595163, 0.00073317235020506322,
    -0.00097737277931345037, 0.00118009349357532760, -0.00131510562601075520, 0.00135538037357331550,
    -0.00127500043910553090, 0.00105123825819202600, -0.00066665041677132217, 0.00011104639797060849,
    0.00061679932385032938, -0.00150780221571094790, 0.00254166657955889630, -0.00368645579334862150,
    0.00489859007829948190, -0.00612312239712655120, 0.00729397293927114360, -0.00833349334762865240,
    0.00915008850822176390, -0.00963109606102358420, 0.00962397381607586060, -0.00888569409203613170,
    0.00692925189763634630, -0.00242922862603297970, -0.01055445373187310400, 0.22768231662660360000,
    0.04986002935449357400, -0.03013132631750860600, 0.02307428102189461700, -0.01874255508564667800,
    0.01539915877460620900, -0.01253265629179869600, 0.00996561045112854400, -0.00763998795040821650,
    0.00554578705976306860, -0.00369167469794095570,  0.00209096873067570130, -0.00075432069698170971,
    -0.00031415360657057714, 0.00111868364406304500, -0.00167227201341311420, 0.00199613334012931730,
    -0.00211841938609982860, 0.00207245646869515550, -0.00189469852631972690, 0.00162258415863347420,
    -0.00129247065827479230, 0.00093780019181280015, -0.00058763442732559881, 0.00026568139742544350,
    0.00001005360016564356, -0.00022679507241807901, 0.00037528520356352324, -0.00044576332941102142,
    0.00041857898393486150, -0.00023999831743559219, -0.00023989303973949967, -0.00000550193565312270
};

static const double polyphase4xfilter3D[64] =
{
    -0.00014500093992842017, -0.00035938108509561283, 0.00036568444376011244, -0.00022811348535030983,
    0.00002636709946986627, 0.00021240498037212825, -0.00047134868059595163, 0.00073317235020506322,
    -0.00097737277931345037, 0.00118009349357532760, -0.00131510562601075520, 0.00135538037357331550,
    -0.00127500043910553090, 0.00105123825819202600, -0.00066665041677132217, 0.00011104639797060849,
    0.00061679932385032938, -0.00150780221571094790, 0.00254166657955889630, -0.00368645579334862150,
    0.00489859007829948190, -0.00612312239712655120, 0.00729397293927114360, -0.00833349334762865240,
    0.00915008850822176390, -0.00963109606102358420, 0.00962397381607586060, -0.00888569409203613170,
    0.00692925189763634630, -0.00242922862603297970, -0.01055445373187310400, 0.22768231662660360000,
    0.04986002935449357400, -0.03013132631750860600, 0.02307428102189461700, -0.01874255508564667800,
    0.01539915877460620900, -0.01253265629179869600, 0.00996561045112854400, -0.00763998795040821650,
    0.00554578705976306860, -0.00369167469794095570,  0.00209096873067570130, -0.00075432069698170971,
    -0.00031415360657057714, 0.00111868364406304500, -0.00167227201341311420, 0.00199613334012931730,
    -0.00211841938609982860, 0.00207245646869515550, -0.00189469852631972690, 0.00162258415863347420,
    -0.00129247065827479230, 0.00093780019181280015, -0.00058763442732559881, 0.00026568139742544350,
    0.00001005360016564356, -0.00022679507241807901, 0.00037528520356352324, -0.00044576332941102142,
    0.00041857898393486150, -0.00023999831743559219, -0.00023989303973949967, -0.00000550193565312270
};


static const float* polyphase4x[4] =
{
    polyphase4xfilter0,
    polyphase4xfilter1,
    polyphase4xfilter2,
    polyphase4xfilter3
};

static const double* polyphase4xD[4] =
{
    polyphase4xfilter0D,
    polyphase4xfilter1D,
    polyphase4xfilter2D,
    polyphase4xfilter3D
};

/******************************************************************************
 * 8x Polyphase Coefficients
 *****************************************************************************/

static const float polyphase8xfilter0[64] =
{
    -0.00000065139115950063, -0.00004102501562491190, 0.00007454552789152059, 0.00067572747586549999,
    -0.00006343701412881912, -0.00060210193780960727, 0.00085677626062353681, -0.00068374678593982705,
    0.00022508973131157744, 0.00034416140819020045, -0.00086697375710736936, 0.00122392123868685560,
    -0.00133970651621087830, 0.00118386211690219210, -0.00076830302548444487, 0.00014256121428106461,
    0.00061301956994263369, -0.00139595818775825930, 0.00209182396262220810, -0.00258537407535367780,
    0.00277185772937339320, -0.00256754278849754100, 0.00191845357598483200, -0.00080635046129673850,
    -0.00074884899307910678, 0.00269106189201281710, -0.00493523685432864770, 0.00738159244609435180,
    -0.00994256584967902380, 0.01260511444709402700, -0.01564820095328426800, 0.02116771695102326000,
    0.11058686170616823000, 0.00645758092866270850, -0.00957449662395692440, 0.00984001261556627540,
    -0.00904206471599113170, 0.00765208391982513670, -0.00592783312693009010, 0.00406863057583288240,
    -0.00224523360091431110, 0.00060189327685675881, 0.00074853435796320012, -0.00172997883882847830,
    0.00230658176248078620, -0.00248303216510402730, 0.00230162693883983070, -0.00183622557500364680,
    0.00118365975641694390, -0.00045337424289750575, -0.00024385296497564968, 0.00080732204284516727,
    -0.00115732657922609960, 0.00124514852293933480, -0.00106196556316837370, 0.00064655676276784320,
    -0.00009140135720076620, -0.00045476803926232661, 0.00079448278839194797, -0.00071843423633693635,
    0.00012637883571250220, 0.00062442144337030445, 0.00002569631187128541, -0.00003279706766106838
};

static const float polyphase8xfilter1[64] =
{
    -0.00000182984414357783, -0.00004831055001691508, 0.00013626225291498421, 0.00070065829999757436,
    -0.00025547582772182008, -0.00043265212159530387, 0.00084162909240903571, -0.00085050570414972073,
    0.00052522292321950030, -0.00000143693453915044, -0.00057289470445241102, 0.00106555149724676140,
    -0.00137406123500243830, 0.00143269136290550630, -0.00121536388445158540, 0.00073554851430502580,
    -0.00004348492312229754, -0.00077915084816193418, 0.00162691111436317150, -0.00238015291908662140,
    0.00291584964226445760, -0.00311871005955389980, 0.00289147118196876480, -0.00216307979887024520,
    0.00089338246990385678, 0.00092760292983492458, -0.00328877553092718580, 0.00618233379133788660,
    -0.00966417201297737480, 0.01403037645222736200, -0.02054581611003392100, 0.03746277721486017400,
    0.10613402761530316000, -0.00584676809795404460, -0.00311138376108936220, 0.00615580584287418050,
    -0.00714024403544936260, 0.00701012531795208030, -0.00617923068231396570, 0.00491363118474887150,
    -0.00342501037656084200, 0.00189442023648568910, -0.00047490028164087116, -0.00071241642862137494,
    0.00158317444638181860, -0.00209304708081012490, 0.00223830512847985300, -0.00205327074145764380,
    0.00160476729027759600, -0.00098402602826872031, 0.00029670356071559174, 0.00034829156389369192,
    -0.00085055664694030269, 0.00113011259388196740, -0.00113937588749536430, 0.00087531359141694673,
    -0.00039197851653214976, -0.00018732750967400416, 0.00066293773127700027, -0.00077461306480611420,
    0.00030143757577000946, 0.00055349978393511490, -0.00001032580058558277, -0.00002472111034662811
};

static const float polyphase8xfilter2[64] =
{
    -0.00000388695825192831, -0.00005318884282459066, 0.00020974856947313336, 0.00069349157210732459,
    -0.00043577094717956485, -0.00022221919156823275, 0.00074782136190461302, -0.00093664861877437256,
    0.00077740863525184738, -0.00035428397707766530, -0.00020962086499279835, 0.00078450145550032285,
    -0.00125298740630276000, 0.00152215495203344830, -0.00153121928051926810, 0.00125599179151124560,
    -0.00071043859554279257, -0.00005464802619675554, 0.00095581822472680811, -0.00188376295526464770,
    0.00271226205959396690, -0.00330848547687067520, 0.00354352391311250860, -0.00330168644834485760,
    0.00248673899386953140, -0.00102220275357614380, -0.00115930523989435050, 0.00414981937115260680,
    -0.00815969289584616080, 0.01381682701592122300, -0.02350820651622280100,0.05434387608734767000,
    0.09755145975376086600, -0.01516237260357852300, 0.00300145106345972670, 0.00204115962512772160,
    -0.00451741705979632340, 0.00557638253007973400, -0.00569207291393189830, 0.00514726404042791740,
    -0.00416047157856245470, 0.00292295670535206620, -0.00160651004059010120, 0.00036069553428501186,
    0.00069324910330136524, -0.00146835891655566310, 0.00191604280034704920, -0.00202740273577447240,
    0.00183144503810963020, -0.00139028769170792730, 0.00079175249424842480, -0.00013986723178465493,
    -0.00045621477010817133, 0.00089451555246586410, -0.00109375535771801020, 0.00100851752563541820,
    -0.00064697362622190139, 0.00009197902579005798, 0.00047618788828089582, -0.00076894092872930046,
    0.00045149583778780219, 0.00047008811115116311, -0.00003446958495946170, -0.00001752354118849114
};

static const float polyphase8xfilter3[64] =
{
    -0.00000707856005788827, -0.00005386828034834762, 0.00029260597635645714, 0.00065022085359971421,
    -0.00059001569340461763, 0.00001271442624327296, 0.00058158960488687678, -0.00093123666109248718,
    0.00095408718371994903, -0.00067661601603670300, 0.00018452036040974686, 0.00040938095642024030,
    -0.00098659778327271500, 0.00143817315098353240, -0.00167597953009304670, 0.00164106173258977890,
    -0.00130928010530110450, 0.00069413459212034961, 0.00015338377641278773, -0.00114795228773180200,
    0.00217587137805511180, -0.00310325220590675600, 0.00378546088727036750, -0.00407634981728579190,
    0.00383521477176884070, -0.00292792329706461190, 0.00121555318378888090, 0.00148596888578445710,
    -0.00552647396504433920, 0.01183152935877091100, -0.02388014901113717800, 0.07071205074526545900,
    0.08545741447609654700, -0.02118099310325147900, 0.00813695469487032970, -0.00200465909461610940,
    -0.00151921945966949170, 0.00355693608056953160, -0.00455482533831069740, 0.00476518663684880220,
    -0.00438038015325593890, 0.00357492186862908880 -0.00251606896499480740, 0.00136193071021434250,
    -0.00025460237565429417, -0.00068788009050978823, 0.00137873756292502610, -0.00176750548092106620,
    0.00184225028480482200, -0.00162869371539639560, 0.00118633130281230620, -0.00060183978921851905,
    -0.00001987860045449669, 0.00056676217416478942, -0.00093289737900114494, 0.00103478202101450390,
    -0.00083221312461457934, 0.00035639915780696013, 0.00025254797270647786, -0.00070480811654608476,
    0.00056917349495288030, 0.00038107279953951152, -0.00004832767552106640, -0.00001160047141568629
};

static const float polyphase8xfilter4[64] =
{
    -0.00001160047141568629, -0.00004832767552106640, 0.00038107279953951152, 0.00056917349495288030,
    -0.00070480811654608476, 0.00025254797270647786, 0.00035639915780696013, -0.00083221312461457934,
    0.00103478202101450390, -0.00093289737900114494, 0.00056676217416478942, -0.00001987860045449669,
    -0.00060183978921851905, 0.00118633130281230620, -0.00162869371539639560, 0.00184225028480482200,
    -0.00176750548092106620, 0.00137873756292502610, -0.00068788009050978823, -0.00025460237565429417,
    0.00136193071021434250, -0.00251606896499480740, 0.00357492186862908880, -0.00438038015325593890,
    0.00476518663684880220, -0.00455482533831069740, 0.00355693608056953160, -0.00151921945966949170,
    -0.00200465909461610940, 0.00813695469487032970, -0.02118099310325147900, 0.08545741447609654700,
    0.07071205074526545900, -0.02388014901113717800, 0.01183152935877091100, -0.00552647396504433920,
    0.00148596888578445710, 0.00121555318378888090, -0.00292792329706461190, 0.00383521477176884070,
    -0.00407634981728579190, 0.00378546088727036750, -0.00310325220590675600, 0.00217587137805511180,
    -0.00114795228773180200, 0.00015338377641278773, 0.00069413459212034961, -0.00130928010530110450,
    0.00164106173258977890, -0.00167597953009304670, 0.00143817315098353240, -0.00098659778327271500,
    0.00040938095642024030, 0.00018452036040974686, -0.00067661601603670300, 0.00095408718371994903,
    -0.00093123666109248718, 0.00058158960488687678, 0.00001271442624327296, -0.00059001569340461763,
    0.00065022085359971421, 0.00029260597635645714, -0.00005386828034834762, -0.00000707856005788827
};

static const float polyphase8xfilter5[64] =
{
    -0.00001752354118849114, -0.00003446958495946170, 0.00047008811115116311, 0.00045149583778780219,
    -0.00076894092872930046, 0.00047618788828089582, 0.00009197902579005798, -0.00064697362622190139,
    0.00100851752563541820, -0.00109375535771801020, 0.00089451555246586410, -0.00045621477010817133,
    -0.00013986723178465493, 0.00079175249424842480, -0.00139028769170792730, 0.00183144503810963020,
    -0.00202740273577447240, 0.00191604280034704920, -0.00146835891655566310, 0.00069324910330136524,
    0.00036069553428501186, -0.00160651004059010120, 0.00292295670535206620, -0.00416047157856245470,
    0.00514726404042791740, -0.00569207291393189830, 0.00557638253007973400, -0.00451741705979632340,
    0.00204115962512772160, 0.00300145106345972670, -0.01516237260357852300, 0.09755145975376086600,
    0.05434387608734767000, -0.02350820651622280100, 0.01381682701592122300, -0.00815969289584616080,
    0.00414981937115260680, -0.00115930523989435050, -0.00102220275357614380, 0.00248673899386953140,
    -0.00330168644834485760, 0.00354352391311250860, -0.00330848547687067520, 0.00271226205959396690,
    -0.00188376295526464770, 0.00095581822472680811, -0.00005464802619675554, -0.00071043859554279257,
    0.00125599179151124560, -0.00153121928051926810, 0.00152215495203344830, -0.00125298740630276000,
    0.00078450145550032285, -0.00020962086499279835, -0.00035428397707766530, 0.00077740863525184738,
    -0.00093664861877437256, 0.00074782136190461302, -0.00022221919156823275, -0.00043577094717956485,
    0.00069349157210732459, 0.00020974856947313336, -0.00005318884282459066, -0.00000388695825192831
};

static const float polyphase8xfilter6[64] =
{
    -0.00002472111034662811, -0.00001032580058558277, 0.00055349978393511490, 0.00030143757577000946,
    -0.00077461306480611420, 0.00066293773127700027, -0.00018732750967400416, -0.00039197851653214976,
    0.00087531359141694673, -0.00113937588749536430, 0.00113011259388196740, -0.00085055664694030269,
    0.00034829156389369192, 0.00029670356071559174, -0.00098402602826872031, 0.00160476729027759600,
    -0.00205327074145764380, 0.00223830512847985300, -0.00209304708081012490, 0.00158317444638181860,
    -0.00071241642862137494, -0.00047490028164087116, 0.00189442023648568910, -0.00342501037656084200,
    0.00491363118474887150, -0.00617923068231396570, 0.00701012531795208030, -0.00714024403544936260,
    0.00615580584287418050, -0.00311138376108936220, -0.00584676809795404460, 0.10613402761530316000,
    0.03746277721486017400, -0.02054581611003392100, 0.01403037645222736200, -0.00966417201297737480,
    0.00618233379133788660, -0.00328877553092718580, 0.00092760292983492458, 0.00089338246990385678,
    -0.00216307979887024520, 0.00289147118196876480, -0.00311871005955389980, 0.00291584964226445760,
    -0.00238015291908662140, 0.00162691111436317150, -0.00077915084816193418, -0.00004348492312229754,
    0.00073554851430502580, -0.00121536388445158540, 0.00143269136290550630, -0.00137406123500243830,
    0.00106555149724676140, -0.00057289470445241102, -0.00000143693453915044, 0.00052522292321950030,
    -0.00085050570414972073, 0.00084162909240903571, -0.00043265212159530387, -0.00025547582772182008,
    0.00070065829999757436, 0.00013626225291498421, -0.00004831055001691508, -0.00000182984414357783
};

static const float polyphase8xfilter7[64] =
{
    -0.00003279706766106838, 0.00002569631187128541, 0.00062442144337030445, 0.00012637883571250220,
    -0.00071843423633693635, 0.00079448278839194797, -0.00045476803926232661, -0.00009140135720076620,
    0.00064655676276784320, -0.00106196556316837370, 0.00124514852293933480, -0.00115732657922609960,
    0.00080732204284516727, -0.00024385296497564968, -0.00045337424289750575, 0.00118365975641694390,
    -0.00183622557500364680, 0.00230162693883983070, -0.00248303216510402730, 0.00230658176248078620,
    -0.00172997883882847830, 0.00074853435796320012, 0.00060189327685675881, -0.00224523360091431110,
    0.00406863057583288240, -0.00592783312693009010, 0.00765208391982513670, -0.00904206471599113170,
    0.00984001261556627540, -0.00957449662395692440, 0.00645758092866270850, 0.11058686170616823000,
    0.02116771695102326000, -0.01564820095328426800, 0.01260511444709402700, -0.00994256584967902380,
    0.00738159244609435180, -0.00493523685432864770, 0.00269106189201281710, -0.00074884899307910678,
    -0.00080635046129673850, 0.00191845357598483200, -0.00256754278849754100, 0.00277185772937339320,
    -0.00258537407535367780, 0.00209182396262220810, -0.00139595818775825930, 0.00061301956994263369,
    0.00014256121428106461, -0.00076830302548444487, 0.00118386211690219210, -0.00133970651621087830,
    0.00122392123868685560, -0.00086697375710736936, 0.00034416140819020045, 0.00022508973131157744,
    -0.00068374678593982705, 0.00085677626062353681, -0.00060210193780960727, -0.00006343701412881912,
    0.00067572747586549999, 0.00007454552789152059, -0.00004102501562491190, -0.00000065139115950063
};

static const double polyphase8xfilter0D[64] =
{
    -0.00000065139115950063, -0.00004102501562491190, 0.00007454552789152059, 0.00067572747586549999,
    -0.00006343701412881912, -0.00060210193780960727, 0.00085677626062353681, -0.00068374678593982705,
    0.00022508973131157744, 0.00034416140819020045, -0.00086697375710736936, 0.00122392123868685560,
    -0.00133970651621087830, 0.00118386211690219210, -0.00076830302548444487, 0.00014256121428106461,
    0.00061301956994263369, -0.00139595818775825930, 0.00209182396262220810, -0.00258537407535367780,
    0.00277185772937339320, -0.00256754278849754100, 0.00191845357598483200, -0.00080635046129673850,
    -0.00074884899307910678, 0.00269106189201281710, -0.00493523685432864770, 0.00738159244609435180,
    -0.00994256584967902380, 0.01260511444709402700, -0.01564820095328426800, 0.02116771695102326000,
    0.11058686170616823000, 0.00645758092866270850, -0.00957449662395692440, 0.00984001261556627540,
    -0.00904206471599113170, 0.00765208391982513670, -0.00592783312693009010, 0.00406863057583288240,
    -0.00224523360091431110, 0.00060189327685675881, 0.00074853435796320012, -0.00172997883882847830,
    0.00230658176248078620, -0.00248303216510402730, 0.00230162693883983070, -0.00183622557500364680,
    0.00118365975641694390, -0.00045337424289750575, -0.00024385296497564968, 0.00080732204284516727,
    -0.00115732657922609960, 0.00124514852293933480, -0.00106196556316837370, 0.00064655676276784320,
    -0.00009140135720076620, -0.00045476803926232661, 0.00079448278839194797, -0.00071843423633693635,
    0.00012637883571250220, 0.00062442144337030445, 0.00002569631187128541, -0.00003279706766106838
};

static const double polyphase8xfilter1D[64] =
{
    -0.00000182984414357783, -0.00004831055001691508, 0.00013626225291498421, 0.00070065829999757436,
    -0.00025547582772182008, -0.00043265212159530387, 0.00084162909240903571, -0.00085050570414972073,
    0.00052522292321950030, -0.00000143693453915044, -0.00057289470445241102, 0.00106555149724676140,
    -0.00137406123500243830, 0.00143269136290550630, -0.00121536388445158540, 0.00073554851430502580,
    -0.00004348492312229754, -0.00077915084816193418, 0.00162691111436317150, -0.00238015291908662140,
    0.00291584964226445760, -0.00311871005955389980, 0.00289147118196876480, -0.00216307979887024520,
    0.00089338246990385678, 0.00092760292983492458, -0.00328877553092718580, 0.00618233379133788660,
    -0.00966417201297737480, 0.01403037645222736200, -0.02054581611003392100, 0.03746277721486017400,
    0.10613402761530316000, -0.00584676809795404460, -0.00311138376108936220, 0.00615580584287418050,
    -0.00714024403544936260, 0.00701012531795208030, -0.00617923068231396570, 0.00491363118474887150,
    -0.00342501037656084200, 0.00189442023648568910, -0.00047490028164087116, -0.00071241642862137494,
    0.00158317444638181860, -0.00209304708081012490, 0.00223830512847985300, -0.00205327074145764380,
    0.00160476729027759600, -0.00098402602826872031, 0.00029670356071559174, 0.00034829156389369192,
    -0.00085055664694030269, 0.00113011259388196740, -0.00113937588749536430, 0.00087531359141694673,
    -0.00039197851653214976, -0.00018732750967400416, 0.00066293773127700027, -0.00077461306480611420,
    0.00030143757577000946, 0.00055349978393511490, -0.00001032580058558277, -0.00002472111034662811
};

static const double polyphase8xfilter2D[64] =
{
    -0.00000388695825192831, -0.00005318884282459066, 0.00020974856947313336, 0.00069349157210732459,
    -0.00043577094717956485, -0.00022221919156823275, 0.00074782136190461302, -0.00093664861877437256,
    0.00077740863525184738, -0.00035428397707766530, -0.00020962086499279835, 0.00078450145550032285,
    -0.00125298740630276000, 0.00152215495203344830, -0.00153121928051926810, 0.00125599179151124560,
    -0.00071043859554279257, -0.00005464802619675554, 0.00095581822472680811, -0.00188376295526464770,
    0.00271226205959396690, -0.00330848547687067520, 0.00354352391311250860, -0.00330168644834485760,
    0.00248673899386953140, -0.00102220275357614380, -0.00115930523989435050, 0.00414981937115260680,
    -0.00815969289584616080, 0.01381682701592122300, -0.02350820651622280100,0.05434387608734767000,
    0.09755145975376086600, -0.01516237260357852300, 0.00300145106345972670, 0.00204115962512772160,
    -0.00451741705979632340, 0.00557638253007973400, -0.00569207291393189830, 0.00514726404042791740,
    -0.00416047157856245470, 0.00292295670535206620, -0.00160651004059010120, 0.00036069553428501186,
    0.00069324910330136524, -0.00146835891655566310, 0.00191604280034704920, -0.00202740273577447240,
    0.00183144503810963020, -0.00139028769170792730, 0.00079175249424842480, -0.00013986723178465493,
    -0.00045621477010817133, 0.00089451555246586410, -0.00109375535771801020, 0.00100851752563541820,
    -0.00064697362622190139, 0.00009197902579005798, 0.00047618788828089582, -0.00076894092872930046,
    0.00045149583778780219, 0.00047008811115116311, -0.00003446958495946170, -0.00001752354118849114
};

static const double polyphase8xfilter3D[64] =
{
    -0.00000707856005788827, -0.00005386828034834762, 0.00029260597635645714, 0.00065022085359971421,
    -0.00059001569340461763, 0.00001271442624327296, 0.00058158960488687678, -0.00093123666109248718,
    0.00095408718371994903, -0.00067661601603670300, 0.00018452036040974686, 0.00040938095642024030,
    -0.00098659778327271500, 0.00143817315098353240, -0.00167597953009304670, 0.00164106173258977890,
    -0.00130928010530110450, 0.00069413459212034961, 0.00015338377641278773, -0.00114795228773180200,
    0.00217587137805511180, -0.00310325220590675600, 0.00378546088727036750, -0.00407634981728579190,
    0.00383521477176884070, -0.00292792329706461190, 0.00121555318378888090, 0.00148596888578445710,
    -0.00552647396504433920, 0.01183152935877091100, -0.02388014901113717800, 0.07071205074526545900,
    0.08545741447609654700, -0.02118099310325147900, 0.00813695469487032970, -0.00200465909461610940,
    -0.00151921945966949170, 0.00355693608056953160, -0.00455482533831069740, 0.00476518663684880220,
    -0.00438038015325593890, 0.00357492186862908880 -0.00251606896499480740, 0.00136193071021434250,
    -0.00025460237565429417, -0.00068788009050978823, 0.00137873756292502610, -0.00176750548092106620,
    0.00184225028480482200, -0.00162869371539639560, 0.00118633130281230620, -0.00060183978921851905,
    -0.00001987860045449669, 0.00056676217416478942, -0.00093289737900114494, 0.00103478202101450390,
    -0.00083221312461457934, 0.00035639915780696013, 0.00025254797270647786, -0.00070480811654608476,
    0.00056917349495288030, 0.00038107279953951152, -0.00004832767552106640, -0.00001160047141568629
};

static const double polyphase8xfilter4D[64] =
{
    -0.00001160047141568629, -0.00004832767552106640, 0.00038107279953951152, 0.00056917349495288030,
    -0.00070480811654608476, 0.00025254797270647786, 0.00035639915780696013, -0.00083221312461457934,
    0.00103478202101450390, -0.00093289737900114494, 0.00056676217416478942, -0.00001987860045449669,
    -0.00060183978921851905, 0.00118633130281230620, -0.00162869371539639560, 0.00184225028480482200,
    -0.00176750548092106620, 0.00137873756292502610, -0.00068788009050978823, -0.00025460237565429417,
    0.00136193071021434250, -0.00251606896499480740, 0.00357492186862908880, -0.00438038015325593890,
    0.00476518663684880220, -0.00455482533831069740, 0.00355693608056953160, -0.00151921945966949170,
    -0.00200465909461610940, 0.00813695469487032970, -0.02118099310325147900, 0.08545741447609654700,
    0.07071205074526545900, -0.02388014901113717800, 0.01183152935877091100, -0.00552647396504433920,
    0.00148596888578445710, 0.00121555318378888090, -0.00292792329706461190, 0.00383521477176884070,
    -0.00407634981728579190, 0.00378546088727036750, -0.00310325220590675600, 0.00217587137805511180,
    -0.00114795228773180200, 0.00015338377641278773, 0.00069413459212034961, -0.00130928010530110450,
    0.00164106173258977890, -0.00167597953009304670, 0.00143817315098353240, -0.00098659778327271500,
    0.00040938095642024030, 0.00018452036040974686, -0.00067661601603670300, 0.00095408718371994903,
    -0.00093123666109248718, 0.00058158960488687678, 0.00001271442624327296, -0.00059001569340461763,
    0.00065022085359971421, 0.00029260597635645714, -0.00005386828034834762, -0.00000707856005788827
};

static const double polyphase8xfilter5D[64] =
{
    -0.00001752354118849114, -0.00003446958495946170, 0.00047008811115116311, 0.00045149583778780219,
    -0.00076894092872930046, 0.00047618788828089582, 0.00009197902579005798, -0.00064697362622190139,
    0.00100851752563541820, -0.00109375535771801020, 0.00089451555246586410, -0.00045621477010817133,
    -0.00013986723178465493, 0.00079175249424842480, -0.00139028769170792730, 0.00183144503810963020,
    -0.00202740273577447240, 0.00191604280034704920, -0.00146835891655566310, 0.00069324910330136524,
    0.00036069553428501186, -0.00160651004059010120, 0.00292295670535206620, -0.00416047157856245470,
    0.00514726404042791740, -0.00569207291393189830, 0.00557638253007973400, -0.00451741705979632340,
    0.00204115962512772160, 0.00300145106345972670, -0.01516237260357852300, 0.09755145975376086600,
    0.05434387608734767000, -0.02350820651622280100, 0.01381682701592122300, -0.00815969289584616080,
    0.00414981937115260680, -0.00115930523989435050, -0.00102220275357614380, 0.00248673899386953140,
    -0.00330168644834485760, 0.00354352391311250860, -0.00330848547687067520, 0.00271226205959396690,
    -0.00188376295526464770, 0.00095581822472680811, -0.00005464802619675554, -0.00071043859554279257,
    0.00125599179151124560, -0.00153121928051926810, 0.00152215495203344830, -0.00125298740630276000,
    0.00078450145550032285, -0.00020962086499279835, -0.00035428397707766530, 0.00077740863525184738,
    -0.00093664861877437256, 0.00074782136190461302, -0.00022221919156823275, -0.00043577094717956485,
    0.00069349157210732459, 0.00020974856947313336, -0.00005318884282459066, -0.00000388695825192831
};

static const double polyphase8xfilter6D[64] =
{
    -0.00002472111034662811, -0.00001032580058558277, 0.00055349978393511490, 0.00030143757577000946,
    -0.00077461306480611420, 0.00066293773127700027, -0.00018732750967400416, -0.00039197851653214976,
    0.00087531359141694673, -0.00113937588749536430, 0.00113011259388196740, -0.00085055664694030269,
    0.00034829156389369192, 0.00029670356071559174, -0.00098402602826872031, 0.00160476729027759600,
    -0.00205327074145764380, 0.00223830512847985300, -0.00209304708081012490, 0.00158317444638181860,
    -0.00071241642862137494, -0.00047490028164087116, 0.00189442023648568910, -0.00342501037656084200,
    0.00491363118474887150, -0.00617923068231396570, 0.00701012531795208030, -0.00714024403544936260,
    0.00615580584287418050, -0.00311138376108936220, -0.00584676809795404460, 0.10613402761530316000,
    0.03746277721486017400, -0.02054581611003392100, 0.01403037645222736200, -0.00966417201297737480,
    0.00618233379133788660, -0.00328877553092718580, 0.00092760292983492458, 0.00089338246990385678,
    -0.00216307979887024520, 0.00289147118196876480, -0.00311871005955389980, 0.00291584964226445760,
    -0.00238015291908662140, 0.00162691111436317150, -0.00077915084816193418, -0.00004348492312229754,
    0.00073554851430502580, -0.00121536388445158540, 0.00143269136290550630, -0.00137406123500243830,
    0.00106555149724676140, -0.00057289470445241102, -0.00000143693453915044, 0.00052522292321950030,
    -0.00085050570414972073, 0.00084162909240903571, -0.00043265212159530387, -0.00025547582772182008,
    0.00070065829999757436, 0.00013626225291498421, -0.00004831055001691508, -0.00000182984414357783
};

static const double polyphase8xfilter7D[64] =
{
    -0.00003279706766106838, 0.00002569631187128541, 0.00062442144337030445, 0.00012637883571250220,
    -0.00071843423633693635, 0.00079448278839194797, -0.00045476803926232661, -0.00009140135720076620,
    0.00064655676276784320, -0.00106196556316837370, 0.00124514852293933480, -0.00115732657922609960,
    0.00080732204284516727, -0.00024385296497564968, -0.00045337424289750575, 0.00118365975641694390,
    -0.00183622557500364680, 0.00230162693883983070, -0.00248303216510402730, 0.00230658176248078620,
    -0.00172997883882847830, 0.00074853435796320012, 0.00060189327685675881, -0.00224523360091431110,
    0.00406863057583288240, -0.00592783312693009010, 0.00765208391982513670, -0.00904206471599113170,
    0.00984001261556627540, -0.00957449662395692440, 0.00645758092866270850, 0.11058686170616823000,
    0.02116771695102326000, -0.01564820095328426800, 0.01260511444709402700, -0.00994256584967902380,
    0.00738159244609435180, -0.00493523685432864770, 0.00269106189201281710, -0.00074884899307910678,
    -0.00080635046129673850, 0.00191845357598483200, -0.00256754278849754100, 0.00277185772937339320,
    -0.00258537407535367780, 0.00209182396262220810, -0.00139595818775825930, 0.00061301956994263369,
    0.00014256121428106461, -0.00076830302548444487, 0.00118386211690219210, -0.00133970651621087830,
    0.00122392123868685560, -0.00086697375710736936, 0.00034416140819020045, 0.00022508973131157744,
    -0.00068374678593982705, 0.00085677626062353681, -0.00060210193780960727, -0.00006343701412881912,
    0.00067572747586549999, 0.00007454552789152059, -0.00004102501562491190, -0.00000065139115950063
};


static const float* polyphase8x[8] =
{
    polyphase8xfilter0,
    polyphase8xfilter1,
    polyphase8xfilter2,
    polyphase8xfilter3,
    polyphase8xfilter4,
    polyphase8xfilter5,
    polyphase8xfilter6,
    polyphase8xfilter7
};

static const double* polyphase8xD[8] =
{
    polyphase8xfilter0D,
    polyphase8xfilter1D,
    polyphase8xfilter2D,
    polyphase8xfilter3D,
    polyphase8xfilter4D,
    polyphase8xfilter5D,
    polyphase8xfilter6D,
    polyphase8xfilter7D
};

/******************************************************************************
 * 16x Polyphase Coefficients
 *****************************************************************************/

/*
 float* polyphase16x[16] =
 {
 polyphase16xfilter0,
 polyphase16xfilter1,
 polyphase16xfilter2,
 polyphase16xfilter3,
 polyphase16xfilter4,
 polyphase16xfilter5,
 polyphase16xfilter6,
 polyphase16xfilter7,
 polyphase16xfilter8,
 polyphase16xfilter9,
 polyphase16xfilter10,
 polyphase16xfilter11,
 polyphase16xfilter12,
 polyphase16xfilter13,
 polyphase16xfilter14,
 polyphase16xfilter15
 };
 
 */


/******************************************************************************
 * Polyphase Coefficient matrix
 *****************************************************************************/

const float** PolyphaseCoeffs[N_FACTORS] =
{
    polyphase2x,
    polyphase4x,
    polyphase8x
    /*polyphase16x*/
};

const double** PolyphaseCoeffsD[N_FACTORS] =
{
    polyphase2xD,
    polyphase4xD,
    polyphase8xD,
};


static const float ref[] = {1.0, (float)K12_REF, (float)K14_REF, (float)K20_REF};
static const double refD[] = {1.0, K12_REF, K14_REF, K20_REF};

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
    const unsigned end = 4 * (length / 4);
    #pragma omp simd aligned(dest,src)
    for (i = 0; i < end; i+=4)
    {
        dest[i] = int16ToFloat(*src++);
        dest[i + 1] = int16ToFloat(*src++);
        dest[i + 2] = int16ToFloat(*src++);
        dest[i + 3] = int16ToFloat(*src++);
    }
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
    #pragma omp simd aligned(dest,src)
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
Error_t
FillBuffer(float *dest, unsigned length, float value)
{
#ifdef __APPLE__
    // Use the Accelerate framework if we have it
    vDSP_vfill(&value, dest, 1, length);

#else
    // Otherwise do it manually
    unsigned i = 0;
    const unsigned end = 4 * (length / 4);
    #pragma omp simd aligned(dest,src)
    for (i = 0; i < end; i+=4)
    {
        dest[i] = value;
        dest[i + 1] = value;
        dest[i + 2] = value;
        dest[i + 3] = value;
    }
    for (i = end; i < length; ++i)
    {
        dest[i] = value;
    }
#endif
    return NOERR;
}

/*******************************************************************************
 FillBufferD */
Error_t
FillBufferD(double *dest, unsigned length, double value)
{
#ifdef __APPLE__
    // Use the Accelerate framework if we have it
    vDSP_vfillD(&value, dest, 1, length);

#else
    // Otherwise do it manually
    unsigned i = 0;
    const unsigned end = 4 * (length / 4);
    #pragma omp simd aligned(dest)
    for (i = 0; i < end; i+=4)
    {
        dest[i] = value;
        dest[i + 1] = value;
        dest[i + 2] = value;
        dest[i + 3] = value;
    }
    for (i = end; i < length; ++i)
    {
        dest[i] = value;
    }
#endif
    return NOERR;
}

/*******************************************************************************
 ClearBuffer */
Error_t
ClearBuffer(float *dest, unsigned length)
{
#ifdef __APPLE__
    // Use the Accelerate framework if we have it
    vDSP_vclr(dest, 1, length);

#else
    // Otherwise do it manually. Yes this works for floats
    memset(dest, 0, length * sizeof(float));
#endif
    return NOERR;
}

/*******************************************************************************
 ClearBufferD */
Error_t
ClearBufferD(double *dest, unsigned length)
{
#ifdef __APPLE__
    // Use the Accelerate framework if we have it
    vDSP_vclrD(dest, 1, length);

#else
    // Otherwise do it manually. Yes this works for doubles
    memset(dest, 0, length * sizeof(double));
#endif
    return NOERR;
}


/*******************************************************************************
 CopyBuffer */
Error_t
CopyBuffer(float* dest, const float* src, unsigned length)
{
    if (src != dest)
    {
#if defined(__APPLE__) || defined(USE_BLAS)
        // Use the Accelerate framework if we have it
        cblas_scopy(length, src, 1, dest, 1);
#else
        // Do it the boring way. Prefer memmove to memcpy in case src and dest
        // overlap
        memmove(dest, src, length * sizeof(float));
#endif
    }
    return NOERR;
}


/*******************************************************************************
 CopyBufferD */
Error_t
CopyBufferD(double* dest, const double* src, unsigned length)
{
    if (src != dest)
    {
#if defined(__APPLE__) || defined(USE_BLAS)
        // Use the Accelerate framework if we have it
        cblas_dcopy(length, src, 1, dest, 1);
#else
        // Do it the boring way. Prefer memmove to memcpy in case src and dest
        // overlap
        memmove(dest, src, length * sizeof(double));
#endif
    }
    return NOERR;
}


/*******************************************************************************
 CopyBufferStride */
Error_t
CopyBufferStride(float*         dest,
                 unsigned       dest_stride,
                 const float*   src,
                 unsigned       src_stride,
                 unsigned       length)
{
#if defined(__APPLE__) || defined(USE_BLAS)
    // Use the Accelerate framework if we have it
    cblas_scopy(length, src, src_stride, dest, dest_stride);
#else
	#pragma omp simd aligned(dest,src)
    for (unsigned i = 0; i < length; ++i)
    {
        dest[i * dest_stride] = src[i * src_stride];
    }
#endif
    return NOERR;
}


/*******************************************************************************
 CopyBufferStrideD */
Error_t
CopyBufferStrideD(double*       dest,
                  unsigned      dest_stride,
                  const double* src,
                  unsigned      src_stride,
                  unsigned      length)
{
#if defined(__APPLE__) || defined(USE_BLAS)
    // Use the Accelerate framework if we have it
    cblas_dcopy(length, src, src_stride, dest, dest_stride);
#else
    for (unsigned i = 0; i < length; ++i)
    {
        dest[i * dest_stride] = src[i * src_stride];
    }
#endif
    return NOERR;
}


/*******************************************************************************
 SplitToInterleaved */
Error_t
SplitToInterleaved(float* dest, const float* real, const float* imag, unsigned length)
{
#if defined(__APPLE__)
    DSPSplitComplex in = {.realp = (float*)real, .imagp = (float*)imag};
    vDSP_ztoc(&in, 1, (DSPComplex*)dest, 2, length);

#elif defined(USE_BLAS)
    cblas_scopy(length, real, 1, dest, 2);
    cblas_scopy(length, imag, 1, dest + 1, 2);
#else
    unsigned i;
    unsigned i2;
    const unsigned end = 4 * (length / 4);
    #pragma omp simd aligned(dest,real,imag)
    for (i = 0; i < end; i+=4)
    {
        i2 = i * 2;
        dest[i2] = real[i];
        dest[i2 + 2] = real[i + 1];
        dest[i2 + 4] = real[i + 2];
        dest[i2 + 6] = real[i + 3];

        dest[i2 + 1] = imag[i];
        dest[i2 + 3] = imag[i + 1];
        dest[i2 + 5] = imag[i + 2];
        dest[i2 + 7] = imag[i + 3];
    }
    for (i = end; i < length; ++i)
    {
        i2 = i * 2;
        dest[i2] = real[i];
        dest[i2 + 1] = imag[i];
    }
#endif
    return NOERR;
}


Error_t
SplitToInterleavedD(double* dest, const double* real, const double* imag, unsigned length)
{
#if defined(__APPLE__)
    DSPDoubleSplitComplex in = {.realp = (double*)real, .imagp = (double*)imag};
    vDSP_ztocD(&in, 1, (DSPDoubleComplex*)dest, 2, length);

#elif defined(USE_BLAS)
    cblas_dcopy(length, real, 1, dest, 2);
    cblas_dcopy(length, imag, 1, dest + 1, 2);
#else
    unsigned i;
    unsigned i2;
    const unsigned end = 4 * (length / 4);
    #pragma omp simd aligned(dest,real,imag)
    for (i = 0; i < end; i+=4)
    {
        i2 = i * 2;
        dest[i2] = real[i];
        dest[i2 + 2] = real[i + 1];
        dest[i2 + 4] = real[i + 2];
        dest[i2 + 6] = real[i + 3];

        dest[i2 + 1] = imag[i];
        dest[i2 + 3] = imag[i + 1];
        dest[i2 + 5] = imag[i + 2];
        dest[i2 + 7] = imag[i + 3];
    }
    for (i = end; i < length; ++i)
    {
        i2 = i * 2;
        dest[i2] = real[i];
        dest[i2 + 1] = imag[i];
    }
#endif
    return NOERR;
}


/*******************************************************************************
 InterleavedToSplit */
Error_t
InterleavedToSplit(float*       real,
                   float*       imag,
                   const float* input,
                   unsigned     length)
{
#if defined(__APPLE__)
    DSPSplitComplex out = {.realp = real, .imagp = imag};
    vDSP_ctoz((const DSPComplex*)input, 2, &out, 1, length);

#elif defined(USE_BLAS)
    cblas_scopy(length, input, 2, real, 1);
    cblas_scopy(length, input + 1, 2, imag, 1);
#else
    unsigned i;
    unsigned i2;
    const unsigned end = 4 * (length / 4);
    #pragma omp simd aligned(real,imag,input)
    for (i = 0; i < end; i+=4)
    {
        i2 = i * 2;
        real[i] = input[i2];
        real[i + 1] = input[i2 + 2];
        real[i + 2] = input[i2 + 4];
        real[i + 3] = input[i2 + 6];

        imag[i] = input[i2 + 1];
        imag[i + 1] = input[i2 + 3];
        imag[i + 2] = input[i2 + 5];
        imag[i + 3] = input[i2 + 7];
    }
    for (i = end; i < length; ++i)
    {
        i2 = i * 2;
        real[i] = input[i2];
        imag[i] = input[i2 + 1];
    }

#endif
    return NOERR;
}


Error_t
InterleavedToSplitD(double*         real,
                    double*         imag,
                    const double*   input,
                    unsigned        length)
{
#if defined(__APPLE__)
    DSPDoubleSplitComplex out = {.realp = real, .imagp = imag};
    vDSP_ctozD((const DSPDoubleComplex*)input, 2, &out, 1, length);

#elif defined(USE_BLAS)
    cblas_dcopy(length, input, 2, real, 1);
    cblas_dcopy(length, input + 1, 2, imag, 1);
#else

    unsigned i;
    unsigned i2;
    const unsigned end = 4 * (length / 4);
    #pragma omp simd aligned(real,imag,input)
    for (i = 0; i < end; i+=4)
    {
        i2 = i * 2;
        real[i] = input[i2];
        real[i + 1] = input[i2 + 2];
        real[i + 2] = input[i2 + 4];
        real[i + 3] = input[i2 + 6];

        imag[i] = input[i2 + 1];
        imag[i + 1] = input[i2 + 3];
        imag[i + 2] = input[i2 + 5];
        imag[i + 3] = input[i2 + 7];
    }
    for (i = end; i < length; ++i)
    {
        i2 = i * 2;
        real[i] = input[i2];
        imag[i] = input[i2 + 1];
    }

#endif
    return NOERR;
}


/*******************************************************************************
 VectorAbs */
Error_t
VectorAbs(float *dest, const float *in, unsigned length)
{
#ifdef __APPLE__
    vDSP_vabs(in, 1, dest, 1, length);
#else
    // Otherwise do it manually
    unsigned i = 0;
    const unsigned end = 4 * (length / 4);
    #pragma omp simd aligned(dest,src)
    for (i = 0; i < end; i+=4)
    {
        dest[i] = fabsf(in[i]);
        dest[i + 1] = fabsf(in[i + 1]);
        dest[i + 2] = fabsf(in[i + 2]);
        dest[i + 3] = fabsf(in[i + 3]);
    }
    #pragma omp simd aligned(dest,src)
    for (unsigned i = end; i < length; ++i)
    {
        dest[i] = fabsf(in[i]);
    }
#endif
    return NOERR;
}


/*******************************************************************************
 VectorAbsD */
Error_t
VectorAbsD(double *dest, const double *in, unsigned length)
{
#ifdef __APPLE__
    vDSP_vabsD(in, 1, dest, 1, length);
#else
    // Otherwise do it manually
    unsigned i = 0;
    const unsigned end = 4 * (length / 4);
    #pragma omp simd aligned(dest,in)
    for (i = 0; i < end; i+=4)
    {
        dest[i] = fabs(in[i]);
        dest[i + 1] = fabs(in[i + 1]);
        dest[i + 2] = fabs(in[i + 2]);
        dest[i + 3] = fabs(in[i + 3]);
    }
    #pragma omp simd aligned(dest,in)
    for (unsigned i = end; i < length; ++i)

    {
        dest[i] = fabs(in[i]);
    }
#endif
    return NOERR;
}


/*******************************************************************************
 VectorNegate */
Error_t
VectorNegate(float          *dest,
             const float    *in,
             unsigned       length)
{
#if defined(__APPLE__)
    // Use the Accelerate framework if we have it
    vDSP_vneg(in, 1, dest, 1, length);
#else
    // Otherwise do it manually
    unsigned i;
    const unsigned end = 4 * (length / 4);
    #pragma omp simd aligned(dest,in)
    for (i = 0; i < end; i+=4)
    {
        dest[i] = -in[i];
        dest[i + 1] = -in[i + 1];
        dest[i + 2] = -in[i + 2];
        dest[i + 3] = -in[i + 3];
    }
    #pragma omp simd aligned(dest,in)
    for (i = end; i < length; ++i)
    {
        dest[i] = -in[i];
    }

#endif
    return NOERR;
}


/*******************************************************************************
 VectorNegateD */
Error_t
VectorNegateD(double          *dest,
              const double    *in,
              unsigned       length)
{
#ifdef __APPLE__
    // Use the Accelerate framework if we have it
    vDSP_vnegD(in, 1, dest, 1, length);

#else
    // Otherwise do it manually
    unsigned i;
    const unsigned end = 4 * (length / 4);
    #pragma omp simd aligned(dest,in)
    for (i = 0; i < end; i+=4)
    {
        dest[i] = -in[i];
        dest[i + 1] = -in[i + 1];
        dest[i + 2] = -in[i + 2];
        dest[i + 3] = -in[i + 3];
    }
    #pragma omp simd aligned(dest,in)
    for (i = end; i < length; ++i)
    {
        dest[i] = -in[i];
    }
#endif
    return NOERR;
}


/*******************************************************************************
 VectorSum */
float
VectorSum(const float* src, unsigned length)
{
    float res = 0.0;
#ifdef __APPLE__
    // Use the Accelerate framework if we have it
    vDSP_sve(src, 1, &res, length);
#else
	#pragma omp simd aligned(src)
    for (unsigned i = 0; i < length; ++i)
    {
        res += src[i];
    }
#endif
    return res;
}

/*******************************************************************************
 VectorSumD */
double
VectorSumD(const double* src, unsigned length)
{
    double res = 0.0;
#ifdef __APPLE__
    // Use the Accelerate framework if we have it
    vDSP_sveD(src, 1, &res, length);
#else
	#pragma omp simd aligned(src)
    for (unsigned i = 0; i < length; ++i)
    {
        res += src[i];
    }
#endif
    return res;
}


/*******************************************************************************
 VectorMax */
float
VectorMax(const float* src, unsigned length)
{
    float res = FLT_MIN;
#ifdef __APPLE__
    // Use the Accelerate framework if we have it
    vDSP_maxv(src, 1, &res, length);
#else
    for (unsigned i = 0; i < length; ++i)
    {
        if (src[i] > res)
        {
            res = src[i];
        }
    }
#endif
    return res;
}


/*******************************************************************************
 VectorMax */
double
VectorMaxD(const double* src, unsigned length)
{
    double res = DBL_MIN;
#ifdef __APPLE__
    // Use the Accelerate framework if we have it
    vDSP_maxvD(src, 1, &res, length);
#else
	#pragma omp simd aligned(src)
    for (unsigned i = 0; i < length; ++i)
    {
        if (src[i] > res)
        {
            res = src[i];
        }
    }
#endif
    return res;
}


/*******************************************************************************
 VectorMaxVI */
Error_t
VectorMaxVI(float* value, unsigned* index, const float* src, unsigned length)
{
#ifdef __APPLE__
    // Use the Accelerate framework if we have it
    vDSP_maxvi(src, 1, value, (unsigned long*)index, length);
#else
    float res = FLT_MIN;
    #pragma omp simd aligned(value,index,src)
    for (unsigned i = 0; i < length; ++i)
    {
        if (src[i] > res)
        {
            *value = res = src[i];
            *index = i;
        }
    }
#endif
    return NOERR;
}


/*******************************************************************************
 VectorMaxVID*/
Error_t
VectorMaxVID(double* value, unsigned* index, const double* src, unsigned length)
{
#ifdef __APPLE__
    // Use the Accelerate framework if we have it
    vDSP_maxviD(src, 1, value, (unsigned long*)index, length);
#else
    double res = src[0];
    #pragma omp simd aligned(value,index,src)
    for (unsigned i = 0; i < length; ++i)
    {
        if (src[i] > res)
        {
            *value = res = src[i];
            *index = i;
        }
    }
#endif
    return NOERR;
}


/*******************************************************************************
 VectorMin */
float
VectorMin(const float* src, unsigned length)
{
    float res = FLT_MAX;
#ifdef __APPLE__
    // Use the Accelerate framework if we have it
    vDSP_minv(src, 1, &res, length);
#else
	#pragma omp simd aligned(src)
    for (unsigned i = 0; i < length; ++i)
    {
        if (src[i] < res)
        {
            res = src[i];
        }
    }
#endif
    return res;
}


/*******************************************************************************
 VectorMinD */
double
VectorMinD(const double* src, unsigned length)
{
    double res = DBL_MAX;
#ifdef __APPLE__
    // Use the Accelerate framework if we have it
    vDSP_minvD(src, 1, &res, length);
#else
	#pragma omp simd aligned(src)
    for (unsigned i = 0; i < length; ++i)
    {
        if (src[i] < res)
        {
            res = src[i];
        }
    }
#endif
    return res;
}


/*******************************************************************************
 VectorMinVI */
Error_t
VectorMinVI(float* value, unsigned* index, const float* src, unsigned length)
{

#ifdef __APPLE__
    // Use the Accelerate framework if we have it
    vDSP_minvi(src, 1, value, (unsigned long*)index, length);
#else
    float res = src[0];
    #pragma omp simd aligned(value,index,src)
    for (unsigned i = 0; i < length; ++i)
    {
        if (src[i] < res)
        {
            *value = res = src[i];
            *index = i;
        }
    }
#endif
    return NOERR;
}

/*******************************************************************************
 VectorMinVID */
Error_t
VectorMinVID(double* value, unsigned* index, const double* src, unsigned length)
{

#ifdef __APPLE__
    // Use the Accelerate framework if we have it
    vDSP_minviD(src, 1, value, (unsigned long*)index, length);
#else
    double res = src[0];
    #pragma omp simd aligned(value,index,src)
    for (unsigned i = 0; i < length; ++i)
    {
        if (src[i] < res)
        {
            *value = res = src[i];
            *index = i;
        }
    }
#endif
    return NOERR;
}

/*******************************************************************************
 VectorVectorAdd */
Error_t
VectorVectorAdd(float         *dest,
                const float   *in1,
                const float   *in2,
                unsigned      length)
{
#ifdef __APPLE__
    // Use the Accelerate framework if we have it
    vDSP_vadd(in1, 1, in2, 1, dest, 1, length);

#else
    // Otherwise do it manually
    unsigned i;    
    const unsigned end = 4 * (length / 4);
    #pragma omp simd aligned(dest,in1,in2)
    for (i = 0; i < end; i+=4)
    {
        dest[i] = in1[i] + in2[i];
        dest[i + 1] = in1[i + 1] + in2[i + 1];
        dest[i + 2] = in1[i + 2] + in2[i + 2];
        dest[i + 3] = in1[i + 3] + in2[i + 3];
    }
    for (i = end; i < length; ++i)
    {
        dest[i] = in1[i] + in2[i];
    }

#endif
    return NOERR;
}

/*******************************************************************************
 VectorVectorAddD */
Error_t
VectorVectorAddD(double           *dest,
                 const double     *in1,
                 const double     *in2,
                 unsigned         length)
{
#ifdef __APPLE__
    // Use the Accelerate framework if we have it
    vDSP_vaddD(in1, 1, in2, 1, dest, 1, length);

#else
    // Otherwise do it manually
    unsigned i;
    const unsigned end = 4 * (length / 4);
    #pragma omp simd aligned(dest,in1,in2)
    for (i = 0; i < end; i+=4)
    {
        dest[i] = in1[i] + in2[i];
        dest[i + 1] = in1[i + 1] + in2[i + 1];
        dest[i + 2] = in1[i + 2] + in2[i + 2];
        dest[i + 3] = in1[i + 3] + in2[i + 3];
    }
    for (i = end; i < length; ++i)
    {
        dest[i] = in1[i] + in2[i];
    }

#endif
    return NOERR;
}

/*******************************************************************************
 VectorVectorSub */
Error_t
VectorVectorSub(float         *dest,
                const float   *in1,
                const float   *in2,
                unsigned      length)
{
#ifdef __APPLE__
    // Use the Accelerate framework if we have it
    vDSP_vsub(in1, 1, in2, 1, dest, 1, length);

#else
    // Otherwise do it manually
    unsigned i;
    const unsigned end = 4 * (length / 4);
    #pragma omp simd aligned(dest,in1,in2)
    for (i = 0; i < end; i+=4)
    {
        dest[i] = in2[i] - in1[i];
        dest[i + 1] = in2[i + 1] - in1[i + 1];
        dest[i + 2] = in2[i + 2] - in1[i + 2];
        dest[i + 3] = in2[i + 3] - in1[i + 3];
    }
    for (i = end; i < length; ++i)
    {
        dest[i] = in2[i] - in1[i];
    }

#endif
    return NOERR;
}

/*******************************************************************************
 VectorVectorSubD */
Error_t
VectorVectorSubD(double           *dest,
                 const double     *in1,
                 const double     *in2,
                 unsigned         length)
{
#ifdef __APPLE__
    // Use the Accelerate framework if we have it
    vDSP_vsubD(in1, 1, in2, 1, dest, 1, length);

#else
    // Otherwise do it manually
    unsigned i;    
    const unsigned end = 4 * (length / 4);
    #pragma omp simd aligned(dest,in1,in2)
    for (i = 0; i < end; i+=4)
    {
        dest[i] = in2[i] - in1[i];
        dest[i + 1] = in2[i + 1] - in1[i + 1];
        dest[i + 2] = in2[i + 2] - in1[i + 2];
        dest[i + 3] = in2[i + 3] - in1[i + 3];
    }
    for (i = end; i < length; ++i)
    {
        dest[i] = in2[i] - in1[i];
    }

#endif
    return NOERR;
}


/*******************************************************************************
 VectorScalarAdd */
Error_t
VectorScalarAdd(float           *dest,
                const float     *in1,
                const float     scalar,
                unsigned        length)
{
#ifdef __APPLE__
    // Use the Accelerate framework if we have it
    vDSP_vsadd(in1, 1, &scalar, dest, 1, length);
#else
    // Otherwise do it manually
    unsigned i;
    const unsigned end = 4 * (length / 4);
    #pragma omp simd aligned(dest,in1)
    for (i = 0; i < end; i+=4)
    {
        dest[i] = in1[i] + scalar;
        dest[i + 1] = in1[i + 1] + scalar;
        dest[i + 2] = in1[i + 2] + scalar;
        dest[i + 3] = in1[i + 3] + scalar;
    }
    #pragma omp simd aligned(dest,in1)
    for (i = end; i < length; ++i)
    {
        dest[i] = in1[i] + scalar;
    }

#endif
    return NOERR;
}


/*******************************************************************************
 VectorScalarAddD */
Error_t
VectorScalarAddD(double         *dest,
                 const double   *in1,
                 const double   scalar,
                 unsigned       length)
{
#ifdef __APPLE__
    // Use the Accelerate framework if we have it
    vDSP_vsaddD(in1, 1, &scalar, dest, 1, length);
#else
    // Otherwise do it manually
    unsigned i;
    const unsigned end = 4 * (length / 4);
    #pragma omp simd aligned(dest,in1)
    for (i = 0; i < end; i+=4)
    {
        dest[i] = in1[i] + scalar;
        dest[i + 1] = in1[i + 1] + scalar;
        dest[i + 2] = in1[i + 2] + scalar;
        dest[i + 3] = in1[i + 3] + scalar;
    }
    #pragma omp simd aligned(dest,in1)
    for (i = end; i < length; ++i)    
    {
        dest[i] = in1[i] + scalar;
    }

#endif
    return NOERR;
}




/*******************************************************************************
 VectorVectorMultiply */
Error_t
VectorVectorMultiply(float          *dest,
                     const float    *in1,
                     const float    *in2,
                     unsigned       length)
{
#ifdef __APPLE__
    // Use the Accelerate framework if we have it
    vDSP_vmul(in1, 1, in2, 1, dest, 1, length);

#else
    // Otherwise do it manually
    unsigned i;
    const unsigned end = 4 * (length / 4);
    #pragma omp simd aligned(dest,in1,in2)
    for (i = 0; i < end; i+=4)
    {
        dest[i] = in1[i] * in2[i];
        dest[i + 1] = in1[i + 1] * in2[i + 1];
        dest[i + 2] = in1[i + 2] * in2[i + 2];
        dest[i + 3] = in1[i + 3] * in2[i + 3];
    }
    #pragma omp simd aligned(dest,in1,in2)
    for (i = end; i < length; ++i)
    {
        dest[i] = in1[i] * in2[i];
    }

#endif
    return NOERR;
}


/*******************************************************************************
 VectorVectorMultiplyD */
Error_t
VectorVectorMultiplyD(double        *dest,
                      const double  *in1,
                      const double  *in2,
                      unsigned      length)
{
#ifdef __APPLE__
    // Use the Accelerate framework if we have it
    vDSP_vmulD(in1, 1, in2, 1, dest, 1, length);

#else
    // Otherwise do it manually
    unsigned i;
    const unsigned end = 4 * (length / 4);
    #pragma omp simd aligned(dest,in1,in2)
    for (i = 0; i < end; i+=4)
    {
        dest[i] = in1[i] * in2[i];
        dest[i + 1] = in1[i + 1] * in2[i + 1];
        dest[i + 2] = in1[i + 2] * in2[i + 2];
        dest[i + 3] = in1[i + 3] * in2[i + 3];
    }
    #pragma omp simd aligned(dest,in1,in2)
    for (i = end; i < length; ++i)
    {
        dest[i] = in1[i] * in2[i];
    }

#endif
    return NOERR;
}


/*******************************************************************************
 VectorScalarMultiply */
Error_t
VectorScalarMultiply(float          *dest,
                     const float    *in1,
                     const float    scalar,
                     unsigned       length)
{
#ifdef __APPLE__
    // Use the Accelerate framework if we have it
    vDSP_vsmul(in1, 1, &scalar,dest, 1, length);

#else
    // Otherwise do it manually
    unsigned i;
    const unsigned end = 4 * (length / 4);
    #pragma omp simd aligned(dest,in1)
    for (i = 0; i < end; i+=4)
    {
        dest[i] = in1[i] * scalar;
        dest[i + 1] = in1[i + 1] * scalar;
        dest[i + 2] = in1[i + 2] * scalar;
        dest[i + 3] = in1[i + 3] * scalar;
    }
    #pragma omp simd aligned(dest,in1)
    for (i = end; i < length; ++i)
    {
        dest[i] = in1[i] * scalar;
    }
#endif
    return NOERR;
}

/*******************************************************************************
 VectorScalarMultiplyD */
Error_t
VectorScalarMultiplyD(double        *dest,
                      const double  *in1,
                      const double  scalar,
                      unsigned      length)
{
#ifdef __APPLE__
    // Use the Accelerate framework if we have it
    vDSP_vsmulD(in1, 1, &scalar,dest, 1, length);

#else
    // Otherwise do it manually
    unsigned i;
    const unsigned end = 4 * (length / 4);
    #pragma omp simd aligned(dest,in1)
    for (i = 0; i < end; i+=4)
    {
        dest[i] = in1[i] * scalar;
        dest[i + 1] = in1[i + 1] * scalar;
        dest[i + 2] = in1[i + 2] * scalar;
        dest[i + 3] = in1[i + 3] * scalar;
    }
    #pragma omp simd aligned(dest,in1)
    for (i = end; i < length; ++i)
    {
        dest[i] = in1[i] * scalar;
    }

#endif
    return NOERR;
}


/*******************************************************************************
 VectorVectorMix */

Error_t
VectorVectorMix(float        *dest,
                const float  *in1,
                const float  *scalar1,
                const float  *in2,
                const float  *scalar2,
                unsigned     length)
{
#ifdef __APPLE__
  // Use the Accelerate framework if we have it
  vDSP_vsmsma(in1, 1, scalar1, in2, 1, scalar2, dest, 1, length);
#else
  unsigned i;
  #pragma omp simd aligned(dest,in1,scalar1,in2,scalar2)
  for (i = 0; i < length; ++i)
  {
    dest[i] = (in1[i] * (*scalar1)) + (in2[i] * (*scalar2));
  }
#endif
  return NOERR;
}


/*******************************************************************************
 VectorVectorMixD */

Error_t
VectorVectorMixD(double        *dest,
                 const double  *in1,
                 const double  *scalar1,
                 const double  *in2,
                 const double  *scalar2,
                 unsigned     length)
{
#ifdef __APPLE__
  // Use the Accelerate framework if we have it
  vDSP_vsmsmaD(in1, 1, scalar1, in2, 1, scalar2, dest, 1, length);
#else
  unsigned i;
  #pragma omp simd aligned(dest,in1,scalar1,in2,scalar2)
  for (i = 0; i < length; ++i)
  {
    dest[i] = (in1[i] * (*scalar1)) + (in2[i] * (*scalar2));
  }
#endif
  return NOERR;
}


/*******************************************************************************
 VectorVectorSumScale */

Error_t
VectorVectorSumScale(float        *dest,
                     const float  *in1,
                     const float  *in2,
                     const float  *scalar,
                     unsigned     length)
{
#ifdef __APPLE__
// Use the Accelerate framework if we have it
vDSP_vasm(in1, 1, in2, 1, scalar, dest, 1, length);

#else
unsigned i;
#pragma omp simd aligned(dest,in1,in2,scalar)
for (i = 0; i < length; ++i)
{
  dest[i] = (in1[i] + in2[i]) * (*scalar);
}
#endif
return NOERR;
}


/*******************************************************************************
 VectorVectorSumScaleD */
Error_t
VectorVectorSumScaleD(double        *dest,
                      const double  *in1,
                      const double  *in2,
                      const double  *scalar,
                      unsigned     length)
{
#ifdef __APPLE__
  // Use the Accelerate framework if we have it
  vDSP_vasmD(in1, 1, in2, 1, scalar, dest, 1, length);

#else
  unsigned i;
  #pragma omp simd aligned(dest,in1,in2,scalar)
  for (i = 0; i < length; ++i)
  {
    dest[i] = (in1[i] + in2[i]) * (*scalar);
  }
#endif
  return NOERR;
}



/*******************************************************************************
 VectorPower */
Error_t
VectorPower(float* dest, const float* in, float power, unsigned length)
{
#ifdef __APPLE__
    if (power == 2.0)
    {
        vDSP_vsq(in, 1, dest, 1, length);
    }
    else
    {
        unsigned i;
        const unsigned end = 4 * (length / 4);
        #pragma omp simd aligned(dest,in)
        for (i = 0; i < end; i+=4)
        {
            dest[i] = powf(in[i], power);
            dest[i + 1] = powf(in[i + 1], power);
            dest[i + 2] = powf(in[i + 2], power);
            dest[i + 3] = powf(in[i + 3], power);
        }
        #pragma omp simd aligned(dest,in)
        for (i = end; i < length; ++i)
        {
            dest[i] = powf(in[i], power);
        }
    }
#else
    unsigned i;
    const unsigned end = 4 * (length / 4);
    #pragma omp simd aligned(dest,in)
    for (i = 0; i < end; i+=4)
    {
        dest[i] = powf(in[i], power);
        dest[i + 1] = powf(in[i + 1], power);
        dest[i + 2] = powf(in[i + 2], power);
        dest[i + 3] = powf(in[i + 3], power);
    }
    #pragma omp simd aligned(dest,in)
    for (i = end; i < length; ++i)
    {
        dest[i] = powf(in[i], power);
    }
#endif
    return NOERR;

}


/*******************************************************************************
 VectorPowerD */
Error_t
VectorPowerD(double* dest, const double* in, double power, unsigned length)
{
#ifdef __APPLE__
    if (power == 2.0)
    {
        vDSP_vsqD(in, 1, dest, 1, length);
    }
    else
    {
        unsigned i;
        const unsigned end = 4 * (length / 4);
        #pragma omp simd aligned(dest,in)
        for (i = 0; i < end; i+=4)
        {
            dest[i] = pow(in[i], power);
            dest[i + 1] = pow(in[i + 1], power);
            dest[i + 2] = pow(in[i + 2], power);
            dest[i + 3] = pow(in[i + 3], power);
        }
        #pragma omp simd aligned(dest,in)
        for (i = end; i < length; ++i)
        {
            dest[i] = pow(in[i], power);
        }
    }
#else
    unsigned i;
    const unsigned end = 4 * (length / 4);
    #pragma omp simd aligned(dest,in)
    for (i = 0; i < end; i+=4)
    {
        dest[i] = pow(in[i], power);
        dest[i + 1] = pow(in[i + 1], power);
        dest[i + 2] = pow(in[i + 2], power);
        dest[i + 3] = pow(in[i + 3], power);
    }
    #pragma omp simd aligned(dest,in)
    for (i = end; i < length; ++i)
    {
        dest[i] = pow(in[i], power);
    }
#endif
    return NOERR;

}


/*******************************************************************************
 Convolve */
Error_t
Convolve(float       *in1,
         unsigned    in1_length,
         float       *in2,
         unsigned    in2_length,
         float       *dest)
{

    unsigned resultLength = in1_length + (in2_length - 1);
#ifdef __APPLE__
    //Use Native vectorized convolution function if available
    float    *in2_end = in2 + (in2_length - 1);
    unsigned signalLength = (in2_length + resultLength);

    float padded[signalLength];

    //float zero = 0.0;
    ClearBuffer(padded, signalLength);

    // Pad the input signal with (filter_length - 1) zeros.
    cblas_scopy(in1_length, in1, 1, (padded + (in2_length - 1)), 1);
    vDSP_conv(padded, 1, in2_end, -1, dest, 1, resultLength, in2_length);

#else
    // Use (boring, slow) canonical implementation
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


#endif
    return NOERR;
}


/*******************************************************************************
 ConvolveD */
Error_t
ConvolveD(double    *in1,
          unsigned  in1_length,
          double    *in2,
          unsigned  in2_length,
          double    *dest)
{

    unsigned resultLength = in1_length + (in2_length - 1);

#ifdef __APPLE__
    //Use Native vectorized convolution function if available
    double    *in2_end = in2 + (in2_length - 1);
    double signalLength = (in2_length + resultLength);

    // So there's some hella weird requirement that the signal input to
    //vDSP_conv has to be larger than (result_length + filter_length - 1),
    // (the output vector length) and it has to be zero-padded. What. The. Fuck!
    double padded[(unsigned)ceil(signalLength)];

    //float zero = 0.0;
    FillBufferD(padded, signalLength, 0.0);

    // Pad the input signal with (filter_length - 1) zeros.
    cblas_dcopy(in1_length, in1, 1, (padded + (in2_length - 1)), 1);
    vDSP_convD(padded, 1, in2_end, -1, dest, 1, resultLength, in2_length);

#else
    // Use (boring, slow) canonical implementation
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


#endif
    return NOERR;
}


/*******************************************************************************
 VectorDbConvert */
Error_t
VectorDbConvert(float* dest,
                const float* in,
                unsigned length)
{
#ifdef __APPLE__
    float one = 1.0;
    vDSP_vdbcon(in, 1, &one, dest, 1, length, 1);
#else
	#pragma omp simd aligned(dest,in)
    for (unsigned i = 0; i < length; ++i)
    {
        dest[i] = AMP_TO_DB(in[i]);
    }
#endif
    return NOERR;
}

/*******************************************************************************
 VectorDbConvertD */
Error_t
VectorDbConvertD(double*        dest,
                 const double*  in,
                 unsigned       length)
{
#ifdef __APPLE__
    double one = 1.0;
    vDSP_vdbconD(in, 1, &one, dest, 1, length, 1);
#else
	#pragma omp simd aligned(dest,in)
    for (unsigned i = 0; i < length; ++i)

    {
        dest[i] = AMP_TO_DBD(in[i]);
    }
#endif
    return NOERR;
}

/*******************************************************************************
 ComplexMultiply */
Error_t
ComplexMultiply(float*          re,
                float*          im,
                const float*    re1,
                const float*    im1,
                const float*    re2,
                const float*    im2,
                unsigned        length)
{
#if defined(__APPLE__)
    DSPSplitComplex in1 = {.realp = (float*)re1, .imagp = (float*)im1};
    DSPSplitComplex in2 = {.realp = (float*)re2, .imagp = (float*)im2};
    DSPSplitComplex out = {.realp = re, .imagp = im};
    vDSP_zvmul(&in1, 1, &in2, 1, &out, 1, length, 1);
#else
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
#endif
    return NOERR;
}


/*******************************************************************************
 ComplexMultiplyD */
Error_t
ComplexMultiplyD(double*        re,
                 double*        im,
                 const double*  re1,
                 const double*  im1,
                 const double*  re2,
                 const double*  im2,
                 unsigned       length)
{
#if defined(__APPLE__)
    DSPDoubleSplitComplex in1 = {.realp = (double*)re1, .imagp = (double*)im1};
    DSPDoubleSplitComplex in2 = {.realp = (double*)re2, .imagp = (double*)im2};
    DSPDoubleSplitComplex out = {.realp = re, .imagp = im};
    vDSP_zvmulD(&in1, 1, &in2, 1, &out, 1, length, 1);
#else
	#pragma omp simd aligned(re,im,re1,im1,re2,im2)
    for (unsigned i = 0; i < length; ++i)
    {
        double ire1 = re1[i];
        double iim1 = im1[i];
        double ire2 = re2[i];
        double iim2 = im2[i];
        re[i] = (ire1 * ire2 - iim1 * iim2);
        im[i] = (iim1 * ire2 + iim2 * ire1);
    }
#endif
    return NOERR;
}


/*******************************************************************************
 VectorRectToPolar */
Error_t
VectorRectToPolar(float*        magnitude,
                  float*        phase,
                  const float*  real,
                  const float*  imaginary,
                  unsigned      length)
{
#ifdef __APPLE__
    float dest[2*length];
    SplitToInterleaved(dest, real, imaginary, length);
    vDSP_polar(dest, 2, dest, 2, length);
    InterleavedToSplit(magnitude, phase, dest, length);
#else
    unsigned i;    
    const unsigned end = 4 * (length / 4);
    #pragma omp simd aligned(magnitude,phase,real,imaginary)
    for (i = 0; i < end; i+=4)
    {
        RectToPolar(real[i], imaginary[i], &magnitude[i], &phase[i]);
        RectToPolar(real[i+1], imaginary[i+1], &magnitude[i+1], &phase[i+1]);
        RectToPolar(real[i+2], imaginary[i+2], &magnitude[i+2], &phase[i+2]);
        RectToPolar(real[i+3], imaginary[i+3], &magnitude[i+3], &phase[i+3]);
    }
    for (i = end; i < length; ++i)
    {
        RectToPolar(real[i], imaginary[i], &magnitude[i], &phase[i]);
    }
#endif
    return NOERR;
}


/*******************************************************************************
 VectorRectToPolarD */
Error_t
VectorRectToPolarD(double*          magnitude,
                   double*          phase,
                   const double*    real,
                   const double*    imaginary,
                   unsigned         length)
{
#ifdef __APPLE__
    double dest[2*length];
    SplitToInterleavedD(dest, real, imaginary, length);
    vDSP_polarD(dest, 2, dest, 2, length);
    InterleavedToSplitD(magnitude, phase, dest, length);
#else
    unsigned i;
    const unsigned end = 4 * (length / 4);
    #pragma omp simd aligned(magnitude,phase,real,imaginary)
    for (i = 0; i < end; i+=4)
    {
        RectToPolarD(real[i], imaginary[i], &magnitude[i], &phase[i]);
        RectToPolarD(real[i+1], imaginary[i+1], &magnitude[i+1], &phase[i+1]);
        RectToPolarD(real[i+2], imaginary[i+2], &magnitude[i+2], &phase[i+2]);
        RectToPolarD(real[i+3], imaginary[i+3], &magnitude[i+3], &phase[i+3]);
    }
    for (i = end; i < length; ++i)
    {
        RectToPolarD(real[i], imaginary[i], &magnitude[i], &phase[i]);
    }
#endif
    return NOERR;
}


/*******************************************************************************
 MeanSquare */
float
MeanSquare(const float* data, unsigned length)
{
    float result = 0.;
#ifdef __APPLE__
    vDSP_measqv(data, 1, &result, length);
#else
    float scratch[length];
    VectorPower(scratch, data, 2, length);
    result = VectorSum(scratch, length) / length;
#endif
    return result;
}

/*******************************************************************************
 MeanSquareD */
double
MeanSquareD(const double* data, unsigned length)
{
    double result = 0.;
#ifdef __APPLE__
    vDSP_measqvD(data, 1, &result, length);
#else
    double scratch[length];
    VectorPowerD(scratch, data, 2, length);
    result = VectorSumD(scratch, length) / length;
#endif
    return result;
}

/*******************************************************************************
 BiquadFilterInit */
BiquadFilter* BiquadFilterInit(const float *bCoeff, const float *aCoeff)
{

    // Allocate Memory
    BiquadFilter* filter = (BiquadFilter*)malloc(sizeof(BiquadFilter));

    if (filter)
    {
        // Initialize Buffers
        CopyBuffer(filter->b, bCoeff, 3);
        CopyBuffer(filter->a, aCoeff, 2);

        ClearBuffer(filter->x, 2);
        ClearBuffer(filter->y, 2);
        ClearBuffer(filter->w, 2);
    }
    return filter;
}

BiquadFilterD* BiquadFilterInitD(const double  *bCoeff, const double  *aCoeff)
{

    // Allocate Memory
    BiquadFilterD* filter = (BiquadFilterD*)malloc(sizeof(BiquadFilterD));

    if (filter)
    {
        // Initialize Buffers
        CopyBufferD(filter->b, bCoeff, 3);
        CopyBufferD(filter->a, aCoeff, 2);

        ClearBufferD(filter->x, 2);
        ClearBufferD(filter->y, 2);
        ClearBufferD(filter->w, 2);
    }
    return filter;
}


/*******************************************************************************
 BiquadFilterFree */
Error_t
BiquadFilterFree(BiquadFilter * filter)
{
    if (filter)
    {
        free(filter);
        filter = NULL;
    }
    return NOERR;
}

Error_t
BiquadFilterFreeD(BiquadFilterD* filter)
{
    if (filter)
    {
        free(filter);
        filter = NULL;
    }
    return NOERR;
}

/*******************************************************************************
 BiquadFilterFlush */
Error_t
BiquadFilterFlush(BiquadFilter* filter)
{
    FillBuffer(filter->x, 2, 0.0);
	FillBuffer(filter->y, 2, 0.0);
    FillBuffer(filter->w, 2, 0.0);
    return NOERR;
}

Error_t
BiquadFilterFlushD(BiquadFilterD* filter)
{
    FillBufferD(filter->x, 2, 0.0);
    FillBufferD(filter->y, 2, 0.0);
    FillBufferD(filter->w, 2, 0.0);
    return NOERR;
}


/*******************************************************************************
 BiquadFilterProcess */
Error_t
BiquadFilterProcess(BiquadFilter    *filter,
                    float           *outBuffer,
                    const float     *inBuffer,
                    unsigned        n_samples)
{

#ifdef __APPLE__
    // Use accelerate if we have it
    float coeffs[5] = {
        filter->b[0], filter->b[1], filter->b[2],
        filter->a[0], filter->a[1]
    };
    float temp_in[n_samples + 2];
    float temp_out[n_samples + 2];


    // Put filter overlaps into beginning of input and output vectors
    cblas_scopy(2, filter->x, 1, temp_in, 1);
    cblas_scopy(2, filter->y, 1, temp_out, 1);
    cblas_scopy(n_samples, inBuffer, 1, (temp_in + 2), 1);

    // Process
    vDSP_deq22(temp_in, 1, coeffs, temp_out, 1, n_samples);

    // Write overlaps to filter x and y arrays
    cblas_scopy(2, (temp_in + n_samples), 1, filter->x, 1);
    cblas_scopy(2, (temp_out + n_samples), 1, filter->y, 1);

    // Write output
    cblas_scopy(n_samples, (temp_out + 2), 1, outBuffer, 1);


#else

    float buffer[n_samples];
    #pragma omp simd aligned(outBuffer,inBuffer)
    for (unsigned buffer_idx = 0; buffer_idx < n_samples; ++buffer_idx)
    {

        // DF-II Implementation
        buffer[buffer_idx] = filter->b[0] * inBuffer[buffer_idx] + filter->w[0];
        filter->w[0] = filter->b[1] * inBuffer[buffer_idx] - filter->a[0] * \
        buffer[buffer_idx] + filter->w[1];
        filter->w[1] = filter->b[2] * inBuffer[buffer_idx] - filter->a[1] * \
        buffer[buffer_idx];

    }

    // Write output
    CopyBuffer(outBuffer, buffer, n_samples);

#endif
    return NOERR;
}


Error_t
BiquadFilterProcessD(BiquadFilterD  *filter,
                     double         *outBuffer,
                     const double   *inBuffer,
                     unsigned       n_samples)
{

#ifdef __APPLE__
    // Use accelerate if we have it
    double coeffs[5] = {
        filter->b[0], filter->b[1], filter->b[2],
        filter->a[0], filter->a[1]
    };
    double temp_in[n_samples + 2];
    double temp_out[n_samples + 2];


    // Put filter overlaps into beginning of input and output vectors
    cblas_dcopy(2, filter->x, 1, temp_in, 1);
    cblas_dcopy(2, filter->y, 1, temp_out, 1);
    cblas_dcopy(n_samples, inBuffer, 1, (temp_in + 2), 1);

    // Process
    vDSP_deq22D(temp_in, 1, coeffs, temp_out, 1, n_samples);

    // Write overlaps to filter x and y arrays
    cblas_dcopy(2, (temp_in + n_samples), 1, filter->x, 1);
    cblas_dcopy(2, (temp_out + n_samples), 1, filter->y, 1);

    // Write output
    cblas_dcopy(n_samples, (temp_out + 2), 1, outBuffer, 1);


#else

    double buffer[n_samples];
    #pragma omp simd aligned(outBuffer,inBuffer)
    for (unsigned buffer_idx = 0; buffer_idx < n_samples; ++buffer_idx)
    {

        // DF-II Implementation
        buffer[buffer_idx] = filter->b[0] * inBuffer[buffer_idx] + filter->w[0];
        filter->w[0] = filter->b[1] * inBuffer[buffer_idx] - filter->a[0] * \
        buffer[buffer_idx] + filter->w[1];
        filter->w[1] = filter->b[2] * inBuffer[buffer_idx] - filter->a[1] * \
        buffer[buffer_idx];

    }

    // Write output
    CopyBufferD(outBuffer, buffer, n_samples);

#endif
    return NOERR;
}


/*******************************************************************************
 BiquadFilterTick */
float
BiquadFilterTick(BiquadFilter* filter, float in_sample)
{
    float out = filter->b[0] * in_sample + filter->w[0];
    filter->w[0] = filter->b[1] * in_sample - filter->a[0] * out + filter->w[1];
    filter->w[1] = filter->b[2] * in_sample - filter->a[1] * out;
    return out;
}

double
BiquadFilterTickD(BiquadFilterD* filter, double in_sample)
{
    double out = filter->b[0] * in_sample + filter->w[0];
    filter->w[0] = filter->b[1] * in_sample - filter->a[0] * out + filter->w[1];
    filter->w[1] = filter->b[2] * in_sample - filter->a[1] * out;
    return out;
}


/*******************************************************************************
 BiquadFilterUpdateKernel */
Error_t
BiquadFilterUpdateKernel(BiquadFilter   *filter,
                         const float    *bCoeff,
                         const float    *aCoeff)
{

    CopyBuffer(filter->b, bCoeff, 3);
    CopyBuffer(filter->a, aCoeff, 2);
    return NOERR;
}

Error_t
BiquadFilterUpdateKernelD(BiquadFilterD *filter,
                          const double  *bCoeff,
                          const double  *aCoeff)
{

    CopyBufferD(filter->b, bCoeff, 3);
    CopyBufferD(filter->a, aCoeff, 2);
    return NOERR;
}

/* Calculate BS.1770 prefilter coefficients for a given sample rate */
static void
calc_prefilter(float* b, float* a, float sample_rate);

static void
calc_prefilterD(double* b, double* a, double sample_rate);

/* Calculate BS.1770 rlb filter coefficients for a given sample rate */
static void
calc_rlbfilter(float* b, float* a, float sample_rate);

static void
calc_rlbfilterD(double* b, double* a, double sample_rate);

KWeightingFilter*
KWeightingFilterInit(float sample_rate)
{
    float b[3] = {0.};
    float a[2] = {0.};

    KWeightingFilter* filter = (KWeightingFilter*)malloc(sizeof(KWeightingFilter));

    if (filter)
    {
        calc_prefilter(b, a, sample_rate);
        filter->pre_filter = BiquadFilterInit(b, a);
        calc_rlbfilter(b, a, sample_rate);
        filter->rlb_filter = BiquadFilterInit(b, a);
    }

    return filter;
}

KWeightingFilterD*
KWeightingFilterInitD(double sample_rate)
{
    double b[3] = {0.};
    double a[2] = {0.};

    KWeightingFilterD* filter = (KWeightingFilterD*)malloc(sizeof(KWeightingFilterD));

    if (filter)
    {
        calc_prefilterD(b, a, sample_rate);
        filter->pre_filter = BiquadFilterInitD(b, a);
        calc_rlbfilterD(b, a, sample_rate);
        filter->rlb_filter = BiquadFilterInitD(b, a);
    }

    return filter;
}


Error_t
KWeightingFilterProcess(KWeightingFilter*   filter,
                        float*              dest,
                        const float*        src,
                        unsigned            length)
{
    float scratch[length];
    BiquadFilterProcess(filter->pre_filter, scratch, src, length);
    BiquadFilterProcess(filter->rlb_filter, dest, (const float*)scratch, length);
    return NOERR;
}


Error_t
KWeightingFilterProcessD(KWeightingFilterD* filter,
                         double*            dest,
                         const double*      src,
                         unsigned           length)
{
    double scratch[length];
    BiquadFilterProcessD(filter->pre_filter, scratch, src, length);
    BiquadFilterProcessD(filter->rlb_filter, dest, (const double*)scratch, length);
    return NOERR;
}

Error_t
KWeightingFilterFlush(KWeightingFilter* filter)
{
    if (filter)
    {
        BiquadFilterFlush(filter->pre_filter);
        BiquadFilterFlush(filter->rlb_filter);
        return NOERR;
    }
    return NULL_PTR_ERROR;
}

Error_t
KWeightingFilterFlushD(KWeightingFilterD* filter)
{
    if (filter)
    {
        BiquadFilterFlushD(filter->pre_filter);
        BiquadFilterFlushD(filter->rlb_filter);
        return NOERR;
    }
    return NULL_PTR_ERROR;
}


Error_t
KWeightingFilterFree(KWeightingFilter* filter)
{
    if (filter)
    {
        BiquadFilterFree(filter->pre_filter);
        BiquadFilterFree(filter->rlb_filter);
        free(filter);
        filter = NULL;
    }
    return NOERR;
}

Error_t
KWeightingFilterFreeD(KWeightingFilterD* filter)
{
    if (filter)
    {
        BiquadFilterFreeD(filter->pre_filter);
        BiquadFilterFreeD(filter->rlb_filter);
        free(filter);
        filter = NULL;
    }
    return NOERR;
}


BS1770Meter*
BS1770MeterInit(unsigned n_channels, float sample_rate)
{
    BS1770Meter* meter = (BS1770Meter*)malloc(sizeof(BS1770Meter));
    KWeightingFilter** filters = (KWeightingFilter**)malloc(n_channels * sizeof(KWeightingFilter*));
    Upsampler** upsamplers = (Upsampler**)malloc(n_channels * sizeof(Upsampler*));
    CircularBuffer** buffers = (CircularBuffer**)malloc(n_channels * sizeof(CircularBuffer*));
    if (meter && filters && upsamplers && buffers)
    {		
        for (unsigned i = 0; i < n_channels; ++i)        
        {
            filters[i] = KWeightingFilterInit(sample_rate);
            upsamplers[i] = UpsamplerInit(X4);
            buffers[i] = CircularBufferInit((unsigned)(2 * GATE_LENGTH_S * sample_rate));
        }

        meter->sample_count = 0;
        meter->n_channels = n_channels;
        meter->gate_len = (unsigned)(GATE_LENGTH_S * sample_rate);
        meter->overlap_len = (unsigned)(GATE_OVERLAP * meter->gate_len);
        meter->filters= filters;
        meter->upsamplers = upsamplers;
        meter->buffers = buffers;
    }
    else
    {
        if (meter)
        {
            free(meter);
        }

        if (filters)
        {
            free(filters);
        }

        if (upsamplers)
        {
            free(upsamplers);
        }
        if (buffers)
        {
            free(buffers);
        }
        return NULL;
    }
    return meter;
}

BS1770MeterD*
BS1770MeterInitD(unsigned n_channels, double sample_rate)
{
    BS1770MeterD* meter = (BS1770MeterD*)malloc(sizeof(BS1770MeterD));
    KWeightingFilterD** filters = (KWeightingFilterD**)malloc(n_channels * sizeof(KWeightingFilterD*));
    UpsamplerD** upsamplers = (UpsamplerD**)malloc(n_channels * sizeof(UpsamplerD*));
    CircularBufferD** buffers = (CircularBufferD**)malloc(n_channels * sizeof(CircularBufferD*));

    if (meter && filters && upsamplers && buffers)
    {
        for (unsigned i = 0; i < n_channels; ++i)
        {
            filters[i] = KWeightingFilterInitD(sample_rate);
            upsamplers[i] = UpsamplerInitD(X4);
            buffers[i] = CircularBufferInitD((unsigned)(2 * GATE_LENGTH_S * sample_rate));
        }

        meter->sample_count = 0;
        meter->n_channels = n_channels;
        meter->gate_len = (unsigned)(GATE_LENGTH_S * sample_rate);
        meter->overlap_len = (unsigned)(GATE_OVERLAP * meter->gate_len);
        meter->filters= filters;
        meter->upsamplers = upsamplers;
        meter->buffers = buffers;

    }
    else
    {
        if (meter)
        {
            free(meter);
        }

        if (filters)
        {
            free(filters);
        }

        if (upsamplers)
        {
            free(upsamplers);
        }
        if (buffers)
        {
            free(buffers);
        }
        return NULL;
    }
    return meter;
}



Error_t
BS1770MeterProcess(BS1770Meter*     meter,
                   float*           loudness,
                   float**          peaks,
                   const float**    samples,
                   unsigned         n_samples)
{
    unsigned os_length = 4 * n_samples;
    float filtered[n_samples];
    float os_sig[os_length];
    float gate[meter->gate_len];
    float sum = 0.0;

    if (meter)
    {
        *loudness = 0.0;
		#pragma omp simd aligned(loudness,peaks,samples)
        for (unsigned i = 0; i < meter->n_channels; ++i)
        {
            // Calculate peak for each channel
            UpsamplerProcess(meter->upsamplers[i], os_sig, samples[i], n_samples);
            VectorAbs(os_sig, (const float*)os_sig, os_length);
            *peaks[i] = AmpToDb(VectorMax(os_sig, os_length));

            KWeightingFilterProcess(meter->filters[i], filtered, samples[i], n_samples);
            CircularBufferWrite(meter->buffers[i], (const float*)filtered, n_samples);

            if (CircularBufferCount(meter->buffers[i]) >= meter->gate_len)
            {
                CircularBufferRead(meter->buffers[i], gate, meter->gate_len);
                CircularBufferRewind(meter->buffers[i], meter->overlap_len);
                sum += CHANNEL_GAIN[i] * MeanSquare(gate, meter->gate_len);
            }
        }

        *loudness = -0.691 + 10 * log10f(sum);
        return NOERR;
    }
    return NULL_PTR_ERROR;
}

Error_t
BS1770MeterProcessD(BS1770MeterD*   meter,
                    double*         loudness,
                    double**        peaks,
                    const double**  samples,
                    unsigned        n_samples)
{
    unsigned os_length = 4 * n_samples;
    double filtered[n_samples];
    double os_sig[os_length];
    double gate[meter->gate_len];
    double sum = 0.0;

    if (meter)
    {
        *loudness = 0.0;
		#pragma omp simd aligned(loudness,peaks,samples)
        for (unsigned i = 0; i < meter->n_channels; ++i)
        {
            // Calculate peak for each channel
            UpsamplerProcessD(meter->upsamplers[i], os_sig, samples[i], n_samples);
            VectorAbsD(os_sig, (const double*)os_sig, os_length);
            *peaks[i] = AmpToDbD(VectorMaxD(os_sig, os_length));

            KWeightingFilterProcessD(meter->filters[i], filtered, samples[i], n_samples);
            CircularBufferWriteD(meter->buffers[i], (const double*)filtered, n_samples);

            if (CircularBufferCountD(meter->buffers[i]) >= meter->gate_len)
            {
                CircularBufferReadD(meter->buffers[i], gate, meter->gate_len);
                CircularBufferRewindD(meter->buffers[i], meter->overlap_len);
                sum += CHANNEL_GAIN[i] * MeanSquareD(gate, meter->gate_len);
            }
        }

        *loudness = -0.691 + 10 * log10(sum);
        return NOERR;
    }
    return NULL_PTR_ERROR;
}


Error_t
BS1770MeterFree(BS1770Meter* meter)
{
    if (meter)
    {
        for (unsigned ch = 0; ch < meter->n_channels; ++ch)
        {
            if (meter->filters)
            {
                KWeightingFilterFree(meter->filters[ch]);
            }
            if (meter->buffers)
            {
                CircularBufferFree(meter->buffers[ch]);
            }
            if (meter->upsamplers)
            {
                UpsamplerFree(meter->upsamplers[ch]);
            }
        }
        free(meter);
        return NOERR;
    }
    return NULL_PTR_ERROR;
}

Error_t
BS1770MeterFreeD(BS1770MeterD* meter)
{
    if (meter)
    {
        for (unsigned ch = 0; ch < meter->n_channels; ++ch)
        {
            if (meter->filters)
            {
                KWeightingFilterFreeD(meter->filters[ch]);
            }
            if (meter->buffers)
            {
                CircularBufferFreeD(meter->buffers[ch]);
            }
            if (meter->upsamplers)
            {
                UpsamplerFreeD(meter->upsamplers[ch]);
            }
        }
        free(meter);
        return NOERR;
    }
    return NULL_PTR_ERROR;
}


static void
calc_prefilter(float* b, float* a, float sample_rate)
{

    if (sample_rate == 48000)
    {
        // Use exact coefficients from BS.1770
        b[0] = 1.53512485958697;
        b[1] = -2.69169618940638;
        b[2] = 1.19839281085285;

        a[0] = -1.69065929318241;
        a[1] = 0.73248077421585;
    }

    else
    {

        // Calculate prefilter as a high shelf using RBJ formula with the
        // following params:
        //    Fc = 1650
        //    Q = 0.807348173842115   (gives close approximation at Fs=48000)
        //    G = 4 (dB)
        float A = sqrtf(powf(10., PREFILTER_GAIN/20.));
        float wc = 2 * M_PI * PREFILTER_FC / sample_rate;
        float wS = sinf(wc);
        float wC = cosf(wc);
        float beta = sqrtf(A) / PREFILTER_Q;

        // Normalize filter by a[0]
        float norm = (A+1) - ((A - 1) * wC) + (beta * wS);

        a[0] = (2 * ((A - 1) - ((A + 1) * wC))) / norm;
        a[1] = ((A+1) - ((A - 1) * wC) - (beta * wS)) / norm;
        b[0] = (A * ((A + 1) + ((A - 1) * wC) + (beta * wS))) / norm;
        b[1] = (-2 * A * ((A - 1) + ((A + 1) * wC))) / norm;
        b[2] = (A * ((A + 1) + ((A - 1) * wC) - (beta * wS))) / norm;

    }

}

static void
calc_prefilterD(double* b, double* a, double sample_rate)
{
    if (sample_rate == 48000)
    {
        // Use exact coefficients from BS.1770
        b[0] = 1.53512485958697;
        b[1] = -2.69169618940638;
        b[2] = 1.19839281085285;

        a[0] = -1.69065929318241;
        a[1] = 0.73248077421585;
    }

    else
    {
        // Calculate prefilter as a high shelf using RBJ formula with the
        // following params:
        //    Fc = 1650
        //    Q = 0.807348173842115   (gives close approximation at Fs=4800)
        //    G = 4 (dB)
        double A = sqrt(pow(10., PREFILTER_GAIN / 20.));
        double wc = 2 * M_PI * PREFILTER_FC / sample_rate;
        double wS = sin(wc);
        double wC = cos(wc);
        double beta = sqrt(A) / PREFILTER_Q;

        // Normalize filter by a[0]
        double norm = (A+1) - ((A - 1) * wC) + (beta * wS);

        a[0] = (2 * ((A - 1) - ((A + 1) * wC))) / norm;
        a[1] = ((A+1) - ((A - 1) * wC) - (beta * wS)) / norm;
        b[0] = (A * ((A + 1) + ((A - 1) * wC) + (beta * wS))) / norm;
        b[1] = (-2 * A * ((A - 1) + ((A + 1) * wC))) / norm;
        b[2] = (A * ((A + 1) + ((A - 1) * wC) - (beta * wS))) / norm;
    }
}

static void
calc_rlbfilter(float* b, float* a, float sample_rate)
{

    if (sample_rate == 48000)
    {
        // Use exact coefficients from BS.1770
        b[0] = 1.0;
        b[1] = -2.0;
        b[2] = 1.0;

        a[0] = -1.99004745483398;
        a[1] = 0.99007225036621;
    }

    else
    {

        // Calculate as highpass using RLB formula
        float wc = 2 * M_PI * RLBFILTER_FC / sample_rate;
        float wS = sinf(wc);
        float wC = cosf(wc);
        float alpha = wS / (2.0 * RLBFILTER_Q);

        float norm = 1 + alpha;

        b[0] = ((1 + wC) / 2.0) / norm;
        b[1] = (-(1 + wC)) / norm;
        b[2] = ((1 + wC) / 2.0) / norm;

        a[0] = (-2 * wC) / norm;
        a[1] = (1 - alpha) / norm;

    }

}

static void
calc_rlbfilterD(double* b, double* a, double sample_rate)
{
    if (sample_rate == 48000)
    {
        // Use exact coefficients from BS.1770
        b[0] = 1.0;
        b[1] = -2.0;
        b[2] = 1.0;

        a[0] = -1.99004745483398;
        a[1] = 0.99007225036621;
    }

    else
    {
        // Calculate as highpass using RLB formula
        double wc = 2 * M_PI * RLBFILTER_FC / sample_rate;
        double wS = sin(wc);
        double wC = cos(wc);
        double alpha = wS / (2.0 * RLBFILTER_Q);

        double norm = 1 + alpha;

        b[0] = ((1 + wC) / 2.0) / norm;
        b[1] = (-(1 + wC)) / norm;
        b[2] = ((1 + wC) / 2.0) / norm;

        a[0] = (-2 * wC) / norm;
        a[1] = (1 - alpha) / norm;
    }
}




/*******************************************************************************
 CircularBufferInit */
CircularBuffer*
CircularBufferInit(unsigned length)
{
    CircularBuffer* cb = (CircularBuffer*)malloc(sizeof(CircularBuffer));
    if (cb)
    {
        // use next power of two so we can do a bitwise wrap
        length = next_pow2(length);
        float* buffer = (float*)malloc(length * sizeof(float));

        ClearBuffer(buffer, length);

        cb->length = length;
        cb->wrap = length - 1;
        cb->buffer = buffer;
        cb->read_index = 0;
        cb->write_index = 0;
        cb->count = 0;
    }
    return cb;
}



CircularBufferD*
CircularBufferInitD(unsigned length)
{
    CircularBufferD* cb = (CircularBufferD*)malloc(sizeof(CircularBufferD));
    if (cb)
    {
        // use next power of two so we can do a bitwise wrap
        length = next_pow2(length);
        double* buffer = (double*)malloc(length * sizeof(double));

        ClearBufferD(buffer, length);

        cb->length = length;
        cb->wrap = length - 1;
        cb->buffer = buffer;
        cb->read_index = 0;
        cb->write_index = 0;
        cb->count = 0;
    }
    return cb;
}


/*******************************************************************************
 CircularBufferFree */
Error_t
CircularBufferFree(CircularBuffer* cb)
{
    if (cb)
    {
        if (cb->buffer)
        {
            free(cb->buffer);
            cb->buffer = NULL;
        }
        free(cb);
        cb = NULL;
    }
    return NOERR;
}


Error_t
CircularBufferFreeD(CircularBufferD* cb)
{
    if (cb)
    {
        if (cb->buffer)
        {
            free(cb->buffer);
            cb->buffer = NULL;
        }
        free(cb);
        cb = NULL;
    }
    return NOERR;
}

/*******************************************************************************
 CircularBufferWrite */
Error_t
CircularBufferWrite(CircularBuffer* cb, const float* src, unsigned n_samples)
{
    for (unsigned i=0; i < n_samples; ++i)
    {
        cb->buffer[++cb->write_index & cb->wrap] = *src++;
    }
    cb->count += n_samples;

    if (cb->count <= cb->length)
    {
        return NOERR;
    }
    else
    {
        return VALUE_ERROR;
    }
}

Error_t
CircularBufferWriteD(CircularBufferD* cb, const double* src, unsigned n_samples)
{
    for (unsigned i=0; i < n_samples; ++i)
    {
        cb->buffer[++cb->write_index & cb->wrap] = *src++;
    }
    cb->count += n_samples;

    if (cb->count <= cb->length)
    {
        return NOERR;
    }
    else
    {
        return VALUE_ERROR;
    }
}


/*******************************************************************************
 CircularBufferRead */
Error_t
CircularBufferRead(CircularBuffer* cb, float* dest, unsigned n_samples)
{
    for (unsigned i=0; i < n_samples; ++i)
    {
        *dest++ = cb->buffer[++cb->read_index & cb->wrap];
    }
    cb->count -= n_samples;

    if (cb->count <= cb->length)
    {
        return NOERR;
    }
    else
    {
        return VALUE_ERROR;
    }
}


Error_t
CircularBufferReadD(CircularBufferD* cb, double* dest, unsigned n_samples)
{
    for (unsigned i=0; i < n_samples; ++i)
    {
        *dest++ = cb->buffer[++cb->read_index & cb->wrap];
    }
    cb->count -= n_samples;

    if (cb->count <= cb->length)
    {
        return NOERR;
    }
    else
    {
        return VALUE_ERROR;
    }
}


/*******************************************************************************
 CircularBufferFlush */
Error_t
CircularBufferFlush(CircularBuffer* cb)
{
    ClearBuffer(cb->buffer, cb->length);
    cb->count = 0;
    return NOERR;
}


Error_t
CircularBufferFlushD(CircularBufferD* cb)
{
    ClearBufferD(cb->buffer, cb->length);
    cb->count = 0;
    return NOERR;
}


/*******************************************************************************
 CircularBufferRewind */

Error_t
CircularBufferRewind(CircularBuffer* cb, unsigned samples)
{
    cb->read_index = ((cb->read_index + cb->length) - samples) % cb->length;
    cb->count += samples;

    if (cb->count <= cb->length)
    {
        return NOERR;
    }
    else
    {
        return VALUE_ERROR;
    }
}

Error_t
CircularBufferRewindD(CircularBufferD* cb, unsigned samples)
{
    cb->read_index = ((cb->read_index + cb->length) - samples) % cb->length;
    cb->count += samples;

    if (cb->count <= cb->length)
    {
        return NOERR;
    }
    else
    {
        return VALUE_ERROR;
    }
}

/*******************************************************************************
 CircularBufferCount */
unsigned
CircularBufferCount(CircularBuffer* cb)
{
    return cb->count;
}

unsigned
CircularBufferCountD(CircularBufferD* cb)
{
    return cb->count;
}


/* DecimatorInit *******************************************************/
Decimator*
DecimatorInit(ResampleFactor_t factor)
{
    unsigned n_filters = 1;
    switch(factor)
    {
        case X2:
            n_filters = 2;
            break;
        case X4:
            n_filters = 4;
            break;
        case X8:
            n_filters = 8;
            break;
      /*  case X16:
            n_filters = 16;
            break; */
        default:
            return NULL;
    }

    // Allocate memory for the upsampler
    Decimator* decimator = (Decimator*)malloc(sizeof(Decimator));

    // Allocate memory for the polyphase array
    FIRFilter** polyphase = (FIRFilter**)malloc(n_filters * sizeof(FIRFilter*));

    if (decimator && polyphase)
    {
        decimator->polyphase = polyphase;

        // Create polyphase filters
        unsigned idx;
        for(idx = 0; idx < n_filters; ++idx)
        {
            decimator->polyphase[idx] = FIRFilterInit(PolyphaseCoeffs[factor][idx], 64, DIRECT);
        }

        // Add factor
        decimator->factor = n_filters;

        return decimator;
    }
    else
    {
        if (polyphase)
        {
            free(polyphase);
        }
        if (decimator)
        {
            free(decimator);
        }
        return NULL;
    }
}

DecimatorD*
DecimatorInitD(ResampleFactor_t factor)
{
    unsigned n_filters = 1;
    switch(factor)
    {
        case X2:
            n_filters = 2;
            break;
        case X4:
            n_filters = 4;
            break;
        case X8:
            n_filters = 8;
            break;
        /*
        case X16:
            n_filters = 16;
            break;
        */
        default:
            return NULL;
    }

    // Allocate memory for the upsampler
    DecimatorD* decimator = (DecimatorD*)malloc(sizeof(DecimatorD));

    // Allocate memory for the polyphase array
    FIRFilterD** polyphase = (FIRFilterD**)malloc(n_filters * sizeof(FIRFilterD*));

    if (decimator && polyphase)
    {
        decimator->polyphase = polyphase;

        // Create polyphase filters
        unsigned idx;
        for(idx = 0; idx < n_filters; ++idx)
        {
            decimator->polyphase[idx] = FIRFilterInitD(PolyphaseCoeffsD[factor][idx], 64, DIRECT);
        }

        // Add factor
        decimator->factor = n_filters;

        return decimator;
    }
    else
    {
        if (polyphase)
        {
            free(polyphase);
        }
        if (decimator)
        {
            free(decimator);
        }
        return NULL;
    }
}


/* *****************************************************************************
 DecimatorFree */
Error_t
DecimatorFree(Decimator* decimator)
{
    if (decimator)
    {
        if (decimator->polyphase)
        {
            for (unsigned i = 0; i < decimator->factor; ++i)
            {
                FIRFilterFree(decimator->polyphase[i]);
            }
            free(decimator->polyphase);
        }
        free(decimator);
    }
    return NOERR;
}

Error_t
DecimatorFreeD(DecimatorD* decimator)
{
    if (decimator)
    {
        if (decimator->polyphase)
        {
            for (unsigned i = 0; i < decimator->factor; ++i)
            {
                FIRFilterFreeD(decimator->polyphase[i]);
            }
            free(decimator->polyphase);
        }
        free(decimator);
    }
    return NOERR;
}


/* *****************************************************************************
 DecimatorFlush */
Error_t
DecimatorFlush(Decimator* decimator)
{
    unsigned idx;
    for (idx = 0; idx < decimator->factor; ++idx)
    {
        FIRFilterFlush(decimator->polyphase[idx]);
    }
    return NOERR;
}

Error_t
DecimatorFlushD(DecimatorD* decimator)
{
    unsigned idx;
    for (idx = 0; idx < decimator->factor; ++idx)
    {
        FIRFilterFlushD(decimator->polyphase[idx]);
    }
    return NOERR;
}

/* *****************************************************************************
 DecimatorProcess */
Error_t
DecimatorProcess(Decimator      *decimator,
                 float          *outBuffer,
                 const float    *inBuffer,
                 unsigned       n_samples)
{
    if (decimator && outBuffer)
    {
        unsigned declen = n_samples / decimator->factor;
        float temp_buf[declen];
        ClearBuffer(outBuffer, declen);
		
        for (unsigned filt = 0; filt < decimator->factor; ++filt)
        {
            CopyBufferStride(temp_buf, 1, inBuffer, decimator->factor, declen);
            FIRFilterProcess(decimator->polyphase[filt], temp_buf, temp_buf, declen);
            VectorVectorAdd(outBuffer, (const float*)outBuffer, temp_buf, declen);
        }
        return NOERR;
    }
    else
    {
        return NULL_PTR_ERROR;
    }
}


Error_t
DecimatorProcessD(DecimatorD*   decimator,
                  double*       outBuffer,
                  const double* inBuffer,
                  unsigned      n_samples)
{
    if (decimator && outBuffer)
    {
        unsigned declen = n_samples / decimator->factor;
        double temp_buf[declen];
        ClearBufferD(outBuffer, declen);

        for (unsigned filt = 0; filt < decimator->factor; ++filt)
        {
            CopyBufferStrideD(temp_buf, 1, inBuffer, decimator->factor, declen);
            FIRFilterProcessD(decimator->polyphase[filt], temp_buf, temp_buf, declen);
            VectorVectorAddD(outBuffer, (const double*)outBuffer, temp_buf, declen);
        }
        return NOERR;
    }
    else
    {
        return NULL_PTR_ERROR;
    }
}



/*******************************************************************************
 DiodeRectifierInit */
DiodeRectifier*
DiodeRectifierInit(bias_t bias, float threshold)
{
    /* Allocate memory for diode struct */
    DiodeRectifier* diode = (DiodeRectifier*)malloc(sizeof(DiodeRectifier));

    if (NULL != diode)
    {
        /* Allocate scratch space */
        float* scratch = (float*)malloc(4096 * sizeof(float));

        if (NULL != scratch)
        {
            // Initialization
            diode->bias = bias;
            diode->threshold = threshold;
            diode->scratch = scratch;
            diode->abs_coeff = (bias == FULL_WAVE) ? 1.0 : 0.0;
            DiodeRectifierSetThreshold(diode, threshold);
        }
        else
        {
            free(diode);
            diode = NULL;
        }
    }
    return diode;
}

DiodeRectifierD*
DiodeRectifierInitD(bias_t bias, double threshold)
{
    /* Allocate memory for diode struct */
    DiodeRectifierD* diode = (DiodeRectifierD*)malloc(sizeof(DiodeRectifierD));

    if (NULL != diode)
    {
        /* Allocate scratch space */
        double* scratch = (double*)malloc(4096 * sizeof(double));

        if (NULL != scratch)
        {

            // Initialization
            diode->bias = bias;
            diode->threshold = threshold;
            diode->scratch = scratch;
            diode->abs_coeff = (bias == FULL_WAVE) ? 1.0 : 0.0;
            DiodeRectifierSetThresholdD(diode, threshold);
        }
        else
        {
            free(diode);
            diode = NULL;
        }
    }
    return diode;
}


/*******************************************************************************
 DiodeRectifierFree */
Error_t
DiodeRectifierFree(DiodeRectifier* diode)
{
    if (NULL != diode)
    {
        if (NULL != diode->scratch)
        {
            free(diode->scratch);
        }
        free(diode);
    }
    diode = NULL;
    return NOERR;
}

Error_t
DiodeRectifierFreeD(DiodeRectifierD* diode)
{
    if (NULL != diode)
    {
        if (NULL != diode->scratch)
        {
            free(diode->scratch);
        }
        free(diode);
    }
    diode = NULL;
    return NOERR;
}


/*******************************************************************************
 DiodeRectifierSetThreshold */
Error_t
DiodeRectifierSetThreshold(DiodeRectifier* diode, float threshold)
{
    float scale = 1.0;
    threshold = LIMIT(fabsf(threshold), 0.01, 0.9);
    scale = (1.0 - threshold);
    if (diode->bias== REVERSE_BIAS)
    {
        scale *= -1.0;
        threshold *= -1.0;
    }
    diode->threshold = threshold;
    diode->vt = -0.1738 * threshold + 0.1735;
    diode->scale = scale/(expf((1.0/diode->vt) - 1.));
    return NOERR;
}

Error_t
DiodeRectifierSetThresholdD(DiodeRectifierD* diode, double threshold)
{
    double scale = 1.0;
    threshold = LIMIT(fabs(threshold), 0.01, 0.9);
    scale = (1.0 - threshold);
    if (diode->bias == REVERSE_BIAS)
    {
        scale *= -1.0;
        threshold *= -1.0;
    }
    diode->threshold = threshold;
    diode->vt = -0.1738 * threshold + 0.1735;
    diode->scale = scale/(exp((1.0/diode->vt) - 1.));
    return NOERR;
}


/*******************************************************************************
 DiodeRectifierProcess */
Error_t
DiodeRectifierProcess(DiodeRectifier*   diode,
                      float*            out_buffer,
                      const float*      in_buffer,
                      unsigned          n_samples)
{
    float inv_vt = 1.0 / diode->vt;
    float scale = diode->scale;
    if (diode->bias == FULL_WAVE)
    {
        VectorAbs(diode->scratch, in_buffer, n_samples);
        VectorScalarMultiply(diode->scratch, diode->scratch, inv_vt, n_samples);
    }
    else
    {
        VectorScalarMultiply(diode->scratch, in_buffer, inv_vt, n_samples);
    }
    VectorScalarAdd(diode->scratch, diode->scratch, -1.0, n_samples);
    for (unsigned i = 0; i < n_samples; ++i)
    {
        out_buffer[i] = expf(diode->scratch[i]) * scale;
    }
    return NOERR;
}

Error_t
DiodeRectifierProcessD(DiodeRectifierD* diode,
                       double*          out_buffer,
                       const double*    in_buffer,
                       unsigned         n_samples)
{
    double inv_vt = 1.0 / diode->vt;
    double scale = diode->scale;
    if (diode->bias == FULL_WAVE)
    {
        VectorAbsD(diode->scratch, in_buffer, n_samples);
        VectorScalarMultiplyD(diode->scratch, diode->scratch, inv_vt, n_samples);
    }
    else
    {
        VectorScalarMultiplyD(diode->scratch, in_buffer, inv_vt, n_samples);
    }
    VectorScalarAddD(diode->scratch, diode->scratch, -1.0, n_samples);
    for (unsigned i = 0; i < n_samples; ++i)
    {
        out_buffer[i] = exp(diode->scratch[i]) * scale;
    }
    return NOERR;
}

/*******************************************************************************
 DiodeRectifierTick */
float
DiodeRectifierTick(DiodeRectifier* diode, float in_sample)
{
    return expf((in_sample/diode->vt)-1) * diode->scale;
}

double
DiodeRectifierTickD(DiodeRectifierD* diode, double in_sample)
{
    return exp((in_sample/diode->vt)-1) * diode->scale;
}




#ifdef __APPLE__
#include <Accelerate/Accelerate.h>
#elif defined(USE_BLAS)
#include <cblas.h>
#endif


/*******************************************************************************
 DiodeInit */
DiodeSaturator*
DiodeSaturatorInit(bias_t bias, float amount)
{
    // Create saturator struct
    DiodeSaturator* saturator = (DiodeSaturator*)malloc(sizeof(DiodeSaturator));
    
    // Initialization
    saturator->bias = bias;
    saturator->amount = amount;
    return saturator;
}

DiodeSaturatorD*
DiodeSaturatorInitD(bias_t bias, double amount)
{
    // Create saturator struct
    DiodeSaturatorD* saturator = (DiodeSaturatorD*)malloc(sizeof(DiodeSaturatorD));
    
    // Initialization
    saturator->bias = bias;
    saturator->amount = amount;
    return saturator;
}


/*******************************************************************************
 DiodeFree */
Error_t
DiodeSaturatorFree(DiodeSaturator* saturator)
{
    if(saturator)
        free(saturator);
    saturator = NULL;
    return NOERR;
}

Error_t
DiodeSaturatorFreeD(DiodeSaturatorD* saturator)
{
    if(saturator)
        free(saturator);
    saturator = NULL;
    return NOERR;
}


/*******************************************************************************
 DiodeSetThreshold */
Error_t
DiodeSaturatorSetAmount(DiodeSaturator* saturator, float amount)
{
    saturator->amount = 0.5 * powf(amount, 0.5);
    return NOERR;
}

Error_t
DiodeSaturatorSetThresholdD(DiodeSaturatorD* saturator, double amount)
{
    saturator->amount = 0.5 * pow(amount, 0.5);
    return NOERR;
}


/*******************************************************************************
 DiodeProcess */
Error_t
DiodeSaturatorProcess(DiodeSaturator*   saturator,
                      float*            out_buffer,
                      const float*      in_buffer,
                      unsigned          n_samples)
{
	#pragma omp simd aligned(out_buffer,in_buffer)
    for (unsigned i = 0; i < n_samples; ++i)
    {
        out_buffer[i] = in_buffer[i] - (saturator->amount * (F_EXP((in_buffer[i]/0.7) - 1.0) + E_INV));
    }
    return NOERR;
}

Error_t
DiodeSaturatorProcessD(DiodeSaturatorD* saturator,
                       double*          out_buffer,
                       const double*    in_buffer,
                       unsigned         n_samples)
{
	#pragma omp simd aligned(out_buffer,in_buffer)
    for (unsigned i = 0; i < n_samples; ++i)
    {
        out_buffer[i] = in_buffer[i] - (saturator->amount * (F_EXP((in_buffer[i]/0.7) - 1.0) + E_INV));
    }
    return NOERR;
}


/*******************************************************************************
 DiodeTick */
float
DiodeSaturatorTick(DiodeSaturator* saturator, float in_sample)
{
    return in_sample - (saturator->amount * (F_EXP((in_sample/0.7) - 1.0) + E_INV));
}

double
DiodeSaturatorTickD(DiodeSaturatorD* saturator, double in_sample)
{
    return in_sample - (saturator->amount * (F_EXP((in_sample/0.7) - 1.0) + E_INV));
}







/******************************************************************************
 STATIC FUNCTION DEFINITIONS */
#pragma mark - Static Function Definitions
#ifdef USE_FFTW_FFT
static inline void
interleave_complex(float*dest, const float* real, const float* imag, unsigned length)
{
#if defined(__APPLE__)
    DSPSplitComplex split = {.realp = (float*)real, .imagp = (float*)imag};
    vDSP_ztoc(&split, 1, (DSPComplex*)dest, 2, length/2);
#elif defined(__HAS_BLAS__)
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
interleave_complexD(double* dest, const double* real, const double* imag, unsigned length)
{
#if defined(__APPLE__)
    DSPDoubleSplitComplex split = {.realp = (double*)real, .imagp = (double*)imag};
    vDSP_ztocD(&split, 1, (DSPDoubleComplex*)dest, 2, length/2);
#elif defined(__HAS_BLAS__)
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
split_complex(float* real, float* imag, const float* data, unsigned length)
{
#if defined(__APPLE__)
    DSPSplitComplex split = {.realp = real, .imagp = imag};
    vDSP_ctoz((DSPComplex*)data, 2, &split, 1, length/2);
#elif defined(__HAS_BLAS__)
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


static inline void
split_complexD(double* real, double* imag, const double* data, unsigned length)
{
#if defined(__APPLE__)
    DSPDoubleSplitComplex split = {.realp = real, .imagp = imag};
    vDSP_ctozD((DSPDoubleComplex*)data, 2, &split, 1, length/2);
#elif defined(__HAS_BLAS__)
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

#endif



FFTConfig*
FFTInit(unsigned length)
{
    FFTConfig* fft = (FFTConfig*)malloc(sizeof(FFTConfig));
    float* split_realp = (float*)malloc(length * sizeof(float));
    float* split2_realp = (float*)malloc(length * sizeof(float));

    if (fft && split_realp && split2_realp)
    {
        fft->length = length;
        fft->scale = 1.0 / (fft->length);
        fft->log2n = log2f(fft->length);

        // Store these consecutively in memory
        fft->split.realp = split_realp;
        fft->split2.realp = split2_realp;
        fft->split.imagp = fft->split.realp + (fft->length / 2);
        fft->split2.imagp = fft->split2.realp + (fft->length / 2);
        fftwf_complex* c = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * length);
        float* r = (float*) fftwf_malloc(sizeof(float) * length);
        fft->setup.forward_plan = fftwf_plan_dft_r2c_1d(length, r, c, FFTW_MEASURE | FFTW_UNALIGNED);
        fft->setup.inverse_plan = fftwf_plan_dft_c2r_1d(length, c, r, FFTW_MEASURE | FFTW_UNALIGNED);
        fftwf_free(r);
        fftwf_free(c);
        ClearBuffer(split_realp, fft->length);
        ClearBuffer(split2_realp, fft->length);
        return fft;
    }
    else
    {
        // Cleanup
        if (fft)
            free(fft);
        if (split_realp)
            free(split_realp);
        if (split2_realp)
            free(split2_realp);
        return NULL;
    }
}



FFTConfigD*
FFTInitD(unsigned length)
{
    FFTConfigD* fft = (FFTConfigD*)malloc(sizeof(FFTConfigD));
    double* split_realp = (double*)malloc(length * sizeof(double));
    double* split2_realp = (double*)malloc(length * sizeof(double));

    if (fft && split_realp && split2_realp)
    {
        fft->length = length;
        fft->scale = 1.0 / (fft->length);
        fft->log2n = log2f(fft->length);

        // Store these consecutively in memory
        fft->split.realp = split_realp;
        fft->split2.realp = split2_realp;
        fft->split.imagp = fft->split.realp + (fft->length / 2);
        fft->split2.imagp = fft->split2.realp + (fft->length / 2);

        fftw_complex* c = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * length);
        double* r = (double*) fftw_malloc(sizeof(double) * length);
        fft->setup.forward_plan = fftw_plan_dft_r2c_1d(length, r, c, FFTW_MEASURE | FFTW_UNALIGNED);
        fft->setup.inverse_plan = fftw_plan_dft_c2r_1d(length, c, r, FFTW_MEASURE | FFTW_UNALIGNED);
        fftw_free(r);
        fftw_free(c);

        ClearBufferD(split_realp, fft->length);
        ClearBufferD(split2_realp, fft->length);

        return fft;
    }

    else
    {
        // Cleanup
        if (fft)
            free(fft);
        if (split_realp)
            free(split_realp);
        if (split2_realp)
            free(split2_realp);
        return NULL;
    }
}


Error_t
FFTFree(FFTConfig* fft)
{
    if (fft)
    {
        if (fft->split.realp)
        {
            free(fft->split.realp);
            fft->split.realp = NULL;
        }
        if (fft->split2.realp)
        {
            free(fft->split2.realp);
            fft->split2.realp = NULL;
        }
        if (fft->setup.forward_plan)
            fftwf_destroy_plan(fft->setup.forward_plan);
        if (fft->setup.inverse_plan)
            fftwf_destroy_plan(fft->setup.inverse_plan);
        free(fft);
        fft = NULL;
    }
    return NOERR;
}

Error_t
FFTFreeD(FFTConfigD* fft)
{
    if (fft)
    {
        if (fft->split.realp)
        {
            free(fft->split.realp);
            fft->split.realp = NULL;
        }
        if (fft->split2.realp)
        {
            free(fft->split2.realp);
            fft->split2.realp = NULL;
        }
        if (fft->setup.forward_plan)
            fftw_destroy_plan(fft->setup.forward_plan);
        if (fft->setup.inverse_plan)
            fftw_destroy_plan(fft->setup.inverse_plan);
        free(fft);
        fft = NULL;
    }
    return NOERR;
}



Error_t
FFT_R2C(FFTConfig*      fft,
        const float*    inBuffer,
        float*          real,
        float*          imag)
{
    FFTComplex temp[fft->length];
    fftwf_execute_dft_r2c(fft->setup.forward_plan, (float*)inBuffer, temp);
    split_complex(real, imag, (const float*)temp, fft->length);
    return NOERR;
}

Error_t
FFT_R2CD(FFTConfigD*    fft,
         const double*  inBuffer,
         double*        real,
         double*        imag)
{


    FFTComplexD temp[fft->length];
    fftw_execute_dft_r2c(fft->setup.forward_plan, (double*)inBuffer, temp);
    split_complexD(real, imag, (const double*)temp, fft->length);
    return NOERR;
}


Error_t
FFT_IR_R2C(FFTConfig*       fft,
           const float*     inBuffer,
           FFTSplitComplex  out)
{
    FFTComplex temp[fft->length];
    fftwf_execute_dft_r2c(fft->setup.forward_plan, (float*)inBuffer, temp);
    split_complex(out.realp, out.imagp, (const float*)temp, fft->length);
    out.imagp[0] = ((float*)temp)[fft->length];
    return NOERR;
}


Error_t
FFT_IR_R2CD(FFTConfigD*         fft,
            const double*       inBuffer,
            FFTSplitComplexD    out)
{
    FFTComplexD temp[fft->length];
    fftw_execute_dft_r2c(fft->setup.forward_plan, (double*)inBuffer, temp);
    split_complexD(out.realp, out.imagp, (const double*)temp, fft->length);
    out.imagp[0] = ((double*)temp)[fft->length];
    return NOERR;
}


#pragma mark - IFFT

Error_t
IFFT_C2R(FFTConfig*    fft,
         const float*  inReal,
         const float*  inImag,
         float*        out)
{
    FFTComplex temp[fft->length/2 + 1];
    interleave_complex((float*)temp, inReal, inImag, fft->length);
    ((float*)temp)[fft->length] = inReal[fft->length / 2 - 1];
    fftwf_execute_dft_c2r(fft->setup.inverse_plan, temp, out);
    VectorScalarMultiply(out, out, fft->scale, fft->length);
    return NOERR;
}

Error_t
IFFT_C2RD(FFTConfigD*   fft,
       const double*    inReal,
       const double*    inImag,
       double*          out)
{
    FFTComplexD temp[fft->length/2 + 1];
    interleave_complexD((double*)temp, inReal, inImag, fft->length);
    ((double*)temp)[fft->length] = inReal[fft->length / 2 - 1];
    fftw_execute_dft_c2r(fft->setup.inverse_plan, temp, out);
    VectorScalarMultiplyD(out, out, fft->scale, fft->length);
    return NOERR;
}

Error_t
FFTConvolve(FFTConfig*  fft,
            float       *in1,
            unsigned    in1_length,
            float       *in2,
            unsigned    in2_length,
            float       *dest)
{

    FFTComplex temp[fft->length];

    // Padded input buffers
    float in1_padded[fft->length];
    float in2_padded[fft->length];
    ClearBuffer(in1_padded, fft->length);
    ClearBuffer(in2_padded, fft->length);
    ClearBuffer(fft->split.realp, fft->length);
    ClearBuffer(fft->split2.realp, fft->length);

    // Zero pad the input buffers to FFT length
    CopyBuffer(in1_padded, in1, in1_length);
    CopyBuffer(in2_padded, in2, in2_length);

    fftwf_execute_dft_r2c(fft->setup.forward_plan, (float*)in1_padded, temp);
    float nyquist1 = ((float*)temp)[fft->length];
    split_complex(fft->split.realp, fft->split.imagp, (const float*)temp, fft->length);

    fftwf_execute_dft_r2c(fft->setup.forward_plan, (float*)in2_padded, temp);
    float nyquist2 = ((float*)temp)[fft->length];
    split_complex(fft->split2.realp, fft->split2.imagp, (const float*)temp, fft->length);


    float nyquist_out = nyquist1 * nyquist2;
    ComplexMultiply(fft->split.realp, fft->split.imagp, fft->split.realp,
                    fft->split.imagp, fft->split2.realp, fft->split2.imagp,
                    fft->length/2);
    interleave_complex((float*)temp, fft->split.realp, fft->split.imagp, fft->length);
    ((float*)temp)[fft->length] = nyquist_out;
    fftwf_execute_dft_c2r(fft->setup.inverse_plan, temp, dest);
    VectorScalarMultiply(dest, dest, fft->scale, fft->length);
    // Do the inverse FFT
    //IFFT_C2R(fft, fft->split.realp, fft->split.imagp, dest);
    return NOERR;
}



Error_t
FFTConvolveD(FFTConfigD*    fft,
             const double*  in1,
             unsigned       in1_length,
             const double*  in2,
             unsigned       in2_length,
             double*        dest)
{

    FFTComplexD temp[fft->length];

    // Padded input buffers
    double in1_padded[fft->length];
    double in2_padded[fft->length];
    ClearBufferD(in1_padded, fft->length);
    ClearBufferD(in2_padded, fft->length);
    ClearBufferD(fft->split.realp, fft->length);
    ClearBufferD(fft->split2.realp, fft->length);

    // Zero pad the input buffers to FFT length
    CopyBufferD(in1_padded, in1, in1_length);
    CopyBufferD(in2_padded, in2, in2_length);

    fftw_execute_dft_r2c(fft->setup.forward_plan, (double*)in1_padded, temp);
    double nyquist1 = ((double*)temp)[fft->length];
    split_complexD(fft->split.realp, fft->split.imagp, (const double*)temp, fft->length);

    fftw_execute_dft_r2c(fft->setup.forward_plan, (double*)in2_padded, temp);
    double nyquist2 = ((double*)temp)[fft->length];
    split_complexD(fft->split2.realp, fft->split2.imagp, (const double*)temp, fft->length);


    double nyquist_out = nyquist1 * nyquist2;
    ComplexMultiplyD(fft->split.realp, fft->split.imagp, fft->split.realp,
                     fft->split.imagp, fft->split2.realp, fft->split2.imagp,
                     fft->length / 2);


    interleave_complexD((double*)temp, fft->split.realp, fft->split.imagp, fft->length);
    ((double*)temp)[fft->length] = nyquist_out;
    fftw_execute_dft_c2r(fft->setup.inverse_plan, temp, dest);
    VectorScalarMultiplyD(dest, dest, fft->scale, fft->length);
    return NOERR;    
}


/* FIRFilterInit *******************************************************/
FIRFilter*
FIRFilterInit(const float*       filter_kernel,
                     unsigned           length,
                     ConvolutionMode_t  convolution_mode)
{

    // Array lengths and sizes
    unsigned kernel_length = length;                    // IN SAMPLES!
    unsigned overlap_length = kernel_length - 1;        // IN SAMPLES!



    // Allocate Memory
    FIRFilter* filter = (FIRFilter*)malloc(sizeof(FIRFilter));
    float* kernel = (float*)malloc(kernel_length * sizeof(float));
    float* overlap = (float*)malloc(overlap_length * sizeof(float));
    if (filter && kernel && overlap)
    {
        // Initialize Buffers
        CopyBuffer(kernel, filter_kernel, kernel_length);
        ClearBuffer(overlap, overlap_length);

        // Set up the struct
        filter->kernel = kernel;
        filter->kernel_end = filter_kernel + (kernel_length - 1);
        filter->overlap = overlap;
        filter->kernel_length = kernel_length;
        filter->overlap_length = overlap_length;
        filter->fft_config = NULL;
        filter->fft_kernel.realp = NULL;
        filter->fft_kernel.imagp = NULL;

        if (((convolution_mode == BEST) &&
             (kernel_length < USE_FFT_CONVOLUTION_LENGTH)) ||
            (convolution_mode == DIRECT))
        {
            filter->conv_mode = DIRECT;
        }

        else
        {
            filter->conv_mode = FFT;
        }

        return filter;
    }

    else
    {
        free(filter);
        free(kernel);
        free(overlap);
        return NULL;
    }
}


FIRFilterD*
FIRFilterInitD(const double*        filter_kernel,
               unsigned             length,
               ConvolutionMode_t    convolution_mode)
{

    // Array lengths and sizes
    unsigned kernel_length = length;                    // IN SAMPLES!
    unsigned overlap_length = kernel_length - 1;        // IN SAMPLES!



    // Allocate Memory
    FIRFilterD* filter = (FIRFilterD*)malloc(sizeof(FIRFilterD));
    double* kernel = (double*)malloc(kernel_length * sizeof(double));
    double* overlap = (double*)malloc(overlap_length * sizeof(double));
    if (filter && kernel && overlap)
    {
        // Initialize Buffers
        CopyBufferD(kernel, filter_kernel, kernel_length);
        ClearBufferD(overlap, overlap_length);

        // Set up the struct
        filter->kernel = kernel;
        filter->kernel_end = filter_kernel + (kernel_length); //- 1);
        filter->overlap = overlap;
        filter->kernel_length = kernel_length;
        filter->overlap_length = overlap_length;
        filter->fft_config = NULL;
        filter->fft_kernel.realp = NULL;
        filter->fft_kernel.imagp = NULL;

        if (((convolution_mode == BEST) &&
             (kernel_length < USE_FFT_CONVOLUTION_LENGTH)) ||
            (convolution_mode == DIRECT))
        {
            filter->conv_mode = DIRECT;
        }

        else
        {
            filter->conv_mode = FFT;
        }

        return filter;
    }

    else
    {
        free(filter);
        free(kernel);
        free(overlap);
        return NULL;
    }
}

/* FIRFilterFree *******************************************************/
Error_t
FIRFilterFree(FIRFilter * filter)
{
    if (filter)
    {
        if (filter->kernel)
        {
            free(filter->kernel);
            filter->kernel = NULL;
        }
        if (filter->overlap)
        {
            free(filter->overlap);
            filter->overlap = NULL;
        }

        if (filter->fft_config)
        {
            //FFTFree(filter->fft_config);
            filter->fft_config = NULL;
        }

        if (filter->fft_kernel.realp)
        {
            free(filter->fft_kernel.realp);
            filter->fft_kernel.realp = NULL;
        }
        free(filter);
        filter = NULL;
    }
    return NOERR;
}

Error_t
FIRFilterFreeD(FIRFilterD * filter)
{
    if (filter)
    {
        if (filter->kernel)
        {
            free(filter->kernel);
            filter->kernel = NULL;
        }
        if (filter->overlap)
        {
            free(filter->overlap);
            filter->overlap = NULL;
        }

        if (filter->fft_config)
        {
            FFTFreeD(filter->fft_config);
            filter->fft_config = NULL;
        }

        if (filter->fft_kernel.realp)
        {
            free(filter->fft_kernel.realp);
            filter->fft_kernel.realp = NULL;
        }

        free(filter);
        filter = NULL;
    }
    return NOERR;
}


/* FIRFilterFlush ******************************************************/
Error_t
FIRFilterFlush(FIRFilter* filter)
{
    // The only stateful part of this is the overlap buffer, so this just
    //zeros it out
    ClearBuffer(filter->overlap, filter->overlap_length);
    return NOERR;
}

Error_t
FIRFilterFlushD(FIRFilterD* filter)
{
    // The only stateful part of this is the overlap buffer, so this just
    //zeros it out
    ClearBufferD(filter->overlap, filter->overlap_length);
    return NOERR;
}


/* FIRFilterProcess ****************************************************/
Error_t
FIRFilterProcess(FIRFilter*     filter,
                 float*         outBuffer,
                 const float*   inBuffer,
                 unsigned       n_samples)
{

    if (filter)
    {
        // Do direct convolution
        if (filter->conv_mode == DIRECT)
        {
            unsigned resultLength = n_samples + (filter->kernel_length - 1);
            // Temporary buffer to store full result of filtering..
            float buffer[resultLength];

            Convolve((float*)inBuffer, n_samples,
                            filter->kernel, filter->kernel_length, buffer);

            // Add in the overlap from the last block
            VectorVectorAdd(buffer, filter->overlap, buffer, filter->overlap_length);
            CopyBuffer(filter->overlap, buffer + n_samples, filter->overlap_length);
            CopyBuffer(outBuffer, buffer, n_samples);
        }

        // Otherwise do FFT Convolution
        else
        {
            // Configure the FFT on the first run, that way we can figure out how
            // long the input blocks are going to be. This makes the filter more
            // complicated internally in order to make the convolution transparent.
            // Calculate length of FFT
            if(filter->fft_config == 0)
            {
                // Calculate FFT Length
                filter->fft_length = next_pow2(n_samples + filter->kernel_length - 1);
                filter->fft_config = FFTInit(filter->fft_length);

                // fft kernel buffers
                float padded_kernel[filter->fft_length];

                // Allocate memory for filter kernel
                filter->fft_kernel.realp = (float*) malloc(filter->fft_length * sizeof(float));
                filter->fft_kernel.imagp = filter->fft_kernel.realp +(filter->fft_length / 2);

                // Write zero padded kernel to buffer
                CopyBuffer(padded_kernel, filter->kernel, filter->kernel_length);
                ClearBuffer((padded_kernel + filter->kernel_length), (filter->fft_length - filter->kernel_length));

                // Calculate FFT of filter kernel
                FFT_IR_R2C(filter->fft_config, padded_kernel, filter->fft_kernel);
            }

            // Buffer for transformed input
            float buffer[filter->fft_length];

            // Convolve
            FFTFilterConvolve(filter->fft_config, (float*)inBuffer, n_samples, filter->fft_kernel, buffer);

            // Add in the overlap from the last block
            VectorVectorAdd(buffer, filter->overlap, buffer, filter->overlap_length);
            CopyBuffer(filter->overlap, buffer + n_samples, filter->overlap_length);
            CopyBuffer(outBuffer, buffer, n_samples);

        }
        return NOERR;
    }

    else
    {
        return ERROR;
    }
}


Error_t
FIRFilterProcessD(FIRFilterD*   filter,
                 double*        outBuffer,
                 const double*  inBuffer,
                 unsigned       n_samples)
{

    if (filter)
    {

        // Do direct convolution
        if (filter->conv_mode == DIRECT)
        {
            unsigned resultLength = n_samples + (filter->kernel_length - 1);

            // Temporary buffer to store full result of filtering..
            double buffer[resultLength];

            ConvolveD((double*)inBuffer, n_samples,
                     filter->kernel, filter->kernel_length, buffer);
            // Add in the overlap from the last block
            VectorVectorAddD(buffer, filter->overlap, buffer, filter->overlap_length);
            CopyBufferD(filter->overlap, buffer + n_samples, filter->overlap_length);
            CopyBufferD(outBuffer, buffer, n_samples);
        }

        // Otherwise do FFT Convolution
        else
        {
            // Configure the FFT on the first run, that way we can figure out how
            // long the input blocks are going to be. This makes the filter more
            // complicated internally in order to make the convolution transparent.
            // Calculate length of FFT
            if(filter->fft_config == 0)
            {
                // Calculate FFT Length
                filter->fft_length = next_pow2(n_samples + filter->kernel_length - 1); //2 * next_pow2(n_samples + filter->kernel_length - 1);
                filter->fft_config = FFTInitD(filter->fft_length);

                // fft kernel buffers
                double padded_kernel[filter->fft_length];

                // Allocate memory for filter kernel
                filter->fft_kernel.realp = (double*) malloc(filter->fft_length * sizeof(double));
                filter->fft_kernel.imagp = filter->fft_kernel.realp +(filter->fft_length / 2);

                // Write zero padded kernel to buffer
                CopyBufferD(padded_kernel, filter->kernel, filter->kernel_length);
                ClearBufferD((padded_kernel + filter->kernel_length), (filter->fft_length - filter->kernel_length));

                // Calculate FFT of filter kernel
                FFT_IR_R2CD(filter->fft_config, padded_kernel, filter->fft_kernel);
            }

            // Buffer for transformed input
            double buffer[filter->fft_length];

            // Convolve
            FFTFilterConvolveD(filter->fft_config, (double*)inBuffer, n_samples, filter->fft_kernel, buffer);

            // Add in the overlap from the last block
            VectorVectorAddD(buffer, filter->overlap, buffer, filter->overlap_length);
            CopyBufferD(filter->overlap, buffer + n_samples, filter->overlap_length);
            CopyBufferD(outBuffer, buffer, n_samples);

        }
        return NOERR;
    }

    else
    {
        return ERROR;
    }
}

/* FIRFilterUpdateKernel ***********************************************/
Error_t
FIRFilterUpdateKernel(FIRFilter*  filter, const float* filter_kernel)
{
    // Copy the new kernel into the filter
    CopyBuffer(filter->kernel, filter_kernel, filter->kernel_length);
    return NOERR;
}

Error_t
FIRFilterUpdateKernelD(FIRFilterD*  filter, const double* filter_kernel)
{
    // Copy the new kernel into the filter
    CopyBufferD(filter->kernel, filter_kernel, filter->kernel_length);
    return NOERR;
}


/* LadderFilterInit ****************************************************/
LadderFilter*
LadderFilterInit(float _sample_rate)
{
    LadderFilter *filter = (LadderFilter*)malloc(sizeof(LadderFilter));
    if (filter)
    {
        ClearBuffer(filter->y, 4);
        ClearBuffer(filter->w, 4);
        filter->Vt = 0.026;
        filter->cutoff = 0;
        filter->resonance = 0;
        filter->sample_rate = _sample_rate;
    }
    return filter;

}
Error_t
LadderFilterSetCutoff(LadderFilter* filter, float c)
{
    if(filter) filter->cutoff = c;
    return NOERR;
}
Error_t
LadderFilterSetResonance(LadderFilter* filter, float r)
{
    if(filter) filter->resonance = r;
    return NOERR;
}

/* LadderFilterFree ****************************************************/
Error_t
LadderFilterFree(LadderFilter* filter)
{
    if (filter)
    {
        free(filter);
        filter = NULL;
    }
    return NOERR;
}


/* LadderFilterFlush ***************************************************/
Error_t
LadderFilterFlush(LadderFilter* filter)
{
    FillBuffer(filter->y, 4, 0.0);
    FillBuffer(filter->w, 4, 0.0);
    filter->cutoff = 0;
    filter->resonance = 0;
    return NOERR;
}


/* LadderFilterProcess *************************************************/
Error_t
LadderFilterProcess(LadderFilter *filter, float *outBuffer, float *inBuffer, unsigned n_samples)
{
    // Pre-calculate Scalars
    float TWO_VT_INV = 1.0 / (2 * filter->Vt);
    float TWO_VT_G = 2 * filter->Vt * (1 - exp(-TWO_PI * filter->cutoff / filter->sample_rate));

    // Filter audio
    #pragma omp simd aligned(outBuffer,inBuffer)
    for (unsigned i = 0; i < n_samples; ++i)
    {

        // Stage 1 output
        filter->y[0] = filter->y[0] + TWO_VT_G * (f_tanh(inBuffer[i] - 4 * \
                       f_tanh(2 * filter->resonance * filter->y[3]) * \
                       TWO_VT_INV) - filter->w[0]);

        filter->w[0] = f_tanh(filter->y[0] * TWO_VT_INV);

        // Stage 2 output
        filter->y[1] = filter->y[1] + TWO_VT_G * (filter->w[0]- filter->w[1]);
        filter->w[1] = f_tanh(filter->y[1] * TWO_VT_INV);

        // Stage 3 output
        filter->y[2] = filter->y[2] + TWO_VT_G * (filter->w[1]- filter->w[2]);
        filter->w[2] = f_tanh(filter->y[2] * TWO_VT_INV);

        // Stage 4 output
        filter->y[3] = filter->y[3] + TWO_VT_G * (filter->w[2]- filter->w[3]);
        filter->w[3] = f_tanh(filter->y[3] * TWO_VT_INV);

        // Write output
        outBuffer[i] = filter->y[3];

    }

    return NOERR;
}

Error_t
LadderFilterSetTemperature(LadderFilter *filter, float tempC)
{
    float T = tempC + 273.15;
    filter->Vt = BOLTZMANS_CONSTANT * T / Q;
    return NOERR;
}


/* LadderFilterInit ****************************************************/
LadderFilterD*
LadderFilterInitD(double _sample_rate)
{
    LadderFilterD *filter = (LadderFilterD*)malloc(sizeof(LadderFilterD));
    if (filter)
    {
        ClearBufferD(filter->y, 4);
        ClearBufferD(filter->w, 4);
        filter->Vt = 0.026;
        filter->cutoff = 0;
        filter->resonance = 0;
        filter->sample_rate = _sample_rate;
    }
    return filter;

}
Error_t
LadderFilterSetCutoffD(LadderFilterD* filter, double c)
{
    if(filter) filter->cutoff = c;
    return NOERR;
}
Error_t
LadderFilterSetResonanceD(LadderFilterD* filter, double r)
{
    if(filter) filter->resonance = r;
    return NOERR;
}

/* LadderFilterFree ****************************************************/
Error_t
LadderFilterFreeD(LadderFilterD* filter)
{
    if (filter)
    {
        free(filter);
        filter = NULL;
    }
    return NOERR;
}


/* LadderFilterFlush ***************************************************/
Error_t
LadderFilterFlushD(LadderFilterD* filter)
{
    FillBufferD(filter->y, 4, 0.0);
    FillBufferD(filter->w, 4, 0.0);
    filter->cutoff = 0;
    filter->resonance = 0;
    return NOERR;
}


/* LadderFilterProcess *************************************************/
Error_t
LadderFilterProcessD(LadderFilterD *filter, double *outBuffer, double *inBuffer, unsigned n_samples)
{
    // Pre-calculate Scalars
    double TWO_VT_INV = 1.0 / (2 * filter->Vt);
    double TWO_VT_G = 2 * filter->Vt * (1 - exp(-TWO_PI * filter->cutoff / filter->sample_rate));

    // Filter audio
    #pragma omp simd aligned(outBuffer,inBuffer)
    for (unsigned i = 0; i < n_samples; ++i)
    {

        // Stage 1 output
        filter->y[0] = filter->y[0] + TWO_VT_G * (f_tanh(inBuffer[i] - 4 * \
                       f_tanh(2 * filter->resonance * filter->y[3]) * \
                       TWO_VT_INV) - filter->w[0]);

        filter->w[0] = f_tanh(filter->y[0] * TWO_VT_INV);

        // Stage 2 output
        filter->y[1] = filter->y[1] + TWO_VT_G * (filter->w[0]- filter->w[1]);
        filter->w[1] = f_tanh(filter->y[1] * TWO_VT_INV);

        // Stage 3 output
        filter->y[2] = filter->y[2] + TWO_VT_G * (filter->w[1]- filter->w[2]);
        filter->w[2] = f_tanh(filter->y[2] * TWO_VT_INV);

        // Stage 4 output
        filter->y[3] = filter->y[3] + TWO_VT_G * (filter->w[2]- filter->w[3]);
        filter->w[3] = f_tanh(filter->y[3] * TWO_VT_INV);

        // Write output
        outBuffer[i] = filter->y[3];

    }

    return NOERR;
}

Error_t
LadderFilterSetTemperatureD(LadderFilterD *filter, double tempC)
{
    double T = tempC + 273.15;
    filter->Vt = BOLTZMANS_CONSTANT * T / Q;
    return NOERR;
}

/* LRFilterInit ***********************************************************/
LRFilter*
LRFilterInit(Filter_t   type,
             float 		cutoff,
             float      Q,
             float      sampleRate)
{
    LRFilter *filter = (LRFilter*) malloc(sizeof(LRFilter));
    filter->type = type;
    filter->cutoff = cutoff;
    filter->Q = Q;
    filter->sampleRate = sampleRate;
    filter->filterA = RBJFilterInit(filter->type, filter->cutoff, filter->sampleRate);
    filter->filterB = RBJFilterInit(filter->type, filter->cutoff, filter->sampleRate);
    RBJFilterSetQ(filter->filterA, filter->Q);
    RBJFilterSetQ(filter->filterB, filter->Q);
    return filter;
}

LRFilterD*
LRFilterInitD(Filter_t  type,
              double 	cutoff,
              double    Q,
              double    sampleRate)
{
    LRFilterD *filter = (LRFilterD*) malloc(sizeof(LRFilterD));
    filter->type = type;
    filter->cutoff = cutoff;
    filter->Q = Q;
    filter->sampleRate = sampleRate;
    filter->filterA = RBJFilterInitD(filter->type, filter->cutoff, filter->sampleRate);
    filter->filterB = RBJFilterInitD(filter->type, filter->cutoff, filter->sampleRate);
    RBJFilterSetQD(filter->filterA, filter->Q);
    RBJFilterSetQD(filter->filterB, filter->Q);
    return filter;
}

/* LRFilterFree ***********************************************************/
Error_t
LRFilterFree(LRFilter* 	filter)
{
    RBJFilterFree(filter->filterA);
    RBJFilterFree(filter->filterB);
    if (filter)
    {
        free(filter);
        filter = NULL;
    }
    
    return NOERR;
}

Error_t
LRFilterFreeD(LRFilterD* filter)
{
    RBJFilterFreeD(filter->filterA);
    RBJFilterFreeD(filter->filterB);
    if (filter)
    {
        free(filter);
        filter = NULL;
    }
    
    return NOERR;
}

/* LRFilterFlush **********************************************************/
Error_t
LRFilterFlush(LRFilter* filter)
{
    RBJFilterFlush(filter->filterA);
    RBJFilterFlush(filter->filterB);
    
    return NOERR;
}

Error_t
LRFilterFlushD(LRFilterD* filter)
{
    RBJFilterFlushD(filter->filterA);
    RBJFilterFlushD(filter->filterB);
    
    return NOERR;
}

/* LRFilterSetParams ******************************************************/
Error_t
LRFilterSetParams(LRFilter* filter,
                      Filter_t  type,
                      float     cutoff,
                      float     Q)
{
    filter->type = type;
    filter->cutoff = cutoff;
    filter->Q = Q;
    RBJFilterSetParams(filter->filterA, filter->type, filter->cutoff, filter->Q);
    RBJFilterSetParams(filter->filterB, filter->type, filter->cutoff, filter->Q);
    
    return NOERR;
}

Error_t
LRFilterSetParamsD(LRFilterD*   filter,
                  Filter_t      type,
                  double        cutoff,
                  double        Q)
{
    filter->type = type;
    filter->cutoff = cutoff;
    filter->Q = Q;
    RBJFilterSetParamsD(filter->filterA, filter->type, filter->cutoff, filter->Q);
    RBJFilterSetParamsD(filter->filterB, filter->type, filter->cutoff, filter->Q);
    
    return NOERR;
}

/* LRFilterProcess ********************************************************/
Error_t
LRFilterProcess(LRFilter*       filter,
                float*          outBuffer,
                const float* 	inBuffer,
                unsigned 		n_samples)

{
    float tempBuffer[n_samples];
    RBJFilterProcess(filter->filterA, tempBuffer, inBuffer, n_samples);
    RBJFilterProcess(filter->filterB, outBuffer, tempBuffer, n_samples);
    return NOERR;
}

Error_t
LRFilterProcessD(LRFilterD*     filter,
                 double*        outBuffer,
                 const double*  inBuffer,
                 unsigned 		n_samples)

{
    double tempBuffer[n_samples];
    RBJFilterProcessD(filter->filterA, tempBuffer, inBuffer, n_samples);
    RBJFilterProcessD(filter->filterB, outBuffer, tempBuffer, n_samples);
    return NOERR;
}




float
phase_correlation(float* left, float* right, unsigned n_samples)
{
    float product = 0.0;
    float lsq = 0.0;
    float rsq = 0.0;
    float denom = 0.0;

    for (unsigned i = 0; i < n_samples; ++i)
    {
        float left_sample = left[i];
        float right_sample = right[i];

        product += left_sample * right_sample;
        lsq += left_sample * left_sample;
        rsq += right_sample * right_sample;
    }

    denom = lsq * rsq;

    if (denom > 0.0)
    {
        return product / sqrtf(denom);
    }
    else
    {
        return 1.0;
    }
}


double
phase_correlationD(double* left, double* right, unsigned n_samples)
{
    double product = 0.0;
    double lsq = 0.0;
    double rsq = 0.0;
    double denom = 0.0;

    for (unsigned i = 0; i < n_samples; ++i)
    {
        double left_sample = left[i];
        double right_sample = right[i];

        product += left_sample * right_sample;
        lsq += left_sample * left_sample;
        rsq += right_sample * right_sample;
    }

    denom = lsq * rsq;

    if (denom > 0.0)
    {
        return product / sqrt(denom);
    }
    else
    {
        return 1.0;
    }
}


float
balance(float* left, float* right, unsigned n_samples)
{
    float r = 0.0;
    float l = 0.0;
    float rbuf[n_samples];
    float lbuf[n_samples];
    VectorPower(rbuf, right, 2.0, n_samples);
    VectorPower(lbuf, left, 2.0, n_samples);
    r = VectorSum(rbuf, n_samples);
    l = VectorSum(lbuf, n_samples);
    return  (r - l) / ((r + l) + FLT_MIN);
}


double
balanceD(double* left, double* right, unsigned n_samples)
{
    double r = 0.0;
    double l = 0.0;
    double rbuf[n_samples];
    double lbuf[n_samples];
    VectorPowerD(rbuf, right, 2.0, n_samples);
    VectorPowerD(lbuf, left, 2.0, n_samples);
    r = VectorSumD(rbuf, n_samples);
    l = VectorSumD(lbuf, n_samples);
    return  (r - l) / ((r + l) + DBL_MIN);
}


float
vu_peak(float* signal, unsigned n_samples, MeterScale scale)
{
    float peak = VectorMax(signal, n_samples);
    return 20.0 * log10f(peak / ref[scale]);
}

double
vu_peakD(double* signal, unsigned n_samples, MeterScale scale)
{
    {
        double peak = VectorMaxD(signal, n_samples);
        return 20 * log10(peak / refD[scale]);
    }
}



/*******************************************************************************
 MultibandFilterInit */

MultibandFilter*
MultibandFilterInit(float   lowCutoff,
                    float   highCutoff,
                    float  sampleRate)
{
    MultibandFilter* filter = (MultibandFilter*) malloc(sizeof(MultibandFilter));
    filter->lowCutoff = lowCutoff;
    filter->highCutoff = highCutoff;
    filter->sampleRate = sampleRate;
    filter->LPA = LRFilterInit(LOWPASS, filter->lowCutoff, FILT_Q, filter->sampleRate);
    filter->HPA = LRFilterInit(HIGHPASS, filter->lowCutoff, FILT_Q, filter->sampleRate);
    filter->LPB = LRFilterInit(LOWPASS, filter->highCutoff, FILT_Q, filter->sampleRate);
    filter->HPB = LRFilterInit(HIGHPASS, filter->highCutoff, FILT_Q, filter->sampleRate);
    filter->APF = RBJFilterInit(ALLPASS, filter->sampleRate/2.0, filter->sampleRate);
    RBJFilterSetQ(filter->APF, 0.5);
    return filter;
}


MultibandFilterD*
MultibandFilterInitD(double lowCutoff,
                     double highCutoff,
                     double sampleRate)
{
    MultibandFilterD* filter = (MultibandFilterD*) malloc(sizeof(MultibandFilterD));
    filter->lowCutoff = lowCutoff;
    filter->highCutoff = highCutoff;
    filter->sampleRate = sampleRate;
    filter->LPA = LRFilterInitD(LOWPASS, filter->lowCutoff, FILT_Q, filter->sampleRate);
    filter->HPA = LRFilterInitD(HIGHPASS, filter->lowCutoff, FILT_Q, filter->sampleRate);
    filter->LPB = LRFilterInitD(LOWPASS, filter->highCutoff, FILT_Q, filter->sampleRate);
    filter->HPB = LRFilterInitD(HIGHPASS, filter->highCutoff, FILT_Q, filter->sampleRate);
    filter->APF = RBJFilterInitD(ALLPASS, filter->sampleRate/2.0, filter->sampleRate);
    RBJFilterSetQD(filter->APF, 0.5);
    return filter;
}


/*******************************************************************************
 MultibandFilterFree */

Error_t
MultibandFilterFree(MultibandFilter* filter)
{
    LRFilterFree(filter->LPA);
    LRFilterFree(filter->LPB);
    LRFilterFree(filter->HPA);
    LRFilterFree(filter->HPB);
    RBJFilterFree(filter->APF);
    if (filter)
    {
        free(filter);
        filter = NULL;
    }

    return NOERR;
}

Error_t
MultibandFilterFreeD(MultibandFilterD* filter)
{
    LRFilterFreeD(filter->LPA);
    LRFilterFreeD(filter->LPB);
    LRFilterFreeD(filter->HPA);
    LRFilterFreeD(filter->HPB);
    RBJFilterFreeD(filter->APF);
    if (filter)
    {
        free(filter);
        filter = NULL;
    }

    return NOERR;
}


/*******************************************************************************
 MultibandFilterFlush */

Error_t
MultibandFilterFlush(MultibandFilter* filter)
{
    LRFilterFlush(filter->LPA);
    LRFilterFlush(filter->LPB);
    LRFilterFlush(filter->HPA);
    LRFilterFlush(filter->HPB);
    RBJFilterFlush(filter->APF);

    return NOERR;
}

Error_t
MultibandFilterFlushD(MultibandFilterD* filter)
{
    LRFilterFlushD(filter->LPA);
    LRFilterFlushD(filter->LPB);
    LRFilterFlushD(filter->HPA);
    LRFilterFlushD(filter->HPB);
    RBJFilterFlushD(filter->APF);

    return NOERR;
}


/*******************************************************************************
 MultibandFilterSetLowCutoff */
Error_t
MultibandFilterSetLowCutoff(MultibandFilter* filter, float lowCutoff)
{
    filter->lowCutoff = lowCutoff;
    LRFilterSetParams(filter->LPA, LOWPASS, lowCutoff, FILT_Q);
    LRFilterSetParams(filter->HPA, HIGHPASS, lowCutoff, FILT_Q);
    return NOERR;
}

Error_t
MultibandFilterSetLowCutoffD(MultibandFilterD* filter, double lowCutoff)
{
    filter->lowCutoff = lowCutoff;
    LRFilterSetParamsD(filter->LPA, LOWPASS, lowCutoff, FILT_Q);
    LRFilterSetParamsD(filter->HPA, HIGHPASS, lowCutoff, FILT_Q);
    return NOERR;
}

/*******************************************************************************
 MultibandFilterSetHighCutoff */

Error_t
MultibandFilterSetHighCutoff(MultibandFilter* filter, float highCutoff)
{
    filter->highCutoff = highCutoff;
    LRFilterSetParams(filter->LPB, LOWPASS, highCutoff, FILT_Q);
    LRFilterSetParams(filter->HPB, HIGHPASS, highCutoff, FILT_Q);
    return NOERR;
}


Error_t
MultibandFilterSetHighCutoffD(MultibandFilterD* filter, double highCutoff)
{
    filter->highCutoff = highCutoff;
    LRFilterSetParamsD(filter->LPB, LOWPASS, highCutoff, FILT_Q);
    LRFilterSetParamsD(filter->HPB, HIGHPASS, highCutoff, FILT_Q);
    return NOERR;
}


/*******************************************************************************
 MultibandFilterUpdate */

Error_t
MultibandFilterUpdate(MultibandFilter*  filter,
                      float             lowCutoff,
                      float             highCutoff)
{
    filter->lowCutoff = lowCutoff;
    filter->highCutoff = highCutoff;
    LRFilterSetParams(filter->LPA, LOWPASS, lowCutoff, FILT_Q);
    LRFilterSetParams(filter->HPA, HIGHPASS, lowCutoff, FILT_Q);
    LRFilterSetParams(filter->LPB, LOWPASS, highCutoff, FILT_Q);
    LRFilterSetParams(filter->HPB, HIGHPASS, highCutoff, FILT_Q);
    return NOERR;
}


Error_t
MultibandFilterUpdateD(MultibandFilterD*    filter,
                       double               lowCutoff,
                       double               highCutoff)
{
    filter->lowCutoff = lowCutoff;
    filter->highCutoff = highCutoff;
    LRFilterSetParamsD(filter->LPA, LOWPASS, lowCutoff, FILT_Q);
    LRFilterSetParamsD(filter->HPA, HIGHPASS, lowCutoff, FILT_Q);
    LRFilterSetParamsD(filter->LPB, LOWPASS, highCutoff, FILT_Q);
    LRFilterSetParamsD(filter->HPB, HIGHPASS, highCutoff, FILT_Q);
    return NOERR;
}



/*******************************************************************************
 MultibandFilterProcess */

Error_t
MultibandFilterProcess(MultibandFilter* filter,
                           float*               lowOut,
                           float*               midOut,
                           float*               highOut,
                           const float*         inBuffer,
                           unsigned             n_samples)
{
    float tempLow[n_samples];
    float tempHi[n_samples];

    LRFilterProcess(filter->LPA, tempLow, inBuffer, n_samples);
    LRFilterProcess(filter->HPA, tempHi, inBuffer, n_samples);

    RBJFilterProcess(filter->APF, lowOut, tempLow, n_samples);
    LRFilterProcess(filter->LPB, midOut, tempHi, n_samples);
    LRFilterProcess(filter->HPB, highOut, tempHi, n_samples);

    return NOERR;
}


Error_t
MultibandFilterProcessD(MultibandFilterD*   filter,
                        double*             lowOut,
                        double*             midOut,
                        double*             highOut,
                        const double*       inBuffer,
                        unsigned            n_samples)
{
    double tempLow[n_samples];
    double tempHi[n_samples];

    LRFilterProcessD(filter->LPA, tempLow, inBuffer, n_samples);
    LRFilterProcessD(filter->HPA, tempHi, inBuffer, n_samples);

    RBJFilterProcessD(filter->APF, lowOut, tempLow, n_samples);
    LRFilterProcessD(filter->LPB, midOut, tempHi, n_samples);
    LRFilterProcessD(filter->HPB, highOut, tempHi, n_samples);

    return NOERR;
}


/* OnePoleFilterInit ***************************************************/
OnePole*
OnePoleInit(float cutoff, float sampleRate, Filter_t type)
{
    OnePole *filter = (OnePole*)malloc(sizeof(OnePole));
    if (filter)
    {
        filter->a0 = 1;
        filter->b1 = 0;
        filter->y1 = 0;
        filter->type = type;
        filter->sampleRate = sampleRate;
        OnePoleSetCutoff(filter, cutoff);
    }

    return filter;
}

OnePoleD*
OnePoleInitD(double cutoff, double sampleRate, Filter_t type)
{
    OnePoleD *filter = (OnePoleD*)malloc(sizeof(OnePoleD));
    if (filter)
    {
        filter->a0 = 1;
        filter->b1 = 0;
        filter->y1 = 0;
        filter->type = type;
        filter->sampleRate = sampleRate;
        OnePoleSetCutoffD(filter, cutoff);
    }

    return filter;
}


OnePole*
OnePoleRawInit(float beta, float alpha)
{
  OnePole *filter = (OnePole*)malloc(sizeof(OnePole));
  if (filter)
  {
    filter->a0 = alpha;
    filter->b1 = beta;
    filter->y1 = 0.0;
    filter->type = LOWPASS;
    filter->sampleRate = 0;
  }
  return filter;
}

OnePoleD*
OnePoleRawInitD(double beta, double alpha)
{
  OnePoleD *filter = (OnePoleD*)malloc(sizeof(OnePoleD));
  if (filter)
  {
    filter->a0 = alpha;
    filter->b1 = beta;
    filter->y1 = 0;
    filter->type = LOWPASS;
    filter->sampleRate = 0;
  }
  return filter;
}



Error_t
OnePoleSetType(OnePole* filter, Filter_t type)
{
    if (filter && (type == LOWPASS || type == HIGHPASS))
    {
        filter->type = type;
        OnePoleSetCutoff(filter, filter->cutoff);
        return NOERR;
    }
    else if (filter)
    {
        return VALUE_ERROR;
    }
    else
    {
        return NULL_PTR_ERROR;
    }
}

Error_t
OnePoleSetTypeD(OnePoleD* filter, Filter_t type)
{
    if (filter && (type == LOWPASS || type == HIGHPASS))
    {
        filter->type = type;
        OnePoleSetCutoffD(filter, filter->cutoff);
        return NOERR;
    }
    else if (filter)
    {
        return VALUE_ERROR;
    }
    else
    {
        return NULL_PTR_ERROR;
    }
}


/* OnePoleFilterSetCutoff **********************************************/
Error_t
OnePoleSetCutoff(OnePole* filter, float cutoff)
{
    filter->cutoff = cutoff;
    if (filter->type == LOWPASS)
    {
        filter->b1 = expf(-2.0 * M_PI * (filter->cutoff / filter->sampleRate));
        filter->a0 = 1.0 - filter->b1;
    }
    else
    {
        filter->b1 = -expf(-2.0 * M_PI * (0.5 - (filter->cutoff / filter->sampleRate)));
        filter->a0 = 1.0 + filter->b1;
    }
    return NOERR;
}


Error_t
OnePoleSetCutoffD(OnePoleD* filter, double cutoff)
{
    filter->cutoff = cutoff;
    if (filter->type == LOWPASS)
    {
        filter->b1 = exp(-2.0 * M_PI * (filter->cutoff / filter->sampleRate));
        filter->a0 = 1.0 - filter->b1;
    }
    else
    {
        filter->b1 = -exp(-2.0 * M_PI * (0.5 - (filter->cutoff / filter->sampleRate)));
        filter->a0 = 1.0 + filter->b1;
    }
    return NOERR;
}

/* OnePoleFilterSetSampleRate **********************************************/
Error_t
OnePoleSetSampleRate(OnePole* filter, float sampleRate)
{
    filter->sampleRate = sampleRate;
    OnePoleSetCutoff(filter, filter->cutoff);
    return NOERR;
}


Error_t
OnePoleSetSampleRateD(OnePoleD* filter, double sampleRate)
{
    filter->sampleRate = sampleRate;
    OnePoleSetCutoffD(filter, filter->cutoff);
    return NOERR;
}

Error_t
OnePoleSetCoefficients(OnePole* filter, float* beta, float* alpha)
{
  filter->b1 = *beta;
  filter->a0 = *alpha;
  return NOERR;
}

Error_t
OnePoleSetCoefficientsD(OnePoleD* filter, double* beta, double* alpha)
{
  filter->b1 = *beta;
  filter->a0 = *alpha;
  return NOERR;
}





/* OnePoleFilterFree ***************************************************/
Error_t
OnePoleFree(OnePole *filter)
{
    if (filter)
    {
        free(filter);
        filter = NULL;
    }
    return NOERR;
}


Error_t
OnePoleFreeD(OnePoleD *filter)
{
    if (filter)
    {
        free(filter);
        filter = NULL;
    }
    return NOERR;
}


Error_t
OnePoleFlush(OnePole *filter)
{
  if (filter)
  {
    filter->y1 = 0.0;
  }
  return NOERR;
}

Error_t
OnePoleFlushD(OnePoleD *filter)
{
  if (filter)
  {
    filter->y1 = 0.0;
  }
  return NOERR;
}

/* OnePoleFilterProcess ************************************************/
Error_t
OnePoleProcess(OnePole*     filter,
                 float*         outBuffer,
                 const float*   inBuffer,
                 unsigned       n_samples)
{
	#pragma omp simd aligned(outBuffer,inBuffer)
    for (unsigned i = 0; i < n_samples; ++i)
    {
        outBuffer[i] = filter->y1 = inBuffer[i] * filter->a0 + filter->y1 * filter->b1;
    }
    return NOERR;
}

Error_t
OnePoleProcessD(OnePoleD*   filter,
                 double*        outBuffer,
                 const double*  inBuffer,
                 unsigned       n_samples)
{
	#pragma omp simd aligned(outBuffer,inBuffer)
    for (unsigned i = 0; i < n_samples; ++i)
    {
        outBuffer[i] = filter->y1 = inBuffer[i] * filter->a0 + filter->y1 * filter->b1;
    }
    return NOERR;
}


/* OnePoleFilterTick ***************************************************/
float
OnePoleTick(OnePole*    filter,
            float       inSample)
{
    return  filter->y1 = inSample * filter->a0 + filter->y1 * filter->b1;
}

double
OnePoleTickD(OnePoleD*  filter,
             double     inSample)
{
    return  filter->y1 = inSample * filter->a0 + filter->y1 * filter->b1;
}


/*******************************************************************************
 OnePoleAlpha */

float
OnePoleAlpha(OnePole* filter)
{
    return filter->a0;
}

double
OnePoleAlphaD(OnePoleD* filter)
{
    return filter->a0;
}

/*******************************************************************************
 OnePoleBeta */

float
OnePoleBeta(OnePole* filter)
{
    return filter->b1;
}

double
OnePoleBetaD(OnePoleD* filter)
{
    return filter->b1;
}


/* Utility Functions **********************************************************/

/* Scale a sample to the output curve given for the given optocoupler type
 */
static inline double
scale_sample(double sample, Opto_t opto_type)
{

    double out = 0.0;

    switch (opto_type)
    {
        case OPTO_LDR:
            out = 3.0e-6 / pow(sample+DBL_MIN, 0.7);
            out = out > 1.0? 1.0 : out;
            break;
        case OPTO_PHOTOTRANSISTOR:
            out = sample;
            break;
        default:
            break;
    }
    return out;
}


/* Calculate the turn-on times [in seconds] for the given optocoupler type with
 the specified delay value
 */
static inline double
calculate_on_time(double delay, Opto_t opto_type)
{
    /* Prevent Denormals */
    double time = DBL_MIN;

    double delay_sq = delay*delay;

    switch (opto_type)
    {
        case OPTO_LDR:
            time = 0.01595 * delay_sq + 0.02795 * delay + 1e-5;
            break;

        case OPTO_PHOTOTRANSISTOR:
            time = 0.01595 * delay_sq + 0.02795 * delay + 1e-5;
            break;

        default:
            break;
    }
    return time;
}


/* Calculate the turn-off times [in seconds] for the given optocoupler type with
 the specified delay value
 */
static inline double
calculate_off_time(double delay, Opto_t opto_type)
{
    /* Prevent Denormals */
    double time = DBL_MIN;

    switch (opto_type)
    {
        case OPTO_LDR:
            time = 1.5*powf(delay+FLT_MIN,3.5);
            break;

        case OPTO_PHOTOTRANSISTOR:
            time = 1.5*powf(delay+FLT_MIN,3.5);
            break;
        default:
            break;
    }
    return time;
}



/* OptoInit ***************************************************************/
Opto*
OptoInit(Opto_t opto_type, float delay, float sample_rate)
{
    // Create opto struct
	Opto* opto = (Opto*)malloc(sizeof(Opto));
    if (opto)
    {
        // Initialization
        opto->type = opto_type;
        opto->sample_rate = sample_rate;
        opto->previous = 0;
        opto->delta_sign = 1;
        OptoSetDelay(opto, delay);
        opto->lp = OnePoleInit(opto->on_cutoff, opto->sample_rate, LOWPASS);
    }
    return opto;
}


OptoD*
OptoInitD(Opto_t opto_type, double delay, double sample_rate)
{
    // Create opto struct
    OptoD* opto = (OptoD*)malloc(sizeof(OptoD));
    if (opto)
    {
        // Initialization
        opto->type = opto_type;
        opto->sample_rate = sample_rate;
        opto->previous = 0;
        opto->delta_sign = 1;
        OptoSetDelayD(opto, delay);
        opto->lp = OnePoleInitD(opto->on_cutoff, opto->sample_rate, LOWPASS);
    }
    return opto;
}

/* OptoFree ***************************************************************/
Error_t
OptoFree(Opto* optocoupler)
{
    if (optocoupler)
    {
        if (optocoupler->lp)
        {
            OnePoleFree(optocoupler->lp);
        }
        free(optocoupler);
    }
     return NOERR;
}


Error_t
OptoFreeD(OptoD* optocoupler)
{
    if (optocoupler)
    {
        if (optocoupler->lp)
        {
            OnePoleFreeD(optocoupler->lp);
        }
        free(optocoupler);
    }
    return NOERR;
}


/* OptoSetDelay ***********************************************************/
Error_t
OptoSetDelay(Opto* optocoupler, float delay)
{
    optocoupler->delay = delay;
    optocoupler->on_cutoff = 1.0/(float)calculate_on_time((double)delay, optocoupler->type);
    optocoupler->off_cutoff = 1.0/(float)calculate_off_time((double)delay, optocoupler->type);
    return NOERR;
}

Error_t
OptoSetDelayD(OptoD* optocoupler, double delay)
{
    optocoupler->delay = delay;
    optocoupler->on_cutoff = 1.0/calculate_on_time(delay, optocoupler->type);
    optocoupler->off_cutoff = 1.0/calculate_off_time(delay, optocoupler->type);
    return NOERR;
}

/* OptoProcess ************************************************************/
Error_t
OptoProcess(Opto*           optocoupler,
            float*          out_buffer,
            const float*    in_buffer,
            unsigned        n_samples)
{
	#pragma omp simd aligned(out_buffer,in_buffer)
    for (unsigned i = 0; i < n_samples; ++i)
    {
        out_buffer[i] = OptoTick(optocoupler, in_buffer[i]);
    }
    return NOERR;
}

Error_t
OptoProcessD(OptoD*         optocoupler,
             double*        out_buffer,
             const double*  in_buffer,
             unsigned       n_samples)
{
	#pragma omp simd aligned(out_buffer,in_buffer)
    for (unsigned i = 0; i < n_samples; ++i)
    {
        out_buffer[i] = OptoTickD(optocoupler, in_buffer[i]);
    }
    return NOERR;
}


/* OptoTick ***************************************************************/
float
OptoTick(Opto* opto, float in_sample)
{
    float out;
    char prev_delta;

    /* Check sign of dv/dt */
    prev_delta = opto->delta_sign;
    opto->delta_sign = (in_sample - opto->previous) >= 0 ? 1 : -1;

    /* Update lopwass if sign changed */
    if (opto->delta_sign != prev_delta)
    {
        if (opto->delta_sign == 1)
        {
            OnePoleSetCutoff(opto->lp, opto->on_cutoff);
        }
        else
        {
            OnePoleSetCutoff(opto->lp, opto->off_cutoff);
        }
    }

    /* Do Delay model */
    out = OnePoleTick(opto->lp, in_sample);
    opto->previous = out;
    out = (float)scale_sample((double)out, opto->type);

    /* spit out sample */
    return out;
}


double
OptoTickD(OptoD* opto, double in_sample)
{
    double out;
    char prev_delta;

    /* Check sign of dv/dt */
    prev_delta = opto->delta_sign;
    opto->delta_sign = (in_sample - opto->previous) >= 0 ? 1 : -1;

    /* Update lopwass if sign changed */
    if (opto->delta_sign != prev_delta)
    {
        if (opto->delta_sign == 1)
        {
            OnePoleSetCutoffD(opto->lp, opto->on_cutoff);
        }
        else
        {
            OnePoleSetCutoffD(opto->lp, opto->off_cutoff);
        }
    }

    /* Do Delay model */
    out = OnePoleTickD(opto->lp, in_sample);
    opto->previous = out;
    out = scale_sample(out, opto->type);

    /* spit out sample */
    return out;
}


Error_t
linear_pan(float control, float *l_gain, float *r_gain)
{
    *l_gain = 1.0f - control;
    *r_gain = control;
    return NOERR;
}

Error_t
linear_panD(double control, double *l_gain, double *r_gain)
{
    *l_gain = 1.0 - control;
    *r_gain = control;
    return NOERR;
}

Error_t
equal_power_3dB_pan(float control, float *l_gain, float *r_gain)
{
    *l_gain = sinf((1.0 - control) * M_PI_2);
    *r_gain = sinf(control * M_PI_2);
    return NOERR;
}

Error_t
equal_power_3dB_panD(double control, double *l_gain, double *r_gain)
{
    *l_gain = sin((1.0 - control) * M_PI_2);
    *r_gain = sin(control * M_PI_2);
    return NOERR;
}

Error_t
equal_power_6dB_pan(float control, float *l_gain, float *r_gain)
{
    *l_gain = powf(sinf((1.0 - control) * M_PI_2), 2.0);
    *r_gain = powf(sinf(control * M_PI_2), 2.0);
    return NOERR;
}

Error_t
equal_power_6dB_panD(double control, double *l_gain, double *r_gain)
{
    *l_gain = pow(sin((1.0 - control) * M_PI_2),2.0);
    *r_gain = pow(sin(control * M_PI_2), 2.0);
    return NOERR;
}




/*******************************************************************************
 PolySaturatorInit */
PolySaturator*
PolySaturatorInit(float n)
{
    PolySaturator* saturator = (PolySaturator*)malloc(sizeof(PolySaturator));
    if (saturator)
    {
        PolySaturatorSetN(saturator, n);
        return saturator;
    }
    else
    {
        return NULL;
    }
}

PolySaturatorD*
PolySaturatorInitD(double n)
{
    PolySaturatorD* saturator = (PolySaturatorD*)malloc(sizeof(PolySaturatorD));
    if (saturator)
    {
        PolySaturatorSetND(saturator, n);
        return saturator;
    }
    else
    {
        return NULL;
    }
}


/*******************************************************************************
 PolySaturatorFree */
Error_t
PolySaturatorFree(PolySaturator* saturator)
{
    if (saturator)
    {
        free(saturator);
    }
    return NOERR;
}

Error_t
PolySaturatorFreeD(PolySaturatorD* saturator)
{
    if (saturator)
    {
        free(saturator);
    }
    return NOERR;
}


/*******************************************************************************
 PolySaturatorSetN */
Error_t
PolySaturatorSetN(PolySaturator* saturator, float n)
{
    if (saturator)
    {
        saturator->a = powf(1./n, 1./n);
        saturator->b = (n + 1) / n;
        saturator->n = n;
        return NOERR;
    }
    else
    {
        return ERROR;
    }
}

Error_t
PolySaturatorSetND(PolySaturatorD* saturator, double n)
{
    if (saturator)
    {
        saturator->a = pow(1./n, 1./n);
        saturator->b = (n + 1) / n;
        saturator->n = n;
        return NOERR;
    }
    else
    {
        return ERROR;
    }
}


/*******************************************************************************
 PolySaturatorProcess */
Error_t
PolySaturatorProcess(PolySaturator*     saturator,
                     float*             out_buffer,
                     const float*       in_buffer,
                     unsigned           n_samples)
{
    float buf[n_samples];
    VectorScalarMultiply(buf, (float*)in_buffer, saturator->a, n_samples);
    VectorAbs(buf, buf, n_samples);
    VectorPower(buf, buf, saturator->n, n_samples);
    VectorScalarAdd(buf, buf, -saturator->b, n_samples);
    VectorNegate(buf, buf, n_samples);
    VectorVectorMultiply(out_buffer, in_buffer, buf, n_samples);
    return NOERR;
}

Error_t
PolySaturatorProcessD(PolySaturatorD*   saturator,
                      double*           out_buffer,
                      const double*     in_buffer,
                      unsigned          n_samples)
{
    double buf[n_samples];
    VectorScalarMultiplyD(buf, (double*)in_buffer, saturator->a, n_samples);
    VectorAbsD(buf, buf, n_samples);
    VectorPowerD(buf, buf, saturator->n, n_samples);
    VectorScalarAddD(buf, buf, -saturator->b, n_samples);
    VectorNegateD(buf, buf, n_samples);
    VectorVectorMultiplyD(out_buffer, in_buffer, buf, n_samples);
    return NOERR;
}

/*******************************************************************************
 PolySaturatorTick */
float
PolySaturatorTick(PolySaturator* saturator, float in_sample)
{
    return -(powf(fabsf(saturator->a * in_sample), saturator->n) - saturator->b) * in_sample;
}

double
PolySaturatorTickD(PolySaturatorD* saturator, double in_sample)
{
    return -(pow(fabs(saturator->a * in_sample), saturator->n) - saturator->b) * in_sample;
}


/* RBJFilterUpdate *****************************************************/

static Error_t
RBJFilterUpdate(RBJFilter* filter)
{
    filter->cosOmega = cos(filter->omega);
    filter->sinOmega = sin(filter->omega);

    switch (filter->type)
    {
    case LOWPASS:
        filter->alpha = filter->sinOmega / (2.0 * filter->Q);
        filter->b[0] = (1 - filter->cosOmega) / 2;
        filter->b[1] = 1 - filter->cosOmega;
        filter->b[2] = filter->b[0];
        filter->a[0] = 1 + filter->alpha;
        filter->a[1] = -2 * filter->cosOmega;
        filter->a[2] = 1 - filter->alpha;
        break;

    case HIGHPASS:
        filter->alpha = filter->sinOmega / (2.0 * filter->Q);
        filter->b[0] = (1 + filter->cosOmega) / 2;
        filter->b[1] = -(1 + filter->cosOmega);
        filter->b[2] = filter->b[0];
        filter->a[0] = 1 + filter->alpha;
        filter->a[1] = -2 * filter->cosOmega;
        filter->a[2] = 1 - filter->alpha;
        break;

    case BANDPASS:
        filter->alpha = filter->sinOmega * sinhf(logf(2.0) / 2.0 * \
            filter->Q * filter->omega/filter->sinOmega);
        filter->b[0] = filter->sinOmega / 2;
        filter->b[1] = 0;
        filter->b[2] = -filter->b[0];
        filter->a[0] = 1 + filter->alpha;
        filter->a[1] = -2 * filter->cosOmega;
        filter->a[2] = 1 - filter->alpha;
        break;

    case ALLPASS:
        filter->alpha = filter->sinOmega / (2.0 * filter->Q);
        filter->b[0] = 1 - filter->alpha;
        filter->b[1] = -2 * filter->cosOmega;
        filter->b[2] = 1 + filter->alpha;
        filter->a[0] = filter->b[2];
        filter->a[1] = filter->b[1];
        filter->a[2] = filter->b[0];
        break;

    case NOTCH:
        filter->alpha = filter->sinOmega * sinhf(logf(2.0) / 2.0 * \
            filter->Q * filter->omega/filter->sinOmega);
        filter->b[0] = 1;
        filter->b[1] = -2 * filter->cosOmega;
        filter->b[2] = 1;
        filter->a[0] = 1 + filter->alpha;
        filter->a[1] = filter->b[1];
        filter->a[2] = 1 - filter->alpha;
        break;

    case PEAK:
        filter->alpha = filter->sinOmega * sinhf(logf(2.0) / 2.0 * \
            filter->Q * filter->omega/filter->sinOmega);
        filter->b[0] = 1 + (filter->alpha * filter->A);
        filter->b[1] = -2 * filter->cosOmega;
        filter->b[2] = 1 - (filter->alpha * filter->A);
        filter->a[0] = 1 + (filter->alpha / filter->A);
        filter->a[1] = filter->b[1];
        filter->a[2] = 1 - (filter->alpha / filter->A);
        break;

    case LOW_SHELF:
        filter->alpha = filter->sinOmega / 2.0 * sqrt( (filter->A + 1.0 / \
            filter->A) * (1.0 / filter->Q - 1.0) + 2.0);
        filter->b[0] = filter->A * ((filter->A + 1) - ((filter->A - 1) *       \
            filter->cosOmega) + (2 * sqrtf(filter->A) * filter->alpha));
        filter->b[1] = 2 * filter->A * ((filter->A - 1) - ((filter->A + 1) *   \
            filter->cosOmega));
        filter->b[2] = filter->A * ((filter->A + 1) - ((filter->A - 1) *       \
            filter->cosOmega) - (2 * sqrtf(filter->A) * filter->alpha));
        filter->a[0] = ((filter->A + 1) + ((filter->A - 1) *                   \
            filter->cosOmega) + (2 * sqrtf(filter->A) * filter->alpha));
        filter->a[1] = -2 * ((filter->A - 1) + ((filter->A + 1) *              \
            filter->cosOmega));
        filter->a[2] = ((filter->A + 1) + ((filter->A - 1) *                   \
            filter->cosOmega) - (2 * sqrtf(filter->A) * filter->alpha));
        break;

    case HIGH_SHELF:
        filter->alpha = filter->sinOmega / 2.0 * sqrt( (filter->A + 1.0 / \
            filter->A) * (1.0 / filter->Q - 1.0) + 2.0);
        filter->b[0] = filter->A * ((filter->A + 1) + ((filter->A - 1) *       \
            filter->cosOmega) + (2 * sqrtf(filter->A) * filter->alpha));
        filter->b[1] = -2 * filter->A * ((filter->A - 1) + ((filter->A + 1) *  \
            filter->cosOmega));
        filter->b[2] = filter->A * ((filter->A + 1) + ((filter->A - 1) *       \
            filter->cosOmega) - (2 * sqrtf(filter->A) * filter->alpha));
        filter->a[0] = ((filter->A + 1) - ((filter->A - 1) *                   \
            filter->cosOmega) + (2 * sqrtf(filter->A) * filter->alpha));
        filter->a[1] = 2 * ((filter->A - 1) - ((filter->A + 1) *               \
            filter->cosOmega));
        filter->a[2] = ((filter->A + 1) - ((filter->A - 1) *                   \
            filter->cosOmega) - (2 * sqrtf(filter->A) * filter->alpha));
        break;

    default:
        return ERROR;
        break;
    }

    // Normalize filter coefficients
    float factor = 1.0 / filter->a[0];
    float norm_a[2];
    float norm_b[3];
    VectorScalarMultiply(norm_a, &filter->a[1], factor, 2);
    VectorScalarMultiply(norm_b, filter->b, factor, 3);
    BiquadFilterUpdateKernel(filter->biquad, norm_b, norm_a);
    return NOERR;
}


static Error_t
RBJFilterUpdateD(RBJFilterD* filter)
{
    filter->cosOmega = cos(filter->omega);
    filter->sinOmega = sin(filter->omega);

    switch (filter->type)
    {
        case LOWPASS:
            filter->alpha = filter->sinOmega / (2.0 * filter->Q);
            filter->b[0] = (1 - filter->cosOmega) / 2;
            filter->b[1] = 1 - filter->cosOmega;
            filter->b[2] = filter->b[0];
            filter->a[0] = 1 + filter->alpha;
            filter->a[1] = -2 * filter->cosOmega;
            filter->a[2] = 1 - filter->alpha;
            break;

        case HIGHPASS:
            filter->alpha = filter->sinOmega / (2.0 * filter->Q);
            filter->b[0] = (1 + filter->cosOmega) / 2;
            filter->b[1] = -(1 + filter->cosOmega);
            filter->b[2] = filter->b[0];
            filter->a[0] = 1 + filter->alpha;
            filter->a[1] = -2 * filter->cosOmega;
            filter->a[2] = 1 - filter->alpha;
            break;

        case BANDPASS:
            filter->alpha = filter->sinOmega * sinh(logf(2.0) / 2.0 * \
                                                     filter->Q * filter->omega/filter->sinOmega);
            filter->b[0] = filter->sinOmega / 2;
            filter->b[1] = 0;
            filter->b[2] = -filter->b[0];
            filter->a[0] = 1 + filter->alpha;
            filter->a[1] = -2 * filter->cosOmega;
            filter->a[2] = 1 - filter->alpha;
            break;

        case ALLPASS:
            filter->alpha = filter->sinOmega / (2.0 * filter->Q);
            filter->b[0] = 1 - filter->alpha;
            filter->b[1] = -2 * filter->cosOmega;
            filter->b[2] = 1 + filter->alpha;
            filter->a[0] = filter->b[2];
            filter->a[1] = filter->b[1];
            filter->a[2] = filter->b[0];
            break;

        case NOTCH:
            filter->alpha = filter->sinOmega * sinh(logf(2.0) / 2.0 * \
                                                     filter->Q * filter->omega/filter->sinOmega);
            filter->b[0] = 1;
            filter->b[1] = -2 * filter->cosOmega;
            filter->b[2] = 1;
            filter->a[0] = 1 + filter->alpha;
            filter->a[1] = filter->b[1];
            filter->a[2] = 1 - filter->alpha;
            break;

        case PEAK:
            filter->alpha = filter->sinOmega * sinh(logf(2.0) / 2.0 * \
                                                     filter->Q * filter->omega/filter->sinOmega);
            filter->b[0] = 1 + (filter->alpha * filter->A);
            filter->b[1] = -2 * filter->cosOmega;
            filter->b[2] = 1 - (filter->alpha * filter->A);
            filter->a[0] = 1 + (filter->alpha / filter->A);
            filter->a[1] = filter->b[1];
            filter->a[2] = 1 - (filter->alpha / filter->A);
            break;

        case LOW_SHELF:
            filter->alpha = filter->sinOmega / 2.0 * sqrt( (filter->A + 1.0 / \
                                                            filter->A) * (1.0 / filter->Q - 1.0) + 2.0);
            filter->b[0] = filter->A * ((filter->A + 1) - ((filter->A - 1) *       \
                                                           filter->cosOmega) + (2 * sqrtf(filter->A) * filter->alpha));
            filter->b[1] = 2 * filter->A * ((filter->A - 1) - ((filter->A + 1) *   \
                                                               filter->cosOmega));
            filter->b[2] = filter->A * ((filter->A + 1) - ((filter->A - 1) *       \
                                                           filter->cosOmega) - (2 * sqrtf(filter->A) * filter->alpha));
            filter->a[0] = ((filter->A + 1) + ((filter->A - 1) *                   \
                                               filter->cosOmega) + (2 * sqrtf(filter->A) * filter->alpha));
            filter->a[1] = -2 * ((filter->A - 1) + ((filter->A + 1) *              \
                                                    filter->cosOmega));
            filter->a[2] = ((filter->A + 1) + ((filter->A - 1) *                   \
                                               filter->cosOmega) - (2 * sqrtf(filter->A) * filter->alpha));
            break;

        case HIGH_SHELF:
            filter->alpha = filter->sinOmega / 2.0 * sqrt( (filter->A + 1.0 / \
                                                            filter->A) * (1.0 / filter->Q - 1.0) + 2.0);
            filter->b[0] = filter->A * ((filter->A + 1) + ((filter->A - 1) *       \
                                                           filter->cosOmega) + (2 * sqrtf(filter->A) * filter->alpha));
            filter->b[1] = -2 * filter->A * ((filter->A - 1) + ((filter->A + 1) *  \
                                                                filter->cosOmega));
            filter->b[2] = filter->A * ((filter->A + 1) + ((filter->A - 1) *       \
                                                           filter->cosOmega) - (2 * sqrtf(filter->A) * filter->alpha));
            filter->a[0] = ((filter->A + 1) - ((filter->A - 1) *                   \
                                               filter->cosOmega) + (2 * sqrtf(filter->A) * filter->alpha));
            filter->a[1] = 2 * ((filter->A - 1) - ((filter->A + 1) *               \
                                                   filter->cosOmega));
            filter->a[2] = ((filter->A + 1) - ((filter->A - 1) *                   \
                                               filter->cosOmega) - (2 * sqrtf(filter->A) * filter->alpha));
            break;

        default:
            return ERROR;
            break;
    }

    // Normalize filter coefficients
    double factor = 1.0 / filter->a[0];
    double norm_a[2];
    double norm_b[3];
    VectorScalarMultiplyD(norm_a, &filter->a[1], factor, 2);
    VectorScalarMultiplyD(norm_b, filter->b, factor, 3);
    BiquadFilterUpdateKernelD(filter->biquad, norm_b, norm_a);
    return NOERR;
}


/* RBJFilterInit **********************************************************/
RBJFilter*
RBJFilterInit(Filter_t type, float cutoff, float sampleRate)
{
    // Create the filter
    RBJFilter* filter = (RBJFilter*)malloc(sizeof(RBJFilter));

    if (filter)
    {
        // Initialization
        filter->type = type;
        filter->omega =  HZ_TO_RAD(cutoff) / sampleRate; //hzToRadians(cutoff, sampleRate);
        filter->Q = 1;
        filter->A = 1;
        filter->dbGain = 0;
        filter->sampleRate = sampleRate;


        // Initialize biquad
        float b[3] = {0, 0, 0};
        float a[2] = {0, 0};
        filter->biquad = BiquadFilterInit(b,a);

        // Calculate coefficients
        RBJFilterUpdate(filter);
    }

    return filter;
}

RBJFilterD*
RBJFilterInitD(Filter_t type, double cutoff, double sampleRate)
{
    // Create the filter
    RBJFilterD* filter = (RBJFilterD*)malloc(sizeof(RBJFilterD));

    if (filter)
    {
        // Initialization
        filter->type = type;
        filter->omega =  HZ_TO_RAD(cutoff) / sampleRate;
        filter->Q = 1;
        filter->A = 1;
        filter->dbGain = 0;
        filter->sampleRate = sampleRate;


        // Initialize biquad
        double b[3] = {0, 0, 0};
        double a[2] = {0, 0};
        filter->biquad = BiquadFilterInitD(b, a);

        // Calculate coefficients
        RBJFilterUpdateD(filter);
    }

    return filter;
}


/* RBJFilterFree *******************************************************/
Error_t
RBJFilterFree(RBJFilter* filter)
{
    BiquadFilterFree(filter->biquad);
    if (filter)
    {
        free(filter);
        filter = NULL;
    }
    return NOERR;
}

Error_t
RBJFilterFreeD(RBJFilterD* filter)
{
    BiquadFilterFreeD(filter->biquad);
    if (filter)
    {
        free(filter);
        filter = NULL;
    }
    return NOERR;
}

/* RBJFilterSetType ****************************************************/
Error_t
RBJFilterSetType(RBJFilter* filter, Filter_t type)
{
    filter->type = type;
    RBJFilterUpdate(filter);
    return NOERR;
}


Error_t
RBJFilterSetTypeD(RBJFilterD* filter, Filter_t type)
{
    filter->type = type;
    RBJFilterUpdateD(filter);
    return NOERR;
}

/* RBJFilterSetCutoff **************************************************/
Error_t
RBJFilterSetCutoff(RBJFilter* filter, float cutoff)
{
    filter->omega = HZ_TO_RAD(cutoff) / filter->sampleRate;
    RBJFilterUpdate(filter);
    return NOERR;
}


Error_t
RBJFilterSetCutoffD(RBJFilterD* filter, double cutoff)
{
    filter->omega = HZ_TO_RAD(cutoff) / filter->sampleRate;
    RBJFilterUpdateD(filter);
    return NOERR;
}


/* RBJFilterSetQ *******************************************************/
Error_t
RBJFilterSetQ(RBJFilter* filter, float Q)
{
    filter->Q = Q;
    RBJFilterUpdate(filter);
    return NOERR;
}

Error_t
RBJFilterSetQD(RBJFilterD* filter, double Q)
{
    filter->Q = Q;
    RBJFilterUpdateD(filter);
    return NOERR;
}

/* RBJFilterSetParams **************************************************/
Error_t
RBJFilterSetParams(RBJFilter*   filter,
                   Filter_t     type,
                   float        cutoff,
                   float        Q)
{
    filter->type = type;
    filter->omega = HZ_TO_RAD(cutoff) / filter->sampleRate;
    filter->Q = Q;
    RBJFilterUpdate(filter);
    return NOERR;
}

Error_t
RBJFilterSetParamsD(RBJFilterD* filter,
                    Filter_t    type,
                    double      cutoff,
                    double      Q)
{
    filter->type = type;
    filter->omega = HZ_TO_RAD(cutoff) / filter->sampleRate;
    filter->Q = Q;
    RBJFilterUpdateD(filter);
    return NOERR;
}

/* RBJFilterProcess ****************************************************/
Error_t
RBJFilterProcess(RBJFilter*     filter,
                        float*              outBuffer,
                        const float*        inBuffer,
                        unsigned            n_samples)
{
    BiquadFilterProcess(filter->biquad,outBuffer,inBuffer,n_samples);
    return NOERR;
}

Error_t
RBJFilterProcessD(RBJFilterD*   filter,
                  double*       outBuffer,
                  const double* inBuffer,
                  unsigned      n_samples)
{
    BiquadFilterProcessD(filter->biquad,outBuffer,inBuffer,n_samples);
    return NOERR;
}

/* RBJFilterFlush ******************************************************/
Error_t
RBJFilterFlush(RBJFilter* filter)
{
    BiquadFilterFlush(filter->biquad);
    return NOERR;
}

Error_t
RBJFilterFlushD(RBJFilterD* filter)
{
    BiquadFilterFlushD(filter->biquad);
    return NOERR;
}



/*******************************************************************************
 RMSEstimatorInit */
RMSEstimator*
RMSEstimatorInit(float avgTime, float sampleRate)
{
    RMSEstimator* rms = (RMSEstimator*) malloc(sizeof(RMSEstimator));
    rms->avgTime = avgTime;
    rms->sampleRate = sampleRate;
    rms->RMS = 1;
    rms->avgCoeff = 0.5 * (1.0 - expf( -1.0 / (rms->sampleRate * rms->avgTime)));

    return rms;
}


RMSEstimatorD*
RMSEstimatorInitD(double avgTime, double sampleRate)
{
    RMSEstimatorD* rms = (RMSEstimatorD*) malloc(sizeof(RMSEstimatorD));
    rms->avgTime = avgTime;
    rms->sampleRate = sampleRate;
    rms->RMS = 1;
    rms->avgCoeff = 0.5 * (1.0 - expf( -1.0 / (rms->sampleRate * rms->avgTime)));

    return rms;
}

/*******************************************************************************
 RMSEstimatorFree */
Error_t
RMSEstimatorFree(RMSEstimator* rms)
{
    if (rms)
    {
        free(rms);
        rms = NULL;
    }
    return NOERR;
}

Error_t
RMSEstimatorFreeD(RMSEstimatorD* rms)
{
    if (rms)
    {
        free(rms);
        rms = NULL;
    }
    return NOERR;
}


/*******************************************************************************
 RMSEstimatorFlush */
Error_t
RMSEstimatorFlush(RMSEstimator* rms)
{
    if (rms)
    {
        rms->RMS = 1.0;
        return NOERR;
    }
    return NULL_PTR_ERROR;
}

Error_t
RMSEstimatorFlushD(RMSEstimatorD* rms)
{
    if (rms)
    {
        rms->RMS = 1.0;
        return NOERR;
    }
    return NULL_PTR_ERROR;
}


/*******************************************************************************
 RMSEstimatorSetAvgTime */
Error_t
RMSEstimatorSetAvgTime(RMSEstimator* rms, float avgTime)
{
    rms->avgTime = avgTime;
    rms->avgCoeff = 0.5 * (1.0 - expf( -1.0 / (rms->sampleRate * rms->avgTime)));
    return NOERR;
}

Error_t
RMSEstimatorSetAvgTimeD(RMSEstimatorD* rms, double avgTime)
{
    rms->avgTime = avgTime;
    rms->avgCoeff = 0.5 * (1.0 - expf( -1.0 / (rms->sampleRate * rms->avgTime)));
    return NOERR;
}

/*******************************************************************************
 RMSEstimatorProcess */
Error_t
RMSEstimatorProcess(RMSEstimator*   rms,
                        float*              outBuffer,
                        const float*        inBuffer,
                        unsigned            n_samples)
{
    for (unsigned i = 0; i < n_samples; ++i)
    {
        rms->RMS += rms->avgCoeff * ((f_abs(inBuffer[i])/rms->RMS) - rms->RMS);
        outBuffer[i] = rms->RMS;
    }
    return NOERR;
}

Error_t
RMSEstimatorProcessD(RMSEstimatorD* rms,
                     double*        outBuffer,
                     const double*  inBuffer,
                     unsigned       n_samples)
{
    for (unsigned i = 0; i < n_samples; ++i)
    {
        rms->RMS += rms->avgCoeff * ((f_abs(inBuffer[i])/rms->RMS) - rms->RMS);
        outBuffer[i] = rms->RMS;
    }
    return NOERR;
}

/*******************************************************************************
 RMSEstimatorTick */
float
RMSEstimatorTick(RMSEstimator*  rms,
                     float              inSample)
{
    rms->RMS += rms->avgCoeff * ((f_abs(inSample/rms->RMS)) - rms->RMS);
    return rms->RMS;
}

double
RMSEstimatorTickD(RMSEstimatorD* rms, double inSample)
{
    rms->RMS += rms->avgCoeff * ((f_abs(inSample/rms->RMS)) - rms->RMS);
    return rms->RMS;
}


SpectrumAnalyzer*
SpectrumAnalyzerInit(unsigned fft_length, float sample_rate)
{
    SpectrumAnalyzer* inst = (SpectrumAnalyzer*) malloc(sizeof(SpectrumAnalyzer));
    if (NULL != inst)
    {
        FFTConfig* fft = FFTInit(fft_length);
        WindowFunction* window = WindowFunctionInit(fft_length, BLACKMAN);
        float* frequencies = (float*)malloc((fft_length / 2) * sizeof(float));
        float* real = (float*)malloc((fft_length / 2) * sizeof(float));
        float* imag = (float*)malloc((fft_length / 2) * sizeof(float));
        float* mag = (float*)malloc((fft_length / 2) * sizeof(float));
        float* phase = (float*)malloc((fft_length / 2) * sizeof(float));
        float* root_moment = (float*)malloc((fft_length / 2) * sizeof(float));
        if ((NULL != window) && (NULL != fft) && (NULL != frequencies) && (NULL != real) \
            && (NULL != imag) && (NULL != mag) && (NULL != phase) && (NULL != root_moment))
        {
            inst->fft_length = fft_length;
            inst->bins = fft_length / 2;
            inst->sample_rate = sample_rate;
            inst->mag_sum = 0.0;
            inst->frequencies = frequencies;
            inst->real = real;
            inst->imag = imag;
            inst->mag = mag;
            inst->phase = phase;
            inst->root_moment = root_moment;
            inst->fft = fft;
            inst->window = window;
            inst->window_type = BLACKMAN;
            calculate_bin_frequencies(inst->frequencies, fft_length, sample_rate);
            *(inst->root_moment) = 0.0;
        }
        else
        {
            WindowFunctionFree(window);
            FFTFree(fft);
            if (NULL != root_moment)
            {
                free(root_moment);
            }
            if (NULL != phase)
            {
                free(phase);
            }
            if (NULL != mag)
            {
                free(mag);
            }
            if (NULL != imag)
            {
                free(imag);
            }
            if (NULL != real)
            {
                free(real);
            }
            if (NULL != frequencies)
            {
                free(real);
            }
            inst = NULL;
        }
    }
    return inst;
}

SpectrumAnalyzerD*
SpectrumAnalyzerInitD(unsigned fft_length, double sample_rate)
{
    SpectrumAnalyzerD* inst = (SpectrumAnalyzerD*) malloc(sizeof(SpectrumAnalyzerD));
    if (NULL != inst)
    {
        FFTConfigD* fft = FFTInitD(fft_length);
        WindowFunctionD* window = WindowFunctionInitD(fft_length, BLACKMAN);
        double* frequencies = (double*)malloc((fft_length / 2) * sizeof(double));
        double* real = (double*)malloc((fft_length / 2) * sizeof(double));
        double* imag = (double*)malloc((fft_length / 2) * sizeof(double));
        double* mag = (double*)malloc((fft_length / 2) * sizeof(double));
        double* phase = (double*)malloc((fft_length / 2) * sizeof(double));
        double* root_moment = (double*)malloc((fft_length / 2) * sizeof(double));
        if ((NULL != window) && (NULL != fft) && (NULL != frequencies) && (NULL != real) \
            && (NULL != imag) && (NULL != mag) && (NULL != phase) && (NULL != root_moment))
        {
            inst->fft_length = fft_length;
            inst->bins = fft_length/2;
            inst->sample_rate = sample_rate;
            inst->mag_sum = 0.0;
            inst->frequencies = frequencies;
            inst->real = real;
            inst->imag = imag;
            inst->mag = mag;
            inst->phase = phase;
            inst->root_moment = root_moment;
            inst->fft = fft;
            inst->window = window;
            inst->window_type = BLACKMAN;
            calculate_bin_frequenciesD(inst->frequencies, fft_length, sample_rate);
            *(inst->root_moment) = 0.0;
        }
        else
        {
            WindowFunctionFreeD(window);
            FFTFreeD(fft);
            if (NULL != root_moment)
            {
                free(root_moment);
            }
            if (NULL != phase)
            {
                free(phase);
            }
            if (NULL != mag)
            {
                free(mag);
            }
            if (NULL != imag)
            {
                free(imag);
            }
            if (NULL != real)
            {
                free(real);
            }
            if (NULL != frequencies)
            {
                free(frequencies);
            }
            inst = NULL;
        }
    }
    return inst;
}


void SpectrumAnalyzerFree(SpectrumAnalyzer * analyzer)
{
    WindowFunctionFree(analyzer->window);
    FFTFree(analyzer->fft);
    if (NULL != analyzer->root_moment)
    {
        free(analyzer->root_moment);
    }
    if (NULL != analyzer->phase)
    {
        free(analyzer->phase);
    }
    if (NULL != analyzer->mag)
    {
        free(analyzer->mag);
    }
    if (NULL != analyzer->imag)
    {
        free(analyzer->imag);
    }
    if (NULL != analyzer->real)
    {
        free(analyzer->real);
    }
    if (NULL != analyzer->frequencies)
    {
        free(analyzer->frequencies);
    }
    free(analyzer);
}

void SpectrumAnalyzerFreeD(SpectrumAnalyzerD * analyzer)
{
    WindowFunctionFreeD(analyzer->window);
    FFTFreeD(analyzer->fft);
    if (NULL != analyzer->root_moment)
    {
        free(analyzer->root_moment);
    }
    if (NULL != analyzer->phase)
    {
        free(analyzer->phase);
    }
    if (NULL != analyzer->mag)
    {
        free(analyzer->mag);
    }
    if (NULL != analyzer->imag)
    {
        free(analyzer->imag);
    }
    if (NULL != analyzer->real)
    {
        free(analyzer->real);
    }
    if (NULL != analyzer->frequencies)
    {
        free(analyzer->frequencies);
    }
    free(analyzer);
}

void
SpectrumAnalyzerAnalyze(SpectrumAnalyzer* analyzer, float* signal)
{
    float scratch[analyzer->fft_length];
    WindowFunctionProcess(analyzer->window, scratch, signal, analyzer->fft_length);
    FFT_R2C(analyzer->fft, scratch, analyzer->real, analyzer->imag);
    VectorRectToPolar(analyzer->mag, analyzer->phase, analyzer->real, analyzer->imag, analyzer->bins);
    analyzer->mag_sum = VectorSum(analyzer->mag, analyzer->bins);
    analyzer->root_moment[0] = 0.0;
}

void
SpectrumAnalyzerAnalyzeD(SpectrumAnalyzerD* analyzer, double* signal)
{
    double scratch[analyzer->fft_length];
    WindowFunctionProcessD(analyzer->window, scratch, signal, analyzer->fft_length);
    FFT_R2CD(analyzer->fft, scratch, analyzer->real, analyzer->imag);
    VectorRectToPolarD(analyzer->mag, analyzer->phase, analyzer->real, analyzer->imag, analyzer->bins);
    analyzer->mag_sum = VectorSumD(analyzer->mag, analyzer->bins);
    analyzer->root_moment[0] = 0.0;
}

float
SpectralCentroid(SpectrumAnalyzer* analyzer)
{
    float num[analyzer->bins];
    VectorVectorMultiply(num, analyzer->mag, analyzer->frequencies, analyzer->bins);
    return VectorSum(num, analyzer->bins) / analyzer->mag_sum;
}

double
SpectralCentroidD(SpectrumAnalyzerD* analyzer)
{
    double num[analyzer->bins];
    VectorVectorMultiplyD(num, analyzer->mag, analyzer->frequencies, analyzer->bins);
    return VectorSumD(num, analyzer->bins) / analyzer->mag_sum;
}

float
SpectralSpread(SpectrumAnalyzer* analyzer)
{
    float mu = SpectralCentroid(analyzer);
    float num[analyzer->bins];
    if (analyzer->root_moment[0] == 0.0)
    {
        VectorScalarAdd(analyzer->root_moment, analyzer->frequencies, -mu, analyzer->bins);
    }
    VectorPower(num, analyzer->root_moment, 2, analyzer->bins);
    return VectorSum(num, analyzer->bins) / analyzer->mag_sum;
}

double
SpectralSpreadD(SpectrumAnalyzerD* analyzer)
{
    double mu = SpectralCentroidD(analyzer);
    double num[analyzer->bins];
    if (analyzer->root_moment[0] == 0.0)
    {
        VectorScalarAddD(analyzer->root_moment, analyzer->frequencies, -mu, analyzer->bins);
    }
    VectorPowerD(num, analyzer->root_moment, 2, analyzer->bins);
    return VectorSumD(num, analyzer->bins) / analyzer->mag_sum;
}

float
SpectralSkewness(SpectrumAnalyzer* analyzer)
{
    float mu = SpectralCentroid(analyzer);
    float num[analyzer->bins];
    if (analyzer->root_moment[0] == 0.0)
    {
        VectorScalarAdd(analyzer->root_moment, analyzer->frequencies, -mu, analyzer->bins);
    }
    VectorPower(num, analyzer->root_moment, 3, analyzer->bins);
    return VectorSum(num, analyzer->bins) / analyzer->mag_sum;
}

double
SpectralSkewnessD(SpectrumAnalyzerD* analyzer)
{
    double mu = SpectralCentroidD(analyzer);
    double num[analyzer->bins];
    if (analyzer->root_moment[0] == 0.0)
    {
        VectorScalarAddD(analyzer->root_moment, analyzer->frequencies, -mu, analyzer->bins);
    }
    VectorPowerD(num, analyzer->root_moment, 3, analyzer->bins);
    return VectorSumD(num, analyzer->bins) / analyzer->mag_sum;
}

float
SpectralKurtosis(SpectrumAnalyzer* analyzer)
{
    float mu = SpectralCentroid(analyzer);
    float num[analyzer->bins];
    if (analyzer->root_moment[0] == 0.0)
    {
        VectorScalarAdd(analyzer->root_moment, analyzer->frequencies, -mu, analyzer->bins);
    }
    VectorPower(num, analyzer->root_moment, 4, analyzer->bins);
    return VectorSum(num, analyzer->bins) / analyzer->mag_sum;
}

double
SpectralKurtosisD(SpectrumAnalyzerD* analyzer)
{
    double mu = SpectralCentroidD(analyzer);
    double num[analyzer->bins];
    if (analyzer->root_moment[0] == 0.0)
    {
        VectorScalarAddD(analyzer->root_moment, analyzer->frequencies, -mu, analyzer->bins);
    }
    VectorPowerD(num, analyzer->root_moment, 4, analyzer->bins);
    return VectorSumD(num, analyzer->bins) / analyzer->mag_sum;
}


/*******************************************************************************
 Static Funcdtion Definitions */
static void
calculate_bin_frequencies(float* dest, unsigned fft_length, float sample_rate)
{
    float freq_step = sample_rate / fft_length;
    for(unsigned bin = 0; bin < fft_length / 2; ++bin)
    {
        dest[bin] = bin * freq_step;
    }
}

static void
calculate_bin_frequenciesD(double* dest, unsigned fft_length, double sample_rate)
{
    double freq_step = sample_rate / fft_length;
    for(unsigned bin = 0; bin < fft_length / 2; ++bin)
    {
        dest[bin] = bin * freq_step;
    }
}



void
StereoToMono(float*        dest,
             const float*  left,
             const float*  right,
             unsigned      length)
{
  float scale = SQRT_TWO_OVER_TWO;
  VectorVectorSumScale(dest, left, right, &scale, length);
}


void
StereoToMonoD(double*       dest,
              const double* left,
              const double* right,
              unsigned      length)
{
  double scale = SQRT_TWO_OVER_TWO;
  VectorVectorSumScaleD(dest, left, right, &scale, length);
}



void
MonoToStereo(float*         left,
             float*         right,
             const float*   mono,
             unsigned       length)
{
  float scale = SQRT_TWO_OVER_TWO;
  VectorScalarMultiply(left, mono, scale, length);
  CopyBuffer(right, left, length);
}

void
MonoToStereoD(double*       left,
              double*       right,
              const double* mono,
              unsigned      length)
{
  double scale = SQRT_TWO_OVER_TWO;
  VectorScalarMultiplyD(left, mono, scale, length);
  CopyBufferD(right, left, length);
}



/*******************************************************************************
TapeInit */
Tape*
TapeInit(TapeSpeed speed, float saturation, float hysteresis, float flutter, float sample_rate)
{
    // Create TapeSaturator Struct
    Tape* tape = (Tape*)malloc(sizeof(Tape));
    PolySaturator* saturator = PolySaturatorInit(1);
    // Allocate for longest period...
    unsigned mod_length = (unsigned)(sample_rate / 0.80);
    float* mod = (float*)malloc(mod_length * sizeof(float));
    if (tape && saturator && mod)
    {
        // Initialization
        tape->polysat = saturator;
        tape->sample_rate = sample_rate;
        tape->pos_peak = 0.0;
        tape->neg_peak = 0.0;
        tape->flutter_mod = mod;
        tape->flutter_mod_length = mod_length;
        
        // Need these initialized here.
        tape->speed = speed;
        tape->saturation = saturation;
        
        // Set up
        TapeSetFlutter(tape, flutter);
        TapeSetSaturation(tape, saturation);
        TapeSetSpeed(tape, speed);
        TapeSetHysteresis(tape, hysteresis);
        return tape;
    }
    else
    {
        free(mod);
        free(tape);
        free(saturator);
        return NULL;
    }
}


/*******************************************************************************
 Tape Free */
Error_t
TapeFree(Tape* tape)
{
    if(tape)
        free(tape);
    tape = NULL;
    return NOERR;
}


/*******************************************************************************
 Set Speed */
Error_t
TapeSetSpeed(Tape* tape, TapeSpeed speed)
{
    
    if (tape)
    {
        // Set speed
        tape->speed = speed;
        
        // Update saturation curve
        PolySaturatorSetN(tape->polysat, calculate_n(tape->saturation, speed));
        
        // Clear old flutter/wow modulation waveform
        ClearBuffer(tape->flutter_mod, tape->flutter_mod_length); // Yes, clear the old length...
        
        // Calculate new modulation waveform length...
        tape->flutter_mod_length = (unsigned)(tape->sample_rate / \
                                              (0.80 * powf(2.0, (float)speed)));
        
        // Generate flutter/wow modulation waveform
        float temp_buffer[tape->flutter_mod_length];
        for (unsigned comp = 0; comp < N_FLUTTER_COMPONENTS; ++comp)
        {
            float phase_step = (2.0 * M_PI * comp * powf(2.0, (float)speed)) / tape->sample_rate;
            ClearBuffer(temp_buffer, tape->flutter_mod_length);
            for (unsigned i = 0; i < tape->flutter_mod_length; ++i)
            {
                temp_buffer[i] = sinf(i * phase_step) / N_FLUTTER_COMPONENTS;
            }
            VectorVectorAdd(tape->flutter_mod, tape->flutter_mod,
                            temp_buffer, tape->flutter_mod_length);
        }
        return NOERR;;
    }
    else
    {
        return NULL_PTR_ERROR;
    }
}

/*******************************************************************************
Set Saturation */
Error_t
TapeSetSaturation(Tape* tape, float saturation)
{

    if (tape)
    {
        float n = calculate_n(saturation, tape->speed);
        tape->saturation = saturation;
        return PolySaturatorSetN(tape->polysat, n);
    }
    else
    {
        return NULL_PTR_ERROR;
    }
}


/*******************************************************************************
 Set Hysteresis */
Error_t
TapeSetHysteresis(Tape* tape, float hysteresis)
{
    if (tape)
    {
        tape->hysteresis = hysteresis;
        return NOERR;
    }
    else
    {
        return NULL_PTR_ERROR;
    }
}


/*******************************************************************************
 Set Flutter */
Error_t
TapeSetFlutter(Tape* tape, float flutter)
{
    if (tape)
    {
        tape->flutter = flutter;
        return NOERR;
    }
    else
    {
        return NULL_PTR_ERROR;
    }
}

/*******************************************************************************
 Get Saturation */
float
TapeGetSaturation(Tape* tape)
{
    if (tape)
    {
        return tape->saturation;
    }
    else
    {
        return -1.0;
    }
}

/*******************************************************************************
 Get Hysteresis */
float
TapeGetHysteresis(Tape* tape)
{
    if (tape)
    {
        return tape->hysteresis;
    }
    else
    {
        return -1.0;
    }
}

/*******************************************************************************
 TapeProcess */
Error_t
TapeProcess(Tape*           tape,
            float*          out_buffer,
            const float*    in_buffer,
            unsigned        n_samples)
{
    return NOERR;
}



/*******************************************************************************
TapeTick */
float
TapeTick(Tape* tape, float in_sample)
{
    
    float hysteresis = tape->hysteresis * 0.05;
    float output = 0.0;
    if (in_sample >= 0)
    {
        tape->neg_peak = 0.0;
        if (in_sample > tape->pos_peak)
        {
            tape->pos_peak = in_sample;
            output = in_sample;
        }
        else if (in_sample > (1 - hysteresis) * tape->pos_peak)
        {
            output = tape->pos_peak;
        }
        else
        {
            output = in_sample + hysteresis * tape->pos_peak;
        }
    }
    else
    {
        tape->pos_peak = 0.0;
        if (in_sample < tape->neg_peak)
        {
            tape->neg_peak = in_sample;
            output = in_sample;
        }
        
        else if (in_sample < (1 - hysteresis) * tape->neg_peak)
        {
            output = tape->neg_peak;
        }
        
        else
        {
            output = in_sample + hysteresis * tape->neg_peak;
        }
    }
    return PolySaturatorTick(tape->polysat, output);
}



/* UpsamplerInit *******************************************************/
Upsampler*
UpsamplerInit(ResampleFactor_t factor)
{
    unsigned n_filters = 1;
    switch(factor)
    {
        case X2:
            n_filters = 2;
            break;
        case X4:
            n_filters = 4;
            break;
        case X8:
            n_filters = 8;
            break;
        /*
        case X16:
            n_filters = 16;
            break;
        */
        default:
            return NULL;
    }

    // Allocate memory for the upsampler
    Upsampler* upsampler = (Upsampler*)malloc(sizeof(Upsampler));

    // Allocate memory for the polyphase array
    FIRFilter** polyphase = (FIRFilter**)malloc(n_filters * sizeof(FIRFilter*));

    if (upsampler && polyphase)
    {
        upsampler->polyphase = polyphase;

        // Create polyphase filters
        unsigned idx;
        for(idx = 0; idx < n_filters; ++idx)
        {
            upsampler->polyphase[idx] = FIRFilterInit(PolyphaseCoeffs[factor][idx], 64, DIRECT);
        }

        // Add factor
        upsampler->factor = n_filters;

        return upsampler;
    }
    else
    {
        if (polyphase)
        {
            free(polyphase);
        }
        if (upsampler)
        {
            free(upsampler);
        }
        return NULL;
    }
}

UpsamplerD*
UpsamplerInitD(ResampleFactor_t factor)
{
    unsigned n_filters = 1;
    switch(factor)
    {
        case X2:
            n_filters = 2;
            break;
        case X4:
            n_filters = 4;
            break;
        case X8:
            n_filters = 8;
            break;
        /*
        case X16:
            n_filters = 16;
            break;
        */
        default:
            return NULL;
    }

    // Allocate memory for the upsampler
    UpsamplerD* upsampler = (UpsamplerD*)malloc(sizeof(UpsamplerD));

    // Allocate memory for the polyphase array
    FIRFilterD** polyphase = (FIRFilterD**)malloc(n_filters * sizeof(FIRFilterD*));

    if (upsampler && polyphase)
    {
        upsampler->polyphase = polyphase;

        // Create polyphase filters
        unsigned idx;
        for(idx = 0; idx < n_filters; ++idx)
        {
            upsampler->polyphase[idx] = FIRFilterInitD(PolyphaseCoeffsD[factor][idx], 64, DIRECT);
        }

        // Add factor
        upsampler->factor = n_filters;

        return upsampler;
    }
    else
    {
        if (polyphase)
        {
            free(polyphase);
        }
        if (upsampler)
        {
            free(upsampler);
        }
        return NULL;
    }
}

/* UpsamplerFree *******************************************************/
Error_t
UpsamplerFree(Upsampler* upsampler)
{
    if (upsampler)
    {
        if (upsampler->polyphase)
        {
            for (unsigned i = 0; i < upsampler->factor; ++i)
            {
                FIRFilterFree(upsampler->polyphase[i]);
            }
            free(upsampler->polyphase);
        }
        free(upsampler);
    }
    return NOERR;
}

Error_t
UpsamplerFreeD(UpsamplerD* upsampler)
{
    if (upsampler)
    {
        if (upsampler->polyphase)
        {
            for (unsigned i = 0; i < upsampler->factor; ++i)
            {
                FIRFilterFreeD(upsampler->polyphase[i]);
            }
            free(upsampler->polyphase);
        }
        free(upsampler);
    }
    return NOERR;
}


/* UpsamplerFlush ****************************************************/
Error_t
UpsamplerFlush(Upsampler* upsampler)
{
    unsigned idx;
    for (idx = 0; idx < upsampler->factor; ++idx)
    {
        FIRFilterFlush(upsampler->polyphase[idx]);
    }
    return NOERR;
}

Error_t
UpsamplerFlushD(UpsamplerD* upsampler)
{
    unsigned idx;
    for (idx = 0; idx < upsampler->factor; ++idx)
    {
        FIRFilterFlushD(upsampler->polyphase[idx]);
    }
    return NOERR;
}


/* UpsamplerProcess ****************************************************/
Error_t
UpsamplerProcess(Upsampler      *upsampler,
                 float          *outBuffer,
                 const float    *inBuffer,
                 unsigned       n_samples)
{
    float tempbuf[n_samples];
    if (upsampler && outBuffer)
    {
        for (unsigned filt = 0; filt < upsampler->factor; ++filt)
        {
            FIRFilterProcess(upsampler->polyphase[filt], tempbuf, inBuffer, n_samples);
            CopyBufferStride(outBuffer+filt, upsampler->factor, tempbuf, 1, n_samples);
        }

        VectorScalarMultiply(outBuffer, (const float*)outBuffer,
                             upsampler->factor, n_samples * upsampler->factor);
        return NOERR;
    }
    else
    {
        return NULL_PTR_ERROR;
    }
}

Error_t
UpsamplerProcessD(UpsamplerD*   upsampler,
                 double*        outBuffer,
                 const double*  inBuffer,
                 unsigned       n_samples)
{
    double tempbuf[n_samples];
    if (upsampler && outBuffer)
    {
        for (unsigned filt = 0; filt < upsampler->factor; ++ filt)
        {
            FIRFilterProcessD(upsampler->polyphase[filt], tempbuf, inBuffer, n_samples);
            CopyBufferStrideD(outBuffer+filt, upsampler->factor, tempbuf, 1, n_samples);
        }

        VectorScalarMultiplyD(outBuffer, (const double*)outBuffer,
                             upsampler->factor, n_samples * upsampler->factor);
        return NOERR;
    }
    else
    {
        return NULL_PTR_ERROR;
    }
}


/* Define log2 and log2f for MSVC */
#ifdef _USE_FXDSP_LOG

double
log2(double n)
{
    return log(n) / M_LN2;
}

float
log2f(float n)
{
    return logf(n) / (float)M_LN2;
}

#endif


/* 32 bit "pointer cast" union */
typedef union
{
    float f;
    int i;
} f_pcast32;



int
next_pow2(int x)
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



inline float
f_abs(float f)
{
    // flip the sign bit...
    int i = ((*(int*)&f) & 0x7fffffff);
    return (*(float*)&i);
}





/* f_max **********************************************************************/
inline float
f_max(float x, float a)
{
    x -= a;
    x += fabs(x);
    x *= 0.5;
    x += a;
    return x;
}



/* f_min **********************************************************************/
inline float
f_min(float x, float b)
{
    x = b - x;
    x += fabs(x);
    x *= 0.5;
    x = b - x;
    return x;
}



/* f_clamp ********************************************************************/
inline float
f_clamp(float x, float a, float b)
{
    const float x1 = fabs(x - a);
    const float x2 = fabs(x - b);
    x = x1 + a + b;
    x -= x2;
    x *= 0.5;
    return x;
}



/* f_pow2 *********************************************************************/
inline float
f_pow2(float x)
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
inline float
f_tanh(float x)
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
inline float
int16ToFloat(signed short sample)
{
    return (float)(sample * INT16_TO_FLOAT_SCALAR);
}

/* floatToInt16 ***************************************************************/
inline signed short
floatToInt16(float sample)
{
    return (signed short)(sample * 32767.0);
}


/* ratioToDb ******************************************************************/
inline float
AmpToDb(float amplitude)
{
    return 20.0*log10f(amplitude);
}

inline double
AmpToDbD(double amplitude)
{
    return 20.0*log10(amplitude);
}

/* dbToRatio ******************************************************************/
inline float
DbToAmp(float dB)
{
    return (dB > -150.0f ? expf(dB * LOG_TEN_OVER_TWENTY): 0.0f);
}

inline double
DbToAmpD(double dB)
{
    return (dB > -150.0 ? exp(dB * LOG_TEN_OVER_TWENTY): 0.0);
}

void
RectToPolar(float real, float imag, float* outMag, float* outPhase)
{
    *outMag = sqrtf(powf(real, 2.0) + powf(imag, 2.0));
    double phase = atanf(imag/real);
    if (phase < 0)
        phase += M_PI;
    *outPhase = phase;
}

void
RectToPolarD(double real, double imag, double* outMag, double* outPhase)
{
    *outMag = sqrt(pow(real, 2.0) + pow(imag, 2.0));
    double phase = atan(imag/real);
    if (phase < 0)
        phase += M_PI;
    *outPhase = phase;
}

void
PolarToRect(float mag, float phase, float* outReal, float* outImag)
{
    *outReal = mag * cosf(phase);
    *outImag = mag * sinf(phase);
}

void
PolarToRectD(double mag, double phase, double* outReal, double* outImag)
{
    *outReal = mag * cos(phase);
    *outImag = mag * sin(phase);
}

/*******************************************************************************
 Hann */
Error_t 
hann(unsigned n, float* dest)
{

#ifdef __APPLE__
    // Use the accelerate version if we have it
    // TODO: FIX THIS!!!!!!
    vDSP_hann_window(dest, n - 1, vDSP_HANN_DENORM);

#else
    // Otherwise do it manually
    unsigned buf_idx;
    for (buf_idx = 0; buf_idx < n; ++buf_idx)
    {
        *dest++ = 0.5 * (1 - cosf((2 * M_PI * buf_idx) / (n - 1)));
    }

#endif // __APPLE__
    return NOERR;
}

Error_t
hannD(unsigned n, double* dest)
{
    
#ifdef __APPLE__
    // Use the accelerate version if we have it
    vDSP_hann_windowD(dest, n - 1, vDSP_HANN_DENORM);
    
#else
    // Otherwise do it manually
    unsigned buf_idx;
    for (buf_idx = 0; buf_idx < n; ++buf_idx)
    {
        *dest++ = 0.5 * (1 - cos((2 * M_PI * buf_idx) / (n - 1)));
    }
    
#endif // __APPLE__
    return NOERR;
}



/*******************************************************************************
 Hamming */
Error_t 
hamming(unsigned n, float* dest)
{
    // do it manually. it seems like vDSP_hamm_window is wrong.
    unsigned buf_idx;
    for (buf_idx = 0; buf_idx < n; ++buf_idx)
    {
        *dest++ = 0.54 - 0.46 * cosf((2 * M_PI * buf_idx) / (n - 1));
    }
    return NOERR;
}


Error_t
hammingD(unsigned n, double* dest)
{
    unsigned buf_idx;
    for (buf_idx = 0; buf_idx < n; ++buf_idx)
    {
        *dest++ = 0.54 - 0.46 * cos((2 * M_PI * buf_idx) / (n - 1));
    }
    return NOERR;
}


/*******************************************************************************
 Blackman */
Error_t 
blackman(unsigned n, float a, float* dest)
{
    if (a > 0.15 && a < 0.17)
    {
        #ifdef __APPLE__
        // Use the builtin version for the specific case it implements, if it is
        // available.
        vDSP_blkman_window(dest, (n - 1), 0);
        return NOERR;
        #endif
    }

    // Otherwise do it manually
    float a0 = (1 - a) / 2;
    float a1 = 0.5;
    float a2 = a / 2;
    
    unsigned buf_idx;
    for (buf_idx = 0; buf_idx < n; ++buf_idx)
    {
        *dest++ = a0 - a1 * cosf((2 * M_PI * buf_idx) / (n - 1)) + a2 * cosf((4 * M_PI * buf_idx) / (n - 1));
    }
    return NOERR;
}



Error_t
blackmanD(unsigned n, double a, double* dest)
{
    if (a > 0.15 && a < 0.17)
    {
#ifdef __APPLE__
        // Use the builtin version for the specific case it implements, if it is
        // available.
        vDSP_blkman_windowD(dest, (n - 1), 0);
        return NOERR;
#endif
    }
    
    // Otherwise do it manually
    double a0 = (1 - a) / 2;
    double a1 = 0.5;
    double a2 = a / 2;
    
    unsigned buf_idx;
    for (buf_idx = 0; buf_idx < n; ++buf_idx)
    {
        *dest++ = a0 - a1 * cos((2 * M_PI * buf_idx) / (n - 1)) + a2 * cos((4 * M_PI * buf_idx) / (n - 1));
    }
    return NOERR;
}





/*******************************************************************************
 Tukey */
Error_t 
tukey(unsigned n, float a, float* dest)
{
    float term = a * (n - 1) / 2;
    unsigned buf_idx;
    for (buf_idx = 0; buf_idx < n; ++buf_idx)
    {
        if (buf_idx <= term)
        {
            dest[buf_idx] = 0.5 * (1 + cosf(M_PI * ((2 * buf_idx) /
                                                    (a * (n - 1)) - 1)));
        }
        else if (term <= buf_idx && buf_idx <= (n - 1) * (1 - a / 2))
        {
            dest[buf_idx] = 1.0;
        }
        else
        {
            dest[buf_idx] = 0.5 * (1 + cosf(M_PI *((2 * buf_idx) /
                                                   (a * (n - 1)) - (2 / a) + 1)));
        }
                
    }
    return NOERR;
}


Error_t
tukeyD(unsigned n, double a, double* dest)
{
    double term = a * (n - 1) / 2;
    unsigned buf_idx;
    for (buf_idx = 0; buf_idx < n; ++buf_idx)
    {
        if (buf_idx <= term)
        {
            dest[buf_idx] = 0.5 * (1 + cos(M_PI * ((2 * buf_idx) /
                                                    (a * (n - 1)) - 1)));
        }
        else if (term <= buf_idx && buf_idx <= (n - 1) * (1 - a / 2))
        {
            dest[buf_idx] = 1.0;
        }
        else
        {
            dest[buf_idx] = 0.5 * (1 + cos(M_PI *((2 * buf_idx) /
                                                   (a * (n - 1)) - (2 / a) + 1)));
        }
        
    }
    return NOERR;
}


/*******************************************************************************
 Cosine */
Error_t 
cosine(unsigned n, float* dest)
{
    float N = n - 1.0;
    unsigned buf_idx;
    for (buf_idx = 0; buf_idx < n; ++buf_idx)
    {
        *dest++ = sinf((M_PI * buf_idx) / N);
    }
    return NOERR;
}


Error_t
cosineD(unsigned n, double* dest)
{
    double N = n - 1.0;
    unsigned buf_idx;
    for (buf_idx = 0; buf_idx < n; ++buf_idx)
    {
        *dest++ = sin((M_PI * buf_idx) / N);
    }
    return NOERR;
}

/*******************************************************************************
 Lanczos */
Error_t 
lanczos(unsigned n, float* dest)
{
    float N = n - 1.0;
    unsigned buf_idx;
    for (buf_idx = 0; buf_idx < n; ++buf_idx)
    {
        float term = M_PI * ((2 * n)/ N) - 1.0;
        *dest++ = sinf(term) / term;
    }
    return NOERR;
}

Error_t
lanczosD(unsigned n, double* dest)
{
    double N = n - 1.0;
    unsigned buf_idx;
    for (buf_idx = 0; buf_idx < n; ++buf_idx)
    {
        double term = M_PI * ((2 * n)/ N) - 1.0;
        *dest++ = sin(term) / term;
    }
    return NOERR;
}

/*******************************************************************************
 Bartlett */
Error_t 
bartlett(unsigned n, float* dest)
{
    unsigned buf_idx;
    for (buf_idx = 0; buf_idx < n; ++buf_idx)
    {
        if (buf_idx <= (n - 1) / 2)
        {
            *dest++ = (float)(2 * buf_idx)/(n - 1);
        }
        else
        {
            *dest ++ = 2.0 - (float)(2 * buf_idx) / (n - 1);
        }
    }
    return NOERR;
}


Error_t
bartlettD(unsigned n, double* dest)
{
    unsigned buf_idx;
    for (buf_idx = 0; buf_idx < n; ++buf_idx)
    {
        if (buf_idx <= (n - 1) / 2)
        {
            *dest++ = (double)(2 * buf_idx)/(n - 1);
        }
        else
        {
            *dest ++ = 2.0 - (double)(2 * buf_idx) / (n - 1);
        }
    }
    return NOERR;
}

/*******************************************************************************
 Gaussian */
Error_t 
gaussian(unsigned n, float sigma, float* dest)
{
    float N = n - 1;
    float L = N / 2.0;
    unsigned buf_idx;
    for (buf_idx = 0; buf_idx < n; ++buf_idx)
    {   
        *dest++ = expf(-0.5 * powf((buf_idx - L)/(sigma * L),2));
    }
    return NOERR;
}

Error_t
gaussianD(unsigned n, double sigma, double* dest)
{
    double N = n - 1;
    double L = N / 2.0;
    unsigned buf_idx;
    for (buf_idx = 0; buf_idx < n; ++buf_idx)
    {
        *dest++ = std::exp(-0.5 * std::pow((buf_idx - L)/(sigma * L),2));
    }
    return NOERR;
}


/*******************************************************************************
 Bartlett-Hann */
Error_t 
bartlett_hann(unsigned n, float* dest)
{
    float N = n - 1.0;
    unsigned buf_idx;
    for (buf_idx = 0; buf_idx < n; ++buf_idx)
    {
        float term = ((buf_idx / N) - 0.5);
        *dest++ = 0.62 - 0.48 * fabsf(term) + 0.38 * cosf(2 * M_PI * term);
    }
    return NOERR;
}


Error_t
bartlett_hannD(unsigned n, double* dest)
{
    double N = n - 1.0;
    unsigned buf_idx;
    for (buf_idx = 0; buf_idx < n; ++buf_idx)
    {
        double term = ((buf_idx / N) - 0.5);
        *dest++ = 0.62 - 0.48 * std::fabs(term) + 0.38 * std::cos(2 * M_PI * term);
    }
    return NOERR;
}

/*******************************************************************************
 Kaiser */
Error_t 
kaiser(unsigned n, float a, float* dest)
{
    // Pre-calc
    float beta = M_PI * a;
    float m_2 = (float)(n - 1.0) / 2.0;
    float denom = modZeroBessel(beta);
    
    unsigned buf_idx;
    for (buf_idx = 0; buf_idx < n; ++buf_idx)
    {
        float val = ((buf_idx) - m_2) / m_2;
        val = 1 - (val * val);
        *dest++ = modZeroBessel(beta * sqrtf(val)) / denom;
    }
    return NOERR;
}

Error_t
kaiserD(unsigned n, double a, double* dest)
{
    // Pre-calc
    double beta = M_PI * a;
    double m_2 = (float)(n - 1.0) / 2.0;
    double denom = modZeroBesselD(beta);
    
    unsigned buf_idx;
    for (buf_idx = 0; buf_idx < n; ++buf_idx)
    {
        double val = ((buf_idx) - m_2) / m_2;
        val = 1 - (val * val);
        *dest++ = modZeroBesselD(beta * std::sqrt(val)) / denom;
    }
    return NOERR;
}


/*******************************************************************************
 Nuttall */
Error_t 
nuttall(unsigned n, float* dest)
{
    float term;
    float N = n - 1.0;
    unsigned buf_idx;
    for (buf_idx = 0; buf_idx < n; ++buf_idx)
    {
        term = 2 * M_PI * (buf_idx / N);
        *dest++ = 0.3635819  - 0.4891775 * cosf(term) + 0.1365995 *
        cosf(2 * term) - 0.0106411 * cosf(3 * term);
    }
    return NOERR;
}

Error_t
nuttallD(unsigned n, double* dest)
{
    double term;
    double N = n - 1.0;
    unsigned buf_idx;
    for (buf_idx = 0; buf_idx < n; ++buf_idx)
    {
        term = 2 * M_PI * (buf_idx / N);
        *dest++ = 0.3635819  - 0.4891775 * std::cos(term) + 0.1365995 *
        std::cos(2 * term) - 0.0106411 * cos(3 * term);
    }
    return NOERR;
}


/*******************************************************************************
 Blackman-Harris */
Error_t 
blackman_harris(unsigned n, float* dest)
{
    float term;
    unsigned buf_idx;
    for (buf_idx = 0; buf_idx < n; ++buf_idx)
    {
        term = (2 * M_PI * buf_idx) / (n - 1);
        *dest++ = 0.35875 - 0.48829 * cosf(term)+ 0.14128 * cosf(2 * term) -
        0.01168 * cosf(3 * term);
    }
    return NOERR;    
}


Error_t
blackman_harrisD(unsigned n, double* dest)
{
    float term;
    unsigned buf_idx;
    for (buf_idx = 0; buf_idx < n; ++buf_idx)
    {
        term = (2 * M_PI * buf_idx) / (n - 1);
        *dest++ = 0.35875 - 0.48829 * std::cos(term)+ 0.14128 * std::cos(2 * term) -
        0.01168 * std::cos(3 * term);
    }
    return NOERR;
}


/*******************************************************************************
 Blackman-Nuttall*/
Error_t
blackman_nuttall(unsigned n, float* dest)
{
    float term;
    unsigned buf_idx;
    for (buf_idx = 0; buf_idx < n; ++buf_idx)
    {
        term = (2 * M_PI * buf_idx) / (n - 1);
        *dest++ = 0.3635819 - 0.4891775 * cosf(term)+ 0.1365995 * cosf(2 * term) - 0.0106411 * cosf(3 * term);
    }
    return NOERR;    
}

Error_t
blackman_nuttallD(unsigned n, double* dest)
{
    double term;
    unsigned buf_idx;
    for (buf_idx = 0; buf_idx < n; ++buf_idx)
    {
        term = (2 * M_PI * buf_idx) / (n - 1);
        *dest++ = 0.3635819 - 0.4891775 * std::cos(term)+ 0.1365995 * std::cos(2 * term) - 0.0106411 * std::cos(3 * term);
    }
    return NOERR;
}

/*******************************************************************************
 Flat-Top */
Error_t
flat_top(unsigned n, float* dest)
{
    float N = n - 1.0;
    float term;
    unsigned buf_idx;
    for (buf_idx = 0; buf_idx < n; ++buf_idx)
    {
        term = (2 * M_PI * buf_idx) / N;
        *dest++ = 0.21557895 - 0.41663158 * cosf(term)+ 0.277263158 *
        cosf(2 * term) - 0.083578947 * cosf(3 * term) + 0.006947368 *
        cosf(4 * term);
    }
    return NOERR;    
}


Error_t
flat_topD(unsigned n, double* dest)
{
    double N = n - 1.0;
    double term;
    unsigned buf_idx;
    for (buf_idx = 0; buf_idx < n; ++buf_idx)
    {
        term = (2 * M_PI * buf_idx) / N;
        *dest++ = 0.21557895 - 0.41663158 * std::cos(term)+ 0.277263158 *
        std::cos(2 * term) - 0.083578947 * std::cos(3 * term) + 0.006947368 *
        std::cos(4 * term);
    }
    return NOERR;
}


/*******************************************************************************
 Poisson */
Error_t 
poisson(unsigned n, float D, float* dest)
{
    float term = (n - 1) / 2;
    float tau_inv = 1. / ((n / 2) * (8.69 / D));

    unsigned buf_idx;
    for (buf_idx = 0; buf_idx < n; ++buf_idx)
    {

        *dest++ = expf(-fabs(buf_idx - term) * tau_inv);
    }
    return NOERR;    
}

Error_t
poissonD(unsigned n, double D, double* dest)
{
    double term = (n - 1) / 2;
    double tau_inv = 1. / ((n / 2) * (8.69 / D));
    
    unsigned buf_idx;
    for (buf_idx = 0; buf_idx < n; ++buf_idx)
    {
        
        *dest++ = exp(-fabs(buf_idx - term) * tau_inv);
    }
    return NOERR;
}


/*******************************************************************************
 Chebyshev */
// TODO: FIX This.
Error_t
chebyshev(unsigned n, float A, float *dest)
{
    float max = 0;
    float N = n - 1.0;
    float M = N / 2;
    float tg = powf(10, A / 20.0);
    float x0 = coshf((1.0 / N) * acoshf(tg));
    
    for(unsigned buf_idx=0; buf_idx<(n/2+1); ++buf_idx)
    {
        float y = buf_idx - M;
        float sum = 0;
        for(unsigned i=1; i<=M; i++){
            sum += chebyshev_poly(N, x0 * cosf(M_PI * i / n)) *
            cosf( 2.0 * y * M_PI * i / n);
        }
        dest[buf_idx] = tg + 2 * sum;
        dest[(unsigned)N - buf_idx] = dest[buf_idx];
        if(dest[buf_idx] > max)
        {
            max = dest[buf_idx];
        }
    }
    
    for(unsigned buf_idx = 0; buf_idx < n; ++buf_idx)
    {
        dest[buf_idx] /= max;
    }
    return NOERR;
}



Error_t
chebyshevD(unsigned n, double A, double *dest)
{
    double max=0;
    double N = n - 1.0;
    double M = N / 2;
    double tg = pow(10, A / 20.0);
    double x0 = cosh((1.0 / N) * acosh(tg));

    for(unsigned buf_idx=0; buf_idx<(n/2+1); ++buf_idx)
    {
        double y = buf_idx - M;
        double sum = 0;
        for(unsigned i=1; i<=M; i++){
            sum += chebyshev_polyD(N, x0 * cos(M_PI * i / n)) *
            cos( 2.0 * y * M_PI * i / n);
        }
        dest[buf_idx] = tg + 2 * sum;
        dest[(unsigned)N - buf_idx] = dest[buf_idx];
        if(dest[buf_idx] > max)
        {
            max = dest[buf_idx];
        }
    }

    for(unsigned buf_idx = 0; buf_idx < n; ++buf_idx)
    {
        dest[buf_idx] /= max; 
    }
    return NOERR;
}

/* Modified Bessel function of the first kind */
static float
modZeroBessel(float x)
{
    float x_2 = x/2;
    float num = 1;
    float fact = 1;
    float result = 1;
    
    unsigned i;
    for (i=1 ; i<20 ; i++) 
    {
        num *= x_2 * x_2;
        fact *= i;
        result += num / (fact * fact);
    }
    return result;
}


static double
modZeroBesselD(double x)
{
    double x_2 = x/2;
    double num = 1;
    double fact = 1;
    double result = 1;
    
    unsigned i;
    for (i=1 ; i<20 ; i++)
    {
        num *= x_2 * x_2;
        fact *= i;
        result += num / (fact * fact);
    }
    return result;
}


/* Chebyshev Polynomial */

static float
chebyshev_poly(int n, float x)
{
    float y;
    if (fabsf(x) <= 1)
    {
        y = cosf(n * acosf(x));
    }
    else
    {
        y = coshf(n * acoshf(x));
    }
    
    return y;
}

static double
chebyshev_polyD(int n, double x)
{
    float y;
    if (fabs(x) <= 1)
    {
        y = cos(n * acos(x));
    }
    else
    {
        y = cosh(n * acosh(x));
    }
    
    return y;
}


/*******************************************************************************
 WindowFunctionInit */
WindowFunction*
WindowFunctionInit(unsigned n, Window_t type)
{
    WindowFunction* window = (WindowFunction*)malloc(sizeof(WindowFunction));
    
    window->length = n;
    window->window = (float*)malloc(n * sizeof(float));
    window->type = type;
    
    switch (type)
    {
        case BOXCAR:
            boxcar(window->length, window->window);
            break;
        case HANN:
            hann(window->length, window->window);
            break;
        case HAMMING:
            hamming(window->length, window->window);
            break;
        case BLACKMAN:
            blackman(window->length, 0.16, window->window);
            break;
        case TUKEY:
            tukey(window->length, 0.5, window->window);
            break;
        case COSINE:
            cosine(window->length, window->window);
            break;
        case LANCZOS:
            lanczos(window->length, window->window);
            break;
        case BARTLETT:
            bartlett(window->length, window->window);
            break;
        case GAUSSIAN:
            gaussian(window->length, 0.4, window->window);
            break;
        case BARTLETT_HANN:
            bartlett_hann(window->length, window->window);
            break;
        case KAISER:
            kaiser(window->length, 0.5, window->window);
            break;
        case NUTTALL:
            nuttall(window->length, window->window);
            break;
        case BLACKMAN_HARRIS:
            blackman_harris(window->length, window->window);
            break;
        case BLACKMAN_NUTTALL:
            blackman_nuttall(window->length, window->window);
            break;
        case FLATTOP:
            flat_top(window->length, window->window);
            break;
        case POISSON:
            poisson(window->length, 8.69, window->window);
            
        default:
            boxcar(window->length, window->window);
            break;
    }
    
    return window;
}

WindowFunctionD*
WindowFunctionInitD(unsigned n, Window_t type)
{
    WindowFunctionD* window = (WindowFunctionD*)malloc(sizeof(WindowFunctionD));
    
    window->length = n;
    window->window = (double*)malloc(n * sizeof(double));
    window->type = type;
    
    switch (type)
    {
        case BOXCAR:
            boxcarD(window->length, window->window);
            break;
        case HANN:
            hannD(window->length, window->window);
            break;
        case HAMMING:
            hammingD(window->length, window->window);
            break;
        case BLACKMAN:
            blackmanD(window->length, 0.16, window->window);
            break;
        case TUKEY:
            tukeyD(window->length, 0.5, window->window);
            break;
        case COSINE:
            cosineD(window->length, window->window);
            break;
        case LANCZOS:
            lanczosD(window->length, window->window);
            break;
        case BARTLETT:
            bartlettD(window->length, window->window);
            break;
        case GAUSSIAN:
            gaussianD(window->length, 0.4, window->window);
            break;
        case BARTLETT_HANN:
            bartlett_hannD(window->length, window->window);
            break;
        case KAISER:
            kaiserD(window->length, 0.5, window->window);
            break;
        case NUTTALL:
            nuttallD(window->length, window->window);
            break;
        case BLACKMAN_HARRIS:
            blackman_harrisD(window->length, window->window);
            break;
        case BLACKMAN_NUTTALL:
            blackman_nuttallD(window->length, window->window);
            break;
        case FLATTOP:
            flat_topD(window->length, window->window);
            break;
        case POISSON:
            poissonD(window->length, 8.69, window->window);
            
        default:
            boxcarD(window->length, window->window);
            break;
    }
    
    return window;
}

/*******************************************************************************
 WindowFunctionFree */
Error_t
WindowFunctionFree(WindowFunction* window)
{
    if (window)
    {
        if (window->window)
        {
            free(window->window);
            window->window = NULL;
        }

        free(window);
        window = NULL;
    }
    return NOERR;
}

Error_t
WindowFunctionFreeD(WindowFunctionD* window)
{
    if (window)
    {
        if (window->window)
        {
            free(window->window);
            window->window = NULL;
        }
        
        free(window);
        window = NULL;
    }
    return NOERR;
}

/*******************************************************************************
 WindowFunctionProcess */
Error_t
WindowFunctionProcess(WindowFunction*   window,
                      float*            outBuffer,
                      const float*      inBuffer,
                      unsigned          n_samples)
{
    VectorVectorMultiply(outBuffer, inBuffer, window->window, n_samples);
    return NOERR;
}


Error_t
WindowFunctionProcessD(WindowFunctionD* window,
                       double*          outBuffer,
                       const double*    inBuffer,
                       unsigned         n_samples)
{
    VectorVectorMultiplyD(outBuffer, inBuffer, window->window, n_samples);
    return NOERR;
}


/*******************************************************************************
 Boxcar */
Error_t 
boxcar(unsigned n, float* dest)
{
    FillBuffer(dest, n, 1.0);
    return NOERR;
}

Error_t
boxcarD(unsigned n, double* dest)
{
    FillBufferD(dest, n, 1.0);
    return NOERR;
}

float
midiNoteToFrequency(unsigned note)
{
    return powf(2.0, ((note - 69.0)/12.)) * 440.0;
}

unsigned
frequencyToMidiNote(float f)
{
    return (unsigned)(69 + (12 * log2f(f / 440.0)));
}

#undef USE_FFTW_FFT
#undef E_INV
#undef INT16_TO_FLOAT_SCALAR
#undef INV_LN2
#undef TWO_PI
#undef PI_OVER_TWO
#undef INVERSE_TWO_PI
#undef LOG_TEN_OVER_TWENTY
#undef TWENTY_OVER_LOG_TEN
#undef SQRT_TWO_OVER_TWO
#undef LIMIT
#undef LIN_INTERP
#undef HZ_TO_RAD
#undef RAD_TO_HZ
#undef F_EXP
#undef DB_TO_AMP
#undef DB_TO_AMPD
#undef AMP_TO_DB
#undef AMP_TO_DBD
#undef SMOOTH_ABS
#undef USE_FFT_CONVOLUTION_LENGTH
#undef FILT_Q
#undef K12_REF
#undef K14_REF
#undef K20_REF
#undef PREFILTER_FC
#undef PREFILTER_GAIN 
#undef PREFILTER_Q    
#undef RLBFILTER_FC   
#undef RLBFILTER_Q    
#undef GATE_LENGTH_S  
#undef GATE_OVERLAP   
#undef SPEED_SATURATION_COEFF
#undef N_FLUTTER_COMPONENTS
