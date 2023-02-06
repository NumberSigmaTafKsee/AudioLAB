#pragma once

#define _GNU_SOURCE 1
#define _USE_GNU 1

/* unlocking some standard math calls. */
#define __USE_ISOC99 1
#define __USE_ISOC9X 1
#define _ISOC99_SOURCE 1
#define _ISOC9X_SOURCE 1

#include <stdlib.h>
#include <string.h>

#include <math.h>
#include <float.h>

#include <assert.h>
#include <stdio.h>

#include <complex>


typedef __int8_t			int8;
typedef __uint8_t			uint8;
typedef __int16_t			int16;
typedef __uint16_t		uint16;
typedef __int32_t			int32;
typedef __uint32_t		uint32;
typedef __int64_t			int64;
typedef __uint64_t		uint64;

#define MIN_GAIN 1e-6 /* -120 dB */
/* smallest non-denormal 32 bit IEEE float is 1.18e-38 */
#define NOISE_FLOOR 1e-20 /* -400 dB */

/* //////////////////////////////////////////////////////////////////////// */

typedef float sample_t;
typedef unsigned int uint;
typedef unsigned long ulong;

/* prototype that takes a sample and yields a sample */
typedef sample_t (*clip_func_t) (sample_t);

#ifndef max
template <class X, class Y> X min (X x, Y y) { return x < (X)y ? x : (X)y; }
template <class X, class Y> X max (X x, Y y) { return x > (X)y ? x : (X)y; }
#endif /* ! max */

template <class T>
T clamp (T value, T lower, T upper)
{
	if (value < lower) return lower;
	if (value > upper) return upper;
	return value;
}

static inline float frandom() { return (float) random() / (float) RAND_MAX; }

/* NB: also true if 0  */
inline bool 
is_denormal (float & f)
{
	int32 i = *((int32 *) &f);
	return ((i & 0x7f800000) == 0);
}

/* not used, check validity before using */
inline bool 
is_denormal (double & f)
{
	int64 i = *((int64 *) &f);
	return ((i & 0x7fe0000000000000ll) == 0);
}

/* lovely algorithm from 
  http://graphics.stanford.edu/~seander/bithacks.html#RoundUpPowerOf2Float
*/
inline uint 
next_power_of_2 (uint n)
{
	assert (n <= 0x40000000);

	--n;
	n |= n >> 1;
	n |= n >> 2;
	n |= n >> 4;
	n |= n >> 8;
	n |= n >> 16;

	return ++n;
}

inline double db2lin (double db) { return pow(10, .05*db); }
inline double lin2db (double lin) { return 20*log10(lin); }

#if defined(__i386__) || defined(__amd64__)
	#define TRAP asm ("int $3;")
#else
	#define TRAP
#endif

namespace DSP {

inline float pow2 (float x) { return x * x; }
inline float pow3 (float x) { return x * pow2(x); }
inline float pow4 (float x) { return pow2 (pow2(x)); }
inline float pow5 (float x) { return x * pow4(x); }
inline float pow6 (float x) { return pow3 (pow2(x)); }
inline float pow7 (float x) { return x * (pow6 (x)); }
inline float pow8 (float x) { return pow2 (pow4 (x)); }

inline float 
sgn (float x)
{
	union { float f; uint32 i; } u;
	u.f = x;
	u.i &= 0x80000000;
	u.i |= 0x3F800000;
	return u.f;
}

inline bool
isprime (int v)
{
	if (v <= 3)
		return true;
	
	if (!(v & 1))
		return false;

	for (int i = 3; i < (int) sqrt (v) + 1; i += 2)
		if ((v % i) == 0)
			return false;

	return true;
}

} /* namespace DSP */


namespace DSP {


class Sine
{
	public:
		int z;
		double y[2];
		double b;

	public:
		Sine()
			{ 
				b = 0;
				y[0] = y[1] = 0;
				z = 0;
			}

		Sine (double f, double fs, double phase)
			{
				set_f (f, fs, phase);
			}

		Sine (double w, double phase = 0.)
			{
				set_f (w, phase);
			}

		inline void set_f (double f, double fs, double phase)
			{
				set_f (f*2*M_PI/fs, phase);
			}

		inline void set_f (double w, double phase)
			{
				b = 2*cos(w);
				y[0] = sin (phase - w);
				y[1] = sin (phase - 2*w);
				z = 0;
			}

		/* advance and return 1 sample */
		inline double get()
			{
				double s = b*y[z]; 
				z ^= 1;
				s -= y[z];
				return y[z] = s;
			}

		double get_phase()
			{
				double x0 = y[z], x1 = b*y[z] - y[z^1];
				double phi = asin(x0);
				
				/* slope is falling: into the 2nd half. */
				return x1 < x0 ? M_PI - phi : phi;
			}
};

/* same as above but including a damping coefficient d */
class DampedSine
: public Sine
{
	public:
		double d;

	public:
		DampedSine()
			{ d = 1; }

		inline double get()
			{
				double s = b * y[z]; 
				z ^= 1;
				s -= d * y[z];
				return y[z] = d * s;
			}
};

inline void
sinc (double omega, sample_t * s, int n)
{
	/* initial phase */
	double phi = (n / 2) * -omega;

	Sine sine (omega, phi);
	
	for (int i = 0; i < n; ++i, phi += omega)
	{
		double sin_phi = sine.get();

		if (fabs (phi) < 0.000000001)
			s[i] = 1.;
		else
			s[i] = sin_phi / phi;
	}
}

/* prototypes for window value application ... */
typedef void (*window_sample_func_t) (sample_t &, sample_t);

/* ... which go as template parameters for the window calculation below */
inline void store_sample (sample_t &d, sample_t s) { d = s; }
inline void apply_window (sample_t &d, sample_t s) { d *= s; }

template <window_sample_func_t F>
void
hann (sample_t * s, int n, double step = 1)
{
	step = M_PI*step/n;
	
	double f = 0;
	for (int i=0; i<n; ++i)
	{
		f = i*step;
		float x = sin(f);
		F(s[i], x*x);
	}
}

/* faster but less accurate version */
template <window_sample_func_t F>
void
hann2 (sample_t *s, int n)
{
	double phi = M_PI/(n-1);
	DSP::Sine sin(phi, 0);
	for (int i=0; i<n; ++i)
	{
		float x = sin.get();
		F(s[i], x*x);
	}
}


template <window_sample_func_t F>
void
hamming (sample_t * s, int n)
{
	float in = 1. / n;
	
	for (int i = 0; i < n; ++i)
	{
		 double f = i*in;
		F (s[i], .54 - .46*cos (2*M_PI*f));
	}
}

template <window_sample_func_t F>
void
blackman (sample_t *s, int n)
{
	float in = 1. / n;

	for (int i = 0; i < n; ++i)
	{
		 float f = (float) i;

		 double b = .42f - 
						.5f*cos (2.f*f*M_PI*in) + 
						.08*cos (4.f*f*M_PI*in);

		F (s[i], b);
	}
}

template <window_sample_func_t F>
void
blackman_harris (sample_t *s, int n)
{
	double w1 = 2.f*M_PI / (n - 1);
	double w2 = 2.f*w1;
	double w3 = 3.f*w1;

	for (int i = 0; i < n; ++i)
	{
		 double f = (double) i;

		 double bh = .35875f - 
				.48829f*cos (w1*f) + 
				.14128f*cos (w2*f) - 
				.01168f*cos (w3*f);

		F (s[i], bh);
	}
}

/* by way of dobson and csound */
inline double 
besseli (double x)
{
	double a = fabs(x);
	if (a < 3.75)     
	{
		double y = x/3.75;
		y *= y;
		return 1. + y*(3.5156229 + y*(3.0899424 + y*(1.2067492 +
					y*(.2659732 + y*(.0360768 + y*.0045813)))));
	}
	double y = 3.75/a;
	return (exp(a)/sqrt(a)) * (.39894228 + y*(.01328592 + 
			y*(.00225319 + y*(-.00157565 + y*(.00916281 + y*(-.02057706 + 
			y*(.02635537 + y*(-.01647633 + y*.00392377))))))));
}

/* step = .5 : window [-n to 0] */
template <window_sample_func_t F>
void
kaiser (sample_t *s, int n, double beta, double step = 1)
{
	double bb = besseli(beta);
	int si = 0;

	for(double i = -n/2.+.5; si < n; ++si, i += step)
	{
		double a = 1 - pow((2*i / (n - 1)), 2);
		double k = besseli((beta*(a < 0 ? 0 : sqrt(a))) / bb);
		F(s[si], k);
	}
}

template <window_sample_func_t F>
void
xfade (sample_t *s, int n, int dir) /* dir [-1,1] */
{
	DSP::Sine cos(.5*M_PI/n, dir>0 ? .5*M_PI : 0);
	for (int i=0; i<n; ++i)
	{
		float c = cos.get();
		F(s[i], c*c);
	}
}

/* 
	Brute-force FIR filter with downsampling method (decimating). 
*/
class FIR
{
	public:
		/* kernel length, history length - 1 */
		uint n, m;
		
		/* coefficients, history */
		sample_t * c, * x;

		/* history index */
		int h; 
		
		FIR() { c = x = 0; }
		~FIR() { free(c); free(x); }
		
		void init (uint N)
			{
				n = N;

				/* keeping history size a power of 2 makes it possible to wrap the
				 * history pointer by & instead of %, saving a few cpu cycles. */
				m = next_power_of_2 (n);

				c = (sample_t *) malloc (n * sizeof (sample_t));
				x = (sample_t *) malloc (m * sizeof (sample_t));

				m -= 1;

				reset();
			}
	
		void reset()
			{
				h = 0;
				memset (x, 0, n * sizeof (sample_t));
			}
		
		inline sample_t process (sample_t s)
			{
				x[h] = s;
				
				s *= c[0];

				for (uint Z = 1, z = h - 1; Z < n; --z, ++Z)
					s += c[Z] * x[z & m];

				h = (h + 1) & m;

				return s;
			}

		/* Z is the time, in samples, since the last non-zero sample.
		 * OVER is the oversampling factor. just here for documentation, use
		 * a FIRUpsampler instead.
		 */
		template <uint Z, uint OVER>
		inline sample_t upsample (sample_t s)
			{
				x[h] = s;
				
				s = 0;

				/* for the interpolation, iterate over the history in z ^ -OVER
				 * steps -- all the samples between are 0.
				 */
				for (uint j = Z, z = h - Z; j < n; --z, j += OVER)
					s += c[j] * x[z & m];

				h = (h + 1) & m;

				return s;
			}

		/* used in downsampling */
		inline void store (sample_t s)
			{
				x[h] = s;
				h = (h + 1) & m;
			}
};

/* FIR upsampler, optimised not to store the 0 samples */
template <int N, int Oversample>
class FIRUpsampler
{
	public:
		uint m; /* history length - 1 */
		int h; /* history index */

		sample_t * c, * x; /* coefficients, history */

		FIRUpsampler()
			{
				c = x = 0;
				init();
			}

		~FIRUpsampler()
			{
				if (c) free (c);
				if (x) free (x);
			}
		
		void init()
			{
				/* FIR kernel length must be a multiple of the oversampling ratio */
				assert (N % Oversample == 0);

				/* like FIR, keep the history buffer a power of 2; additionally,
				 * don't store the 0 samples inbetween. */
				m = next_power_of_2 ((N + Oversample - 1) / Oversample);

				c = (sample_t *) malloc (N * sizeof (sample_t));
				x = (sample_t *) malloc (m * sizeof (sample_t));

				m -= 1;

				reset();
			}
	
		void reset()
			{
				h = 0;
				memset (x, 0, (m + 1) * sizeof (sample_t));
			}
		
		/* upsample the given sample */
		inline sample_t upsample (sample_t s)
			{
				x[h] = s;
				
				s = 0;

				for (uint Z = 0, z = h; Z < N; --z, Z += Oversample)
					s += c[Z] * x[z & m];

				h = (h + 1) & m;

				return s;
			}

		/* upsample a zero sample (interleaving), Z being the time, in samples,
		 * since the last non-0 sample. */
		inline sample_t pad (uint Z)
			{
				sample_t s = 0;

				for (uint z = h-1; Z < N; --z, Z += Oversample)
					s += c[Z] * x[z & m];

				return s;
			}

};

/* templating for kernel length allows g++ to optimise aggressively
 * resulting in appreciable performance gains. */
template <int N>
class FIRn
{
	public:
		/* history length - 1 */
		uint m;
		
		/* coefficients, history */
		sample_t c[N], x[N];

		/* history index */
		int h; 
		
		FIRn()
			{
				/* keeping history size a power of 2 makes it possible to wrap the
				 * history pointer by & instead of %, saving a few cpu cycles. */
				m = next_power_of_2 (N) - 1;

				reset();
			}
	
		void reset()
			{
				h = 0;
				memset (x, 0, N * sizeof (sample_t));
			}
		
		inline sample_t process (sample_t s)
			{
				x[h] = s;
				
				s *= c[0];

				for (uint Z = 1, z = h - 1; Z < N; --z, ++Z)
					s += c[Z] * x[z & m];

				h = (h + 1) & m;

				return s;
			}

		/* used in downsampling */
		inline void store (sample_t s)
			{
				x[h] = s;
				h = (h + 1) & m;
			}
};

/* do-nothing 1:1 oversampler substitute, occasionally needed as a 
 * template parameter */
class NoOversampler 
{
	public:
		enum { Ratio = 1 };
		sample_t downsample (sample_t x) { return x; }
		sample_t upsample (sample_t x) { return x; }
		void downstore (sample_t) { }
		sample_t uppad (uint) { return 0; }
};

template <int Oversample, int FIRSize>
class Oversampler
{
	public:
		enum { Ratio = Oversample };
		/* antialias filters */
		struct {
			DSP::FIRUpsampler<FIRSize, Oversample> up;
			DSP::FIRn<FIRSize> down;
		} fir;

		Oversampler()
			{ init(); }

		void init (float fc = .5) 
			{
				double f = fc * M_PI / Oversample;
				
				/* construct the upsampler filter kernel */
				DSP::sinc (f, fir.up.c, FIRSize);
				DSP::kaiser<DSP::apply_window> (fir.up.c, FIRSize, 6.4);

				/* copy upsampler filter kernel for downsampler, make sum */
				double s = 0;
				for (uint i = 0; i < FIRSize; ++i)
					fir.down.c[i] = fir.up.c[i],
					s += fir.up.c[i];
				
				/* scale downsampler kernel for unity gain */
				s = 1/s;
				for (uint i = 0; i < FIRSize; ++i)
					fir.down.c[i] *= s;

				/* scale upsampler kernel for unity gain */
				s *= Oversample;
				for (uint i = 0; i < FIRSize; ++i)
					fir.up.c[i] *= s;
			}

		void reset() 
			{
				fir.up.reset();
				fir.down.reset();
			}

		inline sample_t upsample(sample_t x)
			{ return fir.up.upsample(x); }
		inline sample_t uppad(uint z)
			{ return fir.up.pad(z); }

		inline sample_t downsample(sample_t x)
			{ return fir.down.process(x); }
		inline void downstore(sample_t x)
			{ fir.down.store(x); }
};


/** The type of filter that the State Variable Filter will output. */
enum SVFType {
    SVFLowpass = 0,
    SVFBandpass,
    SVFHighpass,
    SVFUnitGainBandpass,
    SVFBandShelving,
    SVFNotch,
    SVFAllpass,
    SVFPeak
};

//==============================================================================
class VAStateVariableFilter {
public:
    /** Create and initialize the filter with default values defined in constructor. */
    VAStateVariableFilter();

    //------------------------------------------------------------------------------

    ~VAStateVariableFilter();

    //------------------------------------------------------------------------------

    /**    Sets the type of the filter that processAudioSample() or processAudioBlock() will
        output. This filter can choose between 8 different types using the enums listed
        below or the int given to each.
        0: SVFLowpass
        1: SVFBandpass
        2: SVFHighpass
        3: SVFUnitGainBandpass
        4: SVFBandShelving
        5: SVFNotch
        6: SVFAllpass
        7: SVFPeak
    */
    void setFilterType(int newType);

    //------------------------------------------------------------------------------
    /**    Used for changing the filter's cutoff parameter linearly by frequency (Hz) */
    void setCutoffFreq(double newCutoffFreq);

    //------------------------------------------------------------------------------
    /** Used for setting the resonance amount. This is then converted to a Q
        value, which is used by the filter.
        Range: (0-1)
    */
    void setResonance(double newResonance);

    //------------------------------------------------------------------------------
    /** Used for setting the filter's Q amount. This is then converted to a
        damping parameter called R, which is used in the original filter.
    */
    void setQ(double newQ);

    //------------------------------------------------------------------------------
    /**    Sets the gain of the shelf for the BandShelving filter only. */
    void setShelfGain(double newGain);

    //------------------------------------------------------------------------------
    /**    Statically set the filters parameters. */
    void setFilter(int newType, double newCutoff,
                   double newResonance, double newShelfGain);

    //------------------------------------------------------------------------------
    /**    Set the sample rate used by the host. Needs to be used to accurately
        calculate the coefficients of the filter from the cutoff.
        Note: This is often used in AudioProcessor::prepareToPlay
    */
    void setSampleRate(double newSampleRate);

    //------------------------------------------------------------------------------
    /**    Performs the actual processing.
    */
    void process(float gain, const float *input, float *output, unsigned count);

    //------------------------------------------------------------------------------
    /**    Reset the state variables.
    */
    void clear();

    //------------------------------------------------------------------------------
    /**    Compute the transfer function at given frequency.
    */
    std::complex<double> calcTransfer(double freq) const;

    //------------------------------------------------------------------------------


    double getCutoffFreq() const { return cutoffFreq; }

    double getFilterType() const { return filterType; }

    double getQ() const { return Q; }

    double getShelfGain() const { return shelfGain; }

private:
    //==============================================================================
    //    Calculate the coefficients for the filter based on parameters.
    void calcFilter();

    //
    template <int FilterType>
    void processInternally(float gain, const float *input, float *output, unsigned count);

    //    Parameters:
    int filterType;
    double cutoffFreq;
    double Q;
    double shelfGain;

    double sampleRate;

    //    Coefficients:
    double gCoeff;        // gain element
    double RCoeff;        // feedback damping element
    double KCoeff;        // shelf gain element

    double z1_A, z2_A;        // state variables (z^-1)

private:
    static std::complex<double> calcTransferLowpass(double w, double wc, double r);
    static std::complex<double> calcTransferBandpass(double w, double wc, double r);
    static std::complex<double> calcTransferHighpass(double w, double wc, double r);
    static std::complex<double> calcTransferUnitGainBandpass(double w, double wc, double r);
    static std::complex<double> calcTransferBandShelving(double w, double wc, double r, double k);
    static std::complex<double> calcTransferNotch(double w, double wc, double r);
    static std::complex<double> calcTransferAllpass(double w, double wc, double r);
    static std::complex<double> calcTransferPeak(double w, double wc, double r);
};

#if __cplusplus >= 201703L
# define if_constexpr if constexpr
#else
# define if_constexpr if
#endif

//==============================================================================

static double resonanceToQ(double resonance)
{
    return 1.0 / (2.0 * (1.0 - resonance));
}

//==============================================================================

VAStateVariableFilter::VAStateVariableFilter()
{
    sampleRate = 44100.0;                // default sample rate when constructed
    filterType = SVFLowpass;            // lowpass filter by default

    gCoeff = 1.0;
    RCoeff = 1.0;
    KCoeff = 0.0;

    cutoffFreq = 1000.0;
    Q = resonanceToQ(0.5);
    shelfGain = 1.0;

    z1_A = 0.0;
    z2_A = 0.0;
}

VAStateVariableFilter::~VAStateVariableFilter()
{
}

// Member functions for setting the filter's parameters (and sample rate).
//==============================================================================
void VAStateVariableFilter::setFilterType(int newType)
{
    filterType = newType;
}

void VAStateVariableFilter::setCutoffFreq(double newCutoffFreq)
{
    if (cutoffFreq == newCutoffFreq)
        return;

    cutoffFreq = newCutoffFreq;
    calcFilter();
}

void VAStateVariableFilter::setResonance(double newResonance)
{
    setQ(resonanceToQ(newResonance));
}

void VAStateVariableFilter::setQ(double newQ)
{
    if (Q == newQ)
        return;

    Q = newQ;
    calcFilter();
}

void VAStateVariableFilter::setShelfGain(double newGain)
{
    if (shelfGain == newGain)
        return;

    shelfGain = newGain;
    calcFilter();
}

void VAStateVariableFilter::setFilter(int newType, double newCutoffFreq,
                                      double newResonance, double newShelfGain)
{
    double newQ = resonanceToQ(newResonance);

    if (filterType == newType && cutoffFreq == newCutoffFreq && Q == newQ && shelfGain == newShelfGain)
        return;

    filterType = newType;
    cutoffFreq = newCutoffFreq;
    Q = newQ;
    shelfGain = newShelfGain;
    calcFilter();
}

void VAStateVariableFilter::setSampleRate(double newSampleRate)
{
    if (sampleRate == newSampleRate)
        return;

    sampleRate = newSampleRate;
    calcFilter();
}

//==============================================================================
void VAStateVariableFilter::calcFilter()
{
    // prewarp the cutoff (for bilinear-transform filters)
    double wd = cutoffFreq * (2.0 * M_PI);
    double T = 1.0 / sampleRate;
    double wa = (2.0 / T) * std::tan(wd * T / 2.0);

    // Calculate g (gain element of integrator)
    gCoeff = wa * T / 2.0;            // Calculate g (gain element of integrator)

    // Calculate Zavalishin's R from Q (referred to as damping parameter)
    RCoeff = 1.0 / (2.0 * Q);

    // Gain for BandShelving filter
    KCoeff = shelfGain;
}

static double analogSaturate(double x)
{
    // simple filter analog saturation

    if (x > +1)
        x = 2. / 3.;
    else if (x < -1)
        x = -2. / 3.;
    else
        x = x - (x * x * x) * (1.0 / 3.0);

    return x;
}

template <int FilterType>
void VAStateVariableFilter::processInternally(float gain, const float *input, float *output, unsigned count)
{
    const double gCoeff = this->gCoeff;
    const double RCoeff = this->RCoeff;
    const double KCoeff = this->KCoeff;

    double z1_A = this->z1_A;
    double z2_A = this->z2_A;

    for (unsigned i = 0; i < count; ++i) {
        double in = gain * input[i];

        double HP = (in - ((2.0 * RCoeff + gCoeff) * z1_A) - z2_A)
            * (1.0 / (1.0 + (2.0 * RCoeff * gCoeff) + gCoeff * gCoeff));
        double BP = HP * gCoeff + z1_A;
        double LP = BP * gCoeff + z2_A;

        z1_A = analogSaturate(gCoeff * HP + BP);        // unit delay (state variable)
        z2_A = analogSaturate(gCoeff * BP + LP);        // unit delay (state variable)

        // Selects which filter type this function will output.
        double out = 0.0;
        if_constexpr (FilterType == SVFLowpass)
            out = LP;
        else if_constexpr (FilterType == SVFBandpass)
            out = BP;
        else if_constexpr (FilterType == SVFHighpass)
            out = HP;
        else if_constexpr (FilterType == SVFUnitGainBandpass)
            out = 2.0 * RCoeff * BP;
        else if_constexpr (FilterType == SVFBandShelving)
            out = in + 2.0 * RCoeff * KCoeff * BP;
        else if_constexpr (FilterType == SVFNotch)
            out = in - 2.0 * RCoeff * BP;
        else if_constexpr (FilterType == SVFAllpass)
            out = in - 4.0 * RCoeff * BP;
        else if_constexpr (FilterType == SVFPeak)
            out = LP - HP;

        output[i] = out;
    }

    this->z1_A = z1_A;
    this->z2_A = z2_A;
}

void VAStateVariableFilter::process(float gain, const float *input, float *output, unsigned count)
{
    switch (filterType) {
    case SVFLowpass:
        processInternally<SVFLowpass>(gain, input, output, count);
        break;
    case SVFBandpass:
        processInternally<SVFBandpass>(gain, input, output, count);
        break;
    case SVFHighpass:
        processInternally<SVFHighpass>(gain, input, output, count);
        break;
    case SVFUnitGainBandpass:
        processInternally<SVFUnitGainBandpass>(gain, input, output, count);
        break;
    case SVFBandShelving:
        processInternally<SVFBandShelving>(gain, input, output, count);
        break;
    case SVFNotch:
        processInternally<SVFNotch>(gain, input, output, count);
        break;
    case SVFAllpass:
        processInternally<SVFAllpass>(gain, input, output, count);
        break;
    case SVFPeak:
        processInternally<SVFPeak>(gain, input, output, count);
        break;
    default:
        for (unsigned i = 0; i < count; ++i)
            output[i] = gain * input[i];
    }
}

void VAStateVariableFilter::clear()
{
    z1_A = 0;
    z2_A = 0;
}

std::complex<double> VAStateVariableFilter::calcTransfer(double freq) const
{
    double w = 2 * M_PI * freq;
    double wc = 2 * M_PI * cutoffFreq;

    switch (filterType) {
    case SVFLowpass:
        return calcTransferLowpass(w, wc, RCoeff);
    case SVFBandpass:
        return calcTransferBandpass(w, wc, RCoeff);
    case SVFHighpass:
        return calcTransferHighpass(w, wc, RCoeff);
    case SVFUnitGainBandpass:
        return calcTransferUnitGainBandpass(w, wc, RCoeff);
    case SVFBandShelving:
        return calcTransferBandShelving(w, wc, RCoeff, shelfGain);
    case SVFNotch:
        return calcTransferNotch(w, wc, RCoeff);
    case SVFAllpass:
        return calcTransferAllpass(w, wc, RCoeff);
    case SVFPeak:
        return calcTransferPeak(w, wc, RCoeff);
    default:
        return 0.0;
    }
}

//==============================================================================

std::complex<double> VAStateVariableFilter::calcTransferLowpass(double w, double wc, double r)
{
    std::complex<double> s = w * std::complex<double>(0, 1);
    return (wc * wc) / (s * s + 2.0 * r * wc * s + wc * wc);
}

std::complex<double> VAStateVariableFilter::calcTransferBandpass(double w, double wc, double r)
{
    std::complex<double> s = w * std::complex<double>(0, 1);
    return (wc * s) / (s * s + 2.0 * r * wc * s + wc * wc);
}

std::complex<double> VAStateVariableFilter::calcTransferHighpass(double w, double wc, double r)
{
    std::complex<double> s = w * std::complex<double>(0, 1);
    return (s * s) / (s * s + 2.0 * r * wc * s + wc * wc);
}

std::complex<double> VAStateVariableFilter::calcTransferUnitGainBandpass(double w, double wc, double r)
{
    return 2.0 * r * calcTransferBandpass(w, wc, r);
}

std::complex<double> VAStateVariableFilter::calcTransferBandShelving(double w, double wc, double r, double k)
{
    return 1.0 + k * calcTransferUnitGainBandpass(w, wc, r);
}

std::complex<double> VAStateVariableFilter::calcTransferNotch(double w, double wc, double r)
{
    return calcTransferBandShelving(w, wc, r, -1.0);
}

std::complex<double> VAStateVariableFilter::calcTransferAllpass(double w, double wc, double r)
{
    return calcTransferBandShelving(w, wc, r, -2.0);
}

std::complex<double> VAStateVariableFilter::calcTransferPeak(double w, double wc, double r)
{
    std::complex<double> s = w * std::complex<double>(0, 1);
    return (wc * wc - s * s) / (s * s + 2.0 * r * wc * s + wc * wc);
}

} /* namespace DSP */

