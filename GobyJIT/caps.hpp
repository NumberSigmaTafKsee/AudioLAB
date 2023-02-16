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

} /* namespace DSP */
