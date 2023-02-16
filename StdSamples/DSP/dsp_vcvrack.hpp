#pragma once

#include <cmath>
#include <cstring>
#include <algorithm>
#include <cstdint>
#include <complex>
#include <cassert>

#include "speex_resampler.h"
#include <pffft.h>

namespace simd
{
	float ifelse(bool e, float t1, float t2) {
		return e? t1:t2;
	}
}

/** Useful for storing arrays of samples in ring buffers and casting them to `float*` to be used by interleaved processors, like SampleRateConverter */
template <size_t CHANNELS, typename T = float>
struct Frame {
	T samples[CHANNELS];
};

/** A simple cyclic buffer.
S must be a power of 2.
Thread-safe for single producers and consumers.
*/
template <typename T, size_t S>
struct RingBuffer {
	T data[S];
	size_t start = 0;
	size_t end = 0;

	size_t mask(size_t i) const {
		return i & (S - 1);
	}
	void push(T t) {
		size_t i = mask(end++);
		data[i] = t;
	}
	void pushBuffer(const T* t, int n) {
		size_t i = mask(end);
		size_t e1 = i + n;
		size_t e2 = (e1 < S) ? e1 : S;
		std::memcpy(&data[i], t, sizeof(T) * (e2 - i));
		if (e1 > S) {
			std::memcpy(data, &t[S - i], sizeof(T) * (e1 - S));
		}
		end += n;
	}
	T shift() {
		return data[mask(start++)];
	}
	void shiftBuffer(T* t, size_t n) {
		size_t i = mask(start);
		size_t s1 = i + n;
		size_t s2 = (s1 < S) ? s1 : S;
		std::memcpy(t, &data[i], sizeof(T) * (s2 - i));
		if (s1 > S) {
			std::memcpy(&t[S - i], data, sizeof(T) * (s1 - S));
		}
		start += n;
	}
	void clear() {
		start = end;
	}
	bool empty() const {
		return start == end;
	}
	bool full() const {
		return end - start == S;
	}
	size_t size() const {
		return end - start;
	}
	size_t capacity() const {
		return S - size();
	}
};

/** A cyclic buffer which maintains a valid linear array of size S by keeping a copy of the buffer in adjacent memory.
S must be a power of 2.
Thread-safe for single producers and consumers?
*/
template <typename T, size_t S>
struct DoubleRingBuffer {
	T data[S * 2];
	size_t start = 0;
	size_t end = 0;

	size_t mask(size_t i) const {
		return i & (S - 1);
	}
	void push(T t) {
		size_t i = mask(end++);
		data[i] = t;
		data[i + S] = t;
	}
	T shift() {
		return data[mask(start++)];
	}
	void clear() {
		start = end;
	}
	bool empty() const {
		return start == end;
	}
	bool full() const {
		return end - start == S;
	}
	size_t size() const {
		return end - start;
	}
	size_t capacity() const {
		return S - size();
	}
	/** Returns a pointer to S consecutive elements for appending.
	If any data is appended, you must call endIncr afterwards.
	Pointer is invalidated when any other method is called.
	*/
	T* endData() {
		return &data[mask(end)];
	}
	void endIncr(size_t n) {
		size_t e = mask(end);
		size_t e1 = e + n;
		size_t e2 = (e1 < S) ? e1 : S;
		// Copy data forward
		std::memcpy(&data[S + e], &data[e], sizeof(T) * (e2 - e));

		if (e1 > S) {
			// Copy data backward from the doubled block to the main block
			std::memcpy(data, &data[S], sizeof(T) * (e1 - S));
		}
		end += n;
	}
	/** Returns a pointer to S consecutive elements for consumption
	If any data is consumed, call startIncr afterwards.
	*/
	const T* startData() const {
		return &data[mask(start)];
	}
	void startIncr(size_t n) {
		start += n;
	}
};

/** A cyclic buffer which maintains a valid linear array of size S by sliding along a larger block of size N.
The linear array of S elements are moved back to the start of the block once it outgrows past the end.
This happens every N - S pushes, so the push() time is O(1 + S / (N - S)).
For example, a float buffer of size 64 in a block of size 1024 is nearly as efficient as RingBuffer.
Not thread-safe.
*/
template <typename T, size_t S, size_t N>
struct AppleRingBuffer {
	T data[N];
	size_t start = 0;
	size_t end = 0;

	void returnBuffer() {
		// move end block to beginning
		// may overlap, but memmove handles that correctly
		size_t s = size();
		std::memmove(data, &data[start], sizeof(T) * s);
		start = 0;
		end = s;
	}
	void push(T t) {
		if (end + 1 > N) {
			returnBuffer();
		}
		data[end++] = t;
	}
	T shift() {
		return data[start++];
	}
	bool empty() const {
		return start == end;
	}
	bool full() const {
		return end - start == S;
	}
	size_t size() const {
		return end - start;
	}
	size_t capacity() const {
		return S - size();
	}
	/** Returns a pointer to S consecutive elements for appending, requesting to append n elements.
	*/
	T* endData(size_t n) {
		if (end + n > N) {
			returnBuffer();
		}
		return &data[end];
	}
	/** Actually increments the end position
	Must be called after endData(), and `n` must be at most the `n` passed to endData()
	*/
	void endIncr(size_t n) {
		end += n;
	}
	/** Returns a pointer to S consecutive elements for consumption
	If any data is consumed, call startIncr afterwards.
	*/
	const T* startData() const {
		return &data[start];
	}
	void startIncr(size_t n) {
		// This is valid as long as n < S
		start += n;
	}
};


inline double sinc(double x)
{
	if(x == 0.0) return 1.0f;
	x *= M_PI;
    return std::sin(x) / x;
}

/** Performs a direct sum convolution */
inline float convolveNaive(const float* in, const float* kernel, int len) {
	float y = 0.f;
	for (int i = 0; i < len; i++) {
		y += in[len - 1 - i] * kernel[i];
	}
	return y;
}

/** Computes the impulse response of a boxcar lowpass filter */
inline void boxcarLowpassIR(float* out, int len, float cutoff = 0.5f) {
	for (int i = 0; i < len; i++) {
		float t = i - (len - 1) / 2.f;
		out[i] = 2 * cutoff * sinc(2 * cutoff * t);
	}
}


struct RealTimeConvolver {
	// `kernelBlocks` number of contiguous FFT blocks of size `blockSize`
	// indexed by [i * blockSize*2 + j]
	float* kernelFfts = NULL;
	float* inputFfts = NULL;
	float* outputTail = NULL;
	float* tmpBlock = NULL;
	size_t blockSize;
	size_t kernelBlocks = 0;
	size_t inputPos = 0;
	PFFFT_Setup* pffft;

	/** `blockSize` is the size of each FFT block. It should be >=32 and a power of 2. */
	RealTimeConvolver(size_t blockSize) {
		this->blockSize = blockSize;
		pffft = pffft_new_setup(blockSize * 2, PFFFT_REAL);
		outputTail = new float[blockSize];
		std::memset(outputTail, 0, blockSize * sizeof(float));
		tmpBlock = new float[blockSize * 2];
		std::memset(tmpBlock, 0, blockSize * 2 * sizeof(float));
	}

	~RealTimeConvolver() {
		setKernel(NULL, 0);
		delete[] outputTail;
		delete[] tmpBlock;
		pffft_destroy_setup(pffft);
	}

	void setKernel(const float* kernel, size_t length) {
		// Clear existing kernel
		if (kernelFfts) {
			pffft_aligned_free(kernelFfts);
			kernelFfts = NULL;
		}
		if (inputFfts) {
			pffft_aligned_free(inputFfts);
			inputFfts = NULL;
		}
		kernelBlocks = 0;
		inputPos = 0;

		if (kernel && length > 0) {
			// Round up to the nearest factor of `blockSize`
			kernelBlocks = (length - 1) / blockSize + 1;

			// Allocate blocks
			kernelFfts = (float*) pffft_aligned_malloc(sizeof(float) * blockSize * 2 * kernelBlocks);
			inputFfts = (float*) pffft_aligned_malloc(sizeof(float) * blockSize * 2 * kernelBlocks);
			std::memset(inputFfts, 0, sizeof(float) * blockSize * 2 * kernelBlocks);

			for (size_t i = 0; i < kernelBlocks; i++) {
				// Pad each block with zeros
				std::memset(tmpBlock, 0, sizeof(float) * blockSize * 2);
				size_t len = std::min((int) blockSize, (int)(length - i * blockSize));
				std::memcpy(tmpBlock, &kernel[i * blockSize], sizeof(float)*len);
				// Compute fft
				pffft_transform(pffft, tmpBlock, &kernelFfts[blockSize * 2 * i], NULL, PFFFT_FORWARD);
			}
		}
	}

	/** Applies reverb to input
	input and output must be of size `blockSize`
	*/
	void processBlock(const float* input, float* output) {
		if (kernelBlocks == 0) {
			std::memset(output, 0, sizeof(float) * blockSize);
			return;
		}

		// Step input position
		inputPos = (inputPos + 1) % kernelBlocks;
		// Pad block with zeros
		std::memset(tmpBlock, 0, sizeof(float) * blockSize * 2);
		std::memcpy(tmpBlock, input, sizeof(float) * blockSize);
		// Compute input fft
		pffft_transform(pffft, tmpBlock, &inputFfts[blockSize * 2 * inputPos], NULL, PFFFT_FORWARD);
		// Create output fft
		std::memset(tmpBlock, 0, sizeof(float) * blockSize * 2);
		// convolve input fft by kernel fft
		// Note: This is the CPU bottleneck loop
		for (size_t i = 0; i < kernelBlocks; i++) {
			size_t pos = (inputPos - i + kernelBlocks) % kernelBlocks;
			pffft_zconvolve_accumulate(pffft, &kernelFfts[blockSize * 2 * i], &inputFfts[blockSize * 2 * pos], tmpBlock, 1.f);
		}
		// Compute output
		pffft_transform(pffft, tmpBlock, tmpBlock, NULL, PFFFT_BACKWARD);
		// Add block tail from last output block
		for (size_t i = 0; i < blockSize; i++) {
			tmpBlock[i] += outputTail[i];
		}
		// Copy output block to output
		float scale = 1.f / (blockSize * 2);
		for (size_t i = 0; i < blockSize; i++) {
			// Scale based on FFT
			output[i] = tmpBlock[i] * scale;
		}
		// Set tail
		for (size_t i = 0; i < blockSize; i++) {
			outputTail[i] = tmpBlock[i + blockSize];
		}
	}
};


/** Real-valued FFT context.
Wrapper for [PFFFT](https://bitbucket.org/jpommier/pffft/)
`length` must be a multiple of 32.
Buffers must be aligned to 16-byte boundaries. new[] and malloc() do this for you.
*/
struct RealFFT {
	PFFFT_Setup* setup;
	int length;

	RealFFT(size_t length) {
		this->length = length;
		setup = pffft_new_setup(length, PFFFT_REAL);
	}

	~RealFFT() {
		pffft_destroy_setup(setup);
	}

	/** Performs the real FFT.
	Input and output must be aligned using the above align*() functions.
	Input is `length` elements. Output is `2*length` elements.
	Output is arbitrarily ordered for performance reasons.
	However, this ordering is consistent, so element-wise multiplication with line up with other results, and the inverse FFT will return a correctly ordered result.
	*/
	void rfftUnordered(const float* input, float* output) {
		pffft_transform(setup, input, output, NULL, PFFFT_FORWARD);
	}

	/** Performs the inverse real FFT.
	Input is `2*length` elements. Output is `length` elements.
	Scaling is such that IRFFT(RFFT(x)) = N*x.
	*/
	void irfftUnordered(const float* input, float* output) {
		pffft_transform(setup, input, output, NULL, PFFFT_BACKWARD);
	}

	/** Slower than the above methods, but returns results in the "canonical" FFT order as follows.
		output[0] = F(0)
		output[1] = F(n/2)
		output[2] = real(F(1))
		output[3] = imag(F(1))
		output[4] = real(F(2))
		output[5] = imag(F(2))
		...
		output[length - 2] = real(F(n/2 - 1))
		output[length - 1] = imag(F(n/2 - 1))
	*/
	void rfft(const float* input, float* output) {
		pffft_transform_ordered(setup, input, output, NULL, PFFFT_FORWARD);
	}

	void irfft(const float* input, float* output) {
		pffft_transform_ordered(setup, input, output, NULL, PFFFT_BACKWARD);
	}

	/** Scales the RFFT so that `scale(IFFT(FFT(x))) = x`.
	*/
	void scale(float* x) {
		float a = 1.f / length;
		for (int i = 0; i < length; i++) {
			x[i] *= a;
		}
	}
};


/** Complex-valued FFT context.
`length` must be a multiple of 16.
*/
struct ComplexFFT {
	PFFFT_Setup* setup;
	int length;

	ComplexFFT(size_t length) {
		this->length = length;
		setup = pffft_new_setup(length, PFFFT_COMPLEX);
	}

	~ComplexFFT() {
		pffft_destroy_setup(setup);
	}

	/** Performs the complex FFT.
	Input and output must be aligned using the above align*() functions.
	Input is `2*length` elements. Output is `2*length` elements.
	*/
	void fftUnordered(const float* input, float* output) {
		pffft_transform(setup, input, output, NULL, PFFFT_FORWARD);
	}

	/** Performs the inverse complex FFT.
	Input is `2*length` elements. Output is `2*length` elements.
	Scaling is such that FFT(IFFT(x)) = N*x.
	*/
	void ifftUnordered(const float* input, float* output) {
		pffft_transform(setup, input, output, NULL, PFFFT_BACKWARD);
	}

	void fft(const float* input, float* output) {
		pffft_transform_ordered(setup, input, output, NULL, PFFFT_FORWARD);
	}

	void ifft(const float* input, float* output) {
		pffft_transform_ordered(setup, input, output, NULL, PFFFT_BACKWARD);
	}

	void scale(float* x) {
		float a = 1.f / length;
		for (int i = 0; i < length; i++) {
			x[2 * i + 0] *= a;
			x[2 * i + 1] *= a;
		}
	}
};


/** Detects when a boolean changes from false to true */
struct BooleanTrigger {
	bool state = true;

	void reset() {
		state = true;
	}

	bool process(bool state) {
		bool triggered = (state && !this->state);
		this->state = state;
		return triggered;
	}
};


/** Turns HIGH when value reaches 1.f, turns LOW when value reaches 0.f. */
template <typename T = float>
struct TSchmittTrigger {
	T state;
	TSchmittTrigger() {
		reset();
	}
	void reset() {
		state = T::mask();
	}
	T process(T in) {
		T on = (in >= 1.f);
		T off = (in <= 0.f);
		T triggered = ~state & on;
		state = on | (state & ~off);
		return triggered;
	}
};


template <>
struct TSchmittTrigger<float> {
	bool state = true;

	void reset() {
		state = true;
	}

	/** Updates the state of the Schmitt Trigger given a value.
	Returns true if triggered, i.e. the value increases from 0 to 1.
	If different trigger thresholds are needed, use
		process(rescale(in, low, high, 0.f, 1.f))
	for example.
	*/
	bool process(float in) {
		if (state) {
			// HIGH to LOW
			if (in <= 0.f) {
				state = false;
			}
		}
		else {
			// LOW to HIGH
			if (in >= 1.f) {
				state = true;
				return true;
			}
		}
		return false;
	}

	bool isHigh() {
		return state;
	}
};

typedef TSchmittTrigger<> SchmittTrigger;


/** When triggered, holds a high value for a specified time before going low again */
struct PulseGenerator {
	float remaining = 0.f;

	/** Immediately disables the pulse */
	void reset() {
		remaining = 0.f;
	}

	/** Advances the state by `deltaTime`. Returns whether the pulse is in the HIGH state. */
	bool process(float deltaTime) {
		if (remaining > 0.f) {
			remaining -= deltaTime;
			return true;
		}
		return false;
	}

	/** Begins a trigger with the given `duration`. */
	void trigger(float duration = 1e-3f) {
		// Keep the previous pulse if the existing pulse will be held longer than the currently requested one.
		if (duration > remaining) {
			remaining = duration;
		}
	}
};


struct Timer {
	float time = 0.f;

	void reset() {
		time = 0.f;
	}

	/** Returns the time since last reset or initialization. */
	float process(float deltaTime) {
		time += deltaTime;
		return time;
	}
};


struct ClockDivider {
	uint32_t clock = 0;
	uint32_t division = 1;

	void reset() {
		clock = 0;
	}

	void setDivision(uint32_t division) {
		this->division = division;
	}

	uint32_t getDivision() {
		return division;
	}

	uint32_t getClock() {
		return clock;
	}

	/** Returns true when the clock reaches `division` and resets. */
	bool process() {
		clock++;
		if (clock >= division) {
			clock = 0;
			return true;
		}
		return false;
	}
};



/** The simplest possible analog filter using an Euler solver.
https://en.wikipedia.org/wiki/RC_circuit
Use two RC filters in series for a bandpass filter.
*/
template <typename T = float>
struct TRCFilter {
	T c = 0.f;
	T xstate[1];
	T ystate[1];

	TRCFilter() {
		reset();
	}

	void reset() {
		xstate[0] = 0.f;
		ystate[0] = 0.f;
	}

	/** Sets the cutoff angular frequency in radians.
	*/
	void setCutoff(T r) {
		c = 2.f / r;
	}
	/** Sets the cutoff frequency.
	`f` is the ratio between the cutoff frequency and sample rate, i.e. f = f_c / f_s
	*/
	void setCutoffFreq(T f) {
		setCutoff(2.f * M_PI * f);
	}
	void process(T x) {
		T y = (x + xstate[0] - ystate[0] * (1 - c)) / (1 + c);
		xstate[0] = x;
		ystate[0] = y;
	}
	T lowpass() {
		return ystate[0];
	}
	T highpass() {
		return xstate[0] - ystate[0];
	}
};

typedef TRCFilter<> RCFilter;


/** Applies exponential smoothing to a signal with the ODE
\f$ \frac{dy}{dt} = x \lambda \f$.
*/
template <typename T = float>
struct TExponentialFilter {
	T out = 0.f;
	T lambda = 0.f;

	void reset() {
		out = 0.f;
	}

	void setLambda(T lambda) {
		this->lambda = lambda;
	}

	void setTau(T tau) {
		this->lambda = 1 / tau;
	}

	T process(T deltaTime, T in) {
		T y = out + (in - out) * lambda * deltaTime;
		// If no change was made between the old and new output, assume T granularity is too small and snap output to input
		out = simd::ifelse(out == y, in, y);
		return out;
	}

	T process(T in) {
		return process(1.f, in);
	}
};

typedef TExponentialFilter<> ExponentialFilter;


/** Like ExponentialFilter but jumps immediately to higher values.
*/
template <typename T = float>
struct TPeakFilter {
	T out = 0.f;
	T lambda = 0.f;

	void reset() {
		out = 0.f;
	}

	void setLambda(T lambda) {
		this->lambda = lambda;
	}

	void setTau(T tau) {
		this->lambda = 1 / tau;
	}

	T process(T deltaTime, T in) {
		T y = out + (in - out) * lambda * deltaTime;
		out = std::fmax(y, in);
		return out;
	}
	/** Use the return value of process() instead. */
	T peak() {
		return out;
	}
	/** Use setLambda() instead. */
	void setRate(T r) {
		lambda = 1.f - r;
	}
	T process(T x) {
		return process(1.f, x);
	}
};

typedef TPeakFilter<> PeakFilter;


template <typename T = float>
struct TSlewLimiter {
	T out = 0.f;
	T rise = 0.f;
	T fall = 0.f;

	void reset() {
		out = 0.f;
	}

	void setRiseFall(T rise, T fall) {
		this->rise = rise;
		this->fall = fall;
	}
	T process(T deltaTime, T in) {
		out = std::clamp(in, out - fall * deltaTime, out + rise * deltaTime);
		return out;
	}
	T process(T in) {
		return process(1.f, in);
	}
};

typedef TSlewLimiter<> SlewLimiter;


template <typename T = float>
struct TExponentialSlewLimiter {
	T out = 0.f;
	T riseLambda = 0.f;
	T fallLambda = 0.f;

	void reset() {
		out = 0.f;
	}

	void setRiseFall(T riseLambda, T fallLambda) {
		this->riseLambda = riseLambda;
		this->fallLambda = fallLambda;
	}
	T process(T deltaTime, T in) {
		T lambda = simd::ifelse(in > out, riseLambda, fallLambda);
		T y = out + (in - out) * lambda * deltaTime;
		// If the change from the old out to the new out is too small for floats, set `in` directly.
		out = simd::ifelse(out == y, in, y);
		return out;
	}
	T process(T in) {
		return process(1.f, in);
	}
};

typedef TExponentialSlewLimiter<> ExponentialSlewLimiter;


template <typename T = float>
struct TBiquadFilter {
	/** input state */
	T x[2];
	/** output state */
	T y[2];

	/** transfer function numerator coefficients: b_0, b_1, b_2 */
	float b[3];
	/** transfer function denominator coefficients: a_1, a_2
	a_0 is fixed to 1.
	*/
	float a[2];

	enum Type {
		LOWPASS_1POLE,
		HIGHPASS_1POLE,
		LOWPASS,
		HIGHPASS,
		LOWSHELF,
		HIGHSHELF,
		BANDPASS,
		PEAK,
		NOTCH,
		NUM_TYPES
	};

	TBiquadFilter() {
		reset();
		setParameters(LOWPASS, 0.f, 0.f, 1.f);
	}

	void reset() {
		std::memset(x, 0, sizeof(x));
		std::memset(y, 0, sizeof(y));
	}

	T process(T in) {
		// Advance IIR
		T out = b[0] * in + b[1] * x[0] + b[2] * x[1]
		        - a[0] * y[0] - a[1] * y[1];
		// Push input
		x[1] = x[0];
		x[0] = in;
		// Push output
		y[1] = y[0];
		y[0] = out;
		return out;
	}

	/** Calculates and sets the biquad transfer function coefficients.
	f: normalized frequency (cutoff frequency / sample rate), must be less than 0.5
	Q: quality factor
	V: gain
	*/
	void setParameters(Type type, float f, float Q, float V) {
		float K = std::tan(M_PI * f);
		switch (type) {
			case LOWPASS_1POLE: {
				a[0] = -std::exp(-2.f * M_PI * f);
				a[1] = 0.f;
				b[0] = 1.f + a[0];
				b[1] = 0.f;
				b[2] = 0.f;
			} break;

			case HIGHPASS_1POLE: {
				a[0] = std::exp(-2.f * M_PI * (0.5f - f));
				a[1] = 0.f;
				b[0] = 1.f - a[0];
				b[1] = 0.f;
				b[2] = 0.f;
			} break;

			case LOWPASS: {
				float norm = 1.f / (1.f + K / Q + K * K);
				b[0] = K * K * norm;
				b[1] = 2.f * b[0];
				b[2] = b[0];
				a[0] = 2.f * (K * K - 1.f) * norm;
				a[1] = (1.f - K / Q + K * K) * norm;
			} break;

			case HIGHPASS: {
				float norm = 1.f / (1.f + K / Q + K * K);
				b[0] = norm;
				b[1] = -2.f * b[0];
				b[2] = b[0];
				a[0] = 2.f * (K * K - 1.f) * norm;
				a[1] = (1.f - K / Q + K * K) * norm;

			} break;

			case LOWSHELF: {
				float sqrtV = std::sqrt(V);
				if (V >= 1.f) {
					float norm = 1.f / (1.f + M_SQRT2 * K + K * K);
					b[0] = (1.f + M_SQRT2 * sqrtV * K + V * K * K) * norm;
					b[1] = 2.f * (V * K * K - 1.f) * norm;
					b[2] = (1.f - M_SQRT2 * sqrtV * K + V * K * K) * norm;
					a[0] = 2.f * (K * K - 1.f) * norm;
					a[1] = (1.f - M_SQRT2 * K + K * K) * norm;
				}
				else {
					float norm = 1.f / (1.f + M_SQRT2 / sqrtV * K + K * K / V);
					b[0] = (1.f + M_SQRT2 * K + K * K) * norm;
					b[1] = 2.f * (K * K - 1) * norm;
					b[2] = (1.f - M_SQRT2 * K + K * K) * norm;
					a[0] = 2.f * (K * K / V - 1.f) * norm;
					a[1] = (1.f - M_SQRT2 / sqrtV * K + K * K / V) * norm;
				}
			} break;

			case HIGHSHELF: {
				float sqrtV = std::sqrt(V);
				if (V >= 1.f) {
					float norm = 1.f / (1.f + M_SQRT2 * K + K * K);
					b[0] = (V + M_SQRT2 * sqrtV * K + K * K) * norm;
					b[1] = 2.f * (K * K - V) * norm;
					b[2] = (V - M_SQRT2 * sqrtV * K + K * K) * norm;
					a[0] = 2.f * (K * K - 1.f) * norm;
					a[1] = (1.f - M_SQRT2 * K + K * K) * norm;
				}
				else {
					float norm = 1.f / (1.f / V + M_SQRT2 / sqrtV * K + K * K);
					b[0] = (1.f + M_SQRT2 * K + K * K) * norm;
					b[1] = 2.f * (K * K - 1.f) * norm;
					b[2] = (1.f - M_SQRT2 * K + K * K) * norm;
					a[0] = 2.f * (K * K - 1.f / V) * norm;
					a[1] = (1.f / V - M_SQRT2 / sqrtV * K + K * K) * norm;
				}
			} break;

			case BANDPASS: {
				float norm = 1.f / (1.f + K / Q + K * K);
				b[0] = K / Q * norm;
				b[1] = 0.f;
				b[2] = -b[0];
				a[0] = 2.f * (K * K - 1.f) * norm;
				a[1] = (1.f - K / Q + K * K) * norm;
			} break;

			case PEAK: {
				if (V >= 1.f) {
					float norm = 1.f / (1.f + K / Q + K * K);
					b[0] = (1.f + K / Q * V + K * K) * norm;
					b[1] = 2.f * (K * K - 1.f) * norm;
					b[2] = (1.f - K / Q * V + K * K) * norm;
					a[0] = b[1];
					a[1] = (1.f - K / Q + K * K) * norm;
				}
				else {
					float norm = 1.f / (1.f + K / Q / V + K * K);
					b[0] = (1.f + K / Q + K * K) * norm;
					b[1] = 2.f * (K * K - 1.f) * norm;
					b[2] = (1.f - K / Q + K * K) * norm;
					a[0] = b[1];
					a[1] = (1.f - K / Q / V + K * K) * norm;
				}
			} break;

			case NOTCH: {
				float norm = 1.f / (1.f + K / Q + K * K);
				b[0] = (1.f + K * K) * norm;
				b[1] = 2.f * (K * K - 1.f) * norm;
				b[2] = b[0];
				a[0] = b[1];
				a[1] = (1.f - K / Q + K * K) * norm;
			} break;

			default: break;
		}
	}

	void copyParameters(const TBiquadFilter<T>& from) {
		b[0] = from.b[0];
		b[1] = from.b[1];
		b[2] = from.b[2];
		a[0] = from.a[0];
		a[1] = from.a[1];
	}

	/** Computes the gain of a particular frequency
	f: normalized frequency
	*/
	float getFrequencyResponse(float f) {
		// Compute sum(a_k e^(-i k f))
		std::complex<float> bsum = b[0];
		std::complex<float> asum = 1.f;
		for (int i = 1; i < 3; i++) {
			float p = 2 * M_PI * -i * f;
			std::complex<float> e(std::cos(p), std::sin(p));
			bsum += b[i] * e;
			asum += a[i - 1] * e;
		}
		return std::abs(bsum / asum);
	}
};

typedef TBiquadFilter<> BiquadFilter;




/** Resamples by a fixed rational factor. */
template <int MAX_CHANNELS>
struct SampleRateConverter {
	SpeexResamplerState* st = NULL;
	int channels = MAX_CHANNELS;
	int quality = SPEEX_RESAMPLER_QUALITY_DEFAULT;
	int inRate = 44100;
	int outRate = 44100;

	SampleRateConverter() {
		refreshState();
	}
	~SampleRateConverter() {
		if (st) {
			speex_resampler_destroy(st);
		}
	}

	/** Sets the number of channels to actually process. This can be at most MAX_CHANNELS. */
	void setChannels(int channels) {
		assert(channels <= MAX_CHANNELS);
		if (channels == this->channels)
			return;
		this->channels = channels;
		refreshState();
	}

	/** From 0 (worst, fastest) to 10 (best, slowest) */
	void setQuality(int quality) {
		if (quality == this->quality)
			return;
		this->quality = quality;
		refreshState();
	}

	void setRates(int inRate, int outRate) {
		if (inRate == this->inRate && outRate == this->outRate)
			return;
		this->inRate = inRate;
		this->outRate = outRate;
		refreshState();
	}

	void refreshState() {
		if (st) {
			speex_resampler_destroy(st);
			st = NULL;
		}

		if (channels > 0 && inRate != outRate) {
			int err;
			st = speex_resampler_init(channels, inRate, outRate, quality, &err);
			assert(st);
			assert(err == RESAMPLER_ERR_SUCCESS);

			speex_resampler_set_input_stride(st, MAX_CHANNELS);
			speex_resampler_set_output_stride(st, MAX_CHANNELS);
		}
	}

	/** `in` and `out` are interlaced with the number of channels */
	void process(const Frame<MAX_CHANNELS>* in, int* inFrames, Frame<MAX_CHANNELS>* out, int* outFrames) {
		assert(in);
		assert(inFrames);
		assert(out);
		assert(outFrames);
		if (st) {
			// Resample each channel at a time
			spx_uint32_t inLen;
			spx_uint32_t outLen;
			for (int i = 0; i < channels; i++) {
				inLen = *inFrames;
				outLen = *outFrames;
				int err = speex_resampler_process_float(st, i, ((const float*) in) + i, &inLen, ((float*) out) + i, &outLen);
				assert(err == RESAMPLER_ERR_SUCCESS);
			}
			*inFrames = inLen;
			*outFrames = outLen;
		}
		else {
			// Simply copy the buffer without conversion
			int frames = std::min(*inFrames, *outFrames);
			std::memcpy(out, in, frames * sizeof(Frame<MAX_CHANNELS>));
			*inFrames = frames;
			*outFrames = frames;
		}
	}
};


/** Downsamples by an integer factor. */
template <int OVERSAMPLE, int QUALITY, typename T = float>
struct Decimator {
	T inBuffer[OVERSAMPLE * QUALITY];
	float kernel[OVERSAMPLE * QUALITY];
	int inIndex;

	Decimator(float cutoff = 0.9f) {
		boxcarLowpassIR(kernel, OVERSAMPLE * QUALITY, cutoff * 0.5f / OVERSAMPLE);
		blackmanHarrisWindow(kernel, OVERSAMPLE * QUALITY);
		reset();
	}
	void reset() {
		inIndex = 0;
		std::memset(inBuffer, 0, sizeof(inBuffer));
	}
	/** `in` must be length OVERSAMPLE */
	T process(T* in) {
		// Copy input to buffer
		std::memcpy(&inBuffer[inIndex], in, OVERSAMPLE * sizeof(T));
		// Advance index
		inIndex += OVERSAMPLE;
		inIndex %= OVERSAMPLE * QUALITY;
		// Perform naive convolution
		T out = 0.f;
		for (int i = 0; i < OVERSAMPLE * QUALITY; i++) {
			int index = inIndex - 1 - i;
			index = (index + OVERSAMPLE * QUALITY) % (OVERSAMPLE * QUALITY);
			out += kernel[i] * inBuffer[index];
		}
		return out;
	}
};


/** Upsamples by an integer factor. */
template <int OVERSAMPLE, int QUALITY>
struct Upsampler {
	float inBuffer[QUALITY];
	float kernel[OVERSAMPLE * QUALITY];
	int inIndex;

	Upsampler(float cutoff = 0.9f) {
		boxcarLowpassIR(kernel, OVERSAMPLE * QUALITY, cutoff * 0.5f / OVERSAMPLE);
		blackmanHarrisWindow(kernel, OVERSAMPLE * QUALITY);
		reset();
	}
	void reset() {
		inIndex = 0;
		std::memset(inBuffer, 0, sizeof(inBuffer));
	}
	/** `out` must be length OVERSAMPLE */
	void process(float in, float* out) {
		// Zero-stuff input buffer
		inBuffer[inIndex] = OVERSAMPLE * in;
		// Advance index
		inIndex++;
		inIndex %= QUALITY;
		// Naively convolve each sample
		// TODO replace with polyphase filter hierarchy
		for (int i = 0; i < OVERSAMPLE; i++) {
			float y = 0.f;
			for (int j = 0; j < QUALITY; j++) {
				int index = inIndex - 1 - j;
				index = (index + QUALITY) % QUALITY;
				int kernelIndex = OVERSAMPLE * j + i;
				y += kernel[kernelIndex] * inBuffer[index];
			}
			out[i] = y;
		}
	}
};


/** Hann window function.
p: proportion from [0, 1], usually `i / (len - 1)`
https://en.wikipedia.org/wiki/Window_function#Hann_and_Hamming_windows
*/
template <typename T>
inline T hann(T p) {
	return T(0.5) * (1 - std::cos(2 * T(M_PI) * p));
}

/** Multiplies the Hann window by a signal `x` of length `len` in-place. */
inline void hannWindow(float* x, int len) {
	for (int i = 0; i < len; i++) {
		x[i] *= hann(float(i) / (len - 1));
	}
}

/** Blackman window function.
https://en.wikipedia.org/wiki/Window_function#Blackman_window
A typical alpha value is 0.16.
*/
template <typename T>
inline T blackman(T alpha, T p) {
	return
	  + (1 - alpha) / 2
	  - T(1) / 2 * std::cos(2 * T(M_PI) * p)
	  + alpha / 2 * std::cos(4 * T(M_PI) * p);
}

inline void blackmanWindow(float alpha, float* x, int len) {
	for (int i = 0; i < len; i++) {
		x[i] *= blackman(alpha, float(i) / (len - 1));
	}
}


/** Blackman-Nuttall window function.
https://en.wikipedia.org/wiki/Window_function#Blackman%E2%80%93Nuttall_window
*/
template <typename T>
inline T blackmanNuttall(T p) {
	return
	  + T(0.3635819)
	  - T(0.4891775) * std::cos(2 * T(M_PI) * p)
	  + T(0.1365995) * std::cos(4 * T(M_PI) * p)
	  - T(0.0106411) * std::cos(6 * T(M_PI) * p);
}

inline void blackmanNuttallWindow(float* x, int len) {
	for (int i = 0; i < len; i++) {
		x[i] *= blackmanNuttall(float(i) / (len - 1));
	}
}

/** Blackman-Harris window function.
https://en.wikipedia.org/wiki/Window_function#Blackman%E2%80%93Harris_window
*/
template <typename T>
inline T blackmanHarris(T p) {
	return
	  + T(0.35875)
	  - T(0.48829) * std::cos(2 * T(M_PI) * p)
	  + T(0.14128) * std::cos(4 * T(M_PI) * p)
	  - T(0.01168) * std::cos(6 * T(M_PI) * p);
}

inline void blackmanHarrisWindow(float* x, int len) {
	for (int i = 0; i < len; i++) {
		x[i] *= blackmanHarris(float(i) / (len - 1));
	}
}

