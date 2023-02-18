#pragma once

#include <cmath>
#include <array>
#include <vector>
#include <cstdint>
#include <exception>
#include <array>
#include <random>
#include <cassert>
#include <algorithm>
#include <cstring>


#include "Undenormal.hpp"
#include "ClipFunctions.hpp"

#define MOOG_E         2.71828182845904523536028747135266250
#define MOOG_LOG2E     1.44269504088896340735992468100189214
#define MOOG_LOG10E    0.434294481903251827651128918916605082
#define MOOG_LN2       0.693147180559945309417232121458176568
#define MOOG_LN10      2.30258509299404568401799145468436421
#define MOOG_PI        3.14159265358979323846264338327950288
#define MOOG_PI_2      1.57079632679489661923132169163975144
#define MOOG_PI_4      0.785398163397448309615660845819875721
#define MOOG_1_PI      0.318309886183790671537767526745028724
#define MOOG_2_PI      0.636619772367581343075535053490057448
#define MOOG_2_SQRTPI  1.12837916709551257389615890312154517
#define MOOG_SQRT2     1.41421356237309504880168872420969808
#define MOOG_SQRT1_2   0.707106781186547524400844362104849039
#define MOOG_INV_PI_2  0.159154943091895

#define NO_COPY(C) C(const C &) = delete; C & operator = (const C &) = delete
#define NO_MOVE(C) NO_COPY(C); C(C &&) = delete; C & operator = (const C &&) = delete

#define SNAP_TO_ZERO(n)    if (! (n < -1.0e-8 || n > 1.0e-8)) n = 0;

namespace Analog::Moog
{
	// Linear interpolation, used to crossfade a gain table
	inline DspFloatType moog_lerp(DspFloatType amount, DspFloatType a, DspFloatType b)
	{
		return (1.0f - amount) * a + amount * b;
	}

	inline DspFloatType moog_min(DspFloatType a, DspFloatType b)
	{
		a = b - a;
		a += fabs(a);
		a *= 0.5f;
		a = b - a;
		return a;
	}

	// Clamp without branching
	// If input - _limit < 0, then it really substracts, and the 0.5 to make it half the 2 inputs.
	// If > 0 then they just cancel, and keeps input normal.
	// The easiest way to understand it is check what happends on both cases.
	inline DspFloatType moog_saturate(DspFloatType input)
	{
		DspFloatType x1 = fabs(input + 0.95f);
		DspFloatType x2 = fabs(input - 0.95f);
		return 0.5f * (x1 - x2);
	}

	// Imitate the (tanh) clipping function of a transistor pair.
	// to 4th order, tanh is x - x*x*x/3; this cubic's
	// plateaus are at +/- 1 so clip to 1 and evaluate the cubic.
	// This is pretty coarse - for instance if you clip a sinusoid this way you
	// can sometimes hear the discontinuity in 4th derivative at the clip point
	inline DspFloatType clip(DspFloatType value, DspFloatType saturation, DspFloatType saturationinverse)
	{
		DspFloatType v2 = (value * saturationinverse > 1 ? 1 :
					(value * saturationinverse < -1 ? -1:
					 value * saturationinverse));
		return (saturation * (v2 - (1./3.) * v2 * v2 * v2));
	}

	#define HZ_TO_RAD(f) (MOOG_PI_2 * f)
	#define RAD_TO_HZ(omega) (MOOG_INV_PI_2 * omega)

	#ifdef __GNUC__
		#define ctz(N) __builtin_ctz(N)
	#else
		template<typename T>
		inline int ctz(T x)
		{
			int p, b;
			for (p = 0, b = 1; !(b & x); b <<= 1, ++p)
				;
			return p;
		}
	#endif

	inline DspFloatType fast_tanh(DspFloatType x) 
	{
		DspFloatType x2 = x * x;
		return x * (27.0 + x2) / (27.0 + 9.0 * x2);
	}

	class BiQuadBase : public FilterProcessor
	{
	public:
		
		BiQuadBase() : FilterProcessor()
		{
			bCoef = {{0.0f, 0.0f, 0.0f}};
			aCoef = {{0.0f, 0.0f}};
			w = {{0.0f, 0.0f}};
		}	
		~BiQuadBase()
		{

		}	
		// DF-II impl	
		void ProcessSIMD(uint32_t n, DspFloatType *in, DspFloatType * output)
		{
			Undenormal denormal;
			DspFloatType out = 0;
			#pragma omp simd aligned(in,output)
			for (int s = 0; s < n; ++s)
			{
				out = bCoef[0]  * in[s] + w[0];
				w[0] = bCoef[1] * in[s] - aCoef[0] * out + w[1];
				w[1] = bCoef[2] * in[s] - aCoef[1] * out;
				output[s] = out;
			}
		}
		void ProcessBlock(size_t n, DspFloatType * in, DspFloatType * out) {
			ProcessSIMD(n,in,out);
		}
		void ProcessInplace(size_t n, DspFloatType * out) {
			ProcessSIMD(n,out,out);
		}
		DspFloatType Tick(DspFloatType s, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1)
		{
			Undenormal denormal;
			DspFloatType out = bCoef[0] * s + w[0];
			w[0] = bCoef[1] * s - aCoef[0] * out + w[1];
			w[1] = bCoef[2] * s - aCoef[1] * out;
			return out;
		}
		void SetBiquadCoefs(std::array<DspFloatType, 3> b, std::array<DspFloatType, 2> a)
		{
			bCoef = b;
			aCoef = a;
		}
		
	protected:
		std::array<DspFloatType, 3> bCoef; // b0, b1, b2
		std::array<DspFloatType, 2> aCoef; // a1, a2
		std::array<DspFloatType, 2> w; // delays
	};


	/*! \brief implement a circular buffer of type T
	*/
	template <class T>
	class CRingBuffer
	{
	public:
		explicit CRingBuffer(int iBufferLengthInSamples) :
			m_iBuffLength(iBufferLengthInSamples),
			m_iReadIdx(0),
			m_iWriteIdx(0),
			m_ptBuff(0)
		{
			assert(iBufferLengthInSamples > 0);

			m_ptBuff = new T[m_iBuffLength];
			reset();
		}

		virtual ~CRingBuffer()
		{
			delete[] m_ptBuff;
			m_ptBuff = 0;
		}

		/*! add a new value of type T to write index and increment write index
		\param tNewValue the new value
		\return void
		*/
		void putPostInc(T tNewValue)
		{
			put(tNewValue);
			incIdx(m_iWriteIdx);
		}

		/*! add new values of type T to write index and increment write index
		\param ptNewBuff: new values
		\param iLength: number of values
		\return void
		*/
		void putPostInc(const T* ptNewBuff, int iLength)
		{
			put(ptNewBuff, iLength);
			incIdx(m_iWriteIdx, iLength);
		}

		/*! add a new value of type T to write index
		\param tNewValue the new value
		\return void
		*/
		void put(T tNewValue)
		{
			m_ptBuff[m_iWriteIdx] = tNewValue;
		}

		/*! add new values of type T to write index
		\param ptNewBuff: new values
		\param iLength: number of values
		\return void
		*/
		void put(const T* ptNewBuff, int iLength)
		{
			assert(iLength <= m_iBuffLength && iLength >= 0);

			// copy two parts: to the end of buffer and after wrap around
			int iNumValues2End = std::min(iLength, m_iBuffLength - m_iWriteIdx);

			std::memcpy (&m_ptBuff[m_iWriteIdx], ptNewBuff, sizeof(T)*iNumValues2End);
			if ((iLength - iNumValues2End) > 0)
				std::memcpy (m_ptBuff, &ptNewBuff[iNumValues2End], sizeof(T)*(iLength - iNumValues2End));
		}

		/*! return the value at the current read index and increment the read pointer
		\return DspFloatType the value from the read index
		*/
		T getPostInc(DspFloatType fOffset = 0)
		{
			T tValue = get(fOffset);
			incIdx(m_iReadIdx);
			return tValue;
		}

		/*! return the values starting at the current read index and increment the read pointer
		\param ptBuff: pointer to where the values will be written
		\param iLength: number of values
		\return void
		*/
		void getPostInc(T* ptBuff, int iLength)
		{
			get(ptBuff, iLength);
			incIdx(m_iReadIdx, iLength);
		}

		/*! return the value at the current read index
		\param fOffset: read at offset from read index
		\return DspFloatType the value from the read index
		*/
		T get(DspFloatType fOffset = 0) const
		{
			if (fOffset == 0)
				return m_ptBuff[m_iReadIdx];
			else
			{

				// compute fraction for linear interpolation 
				int     iOffset = static_cast<int>(std::floor(fOffset));
				DspFloatType   fFrac = fOffset - iOffset;
				int     iRead = m_iReadIdx + iOffset;
				while (iRead > m_iBuffLength - 1)
					iRead -= m_iBuffLength;
				while (iRead < 0)
					iRead += m_iBuffLength;

				return (1 - fFrac) * m_ptBuff[iRead] +
					fFrac * m_ptBuff[(iRead + 1) % m_iBuffLength];
			}
		}

		/*! return the values starting at the current read index
		\param ptBuff to where the values will be written
		\param iLength: number of values
		\return void
		*/
		void get(T* ptBuff, int iLength) const
		{
			assert(iLength <= m_iBuffLength && iLength >= 0);

			// copy two parts: to the end of buffer and after wrap around
			int iNumValues2End = std::min(iLength, m_iBuffLength - m_iReadIdx);

			std::memcpy (ptBuff, &m_ptBuff[m_iReadIdx], sizeof(T)*iNumValues2End);
			if ((iLength - iNumValues2End)>0)
				std::memcpy (&ptBuff[iNumValues2End], m_ptBuff, sizeof(T)*(iLength - iNumValues2End));
		}

		T extractPostInc()
		{
			T value = get();
			m_ptBuff[m_iReadIdx] = 0;
			incIdx(m_iReadIdx);
			return value;
		}

		/*! set buffer content and indices to 0
		\return void
		*/
		void reset()
		{
			std::memset (m_ptBuff, 0, sizeof(T)*m_iBuffLength);
			m_iReadIdx  = 0;
			m_iWriteIdx = 0;
		}

		/*! return the current index for writing/put
		\return int
		*/
		int getWriteIdx() const
		{
			return m_iWriteIdx;
		}

		/*! move the write index to a new position
		\param iNewWriteIdx: new position
		\return void
		*/
		void setWriteIdx(int iNewWriteIdx)
		{
			incIdx(m_iWriteIdx, iNewWriteIdx - m_iWriteIdx);
		}

		/*! return the current index for reading/get
		\return int
		*/
		int getReadIdx() const
		{
			return m_iReadIdx;
		}

		/*! move the read index to a new position
		\param iNewReadIdx: new position
		\return void
		*/
		void setReadIdx(int iNewReadIdx)
		{
			incIdx(m_iReadIdx, iNewReadIdx - m_iReadIdx);
		}

		/*! returns the number of values currently buffered (note: 0 could also mean the buffer is full!)
		\return int
		*/
		int getNumValuesInBuffer() const
		{
			return (m_iWriteIdx - m_iReadIdx + m_iBuffLength) % m_iBuffLength;
		}

		/*! returns the length of the internal buffer
		\return int
		*/
		int getLength() const
		{
			return m_iBuffLength;
		}
	private:
		CRingBuffer();
		CRingBuffer(const CRingBuffer& that);

		void incIdx(int& iIdx, int iOffset = 1)
		{
			while ((iIdx + iOffset) < 0)
			{
				// avoid negative buffer indices
				iOffset += m_iBuffLength;
			}
			iIdx = (iIdx + iOffset) % m_iBuffLength;
		};

		int m_iBuffLength = 0,      //!< length of the internal buffer
			m_iReadIdx = 0,         //!< current read index
			m_iWriteIdx = 0;        //!< current write index

		T* m_ptBuff = 0;            //!< data buffer
	};


	///////////////////////////////////////////////////////////////////////////////////////////
	// Noise kaka
	///////////////////////////////////////////////////////////////////////////////////////////
	// +/-0.05dB above 9.2Hz @ 44,100Hz
	class PinkingFilter
	{
		DspFloatType b0, b1, b2, b3, b4, b5, b6;
	public:
		PinkingFilter() : b0(0), b1(0), b2(0), b3(0), b4(0), b5(0), b6(0) {}
		DspFloatType process(const DspFloatType s)
		{
			Undenormal denormal;
			b0 = 0.99886 * b0 + s * 0.0555179;
			b1 = 0.99332 * b1 + s * 0.0750759;
			b2 = 0.96900 * b2 + s * 0.1538520;
			b3 = 0.86650 * b3 + s * 0.3104856;
			b4 = 0.55000 * b4 + s * 0.5329522;
			b5 = -0.7616 * b5 - s * 0.0168980;
			const DspFloatType pink = (b0 + b1 + b2 + b3 + b4 + b5 + b6 + (s * 0.5362)) * 0.11;
			b6 = s * 0.115926;
			return pink;
		}
		void ProcessSIMD(size_t n, DspFloatType * in, DspFloatType * out)
		{
			#pragma omp simd aligned(in,out)
			for(size_t i = 0; i < n; i++) 
			{
				const DspFloatType s = in[i];
				b0 = 0.99886 * b0 + s * 0.0555179;
				b1 = 0.99332 * b1 + s * 0.0750759;
				b2 = 0.96900 * b2 + s * 0.1538520;
				b3 = 0.86650 * b3 + s * 0.3104856;
				b4 = 0.55000 * b4 + s * 0.5329522;
				b5 = -0.7616 * b5 - s * 0.0168980;
				const DspFloatType pink = (b0 + b1 + b2 + b3 + b4 + b5 + b6 + (s * 0.5362)) * 0.11;
				b6 = s * 0.115926;
				out[i] = pink;			
			}
		}
		void ProcessBlock(size_t n, DspFloatType * in, DspFloatType * out) {
			ProcessSIMD(n,in,out);
		}
		void ProcessInplace(size_t n, DspFloatType * out) {
			ProcessSIMD(n,nullptr,out);
		}
	};

	class BrowningFilter
	{
	DspFloatType l;
	public:
		BrowningFilter() : l(0) {}
		DspFloatType process(const DspFloatType s)
		{
			Undenormal denormal;
			DspFloatType brown = (l + (0.02f * s)) / 1.02f;
			l = brown;
			return brown * 3.5f; // compensate for gain
		}
		void ProcessSIMD(size_t n, DspFloatType * in, DspFloatType * out)
		{
			#pragma omp simd aligned(in,out)
			for(size_t i = 0; i < n; i++) 
			{
				const DspFloatType s = in[i];
				DspFloatType brown = (l + (0.02f * s)) / 1.02f;
				l = brown;
				out[i] =  brown * 3.5f; // compensate for gain		
			}
		}
		void ProcessBlock(size_t n, DspFloatType * in, DspFloatType * out) {
			ProcessSIMD(n,in,out);
		}
		void ProcessInplace(size_t n, DspFloatType * out) {
			ProcessSIMD(n,out,out);
		}
	};

	struct WhiteNoiseSource
	{
		WhiteNoiseSource() : dist(-1, 1) {}
		std::mt19937 engine;
		std::uniform_real_distribution<DspFloatType> dist;
	};

	// Full spectrum noise
	struct WhiteNoise : public WhiteNoiseSource
	{
		DspFloatType operator()() { return dist(engine); }
	};

	// Pink noise has a decrease of 3dB/Octave
	struct PinkNoise : public WhiteNoiseSource
	{
		DspFloatType operator()() { return f.process(dist(engine)); }
		PinkingFilter f;
	};

	 // Brown noise has a decrease of 6dB/Octave
	struct BrownNoise : public WhiteNoiseSource
	{
		DspFloatType operator()() { return f.process(dist(engine)); }
		BrowningFilter f;
	};

	// Note! This noise is only valid for 44100 because of the hard-coded filter coefficients
	struct NoiseGenerator
	{
		enum NoiseType
		{
			WHITE,
			PINK,
			BROWN,
		} noise_type = PINK;
		
		std::vector<DspFloatType> produce(NoiseType t, int sampleRate, int channels, DspFloatType seconds)
		{
			int samplesToGenerate = sampleRate * seconds * channels;
			std::vector<DspFloatType> samples;
			samples.resize(samplesToGenerate);
			
			switch (t)
			{
			case NoiseType::WHITE:
			{
				WhiteNoise n;
				for(int s = 0; s < samplesToGenerate; s++) samples[s] = n();
			} break;
			case NoiseType::PINK:
			{
				PinkNoise n;
				for(int s = 0; s < samplesToGenerate; s++) samples[s] = n();
			} break;
			case NoiseType::BROWN:
			{
				BrownNoise n;
				for(int s = 0; s < samplesToGenerate; s++) samples[s] = n();
			} break;
			default: throw std::runtime_error("Invalid noise type");
			}
			return samples;
		}
		
		DspFloatType Tick() {
			switch (noise_type)
			{
			case NoiseType::WHITE:
			{
				WhiteNoise n;
				return n();
			} break;
			case NoiseType::PINK:
			{
				PinkNoise n;
				return n();
			} break;
			case NoiseType::BROWN:
			{
				BrownNoise n;
				return n();
			} break;
			default: throw std::runtime_error("Invalid noise type");
			}
		}
		void ProcessBlock(size_t n, DspFloatType * input, DspFloatType * output) {
			#pragma omp simd aligned(input,output)
			for(size_t i = 0; i < n; i++) output[i] = input[i]*Tick();
		}
		void ProcessInplace(DspFloatType * samples,size_t n) {
			#pragma omp simd aligned(samples)
			for(size_t i = 0; i < n; i++) samples[i] = samples[i]*Tick();
		}	
	};


	///////////////////////////////////////////////////////////////////////////////////////////
	// Moog Ladder
	///////////////////////////////////////////////////////////////////////////////////////////
	class LadderFilterBase
	{
	public:

		LadderFilterBase(DspFloatType sampleRate) : sampleRate(sampleRate) {}
		virtual ~LadderFilterBase() {}

		virtual void ProcessBlock(size_t n, DspFloatType * in, DspFloatType * out) =0;
		virtual void ProcessInplace(size_t n, DspFloatType * out) {
			ProcessBlock(n,out,out);
		}

		virtual void SetResonance(DspFloatType r) = 0;
		virtual void SetCutoff(DspFloatType c) = 0;

		DspFloatType GetResonance() { return resonance; }
		DspFloatType GetCutoff() { return cutoff; }

	protected:

		DspFloatType cutoff;
		DspFloatType resonance;
		DspFloatType sampleRate;
	};

	enum NoiseType
	{
		WHITE,
		PINK,
		BROWN,
	};

	struct NoiseSamples
	{
		NoiseGenerator * noise;
		NoiseType type;
		int sampleRate;
		int channels;

		NoiseSamples(NoiseType type, int sampleRate, int channels)  {
			noise = new NoiseGenerator();
			assert(noise != nullptr);
		}
		~NoiseSamples() {
			if(noise) delete noise;
		}

		std::vector<DspFloatType> produce(DspFloatType seconds) {
			std::vector<DspFloatType> r;
			r = noise->produce((NoiseGenerator::NoiseType)type,sampleRate,channels,seconds);
			return r;
		}
		
	};

}

