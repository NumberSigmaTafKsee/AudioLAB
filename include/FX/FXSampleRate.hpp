#pragma once

#include "FXObjects.hpp"

// --- sample rate conversion
//
// --- supported conversion ratios - you can EASILY add more to this
/**
\enum rateConversionRatio
\ingroup Constants-Enums
\brief
Use this strongly typed enum to easily set up or down sampling ratios.
- enum class rateConversionRatio { k2x, k4x };
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
enum class rateConversionRatio { k2x, k4x };
const unsigned int maxSamplingRatio = 4;

/**
@countForRatio
\ingroup FX-Functions
@brief returns the up or downsample ratio as a numeric value
\param ratio - enum class ratio value
\return the up or downsample ratio as a numeric value
*/
inline unsigned int countForRatio(rateConversionRatio ratio)
{
	if (ratio == rateConversionRatio::k2x)
		return 2;
	else if (ratio == rateConversionRatio::k4x || ratio == rateConversionRatio::k4x)
		return 4;

	return 0;
}

// --- get table pointer for built-in anti-aliasing LPFs
/**
@getFilterIRTable
\ingroup FX-Functions
@brief returns the up or downsample ratio as a numeric value
\param FIRLength - lenght of FIR
\param ratio - the conversinon ratio
\param sampleRate - the sample rate
\return a pointer to the appropriate FIR coefficient table in filters.h or nullptr if not found
*/
inline double* getFilterIRTable(unsigned int FIRLength, rateConversionRatio ratio, unsigned int sampleRate)
{
	// --- we only have built in filters for 44.1 and 48 kHz
	if (sampleRate != 44100 && sampleRate != 48000) return nullptr;

	// --- choose 2xtable
	if (ratio == rateConversionRatio::k2x)
	{
		if (sampleRate == 44100)
		{
			if (FIRLength == 128)
				return &LPF128_882[0];
			else if (FIRLength == 256)
				return &LPF256_882[0];
			else if (FIRLength == 512)
				return &LPF512_882[0];
			else if (FIRLength == 1024)
				return &LPF1024_882[0];
		}
		if (sampleRate == 48000)
		{
			if (FIRLength == 128)
				return &LPF128_96[0];
			else if (FIRLength == 256)
				return &LPF256_96[0];
			else if (FIRLength == 512)
				return &LPF512_96[0];
			else if (FIRLength == 1024)
				return &LPF1024_96[0];
		}
	}

	// --- choose 4xtable
	if (ratio == rateConversionRatio::k4x)
	{
		if (sampleRate == 44100)
		{
			if (FIRLength == 128)
				return &LPF128_1764[0];
			else if (FIRLength == 256)
				return &LPF256_1764[0];
			else if (FIRLength == 512)
				return &LPF512_1764[0];
			else if (FIRLength == 1024)
				return &LPF1024_1764[0];
		}
		if (sampleRate == 48000)
		{
			if (FIRLength == 128)
				return &LPF128_192[0];
			else if (FIRLength == 256)
				return &LPF256_192[0];
			else if (FIRLength == 512)
				return &LPF512_192[0];
			else if (FIRLength == 1024)
				return &LPF1024_192[0];
		}
	}
	return nullptr;
}

// --- get table pointer for built-in anti-aliasing LPFs
/**
@decomposeFilter
\ingroup FX-Functions
@brief performs a polyphase decomposition on a big FIR into a set of sub-band FIRs
\param filterIR - pointer to filter IR array
\param FIRLength - lenght of IR array
\param ratio - up or down sampling ratio
\return a pointer an arry of buffer pointers to the decomposed mini-filters
*/
inline double** decomposeFilter(double* filterIR, unsigned int FIRLength, unsigned int ratio)
{
	unsigned int subBandLength = FIRLength / ratio;
	double ** polyFilterSet = new double*[ratio];
	#pragma omp simd
	for (unsigned int i = 0; i < ratio; i++)
	{
		double* polyFilter = new double[subBandLength];
		polyFilterSet[i] = polyFilter;
	}

	int m = 0;
	#pragma omp simd
	for (int i = 0; i < subBandLength; i++)
	{
		for (int j = ratio - 1; j >= 0; j--)
		{
			double* polyFilter = polyFilterSet[j];
			polyFilter[i] = filterIR[m++];
		}
	}

	return polyFilterSet;
}

/**
\struct InterpolatorOutput
\ingroup FFTW-Objects
\brief
Custom output structure for interpolator; it holds an arry of interpolated output samples.
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
struct InterpolatorOutput
{
	InterpolatorOutput() {}
	double audioData[maxSamplingRatio] = { 0.0, 0.0, 0.0, 0.0 };	///< array of interpolated output samples
	unsigned int count = maxSamplingRatio;			///< number of samples in output array
};


/**
\class Interpolator
\ingroup FFTW-Objects
\brief
The Interpolator object implements a sample rate interpolator. One input sample yields N output samples.
Audio I/O:
- Processes mono input to interpoalted (multi-sample) output.
Control I/F:
- none.
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
class Interpolator
{
public:
	Interpolator() { }		/* C-TOR */
	~Interpolator() { }		/* D-TOR */

	/** setup the sample rate interpolator */
	/**
	\param _FIRLength the interpolator's anti-aliasing filter length
	\param _ratio the conversion ratio (see rateConversionRatio)
	\param _sampleRate the actual sample rate
	\param _polyphase flag to enable polyphase decomposition
	*/
	inline void initialize(unsigned int _FIRLength, rateConversionRatio _ratio, unsigned int _sampleRate, bool _polyphase = true)
	{
		polyphase = _polyphase;
		sampleRate = _sampleRate;
		FIRLength = _FIRLength;
		ratio = _ratio;
		unsigned int count = countForRatio(ratio);
		unsigned int subBandLength = FIRLength / count;

		// --- straight SRC, no polyphase
		convolver.initialize(FIRLength);

		// --- set filterIR from built-in set - user can always override this!
		double* filterTable = getFilterIRTable(FIRLength, ratio, sampleRate);
		if (!filterTable) return;
		convolver.setFilterIR(filterTable);

		if (!polyphase) return;

		// --- decompose filter
		double** polyPhaseFilters = decomposeFilter(filterTable, FIRLength, count);
		if (!polyPhaseFilters)
		{
			polyphase = false;
			return;
		}

		// --- set the individual polyphase filter IRs on the convolvers
		#pragma omp simd
		for (unsigned int i = 0; i < count; i++)
		{
			polyPhaseConvolvers[i].initialize(subBandLength);
			polyPhaseConvolvers[i].setFilterIR(polyPhaseFilters[i]);
			delete[] polyPhaseFilters[i];
		}

		delete[] polyPhaseFilters;
	}

	/** perform the interpolation; the multiple outputs are in an array in the return structure */
	inline InterpolatorOutput interpolateAudio(double xn)
	{
		unsigned int count = countForRatio(ratio);

		// --- setup output
		InterpolatorOutput output;
		output.count = count;

		// --- interpolators need the amp correction
		double ampCorrection = double(count);

		// --- polyphase uses "backwards" indexing for interpolator; see book
		int m = count-1;
		#pragma omp simd
		for (unsigned int i = 0; i < count; i++)
		{
			if (!polyphase)
				output.audioData[i] = i == 0 ? ampCorrection*convolver.processAudioSample(xn) : ampCorrection*convolver.processAudioSample(0.0);
			else
				output.audioData[i] = ampCorrection*polyPhaseConvolvers[m--].processAudioSample(xn);
		}
		return output;
	}

protected:
	// --- for straight, non-polyphase
	FastConvolver convolver; ///< the convolver

	// --- we save these for future expansion, currently only sparsely used
	unsigned int sampleRate = 44100;	///< sample rate
	unsigned int FIRLength = 256;		///< FIR length
	rateConversionRatio ratio = rateConversionRatio::k2x; ///< conversion ration

	// --- polyphase: 4x is max right now
	bool polyphase = true;									///< enable polyphase decomposition
	FastConvolver polyPhaseConvolvers[maxSamplingRatio];	///< a set of sub-band convolvers for polyphase operation
};

/**
\struct DecimatorInput
\ingroup FFTW-Objects
\brief
Custom input structure for DecimatorInput; it holds an arry of input samples that will be decimated down to just one sample.
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
struct DecimatorInput
{
	DecimatorInput() {}
	double audioData[maxSamplingRatio] = { 0.0, 0.0, 0.0, 0.0 };	///< input array of samples to be decimated
	unsigned int count = maxSamplingRatio;			///< count of samples in input array
};

/**
\class Decimator
\ingroup FFTW-Objects
\brief
The Decimator object implements a sample rate decimator. Ana array of M input samples is decimated
to one output sample.
Audio I/O:
- Processes mono input to interpoalted (multi-sample) output.
Control I/F:
- none.
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
class Decimator
{
public:
	Decimator() { }		/* C-TOR */
	~Decimator() { }	/* D-TOR */

	/** setup the sample rate decimator */
	/**
	\param _FIRLength the decimator's anti-aliasing filter length
	\param _ratio the conversion ratio (see rateConversionRatio)
	\param _sampleRate the actual sample rate
	\param _polyphase flag to enable polyphase decomposition
	*/
	inline void initialize(unsigned int _FIRLength, rateConversionRatio _ratio, unsigned int _sampleRate, bool _polyphase = true)
	{
		polyphase = _polyphase;
		sampleRate = _sampleRate;
		FIRLength = _FIRLength;
		ratio = _ratio;
		unsigned int count = countForRatio(ratio);
		unsigned int subBandLength = FIRLength / count;

		// --- straight SRC, no polyphase
		convolver.initialize(FIRLength);

		// --- set filterIR from built-in set - user can always override this!
		double* filterTable = getFilterIRTable(FIRLength, ratio, sampleRate);
		if (!filterTable) return;
		convolver.setFilterIR(filterTable);

		if (!polyphase) return;

		// --- decompose filter
		double** polyPhaseFilters = decomposeFilter(filterTable, FIRLength, count);
		if (!polyPhaseFilters)
		{
			polyphase = false;
			return;
		}

		// --- set the individual polyphase filter IRs on the convolvers
		#pragma omp simd
		for (unsigned int i = 0; i < count; i++)
		{
			polyPhaseConvolvers[i].initialize(subBandLength);
			polyPhaseConvolvers[i].setFilterIR(polyPhaseFilters[i]);
			delete[] polyPhaseFilters[i];
		}

		delete[] polyPhaseFilters;
	}

	/** decimate audio input samples into one outut sample (return value) */
	inline double decimateAudio(DecimatorInput data)
	{
		unsigned int count = countForRatio(ratio);

		// --- setup output
		double output = 0.0;

		// --- polyphase uses "forwards" indexing for decimator; see book
		#pragma omp simd
		for (unsigned int i = 0; i < count; i++)
		{
			if (!polyphase) // overwrites output; only the last output is saved
				output = convolver.processAudioSample(data.audioData[i]);
			else
				output += polyPhaseConvolvers[i].processAudioSample(data.audioData[i]);
		}
		return output;
	}

protected:
	// --- for straight, non-polyphase
	FastConvolver convolver;		 ///< fast convolver

	// --- we save these for future expansion, currently only sparsely used
	unsigned int sampleRate = 44100;	///< sample rate
	unsigned int FIRLength = 256;		///< FIR length
	rateConversionRatio ratio = rateConversionRatio::k2x; ///< conversion ration

	// --- polyphase: 4x is max right now
	bool polyphase = true;									///< enable polyphase decomposition
	FastConvolver polyPhaseConvolvers[maxSamplingRatio];	///< a set of sub-band convolvers for polyphase operation
};
