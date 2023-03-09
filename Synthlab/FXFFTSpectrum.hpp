#pragma once

#include "FXObjects.hpp"

// ------------------------------------------------------------------ //
// --- OBJECTS REQUIRING FFTW --------------------------------------- //
// ------------------------------------------------------------------ //

/**
\enum windowType
\ingroup Constants-Enums
\brief
Use this strongly typed enum to easily set the windowing type for FFT algorithms that use it.
- enum class windowType {kNoWindow, kRectWindow, kHannWindow, kBlackmanHarrisWindow, kHammingWindow };
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
enum class windowType {kNoWindow, kRectWindow, kHannWindow, kBlackmanHarrisWindow, kHammingWindow };

/**
@makeWindow
\ingroup FX-Functions
@brief  creates a new std::unique_ptr<double[]> array for a given window lenght and type.
\param windowLength - length of window array (does NOT need to be power of 2)
\param hopSize - hopSize for vococerf applications, may set to 0 for non vocoder use
\param gainCorrectionValue - return variable that contains the window gain correction value
*/
inline std::unique_ptr<double[]> makeWindow(unsigned int windowLength, unsigned int hopSize, windowType window, double& gainCorrectionValue)
{
	std::unique_ptr<double[]> windowBuffer;
	windowBuffer.reset(new double[windowLength]);

	if (!windowBuffer) return nullptr;

	double overlap = hopSize > 0.0 ? 1.0 - (double)hopSize / (double)windowLength : 0.0;
	gainCorrectionValue = 0.0;

	for (int n = 0; n < windowLength; n++)
	{
		if (window == windowType::kRectWindow)
		{
			if (n >= 1 && n <= windowLength - 1)
				windowBuffer[n] = 1.0;
		}
		else if (window == windowType::kHammingWindow)
		{
			windowBuffer[n] = 0.54 - 0.46*cos((n*2.0*kPi) / (windowLength));
		}
		else if (window == windowType::kHannWindow)
		{
			windowBuffer[n] = 0.5 * (1 - cos((n*2.0*kPi) / (windowLength)));
		}
		else if (window == windowType::kBlackmanHarrisWindow)
		{
			windowBuffer[n] = (0.42323 - (0.49755*cos((n*2.0*kPi) / (windowLength))) + 0.07922*cos((2 * n*2.0*kPi) / (windowLength)));
		}
		else if (window == windowType::kNoWindow)
		{
			windowBuffer[n] = 1.0;
		}

		gainCorrectionValue += windowBuffer[n];
	}

	// --- calculate gain correction factor
	if (window != windowType::kNoWindow)
		gainCorrectionValue = (1.0 - overlap) / gainCorrectionValue;
	else
		gainCorrectionValue = 1.0 / gainCorrectionValue;

	return windowBuffer;
}

// --- FFTW --- to enable, add the statement #define HAVE_FFTW 1 to the top of the file



/**
\class FastFFT
\ingroup FFTW-Objects
\brief
The FastFFT provides a simple wrapper for the FFTW FFT operation - it is ultra-thin and simple to use.
Audio I/O:
- processes mono inputs into FFT outputs.
Control I/F:
- none.
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
class FastFFT
{
public:
	FastFFT() {}		/* C-TOR */
	~FastFFT() {
		if (windowBuffer) delete[] windowBuffer;
		destroyFFTW();
	}	/* D-TOR */

	/** setup the FFT for a given framelength and window type*/
	void initialize(unsigned int _frameLength, windowType _window);

	/** destroy FFTW objects and plans */
	void destroyFFTW();

	/** do the FFT and return real and imaginary arrays */
	fftw_complex* doFFT(double* inputReal, double* inputImag = nullptr);

	/** do the IFFT and return real and imaginary arrays */
	fftw_complex* doInverseFFT(double* inputReal, double* inputImag);

	/** get the current FFT length */
	unsigned int getFrameLength() { return frameLength; }

protected:
	// --- setup FFTW
	fftw_complex*	fft_input = nullptr;		///< array for FFT input
	fftw_complex*	fft_result = nullptr;		///< array for FFT output
	fftw_complex*	ifft_input = nullptr;		///< array for IFFT input
	fftw_complex*	ifft_result = nullptr;		///< array for IFFT output
	fftw_plan       plan_forward = nullptr;		///< FFTW plan for FFT
	fftw_plan		plan_backward = nullptr;	///< FFTW plan for IFFT

	double* windowBuffer = nullptr;				///< buffer for window (naked)
	double windowGainCorrection = 1.0;			///< window gain correction
	windowType window = windowType::kHannWindow; ///< window type
	unsigned int frameLength = 0;				///< current FFT length
};


/**
\brief destroys the FFTW arrays and plans.
*/
void FastFFT::destroyFFTW()
{
#ifdef HAVE_FFTW
	if (plan_forward)
		fftw_destroy_plan(plan_forward);
	if (plan_backward)
		fftw_destroy_plan(plan_backward);

	if (fft_input)
		fftw_free(fft_input);
	if (fft_result)
		fftw_free(fft_result);

	if (ifft_input)
		fftw_free(ifft_input);
	if (ifft_result)
		fftw_free(ifft_result);
#endif
}


/**
\brief initialize the Fast FFT object for operation
- NOTES:<br>
See notes on symmetrical window arrays in comments.<br>
\param _frameLength the FFT length - MUST be a power of 2
\param _window the window type (note: may be set to windowType::kNone)
*/
void FastFFT::initialize(unsigned int _frameLength, windowType _window)
{
	frameLength = _frameLength;
	window = _window;
	windowGainCorrection = 0.0;

	if (windowBuffer)
		delete windowBuffer;

	windowBuffer = new double[frameLength];
	memset(&windowBuffer[0], 0, frameLength * sizeof(double));


	// --- this is from Reiss & McPherson's code
	//     https://code.soundsoftware.ac.uk/projects/audio_effects_textbook_code/repository/entry/effects/pvoc_passthrough/Source/PluginProcessor.cpp
	// NOTE:	"Window functions are typically defined to be symmetrical. This will cause a
	//			problem in the overlap-add process: the windows instead need to be periodic
	//			when arranged end-to-end. As a result we calculate the window of one sample
	//			larger than usual, and drop the last sample. (This works as long as N is even.)
	//			See Julius Smith, "Spectral Audio Signal Processing" for details.
	// --- WP: this is why denominators are (frameLength) rather than (frameLength - 1)
	if (window == windowType::kRectWindow)
	{
		for (int n = 0; n < frameLength - 1; n++)
		{
			windowBuffer[n] = 1.0;
			windowGainCorrection += windowBuffer[n];
		}
	}
	else if (window == windowType::kHammingWindow)
	{
		for (int n = 0; n < frameLength - 1; n++)
		{
			windowBuffer[n] = 0.54 - 0.46*cos((n*2.0*kPi) / (frameLength));
			windowGainCorrection += windowBuffer[n];
		}
	}
	else if (window == windowType::kHannWindow)
	{
		for (int n = 0; n < frameLength; n++)
		{
			windowBuffer[n] = 0.5 * (1 - cos((n*2.0*kPi) / (frameLength)));
			windowGainCorrection += windowBuffer[n];
		}
	}
	else if (window == windowType::kBlackmanHarrisWindow)
	{
		for (int n = 0; n < frameLength; n++)
		{
			windowBuffer[n] = (0.42323 - (0.49755*cos((n*2.0*kPi) / (frameLength))) + 0.07922*cos((2 * n*2.0*kPi) / (frameLength)));
			windowGainCorrection += windowBuffer[n];
		}
	}
	else if (window == windowType::kNoWindow)
	{
		for (int n = 0; n < frameLength; n++)
		{
			windowBuffer[n] = 1.0;
			windowGainCorrection += windowBuffer[n];
		}
	}
	else // --- default to kNoWindow
	{
		for (int n = 0; n < frameLength; n++)
		{
			windowBuffer[n] = 1.0;
			windowGainCorrection += windowBuffer[n];
		}
	}

	// --- calculate gain correction factor
	windowGainCorrection = 1.0 / windowGainCorrection;

	destroyFFTW();
	fft_input = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * frameLength);
	fft_result = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * frameLength);

	ifft_input =  (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * frameLength);
	ifft_result = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * frameLength);

	plan_forward = fftw_plan_dft_1d(frameLength, fft_input, fft_result, FFTW_FORWARD, FFTW_ESTIMATE);
	plan_backward = fftw_plan_dft_1d(frameLength, ifft_input, ifft_result, FFTW_BACKWARD, FFTW_ESTIMATE);
}

/**
\brief perform the FFT operation
- NOTES:<br>
\param inputReal an array of real valued points
\param inputImag an array of imaginary valued points (will be 0 for audio which is real-valued)
\returns a pointer to a fftw_complex array: a 2D array of real (column 0) and imaginary (column 1) parts
*/
fftw_complex* FastFFT::doFFT(double* inputReal, double* inputImag)
{
	// ------ load up the FFT input array
	for (int i = 0; i < frameLength; i++)
	{
		fft_input[i][0] = inputReal[i];		// --- real
		if (inputImag)
			fft_input[i][1] = inputImag[i]; // --- imag
		else
			fft_input[i][1] = 0.0;
	}

	// --- do the FFT
	fftw_execute(plan_forward);

	return fft_result;
}

/**
\brief perform the IFFT operation
- NOTES:<br>
\param inputReal an array of real valued points
\param inputImag an array of imaginary valued points (will be 0 for audio which is real-valued)
\returns a pointer to a fftw_complex array: a 2D array of real (column 0) and imaginary (column 1) parts
*/
fftw_complex* FastFFT::doInverseFFT(double* inputReal, double* inputImag)
{
	// ------ load up the iFFT input array
	for (int i = 0; i < frameLength; i++)
	{
		ifft_input[i][0] = inputReal[i];		// --- real
		if (inputImag)
			ifft_input[i][1] = inputImag[i]; // --- imag
		else
			ifft_input[i][1] = 0.0;
	}

	// --- do the IFFT
	fftw_execute(plan_backward);

	return ifft_result;
}
