#pragma once

#include "FXObjects.hpp"


/**
\class PhaseVocoder
\ingroup FFTW-Objects
\brief
The PhaseVocoder provides a basic phase vocoder that is initialized to N = 4096 and
75% overlap; the de-facto standard for PSM algorithms. The analysis and sythesis
hop sizes are identical.
Audio I/O:
- processes mono input into mono output.
Control I/F:
- none.
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
class PhaseVocoder
{
public:
	PhaseVocoder() {}		/* C-TOR */
	~PhaseVocoder() {
		if (inputBuffer) delete[] inputBuffer;
		if (outputBuffer) delete[] outputBuffer;
		if (windowBuffer) delete[] windowBuffer;
		destroyFFTW();
	}	/* D-TOR */

	/** setup the FFT for a given framelength and window type*/
	void initialize(unsigned int _frameLength, unsigned int _hopSize, windowType _window);

	/** destroy FFTW objects and plans */
	void destroyFFTW();

	/** process audio sample through vocode; check fftReady flag to access FFT output */
	DspFloatType processAudioSample(DspFloatType input, bool& fftReady);

	/** add zero-padding without advancing output read location, for fast convolution */
	bool addZeroPad(unsigned int count);

	/** increment the FFT counter and do the FFT if it is ready */
	bool advanceAndCheckFFT();

	/** get FFT data for manipulation (yes, naked pointer so you can manipulate) */
	fftw_complex* getFFTData() { return fft_result; }

	/** get IFFT data for manipulation (yes, naked pointer so you can manipulate) */
	fftw_complex* getIFFTData() { return ifft_result; }

	/** do the inverse FFT (optional; will be called automatically if not used) */
	void doInverseFFT();

	/** do the overlap-add operation */
	void doOverlapAdd(DspFloatType* outputData = nullptr, int length = 0);

	/** get current FFT length */
	unsigned int getFrameLength() { return frameLength; }

	/** get current hop size ha = hs */
	unsigned int getHopSize() { return hopSize; }

	/** get current overlap as a raw value (75% = 0.75) */
	DspFloatType getOverlap() { return overlap; }

	/** set the vocoder for overlap add only without hop-size */
	// --- for fast convolution and other overlap-add algorithms
	//     that are not hop-size dependent
	void setOverlapAddOnly(bool b){ bool overlapAddOnly = b; }

protected:
	// --- setup FFTW
	fftw_complex*	fft_input = nullptr;		///< array for FFT input
	fftw_complex*	fft_result = nullptr;		///< array for FFT output
	fftw_complex*	ifft_result = nullptr;		///< array for IFFT output
	fftw_plan       plan_forward = nullptr;		///< FFTW plan for FFT
	fftw_plan		plan_backward = nullptr;	///< FFTW plan for IFFT

	// --- linear buffer for window
	DspFloatType*			windowBuffer = nullptr;		///< array for window

	// --- circular buffers for input and output
	DspFloatType*			inputBuffer = nullptr;		///< input timeline (x)
	DspFloatType*			outputBuffer = nullptr;		///< output timeline (y)

	// --- index and wrap masks for input and output buffers
	unsigned int inputWriteIndex = 0;			///< circular buffer index: input write
	unsigned int outputWriteIndex = 0;			///< circular buffer index: output write
	unsigned int inputReadIndex = 0;			///< circular buffer index: input read
	unsigned int outputReadIndex = 0;			///< circular buffer index: output read
	unsigned int wrapMask = 0;					///< input wrap mask
	unsigned int wrapMaskOut = 0;				///< output wrap mask

	// --- amplitude correction factor, aking into account both hopsize (overlap)
	//     and the window power itself
	DspFloatType windowHopCorrection = 1.0;			///< window correction including hop/overlap

	// --- these allow a more robust combination of user interaction
	bool needInverseFFT = false;				///< internal flag to signal IFFT required
	bool needOverlapAdd = false;				///< internal flag to signal overlap/add required

	// --- our window type; you can add more windows if you like
	windowType window = windowType::kHannWindow;///< window type

	// --- counters
	unsigned int frameLength = 0;				///< current FFT length
	unsigned int fftCounter = 0;				///< FFT sample counter

	// --- hop-size and overlap (mathematically related)
	unsigned int hopSize = 0;					///< hop: ha = hs
	DspFloatType overlap = 1.0;						///< overlap as raw value (75% = 0.75)

	// --- flag for overlap-add algorithms that do not involve hop-size, other
	//     than setting the overlap
	bool overlapAddOnly = false;				///< flag for overlap-add-only algorithms

};



// --- PSM Vocoder
const unsigned int PSM_FFT_LEN = 4096;

/**
\struct BinData
\ingroup FFTW-Objects
\brief
Custom structure that holds information about each FFT bin. This includes all information
needed to perform pitch shifting and time stretching with phase locking (optional)
and peak tracking (optional).
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
struct BinData
{
	BinData() {}
	/** all FXObjects parameter objects require overloaded= operator so remember to add new entries if you add new variables. */
	BinData& operator=(const BinData& params)	// need this override for collections to work
	{
		if (this == &params)
			return *this;

		isPeak = params.isPeak;
		magnitude = params.magnitude;
		phi = params.phi;

		psi = params.psi;
		localPeakBin = params.localPeakBin;
		previousPeakBin = params.previousPeakBin;
		updatedPhase = params.updatedPhase;

		return *this;
	}

	/** reset all variables to 0.0 */
	void reset()
	{
		isPeak = false;
		magnitude = 0.0;
		phi = 0.0;

		psi = 0.0;
		localPeakBin = 0;
		previousPeakBin = -1; // -1 is flag
		updatedPhase = 0.0;
	}

	bool isPeak = false;	///< flag for peak bins
	DspFloatType magnitude = 0.0; ///< bin magnitude angle
	DspFloatType phi = 0.0;		///< bin phase angle
	DspFloatType psi = 0.0;		///< bin phase correction
	unsigned int localPeakBin = 0; ///< index of peak-boss
	int previousPeakBin = -1; ///< index of peak bin in previous FFT
	DspFloatType updatedPhase = 0.0; ///< phase update value
};

/**
\struct PSMVocoderParameters
\ingroup FFTW-Objects
\brief
Custom parameter structure for the Biquad object. Default version defines the biquad structure used in the calculation.
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
struct PSMVocoderParameters
{
	PSMVocoderParameters() {}
	/** all FXObjects parameter objects require overloaded= operator so remember to add new entries if you add new variables. */
	PSMVocoderParameters& operator=(const PSMVocoderParameters& params)	// need this override for collections to work
	{
		if (this == &params)
			return *this;

		pitchShiftSemitones = params.pitchShiftSemitones;
		enablePeakPhaseLocking = params.enablePeakPhaseLocking;
		enablePeakTracking = params.enablePeakTracking;

		return *this;
	}

	// --- params
	DspFloatType pitchShiftSemitones = 0.0;	///< pitch shift in half-steps
	bool enablePeakPhaseLocking = false;///< flag to enable phase lock
	bool enablePeakTracking = false;	///< flag to enable peak tracking
};

/**
\class PSMVocoder
\ingroup FFTW-Objects
\brief
The PSMVocoder object implements a phase vocoder pitch shifter. Phase locking and peak tracking
are optional.
Audio I/O:
- Processes mono input to mono output.
Control I/F:
- Use PSMVocoderParameters structure to get/set object params.
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
class PSMVocoder : public IAudioSignalProcessor
{
public:
	PSMVocoder() {
		vocoder.initialize(PSM_FFT_LEN, PSM_FFT_LEN/4, windowType::kHannWindow);  // 75% overlap
	}		/* C-TOR */
	~PSMVocoder() {
		if (windowBuff) delete[] windowBuff;
		if (outputBuff) delete[] outputBuff;

	}	/* D-TOR */

	/** reset members to initialized state */
	virtual bool reset(DspFloatType _sampleRate)
	{
		memset(&phi[0], 0, sizeof(DspFloatType)*PSM_FFT_LEN);
		memset(&psi[0], 0, sizeof(DspFloatType)* PSM_FFT_LEN);
		if(outputBuff)
			memset(outputBuff, 0, sizeof(DspFloatType)*outputBufferLength);

		for (int i = 0; i < PSM_FFT_LEN; i++)
		{
			binData[i].reset();
			binDataPrevious[i].reset();

			peakBins[i] = -1;
			peakBinsPrevious[i] = -1;
		}

		return true;
	}

	/** return false: this object only processes samples */
	virtual bool canProcessAudioFrame() { return false; }

	/** set the pitch shift in semitones (note that this can be fractional too)*/
	void setPitchShift(DspFloatType semitones)
	{
		// --- this is costly so only update when things changed
		DspFloatType newAlpha = pow(2.0, semitones / 12.0);
		DspFloatType newOutputBufferLength = round((1.0/newAlpha)*(DspFloatType)PSM_FFT_LEN);

		// --- check for change
		if (newOutputBufferLength == outputBufferLength)
			return;

		// --- new stuff
		alphaStretchRatio = newAlpha;
		ha = hs / alphaStretchRatio;

		// --- set output resample buffer
		outputBufferLength = newOutputBufferLength;

		// --- create Hann window
		if (windowBuff) delete[] windowBuff;
		windowBuff = new DspFloatType[outputBufferLength];
		windowCorrection = 0.0;
		for (unsigned int i = 0; i < outputBufferLength; i++)
		{
			windowBuff[i] = 0.5 * (1.0 - cos((i*2.0*kPi) / (outputBufferLength)));
			windowCorrection += windowBuff[i];
		}
		windowCorrection = 1.0 / windowCorrection;

		// --- create output buffer
		if (outputBuff) delete[] outputBuff;
		outputBuff = new DspFloatType[outputBufferLength];
		memset(outputBuff, 0, sizeof(DspFloatType)*outputBufferLength);
	}

	/** find bin index of nearest peak bin in previous FFT frame */
	int findPreviousNearestPeak(int peakIndex)
	{
		if (peakBinsPrevious[0] == -1) // first run, there is no peak
			return -1;

		int delta = -1;
		int previousPeak = -1;
		for (int i = 0; i < PSM_FFT_LEN; i++)
		{
			if (peakBinsPrevious[i] < 0)
				break;

			int dist = abs(peakIndex - peakBinsPrevious[i]);
			if (dist > PSM_FFT_LEN/4)
				break;

			if (i == 0)
			{
				previousPeak = i;
				delta = dist;
			}
			else if (dist < delta)
			{
				previousPeak = i;
				delta = dist;
			}
		}

		return previousPeak;
	}

	/** identify peak bins and tag their respective regions of influence */
	void findPeaksAndRegionsOfInfluence()
	{
		// --- FIND PEAKS --- //
		//
		// --- find local maxima in 4-sample window
		DspFloatType localWindow[4] = { 0.0, 0.0, 0.0, 0.0 };
		int m = 0;
		for (int i = 0; i < PSM_FFT_LEN; i++)
		{
			if (i == 0)
			{
				localWindow[0] = 0.0;
				localWindow[1] = 0.0;
				localWindow[2] = binData[i + 1].magnitude;
				localWindow[3] = binData[i + 2].magnitude;
			}
			else  if (i == 1)
			{
				localWindow[0] = 0.0;
				localWindow[1] = binData[i - 1].magnitude;
				localWindow[2] = binData[i + 1].magnitude;
				localWindow[3] = binData[i + 2].magnitude;
			}
			else  if (i == PSM_FFT_LEN - 1)
			{
				localWindow[0] = binData[i - 2].magnitude;
				localWindow[1] = binData[i - 1].magnitude;
				localWindow[2] = 0.0;
				localWindow[3] = 0.0;
			}
			else  if (i == PSM_FFT_LEN - 2)
			{
				localWindow[0] = binData[i - 2].magnitude;
				localWindow[1] = binData[i - 1].magnitude;
				localWindow[2] = binData[i + 1].magnitude;
				localWindow[3] = 0.0;
			}
			else
			{
				localWindow[0] = binData[i - 2].magnitude;
				localWindow[1] = binData[i - 1].magnitude;
				localWindow[2] = binData[i + 1].magnitude;
				localWindow[3] = binData[i + 2].magnitude;
			}

			// --- found peak bin!
			if (binData[i].magnitude > 0.00001 &&
				binData[i].magnitude > localWindow[0]
				&& binData[i].magnitude > localWindow[1]
				&& binData[i].magnitude > localWindow[2]
				&& binData[i].magnitude > localWindow[3])
			{
				binData[i].isPeak = true;
				peakBins[m++] = i;

				// --- for peak bins, assume that it is part of a previous, moving peak
				if (parameters.enablePeakTracking)
					binData[i].previousPeakBin = findPreviousNearestPeak(i);
				else
					binData[i].previousPeakBin = -1;
			}
		}

		// --- assign peak bosses
		if (m > 0)
		{
			int n = 0;
			int bossPeakBin = peakBins[n];
			DspFloatType nextPeak = peakBins[++n];
			int midBoundary = (nextPeak - (DspFloatType)bossPeakBin) / 2.0 + bossPeakBin;

			if (nextPeak >= 0)
			{
				for (int i = 0; i < PSM_FFT_LEN; i++)
				{
					if (i <= bossPeakBin)
					{
						binData[i].localPeakBin = bossPeakBin;
					}
					else if (i < midBoundary)
					{
						binData[i].localPeakBin = bossPeakBin;
					}
					else // > boundary, calc next set
					{
						bossPeakBin = nextPeak;
						nextPeak = peakBins[++n];
						if (nextPeak > bossPeakBin)
							midBoundary = (nextPeak - (DspFloatType)bossPeakBin) / 2.0 + bossPeakBin;
						else // nextPeak == -1
							midBoundary = PSM_FFT_LEN;

						binData[i].localPeakBin = bossPeakBin;
					}
				}
			}
		}
	}

	/** process input sample through PSM vocoder */
	/**
	\param xn input
	\return the processed sample
	*/
	virtual DspFloatType processAudioSample(DspFloatType input)
	{
			bool fftReady = false;
			DspFloatType output = 0.0;

			// --- normal processing
			output = vocoder.processAudioSample(input, fftReady);

			// --- if FFT is here, GO!
			if (fftReady)
			{
			// --- get the FFT data
			fftw_complex* fftData = vocoder.getFFTData();

			if (parameters.enablePeakPhaseLocking)
			{
				// --- get the magnitudes for searching
				for (int i = 0; i < PSM_FFT_LEN; i++)
				{
					binData[i].reset();
					peakBins[i] = -1;

					// --- store mag and phase
					binData[i].magnitude = getMagnitude(fftData[i][0], fftData[i][1]);
					binData[i].phi = getPhase(fftData[i][0], fftData[i][1]);
				}

				findPeaksAndRegionsOfInfluence();

				// --- each bin data should now know its local boss-peak
				//
				// --- now propagate phases accordingly
				//
				//     FIRST: set PSI angles of bosses
				for (int i = 0; i < PSM_FFT_LEN; i++)
				{
					DspFloatType mag_k = binData[i].magnitude;
					DspFloatType phi_k = binData[i].phi;

					// --- horizontal phase propagation
					//
					// --- omega_k = bin frequency(k)
					DspFloatType omega_k = kTwoPi*i / PSM_FFT_LEN;

					// --- phase deviation is actual - expected phase
					//     = phi_k -(phi(last frame) + wk*ha
					DspFloatType phaseDev = phi_k - phi[i] - omega_k*ha;

					// --- unwrapped phase increment
					DspFloatType deltaPhi = omega_k*ha + principalArg(phaseDev);

					// --- save for next frame
					phi[i] = phi_k;

					// --- if peak, assume it could have hopped from a different bin
					if (binData[i].isPeak)
					{
						// --- calculate new phase based on stretch factor; save phase for next time
						if(binData[i].previousPeakBin < 0)
							psi[i] = principalArg(psi[i] + deltaPhi * alphaStretchRatio);
						else
							psi[i] = principalArg(psi[binDataPrevious[i].previousPeakBin] + deltaPhi * alphaStretchRatio);
					}

					// --- save peak PSI (new angle)
					binData[i].psi = psi[i];

					// --- for IFFT
					binData[i].updatedPhase = binData[i].psi;
				}

				// --- now set non-peaks
				for (int i = 0; i < PSM_FFT_LEN; i++)
				{
					if (!binData[i].isPeak)
					{
						int myPeakBin = binData[i].localPeakBin;

						DspFloatType PSI_kp = binData[myPeakBin].psi;
						DspFloatType phi_kp = binData[myPeakBin].phi;

						// --- save for next frame
						// phi[i] = binData[myPeakBin].phi;

						// --- calculate new phase, locked to boss peak
						psi[i] = principalArg(PSI_kp - phi_kp - binData[i].phi);
						binData[i].updatedPhase = psi[i];// principalArg(PSI_kp - phi_kp - binData[i].phi);
					}
				}

				for (int i = 0; i < PSM_FFT_LEN; i++)
				{
					DspFloatType mag_k = binData[i].magnitude;

					// --- convert back
					fftData[i][0] = mag_k*cos(binData[i].updatedPhase);
					fftData[i][1] = mag_k*sin(binData[i].updatedPhase);

					// --- save for next frame
					binDataPrevious[i] = binData[i];
					peakBinsPrevious[i] = peakBins[i];

				}
			}// end if peak locking

			else // ---> old school
			{
				for (int i = 0; i < PSM_FFT_LEN; i++)
				{
					DspFloatType mag_k = getMagnitude(fftData[i][0], fftData[i][1]);
					DspFloatType phi_k = getPhase(fftData[i][0], fftData[i][1]);

					// --- horizontal phase propagation
					//
					// --- omega_k = bin frequency(k)
					DspFloatType omega_k = kTwoPi*i / PSM_FFT_LEN;

					// --- phase deviation is actual - expected phase
					//     = phi_k -(phi(last frame) + wk*ha
					DspFloatType phaseDev = phi_k - phi[i] - omega_k*ha;

					// --- unwrapped phase increment
					DspFloatType deltaPhi = omega_k*ha + principalArg(phaseDev);

					// --- save for next frame
					phi[i] = phi_k;

					// --- calculate new phase based on stretch factor; save phase for next time
					psi[i] = principalArg(psi[i] + deltaPhi * alphaStretchRatio);

					// --- convert back
					fftData[i][0] = mag_k*cos(psi[i]);
					fftData[i][1] = mag_k*sin(psi[i]);
				}
			}


			// --- manually so the IFFT (OPTIONAL)
			vocoder.doInverseFFT();

			// --- can get the iFFT buffers
			fftw_complex* inv_fftData = vocoder.getIFFTData();

			// --- make copy (can speed this up)
			DspFloatType ifft[PSM_FFT_LEN] = { 0.0 };
			for (int i = 0; i < PSM_FFT_LEN; i++)
				ifft[i] = inv_fftData[i][0];

			// --- resample the audio as if it were stretched
			resample(&ifft[0], outputBuff, PSM_FFT_LEN, outputBufferLength, interpolation::kLinear, windowCorrection, windowBuff);

			// --- overlap-add the interpolated buffer to complete the operation
			vocoder.doOverlapAdd(&outputBuff[0], outputBufferLength);
			}

			return output;
	}

	/** get parameters: note use of custom structure for passing param data */
	/**
	\return PSMVocoderParameters custom data structure
	*/
	PSMVocoderParameters getParameters()
	{
		return parameters;
	}

	/** set parameters: note use of custom structure for passing param data */
	/**
	\param PSMVocoderParameters custom data structure
	*/
	void setParameters(const PSMVocoderParameters& params)
	{
		if (params.pitchShiftSemitones != parameters.pitchShiftSemitones)
		{
			setPitchShift(params.pitchShiftSemitones);
		}

		// --- save
		parameters = params;
	}

protected:
	PSMVocoderParameters parameters;	///< object parameters
	PhaseVocoder vocoder;				///< vocoder to perform PSM
	DspFloatType alphaStretchRatio = 1.0;		///< alpha stretch ratio = hs/ha

	// --- FFT is 4096 with 75% overlap
	const DspFloatType hs = PSM_FFT_LEN / 4;	///< hs = N/4 --- 75% overlap
	DspFloatType ha = PSM_FFT_LEN / 4;		///< ha = N/4 --- 75% overlap
	DspFloatType phi[PSM_FFT_LEN] = { 0.0 };	///< array of phase values for classic algorithm
	DspFloatType psi[PSM_FFT_LEN] = { 0.0 };	///< array of phase correction values for classic algorithm

	// --- for peak-locking
	BinData binData[PSM_FFT_LEN];			///< array of BinData structures for current FFT frame
	BinData binDataPrevious[PSM_FFT_LEN];	///< array of BinData structures for previous FFT frame

	int peakBins[PSM_FFT_LEN] = { -1 };		///< array of current peak bin index values (-1 = not peak)
	int peakBinsPrevious[PSM_FFT_LEN] = { -1 }; ///< array of previous peak bin index values (-1 = not peak)

	DspFloatType* windowBuff = nullptr;			///< buffer for window
	DspFloatType* outputBuff = nullptr;			///< buffer for resampled output
	DspFloatType windowCorrection = 0.0;			///< window correction value
	unsigned int outputBufferLength = 0;	///< lenght of resampled output array
};

/**
\brief destroys the FFTW arrays and plans.
*/
void PhaseVocoder::destroyFFTW()
{
	if (plan_forward)
		fftw_destroy_plan(plan_forward);
	if (plan_backward)
		fftw_destroy_plan(plan_backward);

	if (fft_input)
		fftw_free(fft_input);
	if (fft_result)
		fftw_free(fft_result);
	if (ifft_result)
		fftw_free(ifft_result);
}

/**
\brief initialize the Fast FFT object for operation
- NOTES:<br>
See notes on symmetrical window arrays in comments.<br>
\param _frameLength the FFT length - MUST be a power of 2
\param _hopSize the hop size in samples: this object only supports ha = hs (pure real-time operation only)
\param _window the window type (note: may be set to windowType::kNoWindow)
*/
void PhaseVocoder::initialize(unsigned int _frameLength, unsigned int _hopSize, windowType _window)
{
	frameLength = _frameLength;
	wrapMask = frameLength - 1;
	hopSize = _hopSize;
	window = _window;

	// --- this is the overlap as a fraction i.e. 0.75 = 75%
	overlap = hopSize > 0.0 ? 1.0 - (DspFloatType)hopSize / (DspFloatType)frameLength : 0.0;

	// --- gain correction for window + hop size
	windowHopCorrection = 0.0;

	// --- SETUP BUFFERS ---- //
	//     NOTE: input and output buffers are circular, others are linear
	//
	// --- input buffer, for processing the x(n) timeline
	if (inputBuffer)
		delete inputBuffer;

	inputBuffer = new DspFloatType[frameLength];
	memset(&inputBuffer[0], 0, frameLength * sizeof(DspFloatType));

	// --- output buffer, for processing the y(n) timeline and accumulating frames
	if (outputBuffer)
		delete outputBuffer;

	// --- the output buffer is declared as 2x the normal frame size
	//     to accomodate time-stretching/pitch shifting; you can increase the size
	//     here; if so make sure to calculate the wrapMaskOut properly and everything
	//     will work normally you can even dynamically expand and contract the buffer
	//     (not sure why you would do this - and it will surely affect CPU performance)
	//     NOTE: the length of the buffer is only to accomodate accumulations
	//           it does not stretch time or change causality on its own
	outputBuffer = new DspFloatType[frameLength * 4];
	memset(&outputBuffer[0], 0, (frameLength*4.0) * sizeof(DspFloatType));
	wrapMaskOut = (frameLength*4.0) - 1;

	// --- fixed window buffer
	if (windowBuffer)
		delete windowBuffer;

	windowBuffer = new DspFloatType[frameLength];
	memset(&windowBuffer[0], 0, frameLength * sizeof(DspFloatType));

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
			windowHopCorrection += windowBuffer[n];
		}
	}
	else if (window == windowType::kHammingWindow)
	{
		for (int n = 0; n < frameLength - 1; n++)
		{
			windowBuffer[n] = 0.54 - 0.46*cos((n*2.0*kPi) / (frameLength));
			windowHopCorrection += windowBuffer[n];
		}
	}
	else if (window == windowType::kHannWindow)
	{
		for (int n = 0; n < frameLength; n++)
		{
			windowBuffer[n] = 0.5 * (1 - cos((n*2.0*kPi) / (frameLength)));
			windowHopCorrection += windowBuffer[n];
		}
	}
	else if (window == windowType::kBlackmanHarrisWindow)
	{
		for (int n = 0; n < frameLength; n++)
		{
			windowBuffer[n] = (0.42323 - (0.49755*cos((n*2.0*kPi) / (frameLength))) + 0.07922*cos((2 * n*2.0*kPi) / (frameLength)));
			windowHopCorrection += windowBuffer[n];
		}
	}
	else if (window == windowType::kNoWindow)
	{
		for (int n = 0; n < frameLength; n++)
		{
			windowBuffer[n] = 1.0;
			windowHopCorrection += windowBuffer[n];
		}
	}
	else // --- default to kNoWindow
	{
		for (int n = 0; n < frameLength; n++)
		{
			windowBuffer[n] = 1.0;
			windowHopCorrection += windowBuffer[n];
		}
	}

	// --- calculate gain correction factor
	if (window != windowType::kNoWindow)
		windowHopCorrection = (1.0 - overlap) / windowHopCorrection;
	else
		windowHopCorrection = 1.0 / windowHopCorrection;

	// --- set
	inputWriteIndex = 0;
	inputReadIndex = 0;

	outputWriteIndex = 0;
	outputReadIndex = 0;

	fftCounter = 0;

	// --- reset flags
	needInverseFFT = false;
	needOverlapAdd = false;

#ifdef HAVE_FFTW
	destroyFFTW();
	fft_input = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * frameLength);
	fft_result = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * frameLength);
	ifft_result = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * frameLength);

	plan_forward = fftw_plan_dft_1d(frameLength, fft_input, fft_result, FFTW_FORWARD, FFTW_ESTIMATE);
	plan_backward = fftw_plan_dft_1d(frameLength, fft_result, ifft_result, FFTW_BACKWARD, FFTW_ESTIMATE);
#endif
}

/**
\brief zero pad the input timeline
- NOTES:<br>
\param count the number of zero-valued samples to insert
\returns true if the zero-insertion triggered a FFT event, false otherwise
*/
bool PhaseVocoder::addZeroPad(unsigned int count)
{
	bool fftReady = false;
	for (unsigned int i = 0; i < count; i++)
	{
		// --- push into buffer
		inputBuffer[inputWriteIndex++] = 0.0;

		// --- wrap
		inputWriteIndex &= wrapMask;

		// --- check the FFT
		bool didFFT = advanceAndCheckFFT();

		// --- for a zero-padding operation, the last inserted zero
		//     should trigger the FFT; if not something has gone horribly wrong
		if (didFFT && i == count - 1)
			fftReady = true;
	}

	return fftReady;
}

/**
\brief advance the sample counter and check to see if we need to do the FFT.
- NOTES:<br>
\returns true if the advancement triggered a FFT event, false otherwise
*/
bool PhaseVocoder::advanceAndCheckFFT()
{
	// --- inc counter and check count
	fftCounter++;

	if (fftCounter != frameLength)
		return false;

	// --- we have a FFT ready
	// --- load up the input to the FFT
	for (int i = 0; i < frameLength; i++)
	{
		fft_input[i][0] = inputBuffer[inputReadIndex++] * windowBuffer[i];
		fft_input[i][1] = 0.0; // use this if your data is complex valued

		// --- wrap if index > bufferlength - 1
		inputReadIndex &= wrapMask;
	}

	// --- do the FFT
	fftw_execute(plan_forward);

	// --- in case user does not take IFFT, just to prevent zero output
	needInverseFFT = true;
	needOverlapAdd = true;

	// --- fft counter: small hop = more FFTs = less counting before fft
	//
	// --- overlap-add-only algorithms do not involve hop-size in FFT count
	if (overlapAddOnly)
		fftCounter = 0;
	else // normal counter advance
		fftCounter = frameLength - hopSize;

	// --- setup the read index for next time through the loop
	if (!overlapAddOnly)
		inputReadIndex += hopSize;

	// --- wrap if needed
	inputReadIndex &= wrapMask;

	return true;
}

/**
\brief process one input sample throug the vocoder to produce one output sample
- NOTES:<br>
\param input the input sample x(n)
\param fftReady a return flag indicating if the FFT has occurred and FFT data is ready to process
\returns the vocoder output sample y(n)
*/
DspFloatType PhaseVocoder::processAudioSample(DspFloatType input, bool& fftReady)
{
	// --- if user did not manually do fft and overlap, do them here
	//     this allows maximum flexibility in use of the object
	if (needInverseFFT)
		doInverseFFT();
	if(needOverlapAdd)
		doOverlapAdd();

	fftReady = false;

	// --- get the current output sample first
	DspFloatType currentOutput = outputBuffer[outputReadIndex];

	// --- set the buffer to 0.0 in preparation for the next overlap/add process
	outputBuffer[outputReadIndex++] = 0.0;

	// --- wrap
	outputReadIndex &= wrapMaskOut;

	// --- push into buffer
	inputBuffer[inputWriteIndex++] = (DspFloatType)input;

	// --- wrap
	inputWriteIndex &= wrapMask;

	// --- check the FFT
	fftReady = advanceAndCheckFFT();

	return currentOutput;
}

/**
\brief perform the inverse FFT on the processed data
- NOTES:<br>
This function is optional - if you need to sequence the output (synthesis) stage yourself <br>
then you can call this function at the appropriate time - see the PSMVocoder object for an example
*/
void PhaseVocoder::doInverseFFT()
{
	// do the IFFT
	fftw_execute(plan_backward);

	// --- output is now in ifft_result array
	needInverseFFT = false;
}

/**
\brief perform the overlap/add on the IFFT data
- NOTES:<br>
This function is optional - if you need to sequence the output (synthesis) stage yourself <br>
then you can call this function at the appropriate time - see the PSMVocoder object for an example
\param outputData an array of data to overlap/add: if this is NULL then the IFFT data is used
\param length the lenght of the array of data to overlap/add: if this is -1, the normal IFFT length is used
*/
void PhaseVocoder::doOverlapAdd(DspFloatType* outputData, int length)
{
	// --- overlap/add with output buffer
	//     NOTE: this assumes input and output hop sizes are the same!
	outputWriteIndex = outputReadIndex;

	if (outputData)
	{
		for (int i = 0; i < length; i++)
		{
			// --- if you need to window the data, do so prior to this function call
			outputBuffer[outputWriteIndex++] += outputData[i];

			// --- wrap if index > bufferlength - 1
			outputWriteIndex &= wrapMaskOut;
		}
		needOverlapAdd = false;
		return;
	}

	for (int i = 0; i < frameLength; i++)
	{
		// --- accumulate
		outputBuffer[outputWriteIndex++] += windowHopCorrection * ifft_result[i][0];

		// --- wrap if index > bufferlength - 1
		outputWriteIndex &= wrapMaskOut;
	}

	// --- set a flag
	needOverlapAdd = false;
}

/**
\class FastConvolver
\ingroup FFTW-Objects
\brief
The FastConvolver provides a fast convolver - the user supplies the filter IR and the object
snapshots the FFT of that filter IR. Input audio is fast-convovled with the filter FFT using
complex multiplication and zero-padding.
Audio I/O:
- processes mono input into mono output.
Control I/F:
- none.
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
class FastConvolver
{
public:
	FastConvolver() {
		vocoder.setOverlapAddOnly(true);
	}		/* C-TOR */
	~FastConvolver() {
		if (filterIR)
			delete[] filterIR;

		if (filterFFT)
			fftw_free(filterFFT);
	}	/* D-TOR */

	/** setup the FFT for a given IR length */
	/**
	\param _filterImpulseLength the filter IR length, which is 1/2 FFT length due to need for zero-padding (see FX book)
	*/
	void initialize(unsigned int _filterImpulseLength)
	{
		if (filterImpulseLength == _filterImpulseLength)
			return;

		// --- setup a specialized vocoder with 50% hop size
		filterImpulseLength = _filterImpulseLength;
		vocoder.initialize(filterImpulseLength * 2, filterImpulseLength, windowType::kNoWindow);

		// --- initialize the FFT object for capturing the filter FFT
		filterFastFFT.initialize(filterImpulseLength * 2, windowType::kNoWindow);

		// --- array to hold the filter IR; this could be localized to the particular function that uses it
		if (filterIR)
			delete [] filterIR;
		filterIR = new DspFloatType[filterImpulseLength * 2];
		memset(&filterIR[0], 0, filterImpulseLength * 2 * sizeof(DspFloatType));

		// --- allocate the filter FFT arrays
		if(filterFFT)
			fftw_free(filterFFT);

		 filterFFT = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * filterImpulseLength * 2);

		 // --- reset
		 inputCount = 0;
	}

	/** setup the filter IRirBuffer MUST be exactly filterImpulseLength in size, or this will crash! */
	void setFilterIR(DspFloatType* irBuffer)
	{
		if (!irBuffer) return;

		memset(&filterIR[0], 0, filterImpulseLength * 2 * sizeof(DspFloatType));

		// --- copy over first half; filterIR len = filterImpulseLength * 2
		int m = 0;
		for (unsigned int i = 0; i < filterImpulseLength; i++)
		{
			filterIR[i] = irBuffer[i];
		}

		// --- take FFT of the h(n)
		fftw_complex* fftOfFilter = filterFastFFT.doFFT(&filterIR[0]);

		// --- copy the FFT into our local buffer for storage; also
		//     we never want to hold a pointer to a FFT output
		//     for more than one local function's worth
		//     could replace with memcpy( )
		for (int i = 0; i < 2; i++)
		{
			for (unsigned int j = 0; j < filterImpulseLength * 2; j++)
			{
				filterFFT[j][i] = fftOfFilter[j][i];
			}
		}
	}

	/** process an input sample through convolver */
	DspFloatType processAudioSample(DspFloatType input)
	{
		bool fftReady = false;
		DspFloatType output = 0.0;

		if (inputCount == filterImpulseLength)
		{
			fftReady = vocoder.addZeroPad(filterImpulseLength);

			if (fftReady) // should happen on time
			{
				// --- multiply our filter IR with the vocoder FFT
				fftw_complex* signalFFT = vocoder.getFFTData();
				if (signalFFT)
				{
					unsigned int fff = vocoder.getFrameLength();

					// --- complex multiply with FFT of IR
					for (unsigned int i = 0; i < filterImpulseLength * 2; i++)
					{
						// --- get real/imag parts of each FFT
						ComplexNumber signal(signalFFT[i][0], signalFFT[i][1]);
						ComplexNumber filter(filterFFT[i][0], filterFFT[i][1]);

						// --- use complex multiply function; this convolves in the time domain
						ComplexNumber product = complexMultiply(signal, filter);

						// --- overwrite the FFT bins
						signalFFT[i][0] = product.real;
						signalFFT[i][1] = product.imag;
					}
				}
			}

			// --- reset counter
			inputCount = 0;
		}

		// --- process next sample
		output = vocoder.processAudioSample(input, fftReady);
		inputCount++;

		return output;
	}

	/** get current frame length */
	unsigned int getFrameLength() { return vocoder.getFrameLength(); }

	/** get current IR length*/
	unsigned int getFilterIRLength() { return filterImpulseLength; }

protected:
	PhaseVocoder vocoder;				///< vocoder object
	FastFFT filterFastFFT;				///< FastFFT object
	fftw_complex* filterFFT = nullptr;	///< filterFFT output arrays
	DspFloatType* filterIR = nullptr;			///< filter IR
	unsigned int inputCount = 0;		///< input sample counter
	unsigned int filterImpulseLength = 0;///< IR length
};
