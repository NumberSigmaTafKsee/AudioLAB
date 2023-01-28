#pragma once

#include "FXObjects.hpp"

// --- constants
const unsigned int TLD_AUDIO_DETECT_MODE_PEAK = 0;
const unsigned int TLD_AUDIO_DETECT_MODE_MS = 1;
const unsigned int TLD_AUDIO_DETECT_MODE_RMS = 2;
const double TLD_AUDIO_ENVELOPE_ANALOG_TC = -0.99967234081320612357829304641019; // ln(36.7%)

/**
\struct AudioDetectorParameters
\ingroup FX-Objects
\brief
Custom parameter structure for the AudioDetector object. NOTE: this object uses constant defintions:
- const unsigned int TLD_AUDIO_DETECT_MODE_PEAK = 0;
- const unsigned int TLD_AUDIO_DETECT_MODE_MS = 1;
- const unsigned int TLD_AUDIO_DETECT_MODE_RMS = 2;
- const double TLD_AUDIO_ENVELOPE_ANALOG_TC = -0.99967234081320612357829304641019; // ln(36.7%)
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
struct AudioDetectorParameters
{
	AudioDetectorParameters() {}
	/** all FXObjects parameter objects require overloaded= operator so remember to add new entries if you add new variables. */
	AudioDetectorParameters& operator=(const AudioDetectorParameters& params)	// need this override for collections to work
	{
		if (this == &params)
			return *this;
		attackTime_mSec = params.attackTime_mSec;
		releaseTime_mSec = params.releaseTime_mSec;
		detectMode = params.detectMode;
		detect_dB = params.detect_dB;
		clampToUnityMax = params.clampToUnityMax;
		return *this;
	}

	// --- individual parameters
	double attackTime_mSec = 0.0; ///< attack time in milliseconds
	double releaseTime_mSec = 0.0;///< release time in milliseconds
	unsigned int  detectMode = 0;///< detect mode, see TLD_ constants above
	bool detect_dB = false;	///< detect in dB  DEFAULT  = false (linear NOT log)
	bool clampToUnityMax = true;///< clamp output to 1.0 (set false for true log detectors)
};

/**
\class AudioDetector
\ingroup FX-Objects
\brief
The AudioDetector object implements the audio detector defined in the book source below.
NOTE: this detector can receive signals and transmit detection values that are both > 0dBFS
Audio I/O:
- Processes mono input to a detected signal output.
Control I/F:
- Use AudioDetectorParameters structure to get/set object params.
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
class AudioDetector : public IAudioSignalProcessor
{
public:
	AudioDetector() {}	/* C-TOR */
	~AudioDetector() {}	/* D-TOR */

public:
	/** set sample rate dependent time constants and clear last envelope output value */
	virtual bool reset(double _sampleRate)
	{
		setSampleRate(_sampleRate);
		lastEnvelope = 0.0;
		return true;
	}

	/** return false: this object only processes samples */
	virtual bool canProcessAudioFrame() { return false; }

	// --- process audio: detect the log envelope and return it in dB
	/**
	\param xn input
	\return the processed sample
	*/
	virtual double processAudioSample(double xn)
	{
		// --- all modes do Full Wave Rectification
		double input = fabs(xn);

		// --- square it for MS and RMS
		if (audioDetectorParameters.detectMode == TLD_AUDIO_DETECT_MODE_MS ||
			audioDetectorParameters.detectMode == TLD_AUDIO_DETECT_MODE_RMS)
			input *= input;

		// --- to store current
		double currEnvelope = 0.0;

		// --- do the detection with attack or release applied
		if (input > lastEnvelope)
			currEnvelope = attackTime * (lastEnvelope - input) + input;
		else
			currEnvelope = releaseTime * (lastEnvelope - input) + input;

		// --- we are recursive so need to check underflow
		checkFloatUnderflow(currEnvelope);

		// --- bound them; can happen when using pre-detector gains of more than 1.0
		if (audioDetectorParameters.clampToUnityMax)
			currEnvelope = fmin(currEnvelope, 1.0);

		// --- can not be (-)
		currEnvelope = fmax(currEnvelope, 0.0);

		// --- store envelope prior to sqrt for RMS version
		lastEnvelope = currEnvelope;

		// --- if RMS, do the SQRT
		if (audioDetectorParameters.detectMode == TLD_AUDIO_DETECT_MODE_RMS)
			currEnvelope = pow(currEnvelope, 0.5);

		// --- if not dB, we are done
		if (!audioDetectorParameters.detect_dB)
			return currEnvelope;

		// --- setup for log( )
		if (currEnvelope <= 0)
		{
			return -96.0;
		}

		// --- true log output in dB, can go above 0dBFS!
		return 20.0*log10(currEnvelope);
	}

	/** get parameters: note use of custom structure for passing param data */
	/**
	\return AudioDetectorParameters custom data structure
	*/
	AudioDetectorParameters getParameters()
	{
		return audioDetectorParameters;
	}

	/** set parameters: note use of custom structure for passing param data */
	/**
	\param AudioDetectorParameters custom data structure
	*/
	void setParameters(const AudioDetectorParameters& parameters)
	{
		audioDetectorParameters = parameters;

		// --- update structure
		setAttackTime(audioDetectorParameters.attackTime_mSec, true);
		setReleaseTime(audioDetectorParameters.releaseTime_mSec, true);

	}

	/** set sample rate - our time constants depend on it */
	virtual void setSampleRate(double _sampleRate)
	{
		if (sampleRate == _sampleRate)
			return;

		sampleRate = _sampleRate;

		// --- recalculate RC time-constants
		setAttackTime(audioDetectorParameters.attackTime_mSec, true);
		setReleaseTime(audioDetectorParameters.releaseTime_mSec, true);
	}

protected:
	AudioDetectorParameters audioDetectorParameters; ///< parameters for object
	double attackTime = 0.0;	///< attack time coefficient
	double releaseTime = 0.0;	///< release time coefficient
	double sampleRate = 44100;	///< stored sample rate
	double lastEnvelope = 0.0;	///< output register

	/** set our internal atack time coefficients based on times and sample rate */
	void setAttackTime(double attack_in_ms, bool forceCalc = false);

	/** set our internal release time coefficients based on times and sample rate */
	void setReleaseTime(double release_in_ms, bool forceCalc = false);
};


/**
\enum dynamicsProcessorType
\ingroup Constants-Enums
\brief
Use this strongly typed enum to set the dynamics processor type.
- enum class dynamicsProcessorType { kCompressor, kDownwardExpander };
- limiting is the extreme version of kCompressor
- gating is the extreme version of kDownwardExpander
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/

// --- processorType
enum class dynamicsProcessorType { kCompressor, kDownwardExpander };


/**
\struct DynamicsProcessorParameters
\ingroup FX-Objects
\brief
Custom parameter structure for the DynamicsProcessor object. Ths struct includes all information needed from GUI controls.
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
struct DynamicsProcessorParameters
{
	DynamicsProcessorParameters() {}
	/** all FXObjects parameter objects require overloaded= operator so remember to add new entries if you add new variables. */
	DynamicsProcessorParameters& operator=(const DynamicsProcessorParameters& params)	// need this override for collections to work
	{
		if (this == &params)
			return *this;

		ratio = params.ratio;
		threshold_dB = params.threshold_dB;
		kneeWidth_dB = params.kneeWidth_dB;
		hardLimitGate = params.hardLimitGate;
		softKnee = params.softKnee;
		enableSidechain = params.enableSidechain;
		calculation = params.calculation;
		attackTime_mSec = params.attackTime_mSec;
		releaseTime_mSec = params.releaseTime_mSec;
		outputGain_dB = params.outputGain_dB;
		// --- NOTE: do not set outbound variables??
		gainReduction = params.gainReduction;
		gainReduction_dB = params.gainReduction_dB;
		return *this;
	}

	// --- individual parameters
	double ratio = 50.0;				///< processor I/O gain ratio
	double threshold_dB = -10.0;		///< threshold in dB
	double kneeWidth_dB = 10.0;			///< knee width in dB for soft-knee operation
	bool hardLimitGate = false;			///< threshold in dB
	bool softKnee = true;				///< soft knee flag
	bool enableSidechain = false;		///< enable external sidechain input to object
	dynamicsProcessorType calculation = dynamicsProcessorType::kCompressor; ///< processor calculation type
	double attackTime_mSec = 0.0;		///< attack mSec
	double releaseTime_mSec = 0.0;		///< release mSec
	double outputGain_dB = 0.0;			///< make up gain

	// --- outbound values, for owner to use gain-reduction metering
	double gainReduction = 1.0;			///< output value for gain reduction that occurred
	double gainReduction_dB = 0.0;		///< output value for gain reduction that occurred in dB
};

/**
\class DynamicsProcessor
\ingroup FX-Objects
\brief
The DynamicsProcessor object implements a dynamics processor suite: compressor, limiter, downward expander, gate.
Audio I/O:
- Processes mono input to mono output.
Control I/F:
- Use DynamicsProcessorParameters structure to get/set object params.
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
class DynamicsProcessor : public IAudioSignalProcessor
{
public:
	DynamicsProcessor() {}	/* C-TOR */
	~DynamicsProcessor() {}	/* D-TOR */

public:
	/** reset members to initialized state */
	virtual bool reset(double _sampleRate)
	{
		sidechainInputSample = 0.0;
		detector.reset(_sampleRate);
		AudioDetectorParameters detectorParams = detector.getParameters();
		detectorParams.clampToUnityMax = false;
		detectorParams.detect_dB = true;
		detector.setParameters(detectorParams);
		return true;
	}

	/** return false: this object only processes samples */
	virtual bool canProcessAudioFrame() { return false; }

	/** enable sidechaib input */
	virtual void enableAuxInput(bool enableAuxInput){ parameters.enableSidechain = enableAuxInput; }

	/** process the sidechain by saving the value for the upcoming processAudioSample() call */
	virtual double processAuxInputAudioSample(double xn)
	{
		sidechainInputSample = xn;
		return sidechainInputSample;
	}

	/** get parameters: note use of custom structure for passing param data */
	/**
	\return DynamicsProcessorParameters custom data structure
	*/
	DynamicsProcessorParameters getParameters(){ return parameters; }

	/** set parameters: note use of custom structure for passing param data */
	/**
	\param DynamicsProcessorParameters custom data structure
	*/
	void setParameters(const DynamicsProcessorParameters& _parameters)
	{
		parameters = _parameters;

		AudioDetectorParameters detectorParams = detector.getParameters();
		detectorParams.attackTime_mSec = parameters.attackTime_mSec;
		detectorParams.releaseTime_mSec = parameters.releaseTime_mSec;
		detector.setParameters(detectorParams);
	}

	/** process audio using feed-forward dynamics processor flowchart */
	/*
		1. detect input signal
		2. calculate gain
		3. apply to input sample
	*/
	/**
	\param xn input
	\return the processed sample
	*/
	virtual double processAudioSample(double xn)
	{
		// --- detect input
		double detect_dB = 0.0;

		// --- if using the sidechain, process the aux input
		if(parameters.enableSidechain)
			detect_dB = detector.processAudioSample(sidechainInputSample);
		else
			detect_dB = detector.processAudioSample(xn);

		// --- compute gain
		double gr = computeGain(detect_dB);

		// --- makeup gain
		double makeupGain = pow(10.0, parameters.outputGain_dB / 20.0);

		// --- do DCA + makeup gain
		return xn * gr * makeupGain;
	}

protected:
	DynamicsProcessorParameters parameters; ///< object parameters
	AudioDetector detector; ///< the sidechain audio detector

	// --- storage for sidechain audio input (mono only)
	double sidechainInputSample = 0.0; ///< storage for sidechain sample

	/** compute (and save) the current gain value based on detected input (dB) */
	inline double computeGain(double detect_dB)
	{
		double output_dB = 0.0;

		if (parameters.calculation == dynamicsProcessorType::kCompressor)
		{
			// --- hard knee
			if (!parameters.softKnee)
			{
				// --- below threshold, unity
				if (detect_dB <= parameters.threshold_dB)
					output_dB = detect_dB;
				else// --- above threshold, compress
				{
					if (parameters.hardLimitGate) // is limiter?
						output_dB = parameters.threshold_dB;
					else
						output_dB = parameters.threshold_dB + (detect_dB - parameters.threshold_dB) / parameters.ratio;
				}
			}
			else // --- calc gain with knee
			{
				// --- left side of knee, outside of width, unity gain zone
				if (2.0*(detect_dB - parameters.threshold_dB) < -parameters.kneeWidth_dB)
					output_dB = detect_dB;
				// --- else inside the knee,
				else if (2.0*(fabs(detect_dB - parameters.threshold_dB)) <= parameters.kneeWidth_dB)
				{
					if (parameters.hardLimitGate)	// --- is limiter?
						output_dB = detect_dB - pow((detect_dB - parameters.threshold_dB + (parameters.kneeWidth_dB / 2.0)), 2.0) / (2.0*parameters.kneeWidth_dB);
					else // --- 2nd order poly
						output_dB = detect_dB + (((1.0 / parameters.ratio) - 1.0) * pow((detect_dB - parameters.threshold_dB + (parameters.kneeWidth_dB / 2.0)), 2.0)) / (2.0*parameters.kneeWidth_dB);
				}
				// --- right of knee, compression zone
				else if (2.0*(detect_dB - parameters.threshold_dB) > parameters.kneeWidth_dB)
				{
					if (parameters.hardLimitGate) // --- is limiter?
						output_dB = parameters.threshold_dB;
					else
						output_dB = parameters.threshold_dB + (detect_dB - parameters.threshold_dB) / parameters.ratio;
				}
			}
		}
		else if (parameters.calculation == dynamicsProcessorType::kDownwardExpander)
		{
			// --- hard knee
			// --- NOTE: soft knee is not technically possible with a gate because there
			//           is no "left side" of the knee
			if (!parameters.softKnee || parameters.hardLimitGate)
			{
				// --- above threshold, unity gain
				if (detect_dB >= parameters.threshold_dB)
					output_dB = detect_dB;
				else
				{
					if (parameters.hardLimitGate) // --- gate: -inf(dB)
						output_dB = -1.0e34;
					else
						output_dB = parameters.threshold_dB + (detect_dB - parameters.threshold_dB) * parameters.ratio;
				}
			}
			else // --- calc gain with knee
			{
				// --- right side of knee, unity gain zone
				if (2.0*(detect_dB - parameters.threshold_dB) > parameters.kneeWidth_dB)
					output_dB = detect_dB;
				// --- in the knee
				else if (2.0*(fabs(detect_dB - parameters.threshold_dB)) > -parameters.kneeWidth_dB)
					output_dB = ((parameters.ratio - 1.0) * pow((detect_dB - parameters.threshold_dB - (parameters.kneeWidth_dB / 2.0)), 2.0)) / (2.0*parameters.kneeWidth_dB);
				// --- left side of knee, downward expander zone
				else if (2.0*(detect_dB - parameters.threshold_dB) <= -parameters.kneeWidth_dB)
					output_dB = parameters.threshold_dB + (detect_dB - parameters.threshold_dB) * parameters.ratio;
			}
		}

		// --- convert gain; store values for user meters
		parameters.gainReduction_dB = output_dB - detect_dB;
		parameters.gainReduction = pow(10.0, (parameters.gainReduction_dB) / 20.0);

		// --- the current gain coefficient value
		return parameters.gainReduction;
	}
};

/**
\class PeakLimiter
\ingroup FX-Objects
\brief
The PeakLimiter object implements a simple peak limiter; it is really a simplified and hard-wired
versio of the DynamicsProcessor
Audio I/O:
- Processes mono input to mono output.
Control I/F:
- setThreshold_dB(double _threshold_dB) to adjust the limiter threshold
- setMakeUpGain_dB(double _makeUpGain_dB) to adjust the makeup gain
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
class PeakLimiter : public IAudioSignalProcessor
{
public:
	PeakLimiter() { setThreshold_dB(-3.0); }
	~PeakLimiter() {}

	/** reset members to initialized state */
	virtual bool reset(double _sampleRate)
	{
		// --- init; true = analog time-constant
		detector.setSampleRate(_sampleRate);

		AudioDetectorParameters detectorParams = detector.getParameters();
		detectorParams.detect_dB = true;
		detectorParams.attackTime_mSec = 5.0;
		detectorParams.releaseTime_mSec = 25.0;
		detectorParams.clampToUnityMax = false;
		detectorParams.detectMode = ENVELOPE_DETECT_MODE_PEAK;
		detector.setParameters(detectorParams);

		return true;
	}

	/** return false: this object only processes samples */
	virtual bool canProcessAudioFrame() { return false; }

	/** process audio: implement hard limiter */
	/**
	\param xn input
	\return the processed sample
	*/
	virtual double processAudioSample(double xn)
	{
		return dB2Raw(makeUpGain_dB)*xn*computeGain(detector.processAudioSample(xn));
	}

	/** compute the gain reductino value based on detected value in dB */
	double computeGain(double detect_dB)
	{
		double output_dB = 0.0;

		// --- defaults - you can change these here
		bool softknee = true;
		double kneeWidth_dB = 10.0;

		// --- hard knee
		if (!softknee)
		{
			// --- below threshold, unity
			if (detect_dB <= threshold_dB)
				output_dB = detect_dB;
			// --- above threshold, compress
			else
				output_dB = threshold_dB;
		}
		else
		{
			// --- calc gain with knee
			// --- left side of knee, outside of width, unity gain zone
			if (2.0*(detect_dB - threshold_dB) < -kneeWidth_dB)
				output_dB = detect_dB;
			// --- inside the knee,
			else if (2.0*(fabs(detect_dB - threshold_dB)) <= kneeWidth_dB)
				output_dB = detect_dB - pow((detect_dB - threshold_dB + (kneeWidth_dB / 2.0)), 2.0) / (2.0*kneeWidth_dB);
			// --- right of knee, compression zone
			else if (2.0*(detect_dB - threshold_dB) > kneeWidth_dB)
				output_dB = threshold_dB;
		}

		// --- convert difference between threshold and detected to raw
		return  pow(10.0, (output_dB - detect_dB) / 20.0);
	}

	/** adjust threshold in dB */
	void setThreshold_dB(double _threshold_dB) { threshold_dB = _threshold_dB; }

	/** adjust makeup gain in dB*/
	void setMakeUpGain_dB(double _makeUpGain_dB) { makeUpGain_dB = _makeUpGain_dB; }

protected:
	AudioDetector detector;		///< the detector object
	double threshold_dB = 0.0;	///< stored threshold (dB)
	double makeUpGain_dB = 0.0;	///< stored makeup gain (dB)
};



/**
\brief sets the new attack time and re-calculates the time constant
\param attack_in_ms the new attack timme
\param forceCalc flag to force a re-calculation of time constant even if values have not changed.
*/
void AudioDetector::setAttackTime(double attack_in_ms, bool forceCalc)
{
	if (!forceCalc && audioDetectorParameters.attackTime_mSec == attack_in_ms)
		return;

	audioDetectorParameters.attackTime_mSec = attack_in_ms;
	attackTime = exp(TLD_AUDIO_ENVELOPE_ANALOG_TC / (attack_in_ms * sampleRate * 0.001));
}


/**
\brief sets the new release time and re-calculates the time constant
\param release_in_ms the new relase timme
\param forceCalc flag to force a re-calculation of time constant even if values have not changed.
*/
void AudioDetector::setReleaseTime(double release_in_ms, bool forceCalc)
{
	if (!forceCalc && audioDetectorParameters.releaseTime_mSec == release_in_ms)
		return;

	audioDetectorParameters.releaseTime_mSec = release_in_ms;
	releaseTime = exp(TLD_AUDIO_ENVELOPE_ANALOG_TC / (release_in_ms * sampleRate * 0.001));
}
