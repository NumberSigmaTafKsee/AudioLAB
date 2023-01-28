#pragma once

#include "FXObjects.hpp"
/**
\enum delayAlgorithm
\ingroup Constants-Enums
\brief
Use this strongly typed enum to easily set the delay algorithm
- enum class delayAlgorithm { kNormal, kPingPong };
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
enum class delayAlgorithm { kNormal, kPingPong };

/**
\enum delayUpdateType
\ingroup Constants-Enums
\brief
Use this strongly typed enum to easily set the delay update type; this varies depending on the designer's choice
of GUI controls. See the book reference for more details.
- enum class delayUpdateType { kLeftAndRight, kLeftPlusRatio };
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
enum class delayUpdateType { kLeftAndRight, kLeftPlusRatio };


/**
\struct AudioDelayParameters
\ingroup FX-Objects
\brief
Custom parameter structure for the AudioDelay object.
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
struct AudioDelayParameters
{
	AudioDelayParameters() {}
	/** all FXObjects parameter objects require overloaded= operator so remember to add new entries if you add new variables. */
	AudioDelayParameters& operator=(const AudioDelayParameters& params)	// need this override for collections to work
	{
		if (this == &params)
			return *this;

		algorithm = params.algorithm;
		wetLevel_dB = params.wetLevel_dB;
		dryLevel_dB = params.dryLevel_dB;
		feedback_Pct = params.feedback_Pct;

		updateType = params.updateType;
		leftDelay_mSec = params.leftDelay_mSec;
		rightDelay_mSec = params.rightDelay_mSec;
		delayRatio_Pct = params.delayRatio_Pct;

		return *this;
	}

	// --- individual parameters
	delayAlgorithm algorithm = delayAlgorithm::kNormal; ///< delay algorithm
	DspFloatType wetLevel_dB = -3.0;	///< wet output level in dB
	DspFloatType dryLevel_dB = -3.0;	///< dry output level in dB
	DspFloatType feedback_Pct = 0.0;	///< feedback as a % value

	delayUpdateType updateType = delayUpdateType::kLeftAndRight;///< update algorithm
	DspFloatType leftDelay_mSec = 0.0;	///< left delay time
	DspFloatType rightDelay_mSec = 0.0;	///< right delay time
	DspFloatType delayRatio_Pct = 100.0;	///< dela ratio: right length = (delayRatio)*(left length)
};

/**
\class AudioDelay
\ingroup FX-Objects
\brief
The AudioDelay object implements a stereo audio delay with multiple delay algorithms.
Audio I/O:
- Processes mono input to mono output OR stereo output.
Control I/F:
- Use AudioDelayParameters structure to get/set object params.
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
class AudioDelay : public IAudioSignalProcessor
{
public:
	AudioDelay() {}		/* C-TOR */
	~AudioDelay() {}	/* D-TOR */

public:
	/** reset members to initialized state */
	virtual bool reset(DspFloatType _sampleRate)
	{
		// --- if sample rate did not change
		if (sampleRate == _sampleRate)
		{
			// --- just flush buffer and return
			delayBuffer_L.flushBuffer();
			delayBuffer_R.flushBuffer();
			return true;
		}

		// --- create new buffer, will store sample rate and length(mSec)
		createDelayBuffers(_sampleRate, bufferLength_mSec);

		return true;
	}

	/** process MONO audio delay */
	/**
	\param xn input
	\return the processed sample
	*/
	virtual DspFloatType processAudioSample(DspFloatType xn)
	{
		// --- read delay
		DspFloatType yn = delayBuffer_L.readBuffer(delayInSamples_L);

		// --- create input for delay buffer
		DspFloatType dn = xn + (parameters.feedback_Pct / 100.0) * yn;

		// --- write to delay buffer
		delayBuffer_L.writeBuffer(dn);

		// --- form mixture out = dry*xn + wet*yn
		DspFloatType output = dryMix*xn + wetMix*yn;

		return output;
	}

	/** return true: this object can also process frames */
	virtual bool canProcessAudioFrame() { return true; }

	/** process STEREO audio delay in frames */
	virtual bool processAudioFrame(const DspFloatType* inputFrame,		/* ptr to one frame of data: pInputFrame[0] = left, pInputFrame[1] = right, etc...*/
		DspFloatType* outputFrame,
		uint32_t inputChannels,
		uint32_t outputChannels)
	{
		// --- make sure we have input and outputs
		if (inputChannels == 0 || outputChannels == 0)
			return false;

		// --- make sure we support this delay algorithm
		if (parameters.algorithm != delayAlgorithm::kNormal &&
			parameters.algorithm != delayAlgorithm::kPingPong)
			return false;

		// --- if only one output channel, revert to mono operation
		if (outputChannels == 1)
		{
			// --- process left channel only
			outputFrame[0] = processAudioSample(inputFrame[0]);
			return true;
		}

		// --- if we get here we know we have 2 output channels
		//
		// --- pick up inputs
		//
		// --- LEFT channel
		DspFloatType xnL = inputFrame[0];

		// --- RIGHT channel (duplicate left input if mono-in)
		DspFloatType xnR = inputChannels > 1 ? inputFrame[1] : xnL;

		// --- read delay LEFT
		DspFloatType ynL = delayBuffer_L.readBuffer(delayInSamples_L);

		// --- read delay RIGHT
		DspFloatType ynR = delayBuffer_R.readBuffer(delayInSamples_R);

		// --- create input for delay buffer with LEFT channel info
		DspFloatType dnL = xnL + (parameters.feedback_Pct / 100.0) * ynL;

		// --- create input for delay buffer with RIGHT channel info
		DspFloatType dnR = xnR + (parameters.feedback_Pct / 100.0) * ynR;

		// --- decode
		if (parameters.algorithm == delayAlgorithm::kNormal)
		{
			// --- write to LEFT delay buffer with LEFT channel info
			delayBuffer_L.writeBuffer(dnL);

			// --- write to RIGHT delay buffer with RIGHT channel info
			delayBuffer_R.writeBuffer(dnR);
		}
		else if (parameters.algorithm == delayAlgorithm::kPingPong)
		{
			// --- write to LEFT delay buffer with RIGHT channel info
			delayBuffer_L.writeBuffer(dnR);

			// --- write to RIGHT delay buffer with LEFT channel info
			delayBuffer_R.writeBuffer(dnL);
		}

		// --- form mixture out = dry*xn + wet*yn
		DspFloatType outputL = dryMix*xnL + wetMix*ynL;

		// --- form mixture out = dry*xn + wet*yn
		DspFloatType outputR = dryMix*xnR + wetMix*ynR;

		// --- set left channel
		outputFrame[0] = outputL;

		// --- set right channel
		outputFrame[1] = outputR;

		return true;
	}

	void processMono(size_t n, DspFloatType * in, DspFloatType * out)
	{
		for(size_t i = 0; i < n; i++) processAudioFrame(&in[i],&out[i],1,1);
	}
	void processStereo(size_t n, DspFloatType * in, DspFloatType * out)
	{
		DspFloatType frame[2];
		size_t x=0;
		for(size_t i = 0; i < n; i++)
		{
			processAudioFrame(&in[i],frame,2,2);
			out[x++] = frame[0];
			out[x++] = frame[1];
		}
	}
	/** get parameters: note use of custom structure for passing param data */
	/**
	\return AudioDelayParameters custom data structure
	*/
	AudioDelayParameters getParameters() { return parameters; }

	/** set parameters: note use of custom structure for passing param data */
	/**
	\param AudioDelayParameters custom data structure
	*/
	void setParameters(AudioDelayParameters _parameters)
	{
		// --- check mix in dB for calc
		if (_parameters.dryLevel_dB != parameters.dryLevel_dB)
			dryMix = pow(10.0, _parameters.dryLevel_dB / 20.0);
		if (_parameters.wetLevel_dB != parameters.wetLevel_dB)
			wetMix = pow(10.0, _parameters.wetLevel_dB / 20.0);

		// --- save; rest of updates are cheap on CPU
		parameters = _parameters;

		// --- check update type first:
		if (parameters.updateType == delayUpdateType::kLeftAndRight)
		{
			// --- set left and right delay times
			// --- calculate total delay time in samples + fraction
			DspFloatType newDelayInSamples_L = parameters.leftDelay_mSec*(samplesPerMSec);
			DspFloatType newDelayInSamples_R = parameters.rightDelay_mSec*(samplesPerMSec);

			// --- new delay time with fraction
			delayInSamples_L = newDelayInSamples_L;
			delayInSamples_R = newDelayInSamples_R;
		}
		else if (parameters.updateType == delayUpdateType::kLeftPlusRatio)
		{
			// --- get and validate ratio
			DspFloatType delayRatio = parameters.delayRatio_Pct / 100.0;
			boundValue(delayRatio, 0.0, 1.0);

			// --- calculate total delay time in samples + fraction
			DspFloatType newDelayInSamples = parameters.leftDelay_mSec*(samplesPerMSec);

			// --- new delay time with fraction
			delayInSamples_L = newDelayInSamples;
			delayInSamples_R = delayInSamples_L*delayRatio;
		}
	}

	/** creation function */
	void createDelayBuffers(DspFloatType _sampleRate, DspFloatType _bufferLength_mSec)
	{
		// --- store for math
		bufferLength_mSec = _bufferLength_mSec;
		sampleRate = _sampleRate;
		samplesPerMSec = sampleRate / 1000.0;

		// --- total buffer length including fractional part
		bufferLength = (unsigned int)(bufferLength_mSec*(samplesPerMSec)) + 1; // +1 for fractional part

																			   // --- create new buffer
		delayBuffer_L.createCircularBuffer(bufferLength);
		delayBuffer_R.createCircularBuffer(bufferLength);
	}

private:
	AudioDelayParameters parameters; ///< object parameters

	DspFloatType sampleRate = 0.0;		///< current sample rate
	DspFloatType samplesPerMSec = 0.0;	///< samples per millisecond, for easy access calculation
	DspFloatType delayInSamples_L = 0.0;	///< DspFloatType includes fractional part
	DspFloatType delayInSamples_R = 0.0;	///< DspFloatType includes fractional part
	DspFloatType bufferLength_mSec = 0.0;	///< buffer length in mSec
	unsigned int bufferLength = 0;	///< buffer length in samples
	DspFloatType wetMix = 0.707; ///< wet output default = -3dB
	DspFloatType dryMix = 0.707; ///< dry output default = -3dB

	// --- delay buffer of doubles
	CircularBuffer<DspFloatType> delayBuffer_L;	///< LEFT delay buffer of doubles
	CircularBuffer<DspFloatType> delayBuffer_R;	///< RIGHT delay buffer of doubles
};


/**
\enum modDelaylgorithm
\ingroup Constants-Enums
\brief
Use this strongly typed enum to easily set modulated delay algorithm.
- enum class modDelaylgorithm { kFlanger, kChorus, kVibrato };
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
enum class modDelaylgorithm { kFlanger, kChorus, kVibrato };


/**
\struct ModulatedDelayParameters
\ingroup FX-Objects
\brief
Custom parameter structure for the ModulatedDelay object.
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
struct ModulatedDelayParameters
{
	ModulatedDelayParameters() {}
	/** all FXObjects parameter objects require overloaded= operator so remember to add new entries if you add new variables. */
	ModulatedDelayParameters& operator=(const ModulatedDelayParameters& params)	// need this override for collections to work
	{
		if (this == &params)
			return *this;

		algorithm = params.algorithm;
		lfoRate_Hz = params.lfoRate_Hz;
		lfoDepth_Pct = params.lfoDepth_Pct;
		feedback_Pct = params.feedback_Pct;
		return *this;
	}

	// --- individual parameters
	modDelaylgorithm algorithm = modDelaylgorithm::kFlanger; ///< mod delay algorithm
	DspFloatType lfoRate_Hz = 0.0;	///< mod delay LFO rate in Hz
	DspFloatType lfoDepth_Pct = 0.0;	///< mod delay LFO depth in %
	DspFloatType feedback_Pct = 0.0;	///< feedback in %
};

/**
\class ModulatedDelay
\ingroup FX-Objects
\brief
The ModulatedDelay object implements the three basic algorithms: flanger, chorus, vibrato.
Audio I / O :
	-Processes mono input to mono OR stereo output.
Control I / F :
	-Use ModulatedDelayParameters structure to get / set object params.
\author Will Pirkle http ://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed.by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
class ModulatedDelay : public IAudioSignalProcessor
{
public:
	ModulatedDelay() {
	}		/* C-TOR */
	~ModulatedDelay() {}	/* D-TOR */

public:
	/** reset members to initialized state */
	virtual bool reset(DspFloatType _sampleRate)
	{
		// --- create new buffer, 100mSec long
		delay.reset(_sampleRate);
		delay.createDelayBuffers(_sampleRate, 100.0);

		// --- lfo
		lfo.reset(_sampleRate);
		OscillatorParameters params = lfo.getParameters();
		params.waveform = generatorWaveform::kTriangle;
		lfo.setParameters(params);

		return true;
	}

	/** process input sample */
	/**
	\param xn input
	\return the processed sample
	*/
	virtual DspFloatType processAudioSample(DspFloatType xn)
	{
		DspFloatType input = xn;
		DspFloatType output = 0.0;
		processAudioFrame(&input, &output, 1, 1);
		return output;
	}

	/** return true: this object can process frames */
	virtual bool canProcessAudioFrame() { return true; }

	/** process STEREO audio delay of frames */
	virtual bool processAudioFrame(const DspFloatType* inputFrame,		/* ptr to one frame of data: pInputFrame[0] = left, pInputFrame[1] = right, etc...*/
		DspFloatType* outputFrame,
		uint32_t inputChannels,
		uint32_t outputChannels)
	{
		// --- make sure we have input and outputs
		if (inputChannels == 0 || outputChannels == 0)
			return false;

		// --- render LFO
		SignalGenData lfoOutput = lfo.renderAudioOutput();

		// --- setup delay modulation
		AudioDelayParameters params = delay.getParameters();
		DspFloatType minDelay_mSec = 0.0;
		DspFloatType maxDepth_mSec = 0.0;

		// --- set delay times, wet/dry and feedback
		if (parameters.algorithm == modDelaylgorithm::kFlanger)
		{
			minDelay_mSec = 0.1;
			maxDepth_mSec = 7.0;
			params.wetLevel_dB = -3.0;
			params.dryLevel_dB = -3.0;
		}
		if (parameters.algorithm == modDelaylgorithm::kChorus)
		{
			minDelay_mSec = 10.0;
			maxDepth_mSec = 30.0;
			params.wetLevel_dB = -3.0;
			params.dryLevel_dB = -0.0;
			params.feedback_Pct = 0.0;
		}
		if (parameters.algorithm == modDelaylgorithm::kVibrato)
		{
			minDelay_mSec = 0.0;
			maxDepth_mSec = 7.0;
			params.wetLevel_dB = 0.0;
			params.dryLevel_dB = -96.0;
			params.feedback_Pct = 0.0;
		}

		// --- calc modulated delay times
		DspFloatType depth = parameters.lfoDepth_Pct / 100.0;
		DspFloatType modulationMin = minDelay_mSec;
		DspFloatType modulationMax = minDelay_mSec + maxDepth_mSec;

		// --- flanger - unipolar
		if (parameters.algorithm == modDelaylgorithm::kFlanger)
			params.leftDelay_mSec = doUnipolarModulationFromMin(bipolarToUnipolar(depth * lfoOutput.normalOutput),
															     modulationMin, modulationMax);
		else
			params.leftDelay_mSec = doBipolarModulation(depth * lfoOutput.normalOutput, modulationMin, modulationMax);


		// --- set right delay to match (*Hint Homework!)
		params.rightDelay_mSec = params.leftDelay_mSec;

		// --- modulate the delay
		delay.setParameters(params);

		// --- just call the function and pass our info in/out
		return delay.processAudioFrame(inputFrame, outputFrame, inputChannels, outputChannels);
	}

	void processMono(size_t n, DspFloatType * in, DspFloatType * out)
	{
		for(size_t i = 0; i < n; i++) processAudioFrame(&in[i],&out[i],1,1);
	}
	void processStereo(size_t n, DspFloatType * in, DspFloatType * out)
	{
		DspFloatType frame[2];
		size_t x=0;
		for(size_t i = 0; i < n; i++)
		{
			processAudioFrame(&in[i],frame,2,2);
			out[x++] = frame[0];
			out[x++] = frame[1];
		}
	}
	
	/** get parameters: note use of custom structure for passing param data */
	/**
	\return ModulatedDelayParameters custom data structure
	*/
	ModulatedDelayParameters getParameters() { return parameters; }

	/** set parameters: note use of custom structure for passing param data */
	/**
	\param ModulatedDelayParameters custom data structure
	*/
	void setParameters(ModulatedDelayParameters _parameters)
	{
		// --- bulk copy
		parameters = _parameters;

		OscillatorParameters lfoParams = lfo.getParameters();
		lfoParams.frequency_Hz = parameters.lfoRate_Hz;
		if (parameters.algorithm == modDelaylgorithm::kVibrato)
			lfoParams.waveform = generatorWaveform::kSin;
		else
			lfoParams.waveform = generatorWaveform::kTriangle;

		lfo.setParameters(lfoParams);

		AudioDelayParameters adParams = delay.getParameters();
		adParams.feedback_Pct = parameters.feedback_Pct;
		delay.setParameters(adParams);
	}

private:
	ModulatedDelayParameters parameters; ///< object parameters
	AudioDelay delay;	///< the delay to modulate
	LFO lfo;			///< the modulator
};


/**
\struct SimpleDelayParameters
\ingroup FX-Objects
\brief
Custom parameter structure for the SimpleDelay object. Used for reverb algorithms in book.
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
struct SimpleDelayParameters
{
	SimpleDelayParameters() {}
	/** all FXObjects parameter objects require overloaded= operator so remember to add new entries if you add new variables. */
	SimpleDelayParameters& operator=(const SimpleDelayParameters& params)	// need this override for collections to work
	{
		if (this == &params)
			return *this;

		delayTime_mSec = params.delayTime_mSec;
		interpolate = params.interpolate;
		delay_Samples = params.delay_Samples;
		return *this;
	}

	// --- individual parameters
	DspFloatType delayTime_mSec = 0.0;	///< delay tine in mSec
	bool interpolate = false;		///< interpolation flag (diagnostics usually)

	// --- outbound parameters
	DspFloatType delay_Samples = 0.0;		///< current delay in samples; other objects may need to access this information
};

/**
\class SimpleDelay
\ingroup FX-Objects
\brief
The SimpleDelay object implements a basic delay line without feedback.
Audio I/O:
- Processes mono input to mono output.
Control I/F:
- Use SimpleDelayParameters structure to get/set object params.
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
class SimpleDelay : public IAudioSignalProcessor
{
public:
	SimpleDelay(void) {}	/* C-TOR */
	~SimpleDelay(void) {}	/* D-TOR */

public:
	/** reset members to initialized state */
	virtual bool reset(DspFloatType _sampleRate)
	{
		// --- if sample rate did not change
		if (sampleRate == _sampleRate)
		{
			// --- just flush buffer and return
			delayBuffer.flushBuffer();
			return true;
		}

		// --- create new buffer, will store sample rate and length(mSec)
		createDelayBuffer(_sampleRate, bufferLength_mSec);

		return true;
	}

	/** get parameters: note use of custom structure for passing param data */
	/**
	\return SimpleDelayParameters custom data structure
	*/
	SimpleDelayParameters getParameters()
	{
		return simpleDelayParameters;
	}

	/** set parameters: note use of custom structure for passing param data */
	/**
	\param SimpleDelayParameters custom data structure
	*/
	void setParameters(const SimpleDelayParameters& params)
	{
		simpleDelayParameters = params;
		simpleDelayParameters.delay_Samples = simpleDelayParameters.delayTime_mSec*(samplesPerMSec);
		delayBuffer.setInterpolate(simpleDelayParameters.interpolate);
	}

	/** process MONO audio delay */
	/**
	\param xn input
	\return the processed sample
	*/
	virtual DspFloatType processAudioSample(DspFloatType xn)
	{
		// --- read delay
		if (simpleDelayParameters.delay_Samples == 0)
			return xn;

		DspFloatType yn = delayBuffer.readBuffer(simpleDelayParameters.delay_Samples);

		// --- write to delay buffer
		delayBuffer.writeBuffer(xn);

		// --- done
		return yn;
	}

	/** reset members to initialized state */
	virtual bool canProcessAudioFrame() { return false; }

	/** create a new delay buffer */
	void createDelayBuffer(DspFloatType _sampleRate, DspFloatType _bufferLength_mSec)
	{
		// --- store for math
		bufferLength_mSec = _bufferLength_mSec;
		sampleRate = _sampleRate;
		samplesPerMSec = sampleRate / 1000.0;

		// --- total buffer length including fractional part
		bufferLength = (unsigned int)(bufferLength_mSec*(samplesPerMSec)) + 1; // +1 for fractional part

		// --- create new buffer
		delayBuffer.createCircularBuffer(bufferLength);
	}

	/** read delay at current location */
	DspFloatType readDelay()
	{
		// --- simple read
		return delayBuffer.readBuffer(simpleDelayParameters.delay_Samples);
	}

	/** read delay at current location */
	DspFloatType readDelayAtTime_mSec(DspFloatType _delay_mSec)
	{
		// --- calculate total delay time in samples + fraction
		DspFloatType _delay_Samples = _delay_mSec*(samplesPerMSec);

		// --- simple read
		return delayBuffer.readBuffer(_delay_Samples);
	}

	/** read delay at a percentage of total length */
	DspFloatType readDelayAtPercentage(DspFloatType delayPercent)
	{
		// --- simple read
		return delayBuffer.readBuffer((delayPercent / (DspFloatType)100.0)*simpleDelayParameters.delay_Samples);
	}

	/** write a new value into the delay */
	void writeDelay(DspFloatType xn)
	{
		// --- simple write
		delayBuffer.writeBuffer(xn);
	}

private:
	SimpleDelayParameters simpleDelayParameters; ///< object parameters

	DspFloatType sampleRate = 0.0;		///< sample rate
	DspFloatType samplesPerMSec = 0.0;	///< samples per millisecond (for arbitrary access)
	DspFloatType bufferLength_mSec = 0.0; ///< total buffer lenth in mSec
	unsigned int bufferLength = 0;	///< buffer length in samples

	// --- delay buffer of doubles
	CircularBuffer<DspFloatType> delayBuffer; ///< circular buffer for delay
};


/**
\struct CombFilterParameters
\ingroup FX-Objects
\brief
Custom parameter structure for the CombFilter object. Used for reverb algorithms in book.
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
struct CombFilterParameters
{
	CombFilterParameters() {}
	/** all FXObjects parameter objects require overloaded= operator so remember to add new entries if you add new variables. */
	CombFilterParameters& operator=(const CombFilterParameters& params)	// need this override for collections to work
	{
		if (this == &params)
			return *this;

		delayTime_mSec = params.delayTime_mSec;
		RT60Time_mSec = params.RT60Time_mSec;
		enableLPF = params.enableLPF;
		lpf_g = params.lpf_g;
		interpolate = params.interpolate;
		return *this;
	}

	// --- individual parameters
	DspFloatType delayTime_mSec = 0.0;	///< delay time in mSec
	DspFloatType RT60Time_mSec = 0.0;		///< RT 60 time ini mSec
	bool enableLPF = false;			///< enable LPF flag
	DspFloatType lpf_g = 0.0;				///< gain value for LPF (if enabled)
	bool interpolate = false;		///< interpolation flag (diagnostics)
};


/**
\class CombFilter
\ingroup FX-Objects
\brief
The CombFilter object implements a comb filter with optional LPF in feedback loop. Used for reverb algorithms in book.
Audio I/O:
- Processes mono input to mono output.
Control I/F:
- Use CombFilterParameters structure to get/set object params.
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
class CombFilter : public IAudioSignalProcessor
{
public:
	CombFilter(void) {}		/* C-TOR */
	~CombFilter(void) {}	/* D-TOR */

public:
	/** reset members to initialized state */
	virtual bool reset(DspFloatType _sampleRate)
	{
		// --- flush
		lpf_state = 0.0;

		// --- create new buffer, will store sample rate and length(mSec)
		createDelayBuffer(sampleRate, bufferLength_mSec);

		return true;
	}

	/** process CombFilter */
	/**
	\param xn input
	\return the processed sample
	*/
	virtual DspFloatType processAudioSample(DspFloatType xn)
	{
		DspFloatType yn = delay.readDelay();
		DspFloatType input = 0.0;

		// --- form input & write
		if (combFilterParameters.enableLPF)
		{
			// --- apply simple 1st order pole LPF
			DspFloatType g2 = lpf_g*(1.0 - comb_g); // see book for equation 11.27 (old book)
			DspFloatType filteredSignal = yn + g2*lpf_state;
			input = xn + comb_g*(filteredSignal);
			lpf_state = filteredSignal;
		}
		else
		{
			input = xn + comb_g*yn;
		}

		delay.writeDelay(input);

		// --- done
		return yn;
	}

	/** return false: this object only processes samples */
	virtual bool canProcessAudioFrame() { return false; }

	/** get parameters: note use of custom structure for passing param data */
	/**
	\return CombFilterParameters custom data structure
	*/
	CombFilterParameters getParameters()
	{
		return combFilterParameters;
	}

	/** set parameters: note use of custom structure for passing param data */
	/**
	\param CombFilterParameters custom data structure
	*/
	void setParameters(const CombFilterParameters& params)
	{
		combFilterParameters = params;

		// --- update the delay line parameters first
		SimpleDelayParameters delayParams = delay.getParameters();
		delayParams.delayTime_mSec = combFilterParameters.delayTime_mSec;
		delayParams.interpolate = combFilterParameters.interpolate;
		delay.setParameters(delayParams); // this will set the delay time in samples

										  // --- calculate g with RT60 time (requires updated delay above^^)
		DspFloatType exponent = -3.0*delayParams.delay_Samples*(1.0 / sampleRate);
		DspFloatType rt60_mSec = combFilterParameters.RT60Time_mSec / 1000.0; // RT is in mSec!
		comb_g = pow(10.0, exponent / rt60_mSec);

		// --- set LPF g
		lpf_g = combFilterParameters.lpf_g;
	}

	/** create new buffers */
	void createDelayBuffer(DspFloatType _sampleRate, DspFloatType delay_mSec)
	{
		sampleRate = _sampleRate;
		bufferLength_mSec = delay_mSec;

		// --- create new buffer, will store sample rate and length(mSec)
		delay.createDelayBuffer(_sampleRate, delay_mSec);
	}

private:
	CombFilterParameters combFilterParameters; ///< object parameters
	DspFloatType sampleRate = 0.0;	///< sample rate
	DspFloatType comb_g = 0.0;		///< g value for comb filter (feedback, not %)
	DspFloatType bufferLength_mSec = 0.0; ///< total buffer length

	// --- LPF support
	DspFloatType lpf_g = 0.0;		///< LPF g value
	DspFloatType lpf_state = 0.0; ///< state register for LPF (z^-1)

	// --- delay buffer of doubles
	SimpleDelay delay;		///< delay for comb filter
};

/**
\struct DelayAPFParameters
\ingroup FX-Objects
\brief
Custom parameter structure for the DelayAPF object. Used for reverb algorithms in book.
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
struct DelayAPFParameters
{
	DelayAPFParameters() {}
	/** all FXObjects parameter objects require overloaded= operator so remember to add new entries if you add new variables. */
	DelayAPFParameters& operator=(const DelayAPFParameters& params)	// need this override for collections to work
	{
		if (this == &params)
			return *this;

		delayTime_mSec = params.delayTime_mSec;
		apf_g = params.apf_g;
		enableLPF = params.enableLPF;
		lpf_g = params.lpf_g;
		interpolate = params.interpolate;
		enableLFO = params.enableLFO;
		lfoRate_Hz = params.lfoRate_Hz;
		lfoDepth = params.lfoDepth;
		lfoMaxModulation_mSec = params.lfoMaxModulation_mSec;
		return *this;
	}

	// --- individual parameters
	DspFloatType delayTime_mSec = 0.0;	///< APF delay time
	DspFloatType apf_g = 0.0;				///< APF g coefficient
	bool enableLPF = false;			///< flag to enable LPF in structure
	DspFloatType lpf_g = 0.0;				///< LPF g coefficient (if enabled)
	bool interpolate = false;		///< interpolate flag (diagnostics)
	bool enableLFO = false;			///< flag to enable LFO
	DspFloatType lfoRate_Hz = 0.0;		///< LFO rate in Hz, if enabled
	DspFloatType lfoDepth = 0.0;			///< LFO deoth (not in %) if enabled
	DspFloatType lfoMaxModulation_mSec = 0.0;	///< LFO maximum modulation time in mSec

};

/**
\class DelayAPF
\ingroup FX-Objects
\brief
The DelayAPF object implements a delaying APF with optional LPF and optional modulated delay time with LFO.
Audio I/O:
- Processes mono input to mono output.
Control I/F:
- Use DelayAPFParameters structure to get/set object params.
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
class DelayAPF : public IAudioSignalProcessor
{
public:
	DelayAPF(void) {}	/* C-TOR */
	~DelayAPF(void) {}	/* D-TOR */

public:
	/** reset members to initialized state */
	virtual bool reset(DspFloatType _sampleRate)
	{
		// --- reset children
		modLFO.reset(_sampleRate);

		// --- flush
		lpf_state = 0.0;

		// --- create new buffer, will store sample rate and length(mSec)
		createDelayBuffer(sampleRate, bufferLength_mSec);

		return true;
	}

	/** process one input sample through object */
	/**
	\param xn input
	\return the processed sample
	*/
	virtual DspFloatType processAudioSample(DspFloatType xn)
	{
		SimpleDelayParameters delayParams = delay.getParameters();
		if (delayParams.delay_Samples == 0)
			return xn;

		// --- delay line output
		DspFloatType wnD = 0.0;
		DspFloatType apf_g = delayAPFParameters.apf_g;
		DspFloatType lpf_g = delayAPFParameters.lpf_g;
		DspFloatType lfoDepth = delayAPFParameters.lfoDepth;

		// --- for modulated APFs
		if (delayAPFParameters.enableLFO)
		{
			SignalGenData lfoOutput = modLFO.renderAudioOutput();
			DspFloatType maxDelay = delayParams.delayTime_mSec;
			DspFloatType minDelay = maxDelay - delayAPFParameters.lfoMaxModulation_mSec;
			minDelay = fmax(0.0, minDelay); // bound minDelay to 0 as minimum

			// --- calc max-down modulated value with unipolar converted LFO output
			//     NOTE: LFO output is scaled by lfoDepth
			DspFloatType modDelay_mSec = doUnipolarModulationFromMax(bipolarToUnipolar(lfoDepth*lfoOutput.normalOutput),
				minDelay, maxDelay);

			// --- read modulated value to get w(n-D);
			wnD = delay.readDelayAtTime_mSec(modDelay_mSec);
		}
		else
			// --- read the delay line to get w(n-D)
			wnD = delay.readDelay();

		if (delayAPFParameters.enableLPF)
		{
			// --- apply simple 1st order pole LPF, overwrite wnD
			wnD = wnD*(1.0 - lpf_g) + lpf_g*lpf_state;
			lpf_state = wnD;
		}

		// form w(n) = x(n) + gw(n-D)
		DspFloatType wn = xn + apf_g*wnD;

		// form y(n) = -gw(n) + w(n-D)
		DspFloatType yn = -apf_g*wn + wnD;

		// underflow check
		checkFloatUnderflow(yn);

		// write delay line
		delay.writeDelay(wn);

		return yn;
	}

	/** return false: this object only processes samples */
	virtual bool canProcessAudioFrame() { return false; }

	/** get parameters: note use of custom structure for passing param data */
	/**
	\return DelayAPFParameters custom data structure
	*/
	DelayAPFParameters getParameters()
	{
		return delayAPFParameters;
	}

	/** set parameters: note use of custom structure for passing param data */
	/**
	\param DelayAPFParameters custom data structure
	*/
	void setParameters(const DelayAPFParameters& params)
	{
		delayAPFParameters = params;

		// --- update delay line
		SimpleDelayParameters delayParams = delay.getParameters();
		delayParams.delayTime_mSec = delayAPFParameters.delayTime_mSec;
		delay.setParameters(delayParams);
	}

	/** create the delay buffer in mSec */
	void createDelayBuffer(DspFloatType _sampleRate, DspFloatType delay_mSec)
	{
		sampleRate = _sampleRate;
		bufferLength_mSec = delay_mSec;

		// --- create new buffer, will store sample rate and length(mSec)
		delay.createDelayBuffer(_sampleRate, delay_mSec);
	}

protected:
	// --- component parameters
	DelayAPFParameters delayAPFParameters;	///< obeject parameters
	DspFloatType sampleRate = 0.0;				///< current sample rate
	DspFloatType bufferLength_mSec = 0.0;			///< total buffer length in mSec

	// --- delay buffer of doubles
	SimpleDelay delay;						///< delay

	// --- optional LFO
	LFO modLFO;								///< LFO

	// --- LPF support
	DspFloatType lpf_state = 0.0;					///< LPF state register (z^-1)
};


/**
\struct NestedDelayAPFParameters
\ingroup FX-Objects
\brief
Custom parameter structure for the NestedDelayAPF object. Used for reverb algorithms in book.
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
struct NestedDelayAPFParameters
{
	NestedDelayAPFParameters() {}
	/** all FXObjects parameter objects require overloaded= operator so remember to add new entries if you add new variables. */
	NestedDelayAPFParameters& operator=(const NestedDelayAPFParameters& params)	// need this override for collections to work
	{
		if (this == &params)
			return *this;

		outerAPFdelayTime_mSec = params.outerAPFdelayTime_mSec;
		innerAPFdelayTime_mSec = params.innerAPFdelayTime_mSec;
		outerAPF_g = params.outerAPF_g;
		innerAPF_g = params.innerAPF_g;

		// --- outer LFO
		enableLFO = params.enableLFO;
		lfoRate_Hz = params.lfoRate_Hz;
		lfoDepth = params.lfoDepth;
		lfoMaxModulation_mSec = params.lfoMaxModulation_mSec;

		return *this;
	}

	// --- individual parameters
	DspFloatType outerAPFdelayTime_mSec = 0.0;	///< delay time for outer APF
	DspFloatType innerAPFdelayTime_mSec = 0.0;	///< delay time for inner APF
	DspFloatType outerAPF_g = 0.0;				///< g coefficient for outer APF
	DspFloatType innerAPF_g = 0.0;				///< g coefficient for inner APF

	// --- this LFO belongs to the outer APF only
	bool enableLFO = false;					///< flag to enable the modulated delay
	DspFloatType lfoRate_Hz = 0.0;				///< LFO rate in Hz (if enabled)
	DspFloatType lfoDepth = 1.0;					///< LFO depth (not in %) (if enabled)
	DspFloatType lfoMaxModulation_mSec = 0.0;		///< max modulation time if LFO is enabled

};

/**
\class NestedDelayAPF
\ingroup FX-Objects
\brief
The NestedDelayAPF object implements a pair of nested Delaying APF structures. These are labled the
outer and inner APFs. The outer APF's LPF and LFO may be optionally enabled. You might want
to extend this object to enable and use the inner LPF and LFO as well.
Audio I/O:
- Processes mono input to mono output.
Control I/F:
- Use BiquadParameters structure to get/set object params.
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
class NestedDelayAPF : public DelayAPF
{
public:
	NestedDelayAPF(void) { }	/* C-TOR */
	~NestedDelayAPF(void) { }	/* D-TOR */

public:
	/** reset members to initialized state */
	virtual bool reset(DspFloatType _sampleRate)
	{
		// --- call base class reset first
		DelayAPF::reset(_sampleRate);

		// --- then do our stuff
		nestedAPF.reset(_sampleRate);

		return true;
	}

	/** process mono audio input */
	/**
	\param xn input
	\return the processed sample
	*/
	virtual DspFloatType processAudioSample(DspFloatType xn)
	{
		// --- delay line output
		DspFloatType wnD = 0.0;

		SimpleDelayParameters delayParams = delay.getParameters();
		if (delayParams.delay_Samples == 0)
			return xn;

		DspFloatType apf_g = delayAPFParameters.apf_g;
		DspFloatType lpf_g = delayAPFParameters.lpf_g;

		// --- for modulated APFs
		if (delayAPFParameters.enableLFO)
		{
			SignalGenData lfoOutput = modLFO.renderAudioOutput();
			DspFloatType maxDelay = delayParams.delayTime_mSec;
			DspFloatType minDelay = maxDelay - delayAPFParameters.lfoMaxModulation_mSec;
			minDelay = fmax(0.0, minDelay); // bound minDelay to 0 as minimum
			DspFloatType lfoDepth = delayAPFParameters.lfoDepth;

			// --- calc max-down modulated value with unipolar converted LFO output
			//     NOTE: LFO output is scaled by lfoDepth
			DspFloatType modDelay_mSec = doUnipolarModulationFromMax(bipolarToUnipolar(lfoDepth*lfoOutput.normalOutput),
				minDelay, maxDelay);

			// --- read modulated value to get w(n-D);
			wnD = delay.readDelayAtTime_mSec(modDelay_mSec);
		}
		else
			// --- read the delay line to get w(n-D)
			wnD = delay.readDelay();

		if (delayAPFParameters.enableLPF)
		{
			// --- apply simple 1st order pole LPF, overwrite wnD
			wnD = wnD*(1.0 - lpf_g) + lpf_g*lpf_state;
			lpf_state = wnD;
		}

		// --- form w(n) = x(n) + gw(n-D)
		DspFloatType wn = xn + apf_g*wnD;

		// --- process wn through inner APF
		DspFloatType ynInner = nestedAPF.processAudioSample(wn);

		// --- form y(n) = -gw(n) + w(n-D)
		DspFloatType yn = -apf_g*wn + wnD;

		// --- underflow check
		checkFloatUnderflow(yn);

		// --- write delay line
		delay.writeDelay(ynInner);

		return yn;
	}

	/** get parameters: note use of custom structure for passing param data */
	/**
	\return BiquadParameters custom data structure
	*/
	NestedDelayAPFParameters getParameters() { return nestedAPFParameters; }

	/** set parameters: note use of custom structure for passing param data */
	/**
	\param BiquadParameters custom data structure
	*/
	void setParameters(const NestedDelayAPFParameters& params)
	{
		nestedAPFParameters = params;

		DelayAPFParameters outerAPFParameters = DelayAPF::getParameters();
		DelayAPFParameters innerAPFParameters = nestedAPF.getParameters();

		// --- outer APF
		outerAPFParameters.apf_g = nestedAPFParameters.outerAPF_g;
		outerAPFParameters.delayTime_mSec = nestedAPFParameters.outerAPFdelayTime_mSec;

		// --- LFO support
		outerAPFParameters.enableLFO = nestedAPFParameters.enableLFO;
		outerAPFParameters.lfoDepth = nestedAPFParameters.lfoDepth;
		outerAPFParameters.lfoRate_Hz = nestedAPFParameters.lfoRate_Hz;
		outerAPFParameters.lfoMaxModulation_mSec = nestedAPFParameters.lfoMaxModulation_mSec;

		// --- inner APF
		innerAPFParameters.apf_g = nestedAPFParameters.innerAPF_g;
		innerAPFParameters.delayTime_mSec = nestedAPFParameters.innerAPFdelayTime_mSec;

		DelayAPF::setParameters(outerAPFParameters);
		nestedAPF.setParameters(innerAPFParameters);
	}

	/** createDelayBuffers -- note there are two delay times here for inner and outer APFs*/
	void createDelayBuffers(DspFloatType _sampleRate, DspFloatType delay_mSec, DspFloatType nestedAPFDelay_mSec)
	{
		// --- base class
		DelayAPF::createDelayBuffer(_sampleRate, delay_mSec);

		// --- then our stuff
		nestedAPF.createDelayBuffer(_sampleRate, nestedAPFDelay_mSec);
	}

private:
	NestedDelayAPFParameters nestedAPFParameters; ///< object parameters
	DelayAPF nestedAPF;	///< nested APF object
};
