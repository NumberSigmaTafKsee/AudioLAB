#pragma once
#include "FXObjects.hpp"

/**
\class ImpulseConvolver
\ingroup FX-Objects
\brief
The ImpulseConvolver object implements a linear conovlver. NOTE: compile in Release mode or you may experice stuttering,
glitching or other sample-drop activity.
Audio I/O:
- Processes mono input to mono output.
Control I/F:
- none.
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
class ImpulseConvolver : public IAudioSignalProcessor
{
public:
	ImpulseConvolver() {
		init(512);
	}		/* C-TOR */
	~ImpulseConvolver() {}		/* D-TOR */

	/** reset members to initialized state */
	virtual bool reset(double _sampleRate)
	{
		// --- flush signal buffer; IR buffer is static
		signalBuffer.flushBuffer();
		return true;
	}

	/** process one input - note this is CPU intensive as it performs simple linear convolution */
	/**
	\param xn input
	\return the processed sample
	*/
	virtual double processAudioSample(double xn)
	{
		double output = 0.0;

		// --- write buffer; x(n) overwrites oldest value
		//     this is the only time we do not read before write!
		signalBuffer.writeBuffer(xn);

		// --- do the convolution
		for (unsigned int i = 0; i < length; i++)
		{
			// --- y(n) += x(n)h(n)
			//     for signalBuffer.readBuffer(0) -> x(n)
			//		   signalBuffer.readBuffer(n-D)-> x(n-D)
			double signal = signalBuffer.readBuffer((int)i);
			double irrrrr = irBuffer.readBuffer((int)i);
			output += signal*irrrrr;

			//output += signalBuffer.readBuffer((int)i) * irBuffer.readBuffer((int)i);
		}

		return output;
	}

	/** return false: this object only processes samples */
	virtual bool canProcessAudioFrame() { return false; }

	/** create the buffer based on the exact power of 2 */
	void init(unsigned int lengthPowerOfTwo)
	{
		length = lengthPowerOfTwo;
		// --- create (and clear out) the buffers
		signalBuffer.createCircularBufferPowerOfTwo(lengthPowerOfTwo);
		irBuffer.createLinearBuffer(lengthPowerOfTwo);
	}

	/** set the impulse response */
	void setImpulseResponse(double* irArray, unsigned int lengthPowerOfTwo)
	{
		if (lengthPowerOfTwo != length)
		{
			length = lengthPowerOfTwo;
			// --- create (and clear out) the buffers
			signalBuffer.createCircularBufferPowerOfTwo(lengthPowerOfTwo);
			irBuffer.createLinearBuffer(lengthPowerOfTwo);
		}

		// --- load up the IR buffer
		for (unsigned int i = 0; i < length; i++)
		{
			irBuffer.writeBuffer(i, irArray[i]);
		}
	}

protected:
	// --- delay buffer of doubles
	CircularBuffer<double> signalBuffer; ///< circulat buffer for the signal
	LinearBuffer<double> irBuffer;	///< linear buffer for the IR

	unsigned int length = 0;	///< length of convolution (buffer)

};

const unsigned int IR_LEN = 512;
/**
\struct AnalogFIRFilterParameters
\ingroup FX-Objects
\brief
Custom parameter structure for the AnalogFIRFilter object. This is a somewhat silly object that implaments an analog
magnitude response as a FIR filter. NOT DESIGNED to replace virtual analog; rather it is intended to show the
frequency sampling method in an easy (and fun) way.
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
struct AnalogFIRFilterParameters
{
	AnalogFIRFilterParameters() {}
	/** all FXObjects parameter objects require overloaded= operator so remember to add new entries if you add new variables. */
	AnalogFIRFilterParameters& operator=(const AnalogFIRFilterParameters& params)	// need this override for collections to work
	{
		if (this == &params)
			return *this;

		filterType = params.filterType;
		fc = params.fc;
		Q = params.Q;

		return *this;
	}

	// --- individual parameters
	analogFilter filterType = analogFilter::kLPF1; ///< filter type
	double fc = 0.0;	///< filter fc
	double Q = 0.0;		///< filter Q
};

/**
\class AnalogFIRFilter
\ingroup FX-Objects
\brief
The AnalogFIRFilter object implements a somewhat silly algorithm that implaments an analog
magnitude response as a FIR filter. NOT DESIGNED to replace virtual analog; rather it is intended to show the
frequency sampling method in an easy (and fun) way.
Audio I/O:
- Processes mono input to mono output.
Control I/F:
- Use AnalogFIRFilterParameters structure to get/set object params.
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
class AnalogFIRFilter : public IAudioSignalProcessor
{
public:
	AnalogFIRFilter() {}	/* C-TOR */
	~AnalogFIRFilter() {}	/* D-TOR */

public:
	/** reset members to initialized state */
	virtual bool reset(double _sampleRate)
	{
		sampleRate = _sampleRate;
		convolver.reset(_sampleRate);
		convolver.init(IR_LEN);

		memset(&analogMagArray[0], 0, sizeof(double) * IR_LEN);	///< clear
		memset(&irArray[0], 0, sizeof(double) * IR_LEN);	///< clear

		return true;
	}

	/** pefrorm the convolution */
	/**
	\param xn input
	\return the processed sample
	*/
	virtual double processAudioSample(double xn)
	{
		// --- do the linear convolution
		return convolver.processAudioSample(xn);
	}

	/** return false: this object only processes samples */
	virtual bool canProcessAudioFrame() { return false; }

	/** get parameters: note use of custom structure for passing param data */
	/**
	\return AnalogFIRFilterParameters custom data structure
	*/
	AnalogFIRFilterParameters getParameters() { return parameters; }

	/** set parameters: note use of custom structure for passing param data */
	/**
	\param AnalogFIRFilterParameters custom data structure
	*/
	void setParameters(AnalogFIRFilterParameters _parameters)
	{
		if (_parameters.fc != parameters.fc ||
			_parameters.Q != parameters.Q ||
			_parameters.filterType != parameters.filterType)
		{
			// --- set the filter IR for the convolver
			AnalogMagData analogFilterData;
			analogFilterData.sampleRate = sampleRate;
			analogFilterData.magArray = &analogMagArray[0];
			analogFilterData.dftArrayLen = IR_LEN;
			analogFilterData.mirrorMag = false;

			analogFilterData.filterType = _parameters.filterType;
			analogFilterData.fc = _parameters.fc; // 1000.0;
			analogFilterData.Q = _parameters.Q;

			// --- calculate the analog mag array
			calculateAnalogMagArray(analogFilterData);

			// --- frequency sample the mag array
			freqSample(IR_LEN, analogMagArray, irArray, POSITIVE);

			// --- update new frequency response
			convolver.setImpulseResponse(irArray, IR_LEN);
		}

		parameters = _parameters;
	}

private:
	AnalogFIRFilterParameters parameters; ///< object parameters
	ImpulseConvolver convolver; ///< convolver object to perform FIR convolution
	double analogMagArray[IR_LEN]; ///< array for analog magnitude response
	double irArray[IR_LEN]; ///< array to hold calcualted IR
	double sampleRate = 0.0; ///< storage for sample rate
};


