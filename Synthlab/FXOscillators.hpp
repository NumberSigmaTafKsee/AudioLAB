#pragma once

#include "FXObjects.hpp"

/**
\enum generatorWaveform
\ingroup Constants-Enums
\brief
Use this strongly typed enum to easily set the oscillator waveform
- enum  generatorWaveform { kTriangle, kSin, kSaw };
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
enum class generatorWaveform { kTriangle, kSin, kSaw };

/**
\struct OscillatorParameters
\ingroup FX-Objects
\brief
Custom parameter structure for the LFO and DFOscillator objects.
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
struct OscillatorParameters
{
	OscillatorParameters() {}
	/** all FXObjects parameter objects require overloaded= operator so remember to add new entries if you add new variables. */
	OscillatorParameters& operator=(const OscillatorParameters& params)	// need this override for collections to work
	{
		if (this == &params)
			return *this;

		waveform = params.waveform;
		frequency_Hz = params.frequency_Hz;
		return *this;
	}

	// --- individual parameters
	generatorWaveform waveform = generatorWaveform::kTriangle; ///< the current waveform
	double frequency_Hz = 0.0;	///< oscillator frequency
};

/**
\class LFO
\ingroup FX-Objects
\brief
The LFO object implements a mathematically perfect LFO generator for modulation uses only. It should not be used for
audio frequencies except for the sinusoidal output which, though an approximation, has very low TDH.
Audio I/O:
- Output only object: low frequency generator.
Control I/F:
- Use OscillatorParameters structure to get/set object params.
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
class LFO : public IAudioSignalGenerator
{
public:
	LFO() {	srand(time(NULL)); }	/* C-TOR */
	virtual ~LFO() {}				/* D-TOR */

	/** reset members to initialized state */
	virtual bool reset(double _sampleRate)
	{
		sampleRate = _sampleRate;
		phaseInc = lfoParameters.frequency_Hz / sampleRate;

		// --- timebase variables
		modCounter = 0.0;			///< modulo counter [0.0, +1.0]
		modCounterQP = 0.25;		///<Quad Phase modulo counter [0.0, +1.0]

		return true;
	}

	/** get parameters: note use of custom structure for passing param data */
	/**
	\return OscillatorParameters custom data structure
	*/
	OscillatorParameters getParameters(){ return lfoParameters; }

	/** set parameters: note use of custom structure for passing param data */
	/**
	\param OscillatorParameters custom data structure
	*/
	void setParameters(const OscillatorParameters& params)
	{
		if(params.frequency_Hz != lfoParameters.frequency_Hz)
			// --- update phase inc based on osc freq and fs
			phaseInc = params.frequency_Hz / sampleRate;

		lfoParameters = params;
	}

	/** render a new audio output structure */
	virtual const SignalGenData renderAudioOutput();

protected:
	// --- parameters
	OscillatorParameters lfoParameters; ///< obejcgt parameters

	// --- sample rate
	double sampleRate = 0.0;			///< sample rate

	// --- timebase variables
	double modCounter = 0.0;			///< modulo counter [0.0, +1.0]
	double phaseInc = 0.0;				///< phase inc = fo/fs
	double modCounterQP = 0.25;			///<Quad Phase modulo counter [0.0, +1.0]

	/** check the modulo counter and wrap if needed */
	inline bool checkAndWrapModulo(double& moduloCounter, double phaseInc)
	{
		// --- for positive frequencies
		if (phaseInc > 0 && moduloCounter >= 1.0)
		{
			moduloCounter -= 1.0;
			return true;
		}

		// --- for negative frequencies
		if (phaseInc < 0 && moduloCounter <= 0.0)
		{
			moduloCounter += 1.0;
			return true;
		}

		return false;
	}

	/** advanvce the modulo counter, then check the modulo counter and wrap if needed */
	inline bool advanceAndCheckWrapModulo(double& moduloCounter, double phaseInc)
	{
		// --- advance counter
		moduloCounter += phaseInc;

		// --- for positive frequencies
		if (phaseInc > 0 && moduloCounter >= 1.0)
		{
			moduloCounter -= 1.0;
			return true;
		}

		// --- for negative frequencies
		if (phaseInc < 0 && moduloCounter <= 0.0)
		{
			moduloCounter += 1.0;
			return true;
		}

		return false;
	}

	/** advanvce the modulo counter */
	inline void advanceModulo(double& moduloCounter, double phaseInc) { moduloCounter += phaseInc; }

	const double B = 4.0 / kPi;
	const double C = -4.0 / (kPi* kPi);
	const double P = 0.225;
	/** parabolic sinusoidal calcualtion; NOTE: input is -pi to +pi http://devmaster.net/posts/9648/fast-and-accurate-sine-cosine */
	inline double parabolicSine(double angle)
	{
		double y = B * angle + C * angle * fabs(angle);
		y = P * (y * fabs(y) - y) + y;
		return y;
	}
};

/**
\enum DFOscillatorCoeffs
\ingroup Constants-Enums
\brief
Use this non-typed enum to easily access the direct form oscillator coefficients
- enum DFOscillatorCoeffs { df_b1, df_b2, numDFOCoeffs };
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
enum DFOscillatorCoeffs { df_b1, df_b2, numDFOCoeffs };

/**
\enum DFOscillatorStates
\ingroup Constants-Enums
\brief
Use this non-typed enum to easily access the direct form oscillator state registers
- DFOscillatorStates { df_yz1, df_yz2, numDFOStates };
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
enum DFOscillatorStates { df_yz1, df_yz2, numDFOStates };


/**
\class DFOscillator
\ingroup FX-Objects
\brief
The DFOscillator object implements generates a very pure sinusoidal oscillator by placing poles direclty on the unit circle.
Accuracy is excellent even at low frequencies.
Audio I/O:
- Output only object: pitched audio sinusoidal generator.
Control I/F:
- Use OscillatorParameters structure to get/set object params.
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
class DFOscillator : public IAudioSignalGenerator
{
public:
	DFOscillator() { }	/* C-TOR */
	virtual ~DFOscillator() {}				/* D-TOR */

	/** reset members to initialized state */
	virtual bool reset(double _sampleRate)
	{
		sampleRate = _sampleRate;
		memset(&stateArray[0], 0, sizeof(double)*numDFOStates);
		updateDFO();
		return true;
	}

	/** get parameters: note use of custom structure for passing param data */
	/**
	\return OscillatorParameters custom data structure
	*/
	OscillatorParameters getParameters()
	{
		return parameters;
	}

	/** set parameters: note use of custom structure for passing param data */
	/**
	\param OscillatorParameters custom data structure
	*/
	void setParameters(const OscillatorParameters& params)
	{
		if (parameters.frequency_Hz != params.frequency_Hz)
		{
			parameters = params;
			updateDFO();
		}
	}

	/** render the audio signal (pure sinusoid) */
	virtual const SignalGenData renderAudioOutput()
	{
		// --- calculates normal and inverted outputs; quadphase are not used
		SignalGenData output;

		// -- do difference equation y(n) = -b1y(n-2) - b2y(n-2)
		output.normalOutput = (-coeffArray[df_b1]*stateArray[df_yz1] - coeffArray[df_b2]*stateArray[df_yz2]);
		output.invertedOutput = -output.normalOutput;

		// --- update states
		stateArray[df_yz2] = stateArray[df_yz1];
		stateArray[df_yz1] = output.normalOutput;

		return output;
	}

	/** Update the object */
	void updateDFO()
	{
		// --- Oscillation Rate = theta = wT = w/fs
		double wT = (kTwoPi*parameters.frequency_Hz) / sampleRate;

		// --- coefficients to place poles right on unit circle
		coeffArray[df_b1] = -2.0*cos(wT);	// <--- set angle a = -2Rcod(theta)
		coeffArray[df_b2] = 1.0;			// <--- R^2 = 1, so R = 1

		// --- now update states to reflect the new frequency
		//	   re calculate the new initial conditions
		//	   arcsine of y(n-1) gives us wnT
		double wnT1 = asin(stateArray[df_yz1]);

		// find n by dividing wnT by wT
		double n = wnT1 / wT;

		// --- re calculate the new initial conditions
		//	   asin returns values from -pi/2 to +pi/2 where the sinusoid
		//     moves from -1 to +1 -- the leading (rising) edge of the
		//     sinewave. If we are on that leading edge (increasing)
		//     then we use the value 1T behind.
		//
		//     If we are on the falling edge, we use the value 1T ahead
		//     because it mimics the value that would be 1T behind
		if (stateArray[df_yz1] > stateArray[df_yz2])
			n -= 1;
		else
			n += 1;

		// ---  calculate the new (old) sample
		stateArray[df_yz2] = sin((n)*wT);
	}


protected:
	// --- parameters
	OscillatorParameters parameters; ///< object parameters

	// --- implementation of half a biquad - this object is extremely specific
	double stateArray[numDFOStates] = { 0.0, 0.0 };///< array of state registers
	double coeffArray[numDFOCoeffs] = { 0.0, 0.0 };///< array of coefficients

	// --- sample rate
	double sampleRate = 0.0;			///< sample rate
};


/**
\brief generates the oscillator output for one sample interval; note that there are multiple outputs.
*/
const SignalGenData LFO::renderAudioOutput()
{
	// --- always first!
	checkAndWrapModulo(modCounter, phaseInc);

	// --- QP output always follows location of current modulo; first set equal
	modCounterQP = modCounter;

	// --- then, advance modulo by quadPhaseInc = 0.25 = 90 degrees, AND wrap if needed
	advanceAndCheckWrapModulo(modCounterQP, 0.25);

	SignalGenData output;
	generatorWaveform waveform = lfoParameters.waveform;

	// --- calculate the oscillator value
	if (waveform == generatorWaveform::kSin)
	{
		// --- calculate normal angle
		double angle = modCounter*2.0*kPi - kPi;

		// --- norm output with parabolicSine approximation
		output.normalOutput = parabolicSine(-angle);

		// --- calculate QP angle
		angle = modCounterQP*2.0*kPi - kPi;

		// --- calc QP output
		output.quadPhaseOutput_pos = parabolicSine(-angle);
	}
	else if (waveform == generatorWaveform::kTriangle)
	{
		// triv saw
		output.normalOutput = unipolarToBipolar(modCounter);

		// bipolar triagle
		output.normalOutput = 2.0*fabs(output.normalOutput) - 1.0;

		// -- quad phase
		output.quadPhaseOutput_pos = unipolarToBipolar(modCounterQP);

		// bipolar triagle
		output.quadPhaseOutput_pos = 2.0*fabs(output.quadPhaseOutput_pos) - 1.0;
	}
	else if (waveform == generatorWaveform::kSaw)
	{
		output.normalOutput = unipolarToBipolar(modCounter);
		output.quadPhaseOutput_pos = unipolarToBipolar(modCounterQP);
	}

	// --- invert two main outputs to make the opposite versions
	output.quadPhaseOutput_neg = -output.quadPhaseOutput_pos;
	output.invertedOutput = -output.normalOutput;

	// --- setup for next sample period
	advanceModulo(modCounter, phaseInc);

	return output;
}

