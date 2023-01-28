#pragma once

#include "FXObjects.hpp"
#include "FXOscillators.hpp"
#include "FXDynamics.hpp"

// ------------------------------------------------------------------ //
// --- OBJECTS ------------------------------------------------------ //
// ------------------------------------------------------------------ //
/*
Class Declarations :
class name : public IAudioSignalProcessor
	- IAudioSignalProcessor functions
	- member functions that may be called externally
	- mutators & accessors
	- helper functions(may be private / protected if needed)
	- protected member functions
*/

/**
\enum filterCoeff
\ingroup Constants-Enums
\brief
Use this enum to easily access coefficents inside of arrays.
- enum filterCoeff { a0, a1, a2, b1, b2, c0, d0, numCoeffs };
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
enum filterCoeff { a0, a1, a2, b1, b2, c0, d0, numCoeffs };

/**
\enum stateReg
\ingroup Constants-Enums
\brief
Use this enum to easily access z^-1 state values inside of arrays. For some structures, fewer storage units are needed. They are divided as follows:
- Direct Forms: we will allow max of 2 for X (feedforward) and 2 for Y (feedback) data
- Transpose Forms: we will use ONLY the x_z1 and x_z2 registers for the 2 required delays
- enum stateReg { x_z1, x_z2, y_z1, y_z2, numStates };
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/

// --- state array index values
//     z^-1 registers;
//        Direct Forms: we will allow max of 2 for X (feedforward) and 2 for Y (feedback) data
//        Transpose Forms: we will use ONLY the x_z1 and x_z2 registers for the 2 required delays
enum stateReg { x_z1, x_z2, y_z1, y_z2, numStates };

/**
\enum biquadAlgorithm
\ingroup Constants-Enums
\brief
Use this strongly typed enum to easily set the biquad calculation type
- enum class biquadAlgorithm { kDirect, kCanonical, kTransposeDirect, kTransposeCanonical }; //  4 types of biquad calculations, constants (k)
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/

// --- type of calculation (algorithm)
enum class biquadAlgorithm { kDirect, kCanonical, kTransposeDirect, kTransposeCanonical }; //  4 types of biquad calculations, constants (k)


/**
\struct BiquadParameters
\ingroup FX-Objects
\brief
Custom parameter structure for the Biquad object. Default version defines the biquad structure used in the calculation.
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
struct BiquadParameters
{
	BiquadParameters () {}

	/** all FXObjects parameter objects require overloaded= operator so remember to add new entries if you add new variables. */
	BiquadParameters& operator=(const BiquadParameters& params)
	{
		if (this == &params)
			return *this;

		biquadCalcType = params.biquadCalcType;
		return *this;
	}

	biquadAlgorithm biquadCalcType = biquadAlgorithm::kDirect; ///< biquad structure to use
};

/**
\class Biquad
\ingroup FX-Objects
\brief
The Biquad object implements a first or second order H(z) transfer function using one of four standard structures: Direct, Canonical, Transpose Direct, Transpose Canonical.
Audio I/O:
- Processes mono input to mono output.
Control I/F:
- Use BiquadParameters structure to get/set object params.
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
class Biquad : public IAudioSignalProcessor
{
public:
	Biquad() {}		/* C-TOR */
	~Biquad() {}	/* D-TOR */

	// --- IAudioSignalProcessor FUNCTIONS --- //
	//
	/** reset: clear out the state array (flush delays); can safely ignore sampleRate argument - we don't need/use it */
	virtual bool reset(DspFloatType _sampleRate)
	{
		memset(&stateArray[0], 0, sizeof(DspFloatType)*numStates);
		return true;  // handled = true
	}

	/** return false: this object only processes samples */
	virtual bool canProcessAudioFrame() { return false; }

	/** process input x(n) through biquad to produce return value y(n) */
	/**
	\param xn input
	\return the processed sample
	*/
	virtual DspFloatType processAudioSample(DspFloatType xn);

	/** get parameters: note use of custom structure for passing param data */
	/**
	\return BiquadParameters custom data structure
	*/
	BiquadParameters getParameters() { return parameters ; }

	/** set parameters: note use of custom structure for passing param data */
	/**
	\param BiquadParameters custom data structure
	*/
	void setParameters(const BiquadParameters& _parameters){ parameters = _parameters; }

	// --- MUTATORS & ACCESSORS --- //
	/** set the coefficient array NOTE: passing by pointer to array; allows us to use "array notation" with pointers i.e. [ ] */
	void setCoefficients(DspFloatType* coeffs){
		// --- fast block memory copy:
		memcpy(&coeffArray[0], &coeffs[0], sizeof(DspFloatType)*numCoeffs);
	}

	/** get the coefficient array for read/write access to the array (not used in current objects) */
	DspFloatType* getCoefficients()
	{
		// --- read/write access to the array (not used)
		return &coeffArray[0];
	}

	/** get the state array for read/write access to the array (used only in direct form oscillator) */
	DspFloatType* getStateArray()
	{
		// --- read/write access to the array (used only in direct form oscillator)
		return &stateArray[0];
	}

	/** get the structure G (gain) value for Harma filters; see 2nd Ed FX book */
	DspFloatType getG_value() { return coeffArray[a0]; }

	/** get the structure S (storage) value for Harma filters; see 2nd Ed FX book */
	DspFloatType getS_value();// { return storageComponent; }

protected:
	/** array of coefficients */
	DspFloatType coeffArray[numCoeffs] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

	/** array of state (z^-1) registers */
	DspFloatType stateArray[numStates] = { 0.0, 0.0, 0.0, 0.0 };

	/** type of calculation (algorithm  structure) */
	BiquadParameters parameters;

	/** for Harma loop resolution */
	DspFloatType storageComponent = 0.0;
};


/**
\enum filterAlgorithm
\ingroup Constants-Enums
\brief
Use this strongly typed enum to easily set the filter algorithm for the AudioFilter object or any other object that requires filter definitions.
- filterAlgorithm { kLPF1P, kLPF1, kHPF1, kLPF2, kHPF2, kBPF2, kBSF2, kButterLPF2, kButterHPF2, kButterBPF2, kButterBSF2, kMMALPF2, kMMALPF2B, kLowShelf, kHiShelf, kNCQParaEQ, kCQParaEQ, kLWRLPF2, kLWRHPF2, kAPF1, kAPF2, kResonA, kResonB, kMatchLP2A, kMatchLP2B, kMatchBP2A, kMatchBP2B, kImpInvLP1, kImpInvLP2 };
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
enum class filterAlgorithm {
	kLPF1P, kLPF1, kHPF1, kLPF2, kHPF2, kBPF2, kBSF2, kButterLPF2, kButterHPF2, kButterBPF2,
	kButterBSF2, kMMALPF2, kMMALPF2B, kLowShelf, kHiShelf, kNCQParaEQ, kCQParaEQ, kLWRLPF2, kLWRHPF2,
	kAPF1, kAPF2, kResonA, kResonB, kMatchLP2A, kMatchLP2B, kMatchBP2A, kMatchBP2B,
	kImpInvLP1, kImpInvLP2
}; // --- you will add more here...


/**
\struct AudioFilterParameters
\ingroup FX-Objects
\brief
Custom parameter structure for the AudioFilter object.
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
struct AudioFilterParameters
{
	AudioFilterParameters(){}
	/** all FXObjects parameter objects require overloaded= operator so remember to add new entries if you add new variables. */
	AudioFilterParameters& operator=(const AudioFilterParameters& params)	// need this override for collections to work
	{
		if (this == &params)
			return *this;
		algorithm = params.algorithm;
		fc = params.fc;
		Q = params.Q;
		boostCut_dB = params.boostCut_dB;

		return *this;
	}

	// --- individual parameters
	filterAlgorithm algorithm = filterAlgorithm::kLPF1; ///< filter algorithm
	DspFloatType fc = 100.0; ///< filter cutoff or center frequency (Hz)
	DspFloatType Q = 0.707; ///< filter Q
	DspFloatType boostCut_dB = 0.0; ///< filter gain; note not used in all types
};

/**
\class AudioFilter
\ingroup FX-Objects
\brief
The AudioFilter object implements all filters in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
Audio I/O:
- Processes mono input to mono output.
Control I/F:
- Use AudioFilterParameters structure to get/set object params.
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
class AudioFilter : public IAudioSignalProcessor
{
public:
	AudioFilter() {}		/* C-TOR */
	~AudioFilter() {}		/* D-TOR */

	// --- IAudioSignalProcessor
	/** --- set sample rate, then update coeffs */
	virtual bool reset(DspFloatType _sampleRate)
	{
		BiquadParameters bqp = biquad.getParameters();

		// --- you can try both forms - do you hear a difference?
		bqp.biquadCalcType = biquadAlgorithm::kTransposeCanonical; //<- this is the default operation
	//	bqp.biquadCalcType = biquadAlgorithm::kDirect;
		biquad.setParameters(bqp);

		sampleRate = _sampleRate;
		return biquad.reset(_sampleRate);
	}

	/** return false: this object only processes samples */
	virtual bool canProcessAudioFrame() { return false; }

	/** process input x(n) through the filter to produce return value y(n) */
	/**
	\param xn input
	\return the processed sample
	*/
	virtual DspFloatType processAudioSample(DspFloatType xn);

	/** --- sample rate change necessarily requires recalculation */
	virtual void setSampleRate(DspFloatType _sampleRate)
	{
		sampleRate = _sampleRate;
		calculateFilterCoeffs();
	}

	/** --- get parameters */
	AudioFilterParameters getParameters() { return audioFilterParameters; }

	/** --- set parameters */
	void setParameters(const AudioFilterParameters& parameters)
	{
		if (audioFilterParameters.algorithm != parameters.algorithm ||
			audioFilterParameters.boostCut_dB != parameters.boostCut_dB ||
			audioFilterParameters.fc != parameters.fc ||
			audioFilterParameters.Q != parameters.Q)
		{
			// --- save new params
			audioFilterParameters = parameters;
		}
		else
			return;

		// --- don't allow 0 or (-) values for Q
		if (audioFilterParameters.Q <= 0)
			audioFilterParameters.Q = 0.707;

		// --- update coeffs
		calculateFilterCoeffs();
	}

	/** --- helper for Harma filters (phaser) */
	DspFloatType getG_value() { return biquad.getG_value(); }

	/** --- helper for Harma filters (phaser) */
	DspFloatType getS_value() { return biquad.getS_value(); }

protected:
	// --- our calculator
	Biquad biquad; ///< the biquad object

	// --- array to hold coeffs (we need them too)
	DspFloatType coeffArray[numCoeffs] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }; ///< our local copy of biquad coeffs

	// --- object parameters
	AudioFilterParameters audioFilterParameters; ///< parameters
	DspFloatType sampleRate = 44100.0; ///< current sample rate

	/** --- function to recalculate coefficients due to a change in filter parameters */
	bool calculateFilterCoeffs();
};


/**
\struct FilterBankOutput
\ingroup FX-Objects
\brief
Custom output structure for filter bank objects that split the inptu into multiple frequency channels (bands)
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
// --- output for filter bank requires multiple channels (bands)
struct FilterBankOutput
{
	FilterBankOutput() {}

	// --- band-split filter output
	DspFloatType LFOut = 0.0; ///< low frequency output sample
	DspFloatType HFOut = 0.0;	///< high frequency output sample

	// --- add more filter channels here; or use an array[]
};


/**
\struct LRFilterBankParameters
\ingroup FX-Objects
\brief
Custom parameter structure for the LRFilterBank object which splits the input signal into multiple bands.
The stock obejct splits into low and high frequency bands so this structure only requires one split point - add more
split frequencies to support more bands.
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
struct LRFilterBankParameters
{
	LRFilterBankParameters() {}
	/** all FXObjects parameter objects require overloaded= operator so remember to add new entries if you add new variables. */
	LRFilterBankParameters& operator=(const LRFilterBankParameters& params)	// need this override for collections to work
	{
		if (this == &params)
			return *this;
		splitFrequency = params.splitFrequency;
		return *this;
	}

	// --- individual parameters
	DspFloatType splitFrequency = 1000.0; ///< LF/HF split frequency
};


/**
\class LRFilterBank
\ingroup FX-Objects
\brief
The LRFilterBank object implements 2 Linkwitz-Riley Filters in a parallel filter bank to split the signal into two frequency bands.
Note that one channel is inverted (see the FX book below for explanation). You can add more bands here as well.
Audio I/O:
- Processes mono input into a custom FilterBankOutput structure.
NOTE: processAudioSample( ) is inoperable and only returns the input back.
Control I/F:
- Use LRFilterBankParameters structure to get/set object params.
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
class LRFilterBank : public IAudioSignalProcessor
{
public:
	LRFilterBank()		/* C-TOR */
	{
		// --- set filters as Linkwitz-Riley 2nd order
		AudioFilterParameters params = lpFilter.getParameters();
		params.algorithm = filterAlgorithm::kLWRLPF2;
		lpFilter.setParameters(params);

		params = hpFilter.getParameters();
		params.algorithm = filterAlgorithm::kLWRHPF2;
		hpFilter.setParameters(params);
	}

	~LRFilterBank() {}	/* D-TOR */

	/** reset member objects */
	virtual bool reset(DspFloatType _sampleRate)
	{
		lpFilter.reset(_sampleRate);
		hpFilter.reset(_sampleRate);
		return true;
	}

	/** return false: this object only processes samples */
	virtual bool canProcessAudioFrame() { return false; }

	/** this does nothing for this object, see processFilterBank( ) below */
	/**
	\param xn input
	\return the processed sample
	*/
	virtual DspFloatType processAudioSample(DspFloatType xn)
	{
		return xn;
	}

	/** process the filter bank */
	FilterBankOutput processFilterBank(DspFloatType xn)
	{
		FilterBankOutput output;

		// --- process the LPF
		output.LFOut = lpFilter.processAudioSample(xn);

		// --- invert the HP filter output so that recombination will
		//     result in the correct phase and magnitude responses
		output.HFOut = -hpFilter.processAudioSample(xn);

		return output;
	}

	/** get parameters: note use of custom structure for passing param data */
	/**
	\return LRFilterBankParameters custom data structure
	*/
	LRFilterBankParameters getParameters()
	{
		return parameters;
	}

	/** set parameters: note use of custom structure for passing param data */
	/**
	\param LRFilterBankParameters custom data structure
	*/
	void setParameters(const LRFilterBankParameters& _parameters)
	{
		// --- update structure
		parameters = _parameters;

		// --- update member objects
		AudioFilterParameters params = lpFilter.getParameters();
		params.fc = parameters.splitFrequency;
		lpFilter.setParameters(params);

		params = hpFilter.getParameters();
		params.fc = parameters.splitFrequency;
		hpFilter.setParameters(params);
	}

protected:
	AudioFilter lpFilter; ///< low-band filter
	AudioFilter hpFilter; ///< high-band filter

	// --- object parameters
	LRFilterBankParameters parameters; ///< parameters for the object
};


/**
\struct PhaseShifterParameters
\ingroup FX-Objects
\brief
Custom parameter structure for the PhaseShifter object.
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
struct PhaseShifterParameters
{
	PhaseShifterParameters() {}
	/** all FXObjects parameter objects require overloaded= operator so remember to add new entries if you add new variables. */
	PhaseShifterParameters& operator=(const PhaseShifterParameters& params)
	{
		if (this == &params)
			return *this;

		lfoRate_Hz = params.lfoRate_Hz;
		lfoDepth_Pct = params.lfoDepth_Pct;
		intensity_Pct = params.intensity_Pct;
		quadPhaseLFO = params.quadPhaseLFO;
		return *this;
	}

	// --- individual parameters
	DspFloatType lfoRate_Hz = 0.0;	///< phaser LFO rate in Hz
	DspFloatType lfoDepth_Pct = 0.0;	///< phaser LFO depth in %
	DspFloatType intensity_Pct = 0.0;	///< phaser feedback in %
	bool quadPhaseLFO = false;	///< quad phase LFO flag
};

// --- constants for Phaser
const unsigned int PHASER_STAGES = 6;

// --- these are the ideal band definitions
//const DspFloatType apf0_minF = 16.0;
//const DspFloatType apf0_maxF = 1600.0;
//
//const DspFloatType apf1_minF = 33.0;
//const DspFloatType apf1_maxF = 3300.0;
//
//const DspFloatType apf2_minF = 48.0;
//const DspFloatType apf2_maxF = 4800.0;
//
//const DspFloatType apf3_minF = 98.0;
//const DspFloatType apf3_maxF = 9800.0;
//
//const DspFloatType apf4_minF = 160.0;
//const DspFloatType apf4_maxF = 16000.0;
//
//const DspFloatType apf5_minF = 260.0;
//const DspFloatType apf5_maxF = 20480.0;

// --- these are the exact values from the National Semiconductor Phaser design
const DspFloatType apf0_minF = 32.0;
const DspFloatType apf0_maxF = 1500.0;

const DspFloatType apf1_minF = 68.0;
const DspFloatType apf1_maxF = 3400.0;

const DspFloatType apf2_minF = 96.0;
const DspFloatType apf2_maxF = 4800.0;

const DspFloatType apf3_minF = 212.0;
const DspFloatType apf3_maxF = 10000.0;

const DspFloatType apf4_minF = 320.0;
const DspFloatType apf4_maxF = 16000.0;

const DspFloatType apf5_minF = 636.0;
const DspFloatType apf5_maxF = 20480.0;

/**
\class PhaseShifter
\ingroup FX-Objects
\brief
The PhaseShifter object implements a six-stage phaser.
Audio I/O:
- Processes mono input to mono output.
Control I/F:
- Use BiquadParameters structure to get/set object params.
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
class PhaseShifter : public IAudioSignalProcessor
{
public:
	PhaseShifter(void) {
		OscillatorParameters lfoparams = lfo.getParameters();
		lfoparams.waveform = generatorWaveform::kTriangle;	// kTriangle LFO for phaser
	//	lfoparams.waveform = generatorWaveform::kSin;		// kTriangle LFO for phaser
		lfo.setParameters(lfoparams);

		AudioFilterParameters params = apf[0].getParameters();
		params.algorithm = filterAlgorithm::kAPF1; // can also use 2nd order
		// params.Q = 0.001; use low Q if using 2nd order APFs

		for (int i = 0; i < PHASER_STAGES; i++)
		{
			apf[i].setParameters(params);
		}
	}	/* C-TOR */

	~PhaseShifter(void) {}	/* D-TOR */

public:
	/** reset members to initialized state */
	virtual bool reset(DspFloatType _sampleRate)
	{
		// --- reset LFO
		lfo.reset(_sampleRate);

		// --- reset APFs
		for (int i = 0; i < PHASER_STAGES; i++){
			apf[i].reset(_sampleRate);
		}

		return true;
	}

	/** process autio through phaser */
	/**
	\param xn input
	\return the processed sample
	*/
	virtual DspFloatType processAudioSample(DspFloatType xn)
	{
		SignalGenData lfoData = lfo.renderAudioOutput();

		// --- create the bipolar modulator value
		DspFloatType lfoValue = lfoData.normalOutput;
		if (parameters.quadPhaseLFO)
			lfoValue = lfoData.quadPhaseOutput_pos;

		DspFloatType depth = parameters.lfoDepth_Pct / 100.0;
		DspFloatType modulatorValue = lfoValue*depth;

		// --- calculate modulated values for each APF; note they have different ranges
		AudioFilterParameters params = apf[0].getParameters();
		params.fc = doBipolarModulation(modulatorValue, apf0_minF, apf0_maxF);
		apf[0].setParameters(params);

		params = apf[1].getParameters();
		params.fc = doBipolarModulation(modulatorValue, apf1_minF, apf1_maxF);
		apf[1].setParameters(params);

		params = apf[2].getParameters();
		params.fc = doBipolarModulation(modulatorValue, apf2_minF, apf2_maxF);
		apf[2].setParameters(params);

		params = apf[3].getParameters();
		params.fc = doBipolarModulation(modulatorValue, apf3_minF, apf3_maxF);
		apf[3].setParameters(params);

		params = apf[4].getParameters();
		params.fc = doBipolarModulation(modulatorValue, apf4_minF, apf4_maxF);
		apf[4].setParameters(params);

		params = apf[5].getParameters();
		params.fc = doBipolarModulation(modulatorValue, apf5_minF, apf5_maxF);
		apf[5].setParameters(params);

		// --- calculate gamma values
		DspFloatType gamma1 = apf[5].getG_value();
		DspFloatType gamma2 = apf[4].getG_value() * gamma1;
		DspFloatType gamma3 = apf[3].getG_value() * gamma2;
		DspFloatType gamma4 = apf[2].getG_value() * gamma3;
		DspFloatType gamma5 = apf[1].getG_value() * gamma4;
		DspFloatType gamma6 = apf[0].getG_value() * gamma5;

		// --- set the alpha0 value
		DspFloatType K = parameters.intensity_Pct / 100.0;
		DspFloatType alpha0 = 1.0 / (1.0 + K*gamma6);

		// --- create combined feedback
		DspFloatType Sn = gamma5*apf[0].getS_value() + gamma4*apf[1].getS_value() + gamma3*apf[2].getS_value() + gamma2*apf[3].getS_value() + gamma1*apf[4].getS_value() + apf[5].getS_value();

		// --- form input to first APF
		DspFloatType u = alpha0*(xn + K*Sn);

		// --- cascade of APFs (could also nest these in one massive line of code)
		DspFloatType APF1 = apf[0].processAudioSample(u);
		DspFloatType APF2 = apf[1].processAudioSample(APF1);
		DspFloatType APF3 = apf[2].processAudioSample(APF2);
		DspFloatType APF4 = apf[3].processAudioSample(APF3);
		DspFloatType APF5 = apf[4].processAudioSample(APF4);
		DspFloatType APF6 = apf[5].processAudioSample(APF5);

		// --- sum with -3dB coefficients
		//	DspFloatType output = 0.707*xn + 0.707*APF6;

		// --- sum with National Semiconductor design ratio:
		//	   dry = 0.5, wet = 5.0
		// DspFloatType output = 0.5*xn + 5.0*APF6;
		// DspFloatType output = 0.25*xn + 2.5*APF6;
		DspFloatType output = 0.125*xn + 1.25*APF6;

		return output;
	}

	/** return false: this object only processes samples */
	virtual bool canProcessAudioFrame() { return false; }

	/** get parameters: note use of custom structure for passing param data */
	/**
	\return PhaseShifterParameters custom data structure
	*/
	PhaseShifterParameters getParameters() { return parameters; }

	/** set parameters: note use of custom structure for passing param data */
	/**
	\param PhaseShifterParameters custom data structure
	*/
	void setParameters(const PhaseShifterParameters& params)
	{
		// --- update LFO rate
		if (params.lfoRate_Hz != parameters.lfoRate_Hz)
		{
			OscillatorParameters lfoparams = lfo.getParameters();
			lfoparams.frequency_Hz = params.lfoRate_Hz;
			lfo.setParameters(lfoparams);
		}

		// --- save new
		parameters = params;
	}
protected:
	PhaseShifterParameters parameters;  ///< the object parameters
	AudioFilter apf[PHASER_STAGES];		///< six APF objects
	LFO lfo;							///< the one and only LFO
};

/**
\struct SimpleLPFParameters
\ingroup FX-Objects
\brief
Custom parameter structure for the SimpleLPFP object. Used for reverb algorithms in book.
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
struct SimpleLPFParameters
{
	SimpleLPFParameters() {}
	/** all FXObjects parameter objects require overloaded= operator so remember to add new entries if you add new variables. */
	SimpleLPFParameters& operator=(const SimpleLPFParameters& params)	// need this override for collections to work
	{
		if (this == &params)
			return *this;

		g = params.g;
		return *this;
	}

	// --- individual parameters
	DspFloatType g = 0.0; ///< simple LPF g value
};

/**
\class SimpleLPF
\ingroup FX-Objects
\brief
The SimpleLPF object implements a first order one-pole LPF using one coefficient "g" value.
Audio I/O:
- Processes mono input to mono output.
Control I/F:
- Use SimpleLPFParameters structure to get/set object params.
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
class SimpleLPF : public IAudioSignalProcessor
{
public:
	SimpleLPF(void) {}	/* C-TOR */
	~SimpleLPF(void) {}	/* D-TOR */

public:
	/** reset members to initialized state */
	virtual bool reset(DspFloatType _sampleRate)
	{
		state = 0.0;
		return true;
	}

	/** get parameters: note use of custom structure for passing param data */
	/**
	\return SimpleLPFParameters custom data structure
	*/
	SimpleLPFParameters getParameters()
	{
		return simpleLPFParameters;
	}

	/** set parameters: note use of custom structure for passing param data */
	/**
	\param SimpleLPFParameters custom data structure
	*/
	void setParameters(const SimpleLPFParameters& params)
	{
		simpleLPFParameters = params;
	}

	/** process simple one pole FB back filter */
	/**
	\param xn input
	\return the processed sample
	*/
	virtual DspFloatType processAudioSample(DspFloatType xn)
	{
		DspFloatType g = simpleLPFParameters.g;
		DspFloatType yn = (1.0 - g)*xn + g*state;
		state = yn;
		return yn;
	}

	/** return false: this object only processes samples */
	virtual bool canProcessAudioFrame() { return false; }

private:
	SimpleLPFParameters simpleLPFParameters;	///< object parameters
	DspFloatType state = 0.0;							///< single state (z^-1) register
};


/**
\enum vaFilterAlgorithm
\ingroup Constants-Enums
\brief
Use this strongly typed enum to easily set the virtual analog filter algorithm
- enum class vaFilterAlgorithm { kLPF1, kHPF1, kAPF1, kSVF_LP, kSVF_HP, kSVF_BP, kSVF_BS };
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
enum class vaFilterAlgorithm {
	kLPF1, kHPF1, kAPF1, kSVF_LP, kSVF_HP, kSVF_BP, kSVF_BS
}; // --- you will add more here...


/**
\struct ZVAFilterParameters
\ingroup FX-Objects
\brief
Custom parameter structure for the ZVAFilter object.
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
struct ZVAFilterParameters
{
	ZVAFilterParameters() {}
	/** all FXObjects parameter objects require overloaded= operator so remember to add new entries if you add new variables. */
	ZVAFilterParameters& operator=(const ZVAFilterParameters& params)	// need this override for collections to work
	{
		if (this == &params)
			return *this;

		filterAlgorithm = params.filterAlgorithm;
		fc = params.fc;
		Q = params.Q;
		filterOutputGain_dB = params.filterOutputGain_dB;
		enableGainComp = params.enableGainComp;
		matchAnalogNyquistLPF = params.matchAnalogNyquistLPF;
		selfOscillate = params.selfOscillate;
		enableNLP = params.enableNLP;
		return *this;
	}

	// --- individual parameters
	vaFilterAlgorithm filterAlgorithm = vaFilterAlgorithm::kSVF_LP;	///< va filter algorithm
	DspFloatType fc = 1000.0;						///< va filter fc
	DspFloatType Q = 0.707;						///< va filter Q
	DspFloatType filterOutputGain_dB = 0.0;		///< va filter gain (normally unused)
	bool enableGainComp = false;			///< enable gain compensation (see book)
	bool matchAnalogNyquistLPF = false;		///< match analog gain at Nyquist
	bool selfOscillate = false;				///< enable selfOscillation
	bool enableNLP = false;					///< enable non linear processing (use oversampling for best results)
};


/**
\class ZVAFilter
\ingroup FX-Objects
\brief
The ZVAFilter object implements multpile Zavalishin VA Filters.
Audio I/O:
- Processes mono input to mono output.
Control I/F:
- Use BiquadParameters structure to get/set object params.
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
class ZVAFilter : public IAudioSignalProcessor
{
public:
	ZVAFilter() {}		/* C-TOR */
	~ZVAFilter() {}		/* D-TOR */

	/** reset members to initialized state */
	virtual bool reset(DspFloatType _sampleRate)
	{
		sampleRate = _sampleRate;
		integrator_z[0] = 0.0;
		integrator_z[1] = 0.0;

		return true;
	}

	/** get parameters: note use of custom structure for passing param data */
	/**
	\return ZVAFilterParameters custom data structure
	*/
	ZVAFilterParameters getParameters()
	{
		return zvaFilterParameters;
	}

	/** set parameters: note use of custom structure for passing param data */
	/**
	\param ZVAFilterParameters custom data structure
	*/
	void setParameters(const ZVAFilterParameters& params)
	{
		if (params.fc != zvaFilterParameters.fc ||
			params.Q != zvaFilterParameters.Q ||
			params.selfOscillate != zvaFilterParameters.selfOscillate ||
			params.matchAnalogNyquistLPF != zvaFilterParameters.matchAnalogNyquistLPF)
		{
				zvaFilterParameters = params;
				calculateFilterCoeffs();
		}
		else
			zvaFilterParameters = params;
	}

	/** return false: this object only processes samples */
	virtual bool canProcessAudioFrame() { return false; }

	/** process input x(n) through the VA filter to produce return value y(n) */
	/**
	\param xn input
	\return the processed sample
	*/
	virtual DspFloatType processAudioSample(DspFloatType xn)
	{
		// --- with gain comp enabled, we reduce the input by
		//     half the gain in dB at resonant peak
		//     NOTE: you can change that logic here!
		vaFilterAlgorithm filterAlgorithm = zvaFilterParameters.filterAlgorithm;
		bool matchAnalogNyquistLPF = zvaFilterParameters.matchAnalogNyquistLPF;

		if (zvaFilterParameters.enableGainComp)
		{
			DspFloatType peak_dB = dBPeakGainFor_Q(zvaFilterParameters.Q);
			if (peak_dB > 0.0)
			{
				DspFloatType halfPeak_dBGain = dB2Raw(-peak_dB / 2.0);
				xn *= halfPeak_dBGain;
			}
		}

		// --- for 1st order filters:
		if (filterAlgorithm == vaFilterAlgorithm::kLPF1 ||
			filterAlgorithm == vaFilterAlgorithm::kHPF1 ||
			filterAlgorithm == vaFilterAlgorithm::kAPF1)
		{
			// --- create vn node
			DspFloatType vn = (xn - integrator_z[0])*alpha;

			// --- form LP output
			DspFloatType lpf = ((xn - integrator_z[0])*alpha) + integrator_z[0];

			// DspFloatType sn = integrator_z[0];

			// --- update memory
			integrator_z[0] = vn + lpf;

			// --- form the HPF = INPUT = LPF
			DspFloatType hpf = xn - lpf;

			// --- form the APF = LPF - HPF
			DspFloatType apf = lpf - hpf;

			// --- set the outputs
			if (filterAlgorithm == vaFilterAlgorithm::kLPF1)
			{
				// --- this is a very close match as-is at Nyquist!
				if (matchAnalogNyquistLPF)
					return lpf + alpha*hpf;
				else
					return lpf;
			}
			else if (filterAlgorithm == vaFilterAlgorithm::kHPF1)
				return hpf;
			else if (filterAlgorithm == vaFilterAlgorithm::kAPF1)
				return apf;

			// --- unknown filter
			return xn;
		}

		// --- form the HP output first
		DspFloatType hpf = alpha0*(xn - rho*integrator_z[0] - integrator_z[1]);

		// --- BPF Out
		DspFloatType bpf = alpha*hpf + integrator_z[0];
		if (zvaFilterParameters.enableNLP)
			bpf = softClipWaveShaper(bpf, 1.0);

		// --- LPF Out
		DspFloatType lpf = alpha*bpf + integrator_z[1];

		// --- BSF Out
		DspFloatType bsf = hpf + lpf;

		// --- finite gain at Nyquist; slight error at VHF
		DspFloatType sn = integrator_z[0];

		// update memory
		integrator_z[0] = alpha*hpf + bpf;
		integrator_z[1] = alpha*bpf + lpf;

		DspFloatType filterOutputGain = pow(10.0, zvaFilterParameters.filterOutputGain_dB / 20.0);

		// return our selected type
		if (filterAlgorithm == vaFilterAlgorithm::kSVF_LP)
		{
			if (matchAnalogNyquistLPF)
				lpf += analogMatchSigma*(sn);
			return filterOutputGain*lpf;
		}
		else if (filterAlgorithm == vaFilterAlgorithm::kSVF_HP)
			return filterOutputGain*hpf;
		else if (filterAlgorithm == vaFilterAlgorithm::kSVF_BP)
			return filterOutputGain*bpf;
		else if (filterAlgorithm == vaFilterAlgorithm::kSVF_BS)
			return filterOutputGain*bsf;

		// --- unknown filter
		return filterOutputGain*lpf;
	}

	/** recalculate the filter coefficients*/
	void calculateFilterCoeffs()
	{
		DspFloatType fc = zvaFilterParameters.fc;
		DspFloatType Q = zvaFilterParameters.Q;
		vaFilterAlgorithm filterAlgorithm = zvaFilterParameters.filterAlgorithm;

		// --- normal Zavalishin SVF calculations here
		//     prewarp the cutoff- these are bilinear-transform filters
		DspFloatType wd = kTwoPi*fc;
		DspFloatType T = 1.0 / sampleRate;
		DspFloatType wa = (2.0 / T)*tan(wd*T / 2.0);
		DspFloatType g = wa*T / 2.0;

		// --- for 1st order filters:
		if (filterAlgorithm == vaFilterAlgorithm::kLPF1 ||
			filterAlgorithm == vaFilterAlgorithm::kHPF1 ||
			filterAlgorithm == vaFilterAlgorithm::kAPF1)
		{
			// --- calculate alpha
			alpha = g / (1.0 + g);
		}
		else // state variable variety
		{
			// --- note R is the traditional analog damping factor zeta
			DspFloatType R = zvaFilterParameters.selfOscillate ? 0.0 : 1.0 / (2.0*Q);
			alpha0 = 1.0 / (1.0 + 2.0*R*g + g*g);
			alpha = g;
			rho = 2.0*R + g;

			// --- sigma for analog matching version
			DspFloatType f_o = (sampleRate / 2.0) / fc;
			analogMatchSigma = 1.0 / (alpha*f_o*f_o);
		}
	}

	/** set beta value, for filters that aggregate 1st order VA sections*/
	void setBeta(DspFloatType _beta) { beta = _beta; }

	/** get beta value,not used in book projects; for future use*/
	DspFloatType getBeta() { return beta; }

protected:
	ZVAFilterParameters zvaFilterParameters;	///< object parameters
	DspFloatType sampleRate = 44100.0;				///< current sample rate

	// --- state storage
	DspFloatType integrator_z[2];						///< state variables

	// --- filter coefficients
	DspFloatType alpha0 = 0.0;		///< input scalar, correct delay-free loop
	DspFloatType alpha = 0.0;			///< alpha is (wcT/2)
	DspFloatType rho = 0.0;			///< p = 2R + g (feedback)

	DspFloatType beta = 0.0;			///< beta value, not used

	// --- for analog Nyquist matching
	DspFloatType analogMatchSigma = 0.0; ///< analog matching Sigma value (see book)

};

/**
\struct EnvelopeFollowerParameters
\ingroup FX-Objects
\brief
Custom parameter structure for the EnvelopeFollower object.
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
struct EnvelopeFollowerParameters
{
	EnvelopeFollowerParameters() {}
	/** all FXObjects parameter objects require overloaded= operator so remember to add new entries if you add new variables. */
	EnvelopeFollowerParameters& operator=(const EnvelopeFollowerParameters& params)	// need this override for collections to work
	{
		if (this == &params)
			return *this;

		fc = params.fc;
		Q = params.Q;
		attackTime_mSec = params.attackTime_mSec;
		releaseTime_mSec = params.releaseTime_mSec;
		threshold_dB = params.threshold_dB;
		sensitivity = params.sensitivity;

		return *this;
	}

	// --- individual parameters
	DspFloatType fc = 0.0;				///< filter fc
	DspFloatType Q = 0.707;				///< filter Q
	DspFloatType attackTime_mSec = 10.0;	///< detector attack time
	DspFloatType releaseTime_mSec = 10.0;	///< detector release time
	DspFloatType threshold_dB = 0.0;		///< detector threshold in dB
	DspFloatType sensitivity = 1.0;		///< detector sensitivity
};

/**
\class EnvelopeFollower
\ingroup FX-Objects
\brief
The EnvelopeFollower object implements a traditional envelope follower effect modulating a LPR fc value
using the strength of the detected input.
Audio I/O:
- Processes mono input to mono output.
Control I/F:
- Use EnvelopeFollowerParameters structure to get/set object params.
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
class EnvelopeFollower : public IAudioSignalProcessor
{
public:
	EnvelopeFollower() {
		// --- setup the filter
		ZVAFilterParameters filterParams;
		filterParams.filterAlgorithm = vaFilterAlgorithm::kSVF_LP;
		filterParams.fc = 1000.0;
		filterParams.enableGainComp = true;
		filterParams.enableNLP = true;
		filterParams.matchAnalogNyquistLPF = true;
		filter.setParameters(filterParams);

		// --- setup the detector
		AudioDetectorParameters adParams;
		adParams.attackTime_mSec = -1.0;
		adParams.releaseTime_mSec = -1.0;
		adParams.detectMode = TLD_AUDIO_DETECT_MODE_RMS;
		adParams.detect_dB = true;
		adParams.clampToUnityMax = false;
		detector.setParameters(adParams);

	}		/* C-TOR */
	~EnvelopeFollower() {}		/* D-TOR */

	/** reset members to initialized state */
	virtual bool reset(DspFloatType _sampleRate)
	{
		filter.reset(_sampleRate);
		detector.reset(_sampleRate);
		return true;
	}

	/** get parameters: note use of custom structure for passing param data */
	/**
	\return EnvelopeFollowerParameters custom data structure
	*/
	EnvelopeFollowerParameters getParameters() { return parameters; }

	/** set parameters: note use of custom structure for passing param data */
	/**
	\param EnvelopeFollowerParameters custom data structure
	*/
	void setParameters(const EnvelopeFollowerParameters& params)
	{
		ZVAFilterParameters filterParams = filter.getParameters();
		AudioDetectorParameters adParams = detector.getParameters();

		if (params.fc != parameters.fc || params.Q != parameters.Q)
		{
			filterParams.fc = params.fc;
			filterParams.Q = params.Q;
			filter.setParameters(filterParams);
		}
		if (params.attackTime_mSec != parameters.attackTime_mSec ||
			params.releaseTime_mSec != parameters.releaseTime_mSec)
		{
			adParams.attackTime_mSec = params.attackTime_mSec;
			adParams.releaseTime_mSec = params.releaseTime_mSec;
			detector.setParameters(adParams);
		}

		// --- save
		parameters = params;
	}

	/** return false: this object only processes samples */
	virtual bool canProcessAudioFrame() { return false; }

	/** process input x(n) through the envelope follower to produce return value y(n) */
	/**
	\param xn input
	\return the processed sample
	*/
	virtual DspFloatType processAudioSample(DspFloatType xn)
	{
		// --- calc threshold
		DspFloatType threshValue = pow(10.0, parameters.threshold_dB / 20.0);

		// --- detect the signal
		DspFloatType detect_dB = detector.processAudioSample(xn);
		DspFloatType detectValue = pow(10.0, detect_dB / 20.0);
		DspFloatType deltaValue = detectValue - threshValue;

		ZVAFilterParameters filterParams = filter.getParameters();
		filterParams.fc = parameters.fc;

		// --- if above the threshold, modulate the filter fc
		if (deltaValue > 0.0)// || delta_dB > 0.0)
		{
			// --- fc Computer
			DspFloatType modulatorValue = 0.0;

			// --- best results are with linear values when detector is in dB mode
			modulatorValue = (deltaValue * parameters.sensitivity);

			// --- calculate modulated frequency
			filterParams.fc = doUnipolarModulationFromMin(modulatorValue, parameters.fc, kMaxFilterFrequency);
		}

		// --- update with new modulated frequency
		filter.setParameters(filterParams);

		// --- perform the filtering operation
		return filter.processAudioSample(xn);
	}

protected:
	EnvelopeFollowerParameters parameters; ///< object parameters

	// --- 1 filter and 1 detector
	ZVAFilter filter;		///< filter to modulate
	AudioDetector detector; ///< detector to track input signal
};


/**
\brief returns the storage component S(n) for delay-free loop solutions
- NOTES:\n
the storageComponent or "S" value is used for Zavalishin's VA filters as well
as the phaser APFs (Biquad) and is only available on two of the forms: direct
and transposed canonical\n
\returns the storage component of the filter
*/
DspFloatType Biquad::getS_value()
{
	storageComponent = 0.0;
	if (parameters.biquadCalcType == biquadAlgorithm::kDirect)
	{
		// --- 1)  form output y(n) = a0*x(n) + a1*x(n-1) + a2*x(n-2) - b1*y(n-1) - b2*y(n-2)
		storageComponent = coeffArray[a1] * stateArray[x_z1] +
			coeffArray[a2] * stateArray[x_z2] -
			coeffArray[b1] * stateArray[y_z1] -
			coeffArray[b2] * stateArray[y_z2];
	}
	else if (parameters.biquadCalcType == biquadAlgorithm::kTransposeCanonical)
	{
		// --- 1)  form output y(n) = a0*x(n) + stateArray[x_z1]
		storageComponent = stateArray[x_z1];
	}

	return storageComponent;
}

/**
\brief process one sample through the biquad
- RULES:\n
1) do all math required to form the output y(n), reading registers as required - do NOT write registers \n
2) check for underflow, which can happen with feedback structures\n
3) lastly, update the states of the z^-1 registers in the state array just before returning\n
- NOTES:\n
the storageComponent or "S" value is used for Zavalishin's VA filters and is only
available on two of the forms: direct and transposed canonical\n
\param xn the input sample x(n)
\returns the biquad processed output y(n)
*/
DspFloatType Biquad::processAudioSample(DspFloatType xn)
{
	if (parameters.biquadCalcType == biquadAlgorithm::kDirect)
	{
		// --- 1)  form output y(n) = a0*x(n) + a1*x(n-1) + a2*x(n-2) - b1*y(n-1) - b2*y(n-2)
		DspFloatType yn = coeffArray[a0] * xn + 
					coeffArray[a1] * stateArray[x_z1] +
					coeffArray[a2] * stateArray[x_z2] -
					coeffArray[b1] * stateArray[y_z1] -
					coeffArray[b2] * stateArray[y_z2];

		// --- 2) underflow check
		checkFloatUnderflow(yn);

		// --- 3) update states
		stateArray[x_z2] = stateArray[x_z1];
		stateArray[x_z1] = xn;

		stateArray[y_z2] = stateArray[y_z1];
		stateArray[y_z1] = yn;

		// --- return value
		return yn;
	}
	else if (parameters.biquadCalcType == biquadAlgorithm::kCanonical)
	{
		// --- 1)  form output y(n) = a0*w(n) + m_f_a1*stateArray[x_z1] + m_f_a2*stateArray[x_z2][x_z2];
		//
		// --- w(n) = x(n) - b1*stateArray[x_z1] - b2*stateArray[x_z2]
		DspFloatType wn = xn - coeffArray[b1] * stateArray[x_z1] - coeffArray[b2] * stateArray[x_z2];

		// --- y(n):
		DspFloatType yn = coeffArray[a0] * wn + coeffArray[a1] * stateArray[x_z1] + coeffArray[a2] * stateArray[x_z2];

		// --- 2) underflow check
		checkFloatUnderflow(yn);

		// --- 3) update states
		stateArray[x_z2] = stateArray[x_z1];
		stateArray[x_z1] = wn;

		// --- return value
		return yn;
	}
	else if (parameters.biquadCalcType == biquadAlgorithm::kTransposeDirect)
	{
		// --- 1)  form output y(n) = a0*w(n) + stateArray[x_z1]
		//
		// --- w(n) = x(n) + stateArray[y_z1]
		DspFloatType wn = xn + stateArray[y_z1];

		// --- y(n) = a0*w(n) + stateArray[x_z1]
		DspFloatType yn = coeffArray[a0] * wn + stateArray[x_z1];

		// --- 2) underflow check
		checkFloatUnderflow(yn);

		// --- 3) update states
		stateArray[y_z1] = stateArray[y_z2] - coeffArray[b1] * wn;
		stateArray[y_z2] = -coeffArray[b2] * wn;

		stateArray[x_z1] = stateArray[x_z2] + coeffArray[a1] * wn;
		stateArray[x_z2] = coeffArray[a2] * wn;

		// --- return value
		return yn;
	}
	else if (parameters.biquadCalcType == biquadAlgorithm::kTransposeCanonical)
	{
		// --- 1)  form output y(n) = a0*x(n) + stateArray[x_z1]
		DspFloatType yn = coeffArray[a0] * xn + stateArray[x_z1];

		// --- 2) underflow check
		checkFloatUnderflow(yn);

		// --- shuffle/update
		stateArray[x_z1] = coeffArray[a1]*xn - coeffArray[b1]*yn + stateArray[x_z2];
		stateArray[x_z2] = coeffArray[a2]*xn - coeffArray[b2]*yn;

		// --- return value
		return yn;
	}
	return xn; // didn't process anything :(
}

// --- returns true if coeffs were updated
bool AudioFilter::calculateFilterCoeffs()
{
	// --- clear coeff array
	memset(&coeffArray[0], 0, sizeof(DspFloatType)*numCoeffs);

	// --- set default pass-through
	coeffArray[a0] = 1.0;
	coeffArray[c0] = 1.0;
	coeffArray[d0] = 0.0;

	// --- grab these variables, to make calculations look more like the book
	filterAlgorithm algorithm = audioFilterParameters.algorithm;
	DspFloatType fc = audioFilterParameters.fc;
	DspFloatType Q = audioFilterParameters.Q;
	DspFloatType boostCut_dB = audioFilterParameters.boostCut_dB;

	// --- decode filter type and calculate accordingly
	// --- impulse invariabt LPF, matches closely with one-pole version,
	//     but diverges at VHF
	if (algorithm == filterAlgorithm::kImpInvLP1)
	{
		DspFloatType T = 1.0 / sampleRate;
		DspFloatType omega = 2.0*kPi*fc;
		DspFloatType eT = exp(-T*omega);

		coeffArray[a0] = 1.0 - eT; // <--- normalized by 1-e^aT
		coeffArray[a1] = 0.0;
		coeffArray[a2] = 0.0;
		coeffArray[b1] = -eT;
		coeffArray[b2] = 0.0;

		// --- update on calculator
		biquad.setCoefficients(coeffArray);

		// --- we updated
		return true;

	}
	else if (algorithm == filterAlgorithm::kImpInvLP2)
	{
		DspFloatType alpha = 2.0*kPi*fc / sampleRate;
		DspFloatType p_Re = -alpha / (2.0*Q);
		DspFloatType zeta = 1.0 / (2.0 * Q);
		DspFloatType p_Im = alpha*pow((1.0 - (zeta*zeta)), 0.5);
		DspFloatType c_Re = 0.0;
		DspFloatType c_Im = alpha / (2.0*pow((1.0 - (zeta*zeta)), 0.5));

		DspFloatType eP_re = exp(p_Re);
		coeffArray[a0] = c_Re;
		coeffArray[a1] = -2.0*(c_Re*cos(p_Im) + c_Im*sin(p_Im))*exp(p_Re);
		coeffArray[a2] = 0.0;
		coeffArray[b1] = -2.0*eP_re*cos(p_Im);
		coeffArray[b2] = eP_re*eP_re;

		// --- update on calculator
		biquad.setCoefficients(coeffArray);

		// --- we updated
		return true;
	}
	// --- kMatchLP2A = TIGHT fit LPF vicanek algo
	else if (algorithm == filterAlgorithm::kMatchLP2A)
	{
		// http://vicanek.de/articles/BiquadFits.pdf
		DspFloatType theta_c = 2.0*kPi*fc / sampleRate;

		DspFloatType q = 1.0 / (2.0*Q);

		// --- impulse invariant
		DspFloatType b_1 = 0.0;
		DspFloatType b_2 = exp(-2.0*q*theta_c);
		if (q <= 1.0)
		{
			b_1 = -2.0*exp(-q*theta_c)*cos(pow((1.0 - q*q), 0.5)*theta_c);
		}
		else
		{
			b_1 = -2.0*exp(-q*theta_c)*cosh(pow((q*q - 1.0), 0.5)*theta_c);
		}

		// --- TIGHT FIT --- //
		DspFloatType B0 = (1.0 + b_1 + b_2)*(1.0 + b_1 + b_2);
		DspFloatType B1 = (1.0 - b_1 + b_2)*(1.0 - b_1 + b_2);
		DspFloatType B2 = -4.0*b_2;

		DspFloatType phi_0 = 1.0 - sin(theta_c / 2.0)*sin(theta_c / 2.0);
		DspFloatType phi_1 = sin(theta_c / 2.0)*sin(theta_c / 2.0);
		DspFloatType phi_2 = 4.0*phi_0*phi_1;

		DspFloatType R1 = (B0*phi_0 + B1*phi_1 + B2*phi_2)*(Q*Q);
		DspFloatType A0 = B0;
		DspFloatType A1 = (R1 - A0*phi_0) / phi_1;

		if (A0 < 0.0)
			A0 = 0.0;
		if (A1 < 0.0)
			A1 = 0.0;

		DspFloatType a_0 = 0.5*(pow(A0, 0.5) + pow(A1, 0.5));
		DspFloatType a_1 = pow(A0, 0.5) - a_0;
		DspFloatType a_2 = 0.0;

		coeffArray[a0] = a_0;
		coeffArray[a1] = a_1;
		coeffArray[a2] = a_2;
		coeffArray[b1] = b_1;
		coeffArray[b2] = b_2;

		// --- update on calculator
		biquad.setCoefficients(coeffArray);

		// --- we updated
		return true;
	}
	// --- kMatchLP2B = LOOSE fit LPF vicanek algo
	else if (algorithm == filterAlgorithm::kMatchLP2B)
	{
		// http://vicanek.de/articles/BiquadFits.pdf
		DspFloatType theta_c = 2.0*kPi*fc / sampleRate;
		DspFloatType q = 1.0 / (2.0*Q);

		// --- impulse invariant
		DspFloatType b_1 = 0.0;
		DspFloatType b_2 = exp(-2.0*q*theta_c);
		if (q <= 1.0)
		{
			b_1 = -2.0*exp(-q*theta_c)*cos(pow((1.0 - q*q), 0.5)*theta_c);
		}
		else
		{
			b_1 = -2.0*exp(-q*theta_c)*cosh(pow((q*q - 1.0), 0.5)*theta_c);
		}

		// --- LOOSE FIT --- //
		DspFloatType f0 = theta_c / kPi; // note f0 = fraction of pi, so that f0 = 1.0 = pi = Nyquist

		DspFloatType r0 = 1.0 + b_1 + b_2;
		DspFloatType denom = (1.0 - f0*f0)*(1.0 - f0*f0) + (f0*f0) / (Q*Q);
		denom = pow(denom, 0.5);
		DspFloatType r1 = ((1.0 - b_1 + b_2)*f0*f0) / (denom);

		DspFloatType a_0 = (r0 + r1) / 2.0;
		DspFloatType a_1 = r0 - a_0;
		DspFloatType a_2 = 0.0;

		coeffArray[a0] = a_0;
		coeffArray[a1] = a_1;
		coeffArray[a2] = a_2;
		coeffArray[b1] = b_1;
		coeffArray[b2] = b_2;

		// --- update on calculator
		biquad.setCoefficients(coeffArray);

		// --- we updated
		return true;
	}
	// --- kMatchBP2A = TIGHT fit BPF vicanek algo
	else if (algorithm == filterAlgorithm::kMatchBP2A)
	{
		// http://vicanek.de/articles/BiquadFits.pdf
		DspFloatType theta_c = 2.0*kPi*fc / sampleRate;
		DspFloatType q = 1.0 / (2.0*Q);

		// --- impulse invariant
		DspFloatType b_1 = 0.0;
		DspFloatType b_2 = exp(-2.0*q*theta_c);
		if (q <= 1.0)
		{
			b_1 = -2.0*exp(-q*theta_c)*cos(pow((1.0 - q*q), 0.5)*theta_c);
		}
		else
		{
			b_1 = -2.0*exp(-q*theta_c)*cosh(pow((q*q - 1.0), 0.5)*theta_c);
		}

		// --- TIGHT FIT --- //
		DspFloatType B0 = (1.0 + b_1 + b_2)*(1.0 + b_1 + b_2);
		DspFloatType B1 = (1.0 - b_1 + b_2)*(1.0 - b_1 + b_2);
		DspFloatType B2 = -4.0*b_2;

		DspFloatType phi_0 = 1.0 - sin(theta_c / 2.0)*sin(theta_c / 2.0);
		DspFloatType phi_1 = sin(theta_c / 2.0)*sin(theta_c / 2.0);
		DspFloatType phi_2 = 4.0*phi_0*phi_1;

		DspFloatType R1 = B0*phi_0 + B1*phi_1 + B2*phi_2;
		DspFloatType R2 = -B0 + B1 + 4.0*(phi_0 - phi_1)*B2;

		DspFloatType A2 = (R1 - R2*phi_1) / (4.0*phi_1*phi_1);
		DspFloatType A1 = R2 + 4.0*(phi_1 - phi_0)*A2;

		DspFloatType a_1 = -0.5*(pow(A1, 0.5));
		DspFloatType a_0 = 0.5*(pow((A2 + (a_1*a_1)), 0.5) - a_1);
		DspFloatType a_2 = -a_0 - a_1;

		coeffArray[a0] = a_0;
		coeffArray[a1] = a_1;
		coeffArray[a2] = a_2;
		coeffArray[b1] = b_1;
		coeffArray[b2] = b_2;

		// --- update on calculator
		biquad.setCoefficients(coeffArray);

		// --- we updated
		return true;
	}
	// --- kMatchBP2B = LOOSE fit BPF vicanek algo
	else if (algorithm == filterAlgorithm::kMatchBP2B)
	{
		// http://vicanek.de/articles/BiquadFits.pdf
		DspFloatType theta_c = 2.0*kPi*fc / sampleRate;
		DspFloatType q = 1.0 / (2.0*Q);

		// --- impulse invariant
		DspFloatType b_1 = 0.0;
		DspFloatType b_2 = exp(-2.0*q*theta_c);
		if (q <= 1.0)
		{
			b_1 = -2.0*exp(-q*theta_c)*cos(pow((1.0 - q*q), 0.5)*theta_c);
		}
		else
		{
			b_1 = -2.0*exp(-q*theta_c)*cosh(pow((q*q - 1.0), 0.5)*theta_c);
		}

		// --- LOOSE FIT --- //
		DspFloatType f0 = theta_c / kPi; // note f0 = fraction of pi, so that f0 = 1.0 = pi = Nyquist

		DspFloatType r0 = (1.0 + b_1 + b_2) / (kPi*f0*Q);
		DspFloatType denom = (1.0 - f0*f0)*(1.0 - f0*f0) + (f0*f0) / (Q*Q);
		denom = pow(denom, 0.5);

		DspFloatType r1 = ((1.0 - b_1 + b_2)*(f0 / Q)) / (denom);

		DspFloatType a_1 = -r1 / 2.0;
		DspFloatType a_0 = (r0 - a_1) / 2.0;
		DspFloatType a_2 = -a_0 - a_1;

		coeffArray[a0] = a_0;
		coeffArray[a1] = a_1;
		coeffArray[a2] = a_2;
		coeffArray[b1] = b_1;
		coeffArray[b2] = b_2;

		// --- update on calculator
		biquad.setCoefficients(coeffArray);

		// --- we updated
		return true;
	}
	else if (algorithm == filterAlgorithm::kLPF1P)
	{
		// --- see book for formulae
		DspFloatType theta_c = 2.0*kPi*fc / sampleRate;
		DspFloatType gamma = 2.0 - cos(theta_c);

		DspFloatType filter_b1 = pow((gamma*gamma - 1.0), 0.5) - gamma;
		DspFloatType filter_a0 = 1.0 + filter_b1;

		// --- update coeffs
		coeffArray[a0] = filter_a0;
		coeffArray[a1] = 0.0;
		coeffArray[a2] = 0.0;
		coeffArray[b1] = filter_b1;
		coeffArray[b2] = 0.0;

		// --- update on calculator
		biquad.setCoefficients(coeffArray);

		// --- we updated
		return true;
	}
	else if (algorithm == filterAlgorithm::kLPF1)
	{
		// --- see book for formulae
		DspFloatType theta_c = 2.0*kPi*fc / sampleRate;
		DspFloatType gamma = cos(theta_c) / (1.0 + sin(theta_c));

		// --- update coeffs
		coeffArray[a0] = (1.0 - gamma) / 2.0;
		coeffArray[a1] = (1.0 - gamma) / 2.0;
		coeffArray[a2] = 0.0;
		coeffArray[b1] = -gamma;
		coeffArray[b2] = 0.0;

		// --- update on calculator
		biquad.setCoefficients(coeffArray);

		// --- we updated
		return true;
	}
	else if (algorithm == filterAlgorithm::kHPF1)
	{
		// --- see book for formulae
		DspFloatType theta_c = 2.0*kPi*fc / sampleRate;
		DspFloatType gamma = cos(theta_c) / (1.0 + sin(theta_c));

		// --- update coeffs
		coeffArray[a0] = (1.0 + gamma) / 2.0;
		coeffArray[a1] = -(1.0 + gamma) / 2.0;
		coeffArray[a2] = 0.0;
		coeffArray[b1] = -gamma;
		coeffArray[b2] = 0.0;

		// --- update on calculator
		biquad.setCoefficients(coeffArray);

		// --- we updated
		return true;
	}
	else if (algorithm == filterAlgorithm::kLPF2)
	{
		// --- see book for formulae
		DspFloatType theta_c = 2.0*kPi*fc / sampleRate;
		DspFloatType d = 1.0 / Q;
		DspFloatType betaNumerator = 1.0 - ((d / 2.0)*(sin(theta_c)));
		DspFloatType betaDenominator = 1.0 + ((d / 2.0)*(sin(theta_c)));

		DspFloatType beta = 0.5*(betaNumerator / betaDenominator);
		DspFloatType gamma = (0.5 + beta)*(cos(theta_c));
		DspFloatType alpha = (0.5 + beta - gamma) / 2.0;

		// --- update coeffs
		coeffArray[a0] = alpha;
		coeffArray[a1] = 2.0*alpha;
		coeffArray[a2] = alpha;
		coeffArray[b1] = -2.0*gamma;
		coeffArray[b2] = 2.0*beta;

	//	DspFloatType mag = getMagResponse(theta_c, coeffArray[a0], coeffArray[a1], coeffArray[a2], coeffArray[b1], coeffArray[b2]);
		// --- update on calculator
		biquad.setCoefficients(coeffArray);

		// --- we updated
		return true;
	}
	else if (algorithm == filterAlgorithm::kHPF2)
	{
		// --- see book for formulae
		DspFloatType theta_c = 2.0*kPi*fc / sampleRate;
		DspFloatType d = 1.0 / Q;

		DspFloatType betaNumerator = 1.0 - ((d / 2.0)*(sin(theta_c)));
		DspFloatType betaDenominator = 1.0 + ((d / 2.0)*(sin(theta_c)));

		DspFloatType beta = 0.5*(betaNumerator / betaDenominator);
		DspFloatType gamma = (0.5 + beta)*(cos(theta_c));
		DspFloatType alpha = (0.5 + beta + gamma) / 2.0;

		// --- update coeffs
		coeffArray[a0] = alpha;
		coeffArray[a1] = -2.0*alpha;
		coeffArray[a2] = alpha;
		coeffArray[b1] = -2.0*gamma;
		coeffArray[b2] = 2.0*beta;

		// --- update on calculator
		biquad.setCoefficients(coeffArray);

		// --- we updated
		return true;
	}
	else if (algorithm == filterAlgorithm::kBPF2)
	{
		// --- see book for formulae
		DspFloatType K = tan(kPi*fc / sampleRate);
		DspFloatType delta = K*K*Q + K + Q;

		// --- update coeffs
		coeffArray[a0] = K / delta;;
		coeffArray[a1] = 0.0;
		coeffArray[a2] = -K / delta;
		coeffArray[b1] = 2.0*Q*(K*K - 1) / delta;
		coeffArray[b2] = (K*K*Q - K + Q) / delta;

		// --- update on calculator
		biquad.setCoefficients(coeffArray);

		// --- we updated
		return true;
	}
	else if (algorithm == filterAlgorithm::kBSF2)
	{
		// --- see book for formulae
		DspFloatType K = tan(kPi*fc / sampleRate);
		DspFloatType delta = K*K*Q + K + Q;

		// --- update coeffs
		coeffArray[a0] = Q*(1 + K*K) / delta;
		coeffArray[a1] = 2.0*Q*(K*K - 1) / delta;
		coeffArray[a2] = Q*(1 + K*K) / delta;
		coeffArray[b1] = 2.0*Q*(K*K - 1) / delta;
		coeffArray[b2] = (K*K*Q - K + Q) / delta;

		// --- update on calculator
		biquad.setCoefficients(coeffArray);

		// --- we updated
		return true;
	}
	else if (algorithm == filterAlgorithm::kButterLPF2)
	{
		// --- see book for formulae
		DspFloatType theta_c = kPi*fc / sampleRate;
		DspFloatType C = 1.0 / tan(theta_c);

		// --- update coeffs
		coeffArray[a0] = 1.0 / (1.0 + kSqrtTwo*C + C*C);
		coeffArray[a1] = 2.0*coeffArray[a0];
		coeffArray[a2] = coeffArray[a0];
		coeffArray[b1] = 2.0*coeffArray[a0] * (1.0 - C*C);
		coeffArray[b2] = coeffArray[a0] * (1.0 - kSqrtTwo*C + C*C);

		// --- update on calculator
		biquad.setCoefficients(coeffArray);

		// --- we updated
		return true;
	}
	else if (algorithm == filterAlgorithm::kButterHPF2)
	{
		// --- see book for formulae
		DspFloatType theta_c = kPi*fc / sampleRate;
		DspFloatType C = tan(theta_c);

		// --- update coeffs
		coeffArray[a0] = 1.0 / (1.0 + kSqrtTwo*C + C*C);
		coeffArray[a1] = -2.0*coeffArray[a0];
		coeffArray[a2] = coeffArray[a0];
		coeffArray[b1] = 2.0*coeffArray[a0] * (C*C - 1.0);
		coeffArray[b2] = coeffArray[a0] * (1.0 - kSqrtTwo*C + C*C);

		// --- update on calculator
		biquad.setCoefficients(coeffArray);

		// --- we updated
		return true;
	}
	else if (algorithm == filterAlgorithm::kButterBPF2)
	{
		// --- see book for formulae
		DspFloatType theta_c = 2.0*kPi*fc / sampleRate;
		DspFloatType BW = fc / Q;
		DspFloatType delta_c = kPi*BW / sampleRate;
		if (delta_c >= 0.95*kPi / 2.0) delta_c = 0.95*kPi / 2.0;

		DspFloatType C = 1.0 / tan(delta_c);
		DspFloatType D = 2.0*cos(theta_c);

		// --- update coeffs
		coeffArray[a0] = 1.0 / (1.0 + C);
		coeffArray[a1] = 0.0;
		coeffArray[a2] = -coeffArray[a0];
		coeffArray[b1] = -coeffArray[a0] * (C*D);
		coeffArray[b2] = coeffArray[a0] * (C - 1.0);

		// --- update on calculator
		biquad.setCoefficients(coeffArray);

		// --- we updated
		return true;
	}
	else if (algorithm == filterAlgorithm::kButterBSF2)
	{
		// --- see book for formulae
		DspFloatType theta_c = 2.0*kPi*fc / sampleRate;
		DspFloatType BW = fc / Q;
		DspFloatType delta_c = kPi*BW / sampleRate;
		if (delta_c >= 0.95*kPi / 2.0) delta_c = 0.95*kPi / 2.0;

		DspFloatType C = tan(delta_c);
		DspFloatType D = 2.0*cos(theta_c);

		// --- update coeffs
		coeffArray[a0] = 1.0 / (1.0 + C);
		coeffArray[a1] = -coeffArray[a0] * D;
		coeffArray[a2] = coeffArray[a0];
		coeffArray[b1] = -coeffArray[a0] * D;
		coeffArray[b2] = coeffArray[a0] * (1.0 - C);

		// --- update on calculator
		biquad.setCoefficients(coeffArray);

		// --- we updated
		return true;
	}
	else if (algorithm == filterAlgorithm::kMMALPF2 || algorithm == filterAlgorithm::kMMALPF2B)
	{
		// --- see book for formulae
		DspFloatType theta_c = 2.0*kPi*fc / sampleRate;
		DspFloatType resonance_dB = 0;

		if (Q > 0.707)
		{
			DspFloatType peak = Q*Q / pow(Q*Q - 0.25, 0.5);
			resonance_dB = 20.0*log10(peak);
		}

		// --- intermediate vars
		DspFloatType resonance = (cos(theta_c) + (sin(theta_c) * sqrt(pow(10.0, (resonance_dB / 10.0)) - 1))) / ((pow(10.0, (resonance_dB / 20.0)) * sin(theta_c)) + 1);
		DspFloatType g = pow(10.0, (-resonance_dB / 40.0));

		// --- kMMALPF2B disables the GR with increase in Q
		if (algorithm == filterAlgorithm::kMMALPF2B)
			g = 1.0;

		DspFloatType filter_b1 = (-2.0) * resonance * cos(theta_c);
		DspFloatType filter_b2 = resonance * resonance;
		DspFloatType filter_a0 = g * (1 + filter_b1 + filter_b2);

		// --- update coeffs
		coeffArray[a0] = filter_a0;
		coeffArray[a1] = 0.0;
		coeffArray[a2] = 0.0;
		coeffArray[b1] = filter_b1;
		coeffArray[b2] = filter_b2;

		// --- update on calculator
		biquad.setCoefficients(coeffArray);

		// --- we updated
		return true;
	}
	else if (algorithm == filterAlgorithm::kLowShelf)
	{
		// --- see book for formulae
		DspFloatType theta_c = 2.0*kPi*fc / sampleRate;
		DspFloatType mu = pow(10.0, boostCut_dB / 20.0);

		DspFloatType beta = 4.0 / (1.0 + mu);
		DspFloatType delta = beta*tan(theta_c / 2.0);
		DspFloatType gamma = (1.0 - delta) / (1.0 + delta);

		// --- update coeffs
		coeffArray[a0] = (1.0 - gamma) / 2.0;
		coeffArray[a1] = (1.0 - gamma) / 2.0;
		coeffArray[a2] = 0.0;
		coeffArray[b1] = -gamma;
		coeffArray[b2] = 0.0;

		coeffArray[c0] = mu - 1.0;
		coeffArray[d0] = 1.0;

		// --- update on calculator
		biquad.setCoefficients(coeffArray);

		// --- we updated
		return true;
	}
	else if (algorithm == filterAlgorithm::kHiShelf)
	{
		DspFloatType theta_c = 2.0*kPi*fc / sampleRate;
		DspFloatType mu = pow(10.0, boostCut_dB / 20.0);

		DspFloatType beta = (1.0 + mu) / 4.0;
		DspFloatType delta = beta*tan(theta_c / 2.0);
		DspFloatType gamma = (1.0 - delta) / (1.0 + delta);

		coeffArray[a0] = (1.0 + gamma) / 2.0;
		coeffArray[a1] = -coeffArray[a0];
		coeffArray[a2] = 0.0;
		coeffArray[b1] = -gamma;
		coeffArray[b2] = 0.0;

		coeffArray[c0] = mu - 1.0;
		coeffArray[d0] = 1.0;

		// --- update on calculator
		biquad.setCoefficients(coeffArray);

		// --- we updated
		return true;
	}
	else if (algorithm == filterAlgorithm::kCQParaEQ)
	{
		// --- see book for formulae
		DspFloatType K = tan(kPi*fc / sampleRate);
		DspFloatType Vo = pow(10.0, boostCut_dB / 20.0);
		bool bBoost = boostCut_dB >= 0 ? true : false;

		DspFloatType d0 = 1.0 + (1.0 / Q)*K + K*K;
		DspFloatType e0 = 1.0 + (1.0 / (Vo*Q))*K + K*K;
		DspFloatType alpha = 1.0 + (Vo / Q)*K + K*K;
		DspFloatType beta = 2.0*(K*K - 1.0);
		DspFloatType gamma = 1.0 - (Vo / Q)*K + K*K;
		DspFloatType delta = 1.0 - (1.0 / Q)*K + K*K;
		DspFloatType eta = 1.0 - (1.0 / (Vo*Q))*K + K*K;

		// --- update coeffs
		coeffArray[a0] = bBoost ? alpha / d0 : d0 / e0;
		coeffArray[a1] = bBoost ? beta / d0 : beta / e0;
		coeffArray[a2] = bBoost ? gamma / d0 : delta / e0;
		coeffArray[b1] = bBoost ? beta / d0 : beta / e0;
		coeffArray[b2] = bBoost ? delta / d0 : eta / e0;

		// --- update on calculator
		biquad.setCoefficients(coeffArray);

		// --- we updated
		return true;
	}
	else if (algorithm == filterAlgorithm::kNCQParaEQ)
	{
		// --- see book for formulae
		DspFloatType theta_c = 2.0*kPi*fc / sampleRate;
		DspFloatType mu = pow(10.0, boostCut_dB / 20.0);

		// --- clamp to 0.95 pi/2 (you can experiment with this)
		DspFloatType tanArg = theta_c / (2.0 * Q);
		if (tanArg >= 0.95*kPi / 2.0) tanArg = 0.95*kPi / 2.0;

		// --- intermediate variables (you can condense this if you wish)
		DspFloatType zeta = 4.0 / (1.0 + mu);
		DspFloatType betaNumerator = 1.0 - zeta*tan(tanArg);
		DspFloatType betaDenominator = 1.0 + zeta*tan(tanArg);

		DspFloatType beta = 0.5*(betaNumerator / betaDenominator);
		DspFloatType gamma = (0.5 + beta)*(cos(theta_c));
		DspFloatType alpha = (0.5 - beta);

		// --- update coeffs
		coeffArray[a0] = alpha;
		coeffArray[a1] = 0.0;
		coeffArray[a2] = -alpha;
		coeffArray[b1] = -2.0*gamma;
		coeffArray[b2] = 2.0*beta;

		coeffArray[c0] = mu - 1.0;
		coeffArray[d0] = 1.0;

		// --- update on calculator
		biquad.setCoefficients(coeffArray);

		// --- we updated
		return true;
	}
	else if (algorithm == filterAlgorithm::kLWRLPF2)
	{
		// --- see book for formulae
		DspFloatType omega_c = kPi*fc;
		DspFloatType theta_c = kPi*fc / sampleRate;

		DspFloatType k = omega_c / tan(theta_c);
		DspFloatType denominator = k*k + omega_c*omega_c + 2.0*k*omega_c;
		DspFloatType b1_Num = -2.0*k*k + 2.0*omega_c*omega_c;
		DspFloatType b2_Num = -2.0*k*omega_c + k*k + omega_c*omega_c;

		// --- update coeffs
		coeffArray[a0] = omega_c*omega_c / denominator;
		coeffArray[a1] = 2.0*omega_c*omega_c / denominator;
		coeffArray[a2] = coeffArray[a0];
		coeffArray[b1] = b1_Num / denominator;
		coeffArray[b2] = b2_Num / denominator;

		// --- update on calculator
		biquad.setCoefficients(coeffArray);

		// --- we updated
		return true;
	}
	else if (algorithm == filterAlgorithm::kLWRHPF2)
	{
		// --- see book for formulae
		DspFloatType omega_c = kPi*fc;
		DspFloatType theta_c = kPi*fc / sampleRate;

		DspFloatType k = omega_c / tan(theta_c);
		DspFloatType denominator = k*k + omega_c*omega_c + 2.0*k*omega_c;
		DspFloatType b1_Num = -2.0*k*k + 2.0*omega_c*omega_c;
		DspFloatType b2_Num = -2.0*k*omega_c + k*k + omega_c*omega_c;

		// --- update coeffs
		coeffArray[a0] = k*k / denominator;
		coeffArray[a1] = -2.0*k*k / denominator;
		coeffArray[a2] = coeffArray[a0];
		coeffArray[b1] = b1_Num / denominator;
		coeffArray[b2] = b2_Num / denominator;

		// --- update on calculator
		biquad.setCoefficients(coeffArray);

		// --- we updated
		return true;
	}
	else if (algorithm == filterAlgorithm::kAPF1)
	{
		// --- see book for formulae
		DspFloatType alphaNumerator = tan((kPi*fc) / sampleRate) - 1.0;
		DspFloatType alphaDenominator = tan((kPi*fc) / sampleRate) + 1.0;
		DspFloatType alpha = alphaNumerator / alphaDenominator;

		// --- update coeffs
		coeffArray[a0] = alpha;
		coeffArray[a1] = 1.0;
		coeffArray[a2] = 0.0;
		coeffArray[b1] = alpha;
		coeffArray[b2] = 0.0;

		// --- update on calculator
		biquad.setCoefficients(coeffArray);

		// --- we updated
		return true;
	}
	else if (algorithm == filterAlgorithm::kAPF2)
	{
		// --- see book for formulae
		DspFloatType theta_c = 2.0*kPi*fc / sampleRate;
		DspFloatType BW = fc / Q;
		DspFloatType argTan = kPi*BW / sampleRate;
		if (argTan >= 0.95*kPi / 2.0) argTan = 0.95*kPi / 2.0;

		DspFloatType alphaNumerator = tan(argTan) - 1.0;
		DspFloatType alphaDenominator = tan(argTan) + 1.0;
		DspFloatType alpha = alphaNumerator / alphaDenominator;
		DspFloatType beta = -cos(theta_c);

		// --- update coeffs
		coeffArray[a0] = -alpha;
		coeffArray[a1] = beta*(1.0 - alpha);
		coeffArray[a2] = 1.0;
		coeffArray[b1] = beta*(1.0 - alpha);
		coeffArray[b2] = -alpha;

		// --- update on calculator
		biquad.setCoefficients(coeffArray);

		// --- we updated
		return true;
	}
	else if (algorithm == filterAlgorithm::kResonA)
	{
		// --- see book for formulae
		DspFloatType theta_c = 2.0*kPi*fc / sampleRate;
		DspFloatType BW = fc / Q;
		DspFloatType filter_b2 = exp(-2.0*kPi*(BW / sampleRate));
		DspFloatType filter_b1 = ((-4.0*filter_b2) / (1.0 + filter_b2))*cos(theta_c);
		DspFloatType filter_a0 = (1.0 - filter_b2)*pow((1.0 - (filter_b1*filter_b1) / (4.0 * filter_b2)), 0.5);

		// --- update coeffs
		coeffArray[a0] = filter_a0;
		coeffArray[a1] = 0.0;
		coeffArray[a2] = 0.0;
		coeffArray[b1] = filter_b1;
		coeffArray[b2] = filter_b2;

		// --- update on calculator
		biquad.setCoefficients(coeffArray);

		// --- we updated
		return true;
	}
	else if (algorithm == filterAlgorithm::kResonB)
	{
		// --- see book for formulae
		DspFloatType theta_c = 2.0*kPi*fc / sampleRate;
		DspFloatType BW = fc / Q;
		DspFloatType filter_b2 = exp(-2.0*kPi*(BW / sampleRate));
		DspFloatType filter_b1 = ((-4.0*filter_b2) / (1.0 + filter_b2))*cos(theta_c);
		DspFloatType filter_a0 = 1.0 - pow(filter_b2, 0.5); // (1.0 - filter_b2)*pow((1.0 - (filter_b1*filter_b1) / (4.0 * filter_b2)), 0.5);

		// --- update coeffs
		coeffArray[a0] = filter_a0;
		coeffArray[a1] = 0.0;
		coeffArray[a2] = -filter_a0;
		coeffArray[b1] = filter_b1;
		coeffArray[b2] = filter_b2;

		// --- update on calculator
		biquad.setCoefficients(coeffArray);

		// --- we updated
		return true;
	}

	// --- we did n't update :(
	return false;
}

/**
\brief process one sample through the audio filter
- NOTES:\n
Uses the modified biquaqd structure that includes the wet and dry signal coefficients c and d.\n
Here the biquad object does all of the work and we simply combine the wet and dry signals.\n
// return (dry) + (processed): x(n)*d0 + y(n)*c0\n
\param xn the input sample x(n)
\returns the biquad processed output y(n)
*/
DspFloatType AudioFilter::processAudioSample(DspFloatType xn)
{
	// --- let biquad do the grunt-work
	//
	// return (dry) + (processed): x(n)*d0 + y(n)*c0
	return coeffArray[d0] * xn + coeffArray[c0] * biquad.processAudioSample(xn);
}


/**
\struct TwoBandShelvingFilterParameters
\ingroup FX-Objects
\brief
Custom parameter structure for the TwoBandShelvingFilter object. Used for reverb algorithms in book.
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
struct TwoBandShelvingFilterParameters
{
	TwoBandShelvingFilterParameters() {}
	/** all FXObjects parameter objects require overloaded= operator so remember to add new entries if you add new variables. */
	TwoBandShelvingFilterParameters& operator=(const TwoBandShelvingFilterParameters& params)	// need this override for collections to work
	{
		if (this == &params)
			return *this;

		lowShelf_fc = params.lowShelf_fc;
		lowShelfBoostCut_dB = params.lowShelfBoostCut_dB;
		highShelf_fc = params.highShelf_fc;
		highShelfBoostCut_dB = params.highShelfBoostCut_dB;
		return *this;
	}

	// --- individual parameters
	DspFloatType lowShelf_fc = 0.0;			///< fc for low shelf
	DspFloatType lowShelfBoostCut_dB = 0.0;	///< low shelf gain
	DspFloatType highShelf_fc = 0.0;			///< fc for high shelf
	DspFloatType highShelfBoostCut_dB = 0.0;	///< high shelf gain
};

/**
\class TwoBandShelvingFilter
\ingroup FX-Objects
\brief
The TwoBandShelvingFilter object implements two shelving filters in series in the standard "Bass and Treble" configuration.
Audio I/O:
- Processes mono input to mono output.
Control I/F:
- Use TwoBandShelvingFilterParameters structure to get/set object params.
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
class TwoBandShelvingFilter : public IAudioSignalProcessor
{
public:
	TwoBandShelvingFilter()
	{
		AudioFilterParameters params = lowShelfFilter.getParameters();
		params.algorithm = filterAlgorithm::kLowShelf;
		lowShelfFilter.setParameters(params);

		params = highShelfFilter.getParameters();
		params.algorithm = filterAlgorithm::kHiShelf;
		highShelfFilter.setParameters(params);
	}		/* C-TOR */

	~TwoBandShelvingFilter() {}		/* D-TOR */

	/** reset members to initialized state */
	virtual bool reset(DspFloatType _sampleRate)
	{
		lowShelfFilter.reset(_sampleRate);
		highShelfFilter.reset(_sampleRate);
		return true;
	}

	/** return false: this object only processes samples */
	virtual bool canProcessAudioFrame() { return false; }

	/** process a single input through the two filters in series */
	/**
	\param xn input
	\return the processed sample
	*/
	virtual DspFloatType processAudioSample(DspFloatType xn)
	{
		// --- all modes do Full Wave Rectification
		DspFloatType filteredSignal = lowShelfFilter.processAudioSample(xn);
		filteredSignal = highShelfFilter.processAudioSample(filteredSignal);

		return filteredSignal;
	}

	/** get parameters: note use of custom structure for passing param data */
	/**
	\return TwoBandShelvingFilterParameters custom data structure
	*/
	TwoBandShelvingFilterParameters getParameters()
	{
		return parameters;
	}

	/** set parameters: note use of custom structure for passing param data */
	/**
	\param TwoBandShelvingFilterParameters custom data structure
	*/
	void setParameters(const TwoBandShelvingFilterParameters& params)
	{
		parameters = params;
		AudioFilterParameters filterParams = lowShelfFilter.getParameters();
		filterParams.fc = parameters.lowShelf_fc;
		filterParams.boostCut_dB = parameters.lowShelfBoostCut_dB;
		lowShelfFilter.setParameters(filterParams);

		filterParams = highShelfFilter.getParameters();
		filterParams.fc = parameters.highShelf_fc;
		filterParams.boostCut_dB = parameters.highShelfBoostCut_dB;
		highShelfFilter.setParameters(filterParams);
	}

private:
	TwoBandShelvingFilterParameters parameters; ///< object parameters
	AudioFilter lowShelfFilter;					///< filter for low shelf
	AudioFilter highShelfFilter;				///< filter for high shelf
};

