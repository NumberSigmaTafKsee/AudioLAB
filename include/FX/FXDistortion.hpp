#pragma once

#include "FXObjects.hpp"

/**
\enum distortionModel
\ingroup Constants-Enums
\brief
Use this strongly typed enum to easily set the waveshaper model for the Triode objects
- enum class distortionModel { kSoftClip, kArcTan, kFuzzAsym };
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
enum class distortionModel { kSoftClip, kArcTan, kFuzzAsym };

/**
\struct TriodeClassAParameters
\ingroup FX-Objects
\brief
Custom parameter structure for the TriodeClassA object.
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
struct TriodeClassAParameters
{
	TriodeClassAParameters() {}
	/** all FXObjects parameter objects require overloaded= operator so remember to add new entries if you add new variables. */
	TriodeClassAParameters& operator=(const TriodeClassAParameters& params)	// need this override for collections to work
	{
		if (this == &params)
			return *this;

		waveshaper = params.waveshaper;
		saturation = params.saturation;
		asymmetry = params.asymmetry;
		outputGain = params.outputGain;

		invertOutput = params.invertOutput;
		enableHPF = params.enableHPF;
		enableLSF = params.enableLSF;

		hpf_Fc = params.hpf_Fc;
		lsf_Fshelf = params.lsf_Fshelf;
		lsf_BoostCut_dB = params.lsf_BoostCut_dB;

		return *this;
	}

	// --- individual parameters
	distortionModel waveshaper = distortionModel::kSoftClip; ///< waveshaper

	double saturation = 1.0;	///< saturation level
	double asymmetry = 0.0;		///< asymmetry level
	double outputGain = 1.0;	///< outputGain level

	bool invertOutput = true;	///< invertOutput - triodes invert output
	bool enableHPF = true;		///< HPF simulates DC blocking cap on output
	bool enableLSF = false;		///< LSF simulates shelf due to cathode self biasing

	double hpf_Fc = 1.0;		///< fc of DC blocking cap
	double lsf_Fshelf = 80.0;	///< shelf fc from self bias cap
	double lsf_BoostCut_dB = 0.0;///< boost/cut due to cathode self biasing
};

/**
\class TriodeClassA
\ingroup FX-Objects
\brief
The TriodeClassA object simulates a triode in class A configuration. This is a very simple and basic simulation
and a starting point for other designs; it is not intended to be a full-fledged triode simulator.
Audio I/O:
- Processes mono input to mono output.
Control I/F:
- Use TriodeClassAParameters structure to get/set object params.
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
class TriodeClassA : public IAudioSignalProcessor
{
public:
	TriodeClassA() {
		AudioFilterParameters params;
		params.algorithm = filterAlgorithm::kHPF1;
		params.fc = parameters.hpf_Fc;
		outputHPF.setParameters(params);

		params.algorithm = filterAlgorithm::kLowShelf;
		params.fc = parameters.lsf_Fshelf;
		params.boostCut_dB = parameters.lsf_BoostCut_dB;
		outputLSF.setParameters(params);
	}		/* C-TOR */
	~TriodeClassA() {}		/* D-TOR */

	/** reset members to initialized state */
	virtual bool reset(double _sampleRate)
	{
		outputHPF.reset(_sampleRate);
		outputLSF.reset(_sampleRate);

		// ---
		return true;
	}

	/** get parameters: note use of custom structure for passing param data */
	/**
	\return TriodeClassAParameters custom data structure
	*/
	TriodeClassAParameters getParameters() { return parameters; }

	/** set parameters: note use of custom structure for passing param data */
	/**
	\param TriodeClassAParameters custom data structure
	*/
	void setParameters(const TriodeClassAParameters& params)
	{
		parameters = params;

		AudioFilterParameters filterParams;
		filterParams.algorithm = filterAlgorithm::kHPF1;
		filterParams.fc = parameters.hpf_Fc;
		outputHPF.setParameters(filterParams);

		filterParams.algorithm = filterAlgorithm::kLowShelf;
		filterParams.fc = parameters.lsf_Fshelf;
		filterParams.boostCut_dB = parameters.lsf_BoostCut_dB;
		outputLSF.setParameters(filterParams);
	}

	/** return false: this object only processes samples */
	virtual bool canProcessAudioFrame() { return false; }

	/** do the triode simulation to process one input to one output*/
	/**
	\param xn input
	\return the processed sample
	*/
	virtual double processAudioSample(double xn)
	{
		// --- perform waveshaping
		double output = 0.0;

		if (parameters.waveshaper == distortionModel::kSoftClip)
			output = softClipWaveShaper(xn, parameters.saturation);
		else if (parameters.waveshaper == distortionModel::kArcTan)
			output = atanWaveShaper(xn, parameters.saturation);
		else if (parameters.waveshaper == distortionModel::kFuzzAsym)
			output = fuzzExp1WaveShaper(xn, parameters.saturation, parameters.asymmetry);

		// --- inversion, normal for plate of class A triode
		if (parameters.invertOutput)
			output *= -1.0;

		// --- Output (plate) capacitor = HPF, remove DC offset
		if (parameters.enableHPF)
			output = outputHPF.processAudioSample(output);

		// --- if cathode resistor bypass, will create low shelf
		if (parameters.enableLSF)
			output = outputLSF.processAudioSample(output);

		// --- final resistor divider/potentiometer
		output *= parameters.outputGain;

		return output;
	}

protected:
	TriodeClassAParameters parameters;	///< object parameters
	AudioFilter outputHPF;				///< HPF to simulate output DC blocking cap
	AudioFilter outputLSF;				///< LSF to simulate shelf caused by cathode self-biasing cap
};

const unsigned int NUM_TUBES = 4;

/**
\struct ClassATubePreParameters
\ingroup FX-Objects
\brief
Custom parameter structure for the ClassATubePre object.
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
struct ClassATubePreParameters
{
	ClassATubePreParameters() {}
	/** all FXObjects parameter objects require overloaded= operator so remember to add new entries if you add new variables. */
	ClassATubePreParameters& operator=(const ClassATubePreParameters& params)	// need this override for collections to work
	{
		if (this == &params)
			return *this;

		inputLevel_dB = params.inputLevel_dB;
		saturation = params.saturation;
		asymmetry = params.asymmetry;
		outputLevel_dB = params.outputLevel_dB;

		lowShelf_fc = params.lowShelf_fc;
		lowShelfBoostCut_dB = params.lowShelfBoostCut_dB;
		highShelf_fc = params.highShelf_fc;
		highShelfBoostCut_dB = params.highShelfBoostCut_dB;

		return *this;
	}

	// --- individual parameters
	double inputLevel_dB = 0.0;		///< input level in dB
	double saturation = 0.0;		///< input level in dB
	double asymmetry = 0.0;			///< input level in dB
	double outputLevel_dB = 0.0;	///< input level in dB

	// --- shelving filter params
	double lowShelf_fc = 0.0;			///< LSF shelf frequency
	double lowShelfBoostCut_dB = 0.0;	///< LSF shelf gain/cut
	double highShelf_fc = 0.0;			///< HSF shelf frequency
	double highShelfBoostCut_dB = 0.0;	///< HSF shelf frequency

};

/**
\class ClassATubePre
\ingroup FX-Objects
\brief
The ClassATubePre object implements a simple cascade of four (4) triode tube models.
Audio I/O:
- Processes mono input to mono output.
Control I/F:
- Use ClassATubePreParameters structure to get/set object params.
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
class ClassATubePre : public IAudioSignalProcessor
{
public:
	ClassATubePre() {}		/* C-TOR */
	~ClassATubePre() {}		/* D-TOR */

	/** reset members to initialized state */
	virtual bool reset(double _sampleRate)
	{
		TriodeClassAParameters tubeParams = triodes[0].getParameters();
		tubeParams.invertOutput = true;
		tubeParams.enableHPF = true; // remove DC offsets
		tubeParams.outputGain = 1.0;
		tubeParams.saturation = 1.0;
		tubeParams.asymmetry = 0.0;
		tubeParams.enableLSF = true;
		tubeParams.lsf_Fshelf = 88.0;
		tubeParams.lsf_BoostCut_dB = -12.0;
		tubeParams.waveshaper = distortionModel::kFuzzAsym;

		for (int i = 0; i < NUM_TUBES; i++)
		{
			triodes[i].reset(_sampleRate);
			triodes[i].setParameters(tubeParams);
		}

		shelvingFilter.reset(_sampleRate);

		return true;
	}

	/** get parameters: note use of custom structure for passing param data */
	/**
	\return ClassATubePreParameters custom data structure
	*/
	ClassATubePreParameters getParameters() { return parameters; }

	/** set parameters: note use of custom structure for passing param data */
	/**
	\param ClassATubePreParameters custom data structure
	*/
	void setParameters(const ClassATubePreParameters& params)
	{
		// --- check for re-calc
		if (params.inputLevel_dB != parameters.inputLevel_dB)
			inputLevel = pow(10.0, params.inputLevel_dB / 20.0);
		if (params.outputLevel_dB != parameters.outputLevel_dB)
			outputLevel = pow(10.0, params.outputLevel_dB / 20.0);

		// --- store
		parameters = params;

		// --- shelving filter update
		TwoBandShelvingFilterParameters sfParams = shelvingFilter.getParameters();
		sfParams.lowShelf_fc = parameters.lowShelf_fc;
		sfParams.lowShelfBoostCut_dB = parameters.lowShelfBoostCut_dB;
		sfParams.highShelf_fc = parameters.highShelf_fc;
		sfParams.highShelfBoostCut_dB = parameters.highShelfBoostCut_dB;
		shelvingFilter.setParameters(sfParams);

		// --- triode updates
		TriodeClassAParameters tubeParams = triodes[0].getParameters();
		tubeParams.saturation = parameters.saturation;
		tubeParams.asymmetry = parameters.asymmetry;

		for (int i = 0; i < NUM_TUBES; i++)
			triodes[i].setParameters(tubeParams);
	}

	/** return false: this object only processes samples */
	virtual bool canProcessAudioFrame() { return false; }

	/** process the input through the four tube models in series */
	/**
	\param xn input
	\return the processed sample
	*/
	virtual double processAudioSample(double xn)
	{
		double output1 = triodes[0].processAudioSample(xn*inputLevel);
		double output2 = triodes[1].processAudioSample(output1);
		double output3 = triodes[2].processAudioSample(output2);

		// --- filter stage is between 3 and 4
		double outputEQ = shelvingFilter.processAudioSample(output3);
		double output4 = triodes[3].processAudioSample(outputEQ);

		return output4*outputLevel;
	}

protected:
	ClassATubePreParameters parameters;		///< object parameters
	TriodeClassA triodes[NUM_TUBES];		///< array of triode tube objects
	TwoBandShelvingFilter shelvingFilter;	///< shelving filters

	double inputLevel = 1.0;	///< input level (not in dB)
	double outputLevel = 1.0;	///< output level (not in dB)
};


/**
\struct BitCrusherParameters
\ingroup FX-Objects
\brief
Custom parameter structure for the BitCrusher object.
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
struct BitCrusherParameters
{
	BitCrusherParameters() {}
	/** all FXObjects parameter objects require overloaded= operator so remember to add new entries if you add new variables. */
	BitCrusherParameters& operator=(const BitCrusherParameters& params)	// need this override for collections to work
	{
		if (this == &params)
			return *this;

		quantizedBitDepth = params.quantizedBitDepth;

		return *this;
	}

	double quantizedBitDepth = 4.0; ///< bid depth of quantizer
};

/**
\class BitCrusher
\ingroup FX-Objects
\brief
The BitCrusher object implements a quantizing bitcrusher algorithm.
Audio I/O:
- Processes mono input to mono output.
Control I/F:
- Use BitCrusherParameters structure to get/set object params.
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
class BitCrusher : public IAudioSignalProcessor
{
public:
	BitCrusher() {}		/* C-TOR */
	~BitCrusher() {}		/* D-TOR */

	/** reset members to initialized state */
	virtual bool reset(double _sampleRate){ return true; }

	/** get parameters: note use of custom structure for passing param data */
	/**
	\return BitCrusherParameters custom data structure
	*/
	BitCrusherParameters getParameters() { return parameters; }

	/** set parameters: note use of custom structure for passing param data */
	/**
	\param BitCrusherParameters custom data structure
	*/
	void setParameters(const BitCrusherParameters& params)
	{
		// --- calculate and store
		if (params.quantizedBitDepth != parameters.quantizedBitDepth)
			QL = 2.0 / (pow(2.0, params.quantizedBitDepth) - 1.0);

		parameters = params;
	}

	/** return false: this object only processes samples */
	virtual bool canProcessAudioFrame() { return false; }

	/** perform the bitcrushing operation (see FX book for back story and details) */
	/**
	\param xn input
	\return the processed sample
	*/
	virtual double processAudioSample(double xn)
	{
		return QL*(int(xn / QL));
	}

protected:
	BitCrusherParameters parameters; ///< object parameters
	double QL = 1.0;				 ///< the quantization level
};
