#pragma once

#include "FXObjects.hpp"
#include "FXDelays.hpp"
#include "FXFilters.hpp"

// --- constants for reverb tank
const unsigned int NUM_BRANCHES = 4;
const unsigned int NUM_CHANNELS = 2; // stereo

/**
\enum reverbDensity
\ingroup Constants-Enums
\brief
Use this strongly typed enum to easily set the density in the reverb object.
- enum class reverbDensity { kThick, kSparse };
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
enum class reverbDensity { kThick, kSparse };

/**
\struct ReverbTankParameters
\ingroup FX-Objects
\brief
Custom parameter structure for the ReverbTank object.
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
struct ReverbTankParameters
{
	ReverbTankParameters() {}
	/** all FXObjects parameter objects require overloaded= operator so remember to add new entries if you add new variables. */
	ReverbTankParameters& operator=(const ReverbTankParameters& params)	// need this override for collections to work
	{
		if (this == &params)
			return *this;

		density = params.density;

		// --- tweaker variables
		apfDelayMax_mSec = params.apfDelayMax_mSec;
		apfDelayWeight_Pct = params.apfDelayWeight_Pct;
		fixeDelayMax_mSec = params.fixeDelayMax_mSec;
		fixeDelayWeight_Pct = params.fixeDelayWeight_Pct;
		preDelayTime_mSec = params.preDelayTime_mSec;

		lpf_g = params.lpf_g;
		kRT = params.kRT;

		lowShelf_fc = params.lowShelf_fc;
		lowShelfBoostCut_dB = params.lowShelfBoostCut_dB;
		highShelf_fc = params.highShelf_fc;
		highShelfBoostCut_dB = params.highShelfBoostCut_dB;

		wetLevel_dB = params.wetLevel_dB;
		dryLevel_dB = params.dryLevel_dB;
		return *this;
	}

	// --- individual parameters
	reverbDensity density = reverbDensity::kThick;	///< density setting thick or thin

	// --- tweaking parameters - you may not want to expose these
	//     in the final plugin!
	// --- See the book for all the details on how these tweakers work!!
	DspFloatType apfDelayMax_mSec = 5.0;					///< APF max delay time
	DspFloatType apfDelayWeight_Pct = 100.0;				///< APF max delay weighying
	DspFloatType fixeDelayMax_mSec = 50.0;				///< fixed delay max time
	DspFloatType fixeDelayWeight_Pct = 100.0;				///< fixed delay max weighying

	// --- direct control parameters
	DspFloatType preDelayTime_mSec = 0.0;					///< pre-delay time in mSec
	DspFloatType lpf_g = 0.0;								///< LPF g coefficient
	DspFloatType kRT = 0.0;								///< reverb time, 0 to 1

	DspFloatType lowShelf_fc = 0.0;						///< low shelf fc
	DspFloatType lowShelfBoostCut_dB = 0.0;				///< low shelf gain
	DspFloatType highShelf_fc = 0.0;						///< high shelf fc
	DspFloatType highShelfBoostCut_dB = 0.0;				///< high shelf gain

	DspFloatType wetLevel_dB = -3.0;						///< wet output level in dB
	DspFloatType dryLevel_dB = -3.0;						///< dry output level in dB
};


/**
\class ReverbTank
\ingroup FX-Objects
\brief
The ReverbTank object implements the cyclic reverb tank in the FX book listed below.
Audio I/O:
- Processes mono input to mono OR stereo output.
Control I/F:
- Use ReverbTankParameters structure to get/set object params.
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
class ReverbTank : public IAudioSignalProcessor
{
public:
	ReverbTank() {}		/* C-TOR */
	~ReverbTank() {}	/* D-TOR */

	/** reset members to initialized state */
	virtual bool reset(DspFloatType _sampleRate)
	{
		// ---store
		sampleRate = _sampleRate;

		// ---set up preDelay
		preDelay.reset(_sampleRate);
		preDelay.createDelayBuffer(_sampleRate, 100.0);

		for (int i = 0; i < NUM_BRANCHES; i++)
		{
			branchDelays[i].reset(_sampleRate);
			branchDelays[i].createDelayBuffer(_sampleRate, 100.0);

			branchNestedAPFs[i].reset(_sampleRate);
			branchNestedAPFs[i].createDelayBuffers(_sampleRate, 100.0, 100.0);

			branchLPFs[i].reset(_sampleRate);
		}
		for (int i = 0; i < NUM_CHANNELS; i++)
		{
			shelvingFilters[i].reset(_sampleRate);
		}

		return true;
	}

	/** return true: this object can process frames */
	virtual bool canProcessAudioFrame() { return true; }

	/** process mono reverb tank */
	/**
	\param xn input
	\return the processed sample
	*/
	virtual DspFloatType processAudioSample(DspFloatType xn)
	{
		DspFloatType inputs[2] = { 0.0, 0.0 };
		DspFloatType outputs[2] = { 0.0, 0.0 };
		processAudioFrame(inputs, outputs, 1, 1);
		return outputs[0];
	}

	/** process stereo reverb tank */
	virtual bool processAudioFrame(const DspFloatType* inputFrame,
		DspFloatType* outputFrame,
		uint32_t inputChannels,
		uint32_t outputChannels)
	{
		// --- global feedback from delay in last branch
		DspFloatType globFB = branchDelays[NUM_BRANCHES-1].readDelay();

		// --- feedback value
		DspFloatType fb = parameters.kRT*(globFB);

		// --- mono-ized input signal
		DspFloatType xnL = inputFrame[0];
		DspFloatType xnR = inputChannels > 1 ? inputFrame[1] : 0.0;
		DspFloatType monoXn = DspFloatType(1.0 / inputChannels)*xnL + DspFloatType(1.0 / inputChannels)*xnR;

		// --- pre delay output
		DspFloatType preDelayOut = preDelay.processAudioSample(monoXn);

		// --- input to first branch = preDalay + globFB
		DspFloatType input = preDelayOut + fb;
		for (int i = 0; i < NUM_BRANCHES; i++)
		{
			DspFloatType apfOut = branchNestedAPFs[i].processAudioSample(input);
			DspFloatType lpfOut = branchLPFs[i].processAudioSample(apfOut);
			DspFloatType delayOut = parameters.kRT*branchDelays[i].processAudioSample(lpfOut);
			input = delayOut + preDelayOut;
		}
		// --- gather outputs
		/*
		There are 25 prime numbers between 1 and 100.
		They are 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41,
		43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, and 97
		we want 16 of them: 23, 29, 31, 37, 41,
		43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, and 97
		*/

		DspFloatType weight = 0.707;

		DspFloatType outL= 0.0;
		outL += weight*branchDelays[0].readDelayAtPercentage(23.0);
		outL -= weight*branchDelays[1].readDelayAtPercentage(41.0);
		outL += weight*branchDelays[2].readDelayAtPercentage(59.0);
		outL -= weight*branchDelays[3].readDelayAtPercentage(73.0);

		DspFloatType outR = 0.0;
		outR -= weight*branchDelays[0].readDelayAtPercentage(29.0);
		outR += weight*branchDelays[1].readDelayAtPercentage(43.0);
		outR -= weight*branchDelays[2].readDelayAtPercentage(61.0);
		outR += weight*branchDelays[3].readDelayAtPercentage(79.0);

		if (parameters.density == reverbDensity::kThick)
		{
			outL += weight*branchDelays[0].readDelayAtPercentage(31.0);
			outL -= weight*branchDelays[1].readDelayAtPercentage(47.0);
			outL += weight*branchDelays[2].readDelayAtPercentage(67.0);
			outL -= weight*branchDelays[3].readDelayAtPercentage(83.0);

			outR -= weight*branchDelays[0].readDelayAtPercentage(37.0);
			outR += weight*branchDelays[1].readDelayAtPercentage(53.0);
			outR -= weight*branchDelays[2].readDelayAtPercentage(71.0);
			outR += weight*branchDelays[3].readDelayAtPercentage(89.0);
		}

		// ---  filter
		DspFloatType tankOutL = shelvingFilters[0].processAudioSample(outL);
		DspFloatType tankOutR = shelvingFilters[1].processAudioSample(outR);

		// --- sum with dry
		DspFloatType dry = pow(10.0, parameters.dryLevel_dB / 20.0);
		DspFloatType wet = pow(10.0, parameters.wetLevel_dB / 20.0);

		if (outputChannels == 1)
			outputFrame[0] = dry*xnL + wet*(0.5*tankOutL + 0.5*tankOutR);
		else
		{
			outputFrame[0] = dry*xnL + wet*tankOutL;
			outputFrame[1] = dry*xnR + wet*tankOutR;
		}

		return true;
	}

	/** get parameters: note use of custom structure for passing param data */
	/**
	\return ReverbTankParameters custom data structure
	*/
	ReverbTankParameters getParameters() { return parameters; }

	/** set parameters: note use of custom structure for passing param data */
	/**
	\param ReverbTankParameters custom data structure
	*/
	void setParameters(const ReverbTankParameters& params)
	{
		// --- do the updates here, the sub-components will only update themselves if
		//     their parameters changed, so we let those object handle that chore
		TwoBandShelvingFilterParameters filterParams = shelvingFilters[0].getParameters();
		filterParams.highShelf_fc = params.highShelf_fc;
		filterParams.highShelfBoostCut_dB = params.highShelfBoostCut_dB;
		filterParams.lowShelf_fc = params.lowShelf_fc;
		filterParams.lowShelfBoostCut_dB = params.lowShelfBoostCut_dB;

		// --- copy to both channels
		shelvingFilters[0].setParameters(filterParams);
		shelvingFilters[1].setParameters(filterParams);

		SimpleLPFParameters  lpfParams = branchLPFs[0].getParameters();
		lpfParams.g = params.lpf_g;

		for (int i = 0; i < NUM_BRANCHES; i++)
		{
			branchLPFs[i].setParameters(lpfParams);
		}

		// --- update pre delay
		SimpleDelayParameters delayParams = preDelay.getParameters();
		delayParams.delayTime_mSec = params.preDelayTime_mSec;
		preDelay.setParameters(delayParams);

		// --- set apf and delay parameters
		int m = 0;
		NestedDelayAPFParameters apfParams = branchNestedAPFs[0].getParameters();
		delayParams = branchDelays[0].getParameters();

		// --- global max Delay times
		DspFloatType globalAPFMaxDelay = (parameters.apfDelayWeight_Pct / 100.0)*parameters.apfDelayMax_mSec;
		DspFloatType globalFixedMaxDelay = (parameters.fixeDelayWeight_Pct / 100.0)*parameters.fixeDelayMax_mSec;

		// --- lfo
		apfParams.enableLFO = true;
		apfParams.lfoMaxModulation_mSec = 0.3;
		apfParams.lfoDepth = 1.0;

		for (int i = 0; i < NUM_BRANCHES; i++)
		{
			// --- setup APFs
			apfParams.outerAPFdelayTime_mSec = globalAPFMaxDelay*apfDelayWeight[m++];
			apfParams.innerAPFdelayTime_mSec = globalAPFMaxDelay*apfDelayWeight[m++];
			apfParams.innerAPF_g = -0.5;
			apfParams.outerAPF_g = 0.5;
			if (i == 0)
				apfParams.lfoRate_Hz = 0.15;
			else if (i == 1)
				apfParams.lfoRate_Hz = 0.33;
			else if (i == 2)
				apfParams.lfoRate_Hz = 0.57;
			else if (i == 3)
				apfParams.lfoRate_Hz = 0.73;

			branchNestedAPFs[i].setParameters(apfParams);

			// --- fixedDelayWeight
			delayParams.delayTime_mSec = globalFixedMaxDelay*fixedDelayWeight[i];
			branchDelays[i].setParameters(delayParams);
		}

		// --- save our copy
		parameters = params;
	}


private:
	ReverbTankParameters parameters;				///< object parameters

	SimpleDelay  preDelay;							///< pre delay object
	SimpleDelay  branchDelays[NUM_BRANCHES];		///< branch delay objects
	NestedDelayAPF branchNestedAPFs[NUM_BRANCHES];	///< nested APFs for each branch
	SimpleLPF  branchLPFs[NUM_BRANCHES];			///< LPFs in each branch

	TwoBandShelvingFilter shelvingFilters[NUM_CHANNELS]; ///< shelving filters 0 = left; 1 = right

	// --- weighting values to make various and low-correlated APF delay values easily
	DspFloatType apfDelayWeight[NUM_BRANCHES * 2] = { 0.317, 0.873, 0.477, 0.291, 0.993, 0.757, 0.179, 0.575 };///< weighting values to make various and low-correlated APF delay values easily
	DspFloatType fixedDelayWeight[NUM_BRANCHES] = { 1.0, 0.873, 0.707, 0.667 };	///< weighting values to make various and fixed delay values easily
	DspFloatType sampleRate = 0.0;	///< current sample rate
};
