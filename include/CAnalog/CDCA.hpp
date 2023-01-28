#pragma once
#include "synthlite.h"
#include <cstdint>

#include "SoundObject.hpp"
#define AMP_MOD_RANGE -96	// -96dB

struct CDCA : public StereoFXProcessor
{	
	CDCA(void);
	~CDCA(void);

	// --- ATTRIBUTES
	// --- PUBLIC: these variables may be get/set

	DspFloatType m_dGain;

	// --- velocity input from MIDI keyboard
	uint32_t m_uMIDIVelocity; // 0 -> 127

	// --- controls for user GUI (optional)
	DspFloatType m_dAmplitude_dB;		// the user's control setting in dB
	DspFloatType m_dAmplitudeControl;	// the user's control setting, converted from dB

	// --- pan control
	DspFloatType m_dPanControl;	/* -1 to +1 == left to right */

	// --- modulate amplitude
	DspFloatType m_dAmpMod_dB;	/* for tremolo, not true AM */

	// --- input to EGMod is EXPONENTIAL
	DspFloatType m_dEGMod;		 /* modulation input for EG 0 to +1 */

	// --- input to modulate pan control is bipolar
	DspFloatType m_dPanMod;		/* modulation input for -1 to +1 */

	// --- FUNCTIONS: all public
	//
	// --- MIDI controller functions
	inline void setMIDIVelocity(uint32_t u){m_uMIDIVelocity = u;}

	// --- Pan control
	inline void setPanControl(DspFloatType d){m_dPanControl = d;}

	// --- reset mods
	inline void reset()
	{
		m_dEGMod = 1.0;
		m_dAmpMod_dB = 0.0;
	}

	// --- NOTE: -96dB to +24dB
	inline void setAmplitude_dB(DspFloatType d)
	{
		m_dAmplitude_dB = d;
		m_dAmplitudeControl = std::pow((DspFloatType)10.0, m_dAmplitude_dB/(DspFloatType)20.0);
	}

	// --- modulation functions - NOT needed/used if you implement the Modulation Matrix!
	//
	// --- expecting connection from bipolar source (LFO)
	//	   but this component will only be able to attenuate
	//	   so convert to unipolar
	inline void setAmpMod_dB(DspFloatType d){m_dAmpMod_dB = bipolarToUnipolar(d);}

	// --- EG Mod Input Functions
	inline void setEGMod(DspFloatType d){m_dEGMod = d;}

	// --- Pan modulation
	inline void setPanMod(DspFloatType d){m_dPanMod = d;}


	// --- DCA operation functions
	// --- recalculate gain values
	inline void update()
	{
		// --- check polarity
		if(m_dEGMod >= 0)
			m_dGain = m_dEGMod;
		else
			m_dGain = m_dEGMod + 1.0;

		// --- amp mod is attenuation only, in dB
		m_dGain *= std::pow(10.0, m_dAmpMod_dB/(DspFloatType)20.0);

		// --- use MMA MIDI->Atten (convex) transform
		m_dGain *= mmaMIDItoAtten(m_uMIDIVelocity);
	}

	// --- do the DCA: uses pass-by-reference for outputs
	//     For mono-in, just repeat the inputs
	inline void doDCA(DspFloatType dLeftInput, DspFloatType dRightInput, DspFloatType& dLeftOutput, DspFloatType& dRightOutput)
	{
		// total pan value
		DspFloatType dPanTotal = m_dPanControl + m_dPanMod;

		// limit in case pan control is biased
		dPanTotal = fmin(dPanTotal, (DspFloatType)1.0);
		dPanTotal = fmax(dPanTotal, (DspFloatType)-1.0);

		DspFloatType dPanLeft = 0.707;
		DspFloatType dPanRight = 0.707;

		// equal std::power calculation in synthfunction.h
		calculatePanValues(dPanTotal, dPanLeft, dPanRight);

		// form left and right outputs
		dLeftOutput =  dPanLeft*m_dAmplitudeControl*dLeftInput*m_dGain;
		dRightOutput =  dPanRight*m_dAmplitudeControl*dRightInput*m_dGain;
	}

	DspFloatType Tick(DspFloatType iL, DspFloatType iR, DspFloatType &oL, DspFloatType & oR, DspFloatType A=1, DspFloatType X=0, DspFloatType Y=0)
	{            
		doDCA(iL,iR,oL,oR);
		return (oL+oR)*0.5;
	}
	void ProcessBlock(size_t n, DspFloatType ** in, DspFloatType ** out) {
		for(size_t i = 0; i < n; i++) doDCA(in[0][i],in[1][i],out[0][i],out[1][i]);
	}
};

CDCA::CDCA(void) : StereoFXProcessor()
{
	// --- initialize variables
	m_dAmplitudeControl = 1.0;
	m_dAmpMod_dB = 0.0;
	m_dGain = 1.0;
	m_dAmplitude_dB = 0.0;
	m_dEGMod = 1.0;
	m_dPanControl = 0.0;
	m_dPanMod = 0.0;
	m_uMIDIVelocity = 127;
}

// --- destruction
CDCA::~CDCA(void)
{
}
