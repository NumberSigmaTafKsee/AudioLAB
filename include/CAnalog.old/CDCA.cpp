#include "CDCA.h"

namespace Analog::CAnalog
{
	CDCA::CDCA(void) : StereoAmplifierProcessor()
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
}