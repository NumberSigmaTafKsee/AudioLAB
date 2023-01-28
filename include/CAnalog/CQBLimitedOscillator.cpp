#include "CQBLimitedOscillator.h"

CQBLimitedOscillator::CQBLimitedOscillator(void)
{
}

CQBLimitedOscillator::~CQBLimitedOscillator(void)
{
}

void CQBLimitedOscillator::reset()
{
	COscillator::reset();

	if(m_uWaveform == SAW1 || m_uWaveform == SAW2 ||
	   m_uWaveform == SAW3 || m_uWaveform == TRI)
	{
		m_dModulo = 0.5;
	}
}

void CQBLimitedOscillator::startOscillator()
{
	reset();
	m_bNoteOn = true;
}

void CQBLimitedOscillator::stopOscillator()
{
	m_bNoteOn = false;
}