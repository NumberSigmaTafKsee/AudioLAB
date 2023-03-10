
#include "CFilter.h"

// --- construction
CFilter::CFilter(void)
{
	// --- defaults
	m_dSampleRate = 44100;
	m_dQControl = 1.0; // Q is 1 to 10 on GUI
	m_dFc = FILTER_FC_DEFAULT;
	m_dQ = FILTER_Q_DEFAULT;
	m_dFcControl = FILTER_FC_DEFAULT;

	// --- clear
	m_dFcMod = 0.0;
	m_dAuxControl = 0.0; 
	m_uNLP = OFF;
	m_dSaturation = 1.0;
}

CFilter::~CFilter(void)
{
}


// --- flush buffers
void CFilter::reset()
{
	// do nothing
}

// --- optional depending on filter type
void CFilter::setQControl(DspFloatType dQControl)
{
	// do nothing
}