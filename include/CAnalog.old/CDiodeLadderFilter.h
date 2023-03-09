#pragma once
#include "CFilter.h"
#include "synthlite.h"
#include "CVAOnePoleFilter.h"

#include <cstdint>
#include <cmath>

class CDiodeLadderFilter : public CFilter
{
public:
	CDiodeLadderFilter(void);
	~CDiodeLadderFilter(void);

	CVAOnePoleFilter m_LPF1;
	CVAOnePoleFilter m_LPF2;
	CVAOnePoleFilter m_LPF3;
	CVAOnePoleFilter m_LPF4;

	// variables
	DspFloatType m_dK;			// K, set with Q
	DspFloatType m_dGamma;		// needed for final calc and update
	DspFloatType m_dSG1; 
	DspFloatType m_dSG2; 
	DspFloatType m_dSG3; 
	DspFloatType m_dSG4; 

	// -- CFilter Overrides --
	virtual void reset();
	virtual void setQControl(DspFloatType dQControl);
	virtual void update();

	inline virtual DspFloatType doFilter(DspFloatType xn)
	{
		// --- return xn if filter not supported
		if(m_uFilterType != LPF4)
			return xn;

		m_LPF4.setFeedback(0.0);
		m_LPF3.setFeedback(m_LPF4.getFeedbackOutput());
		m_LPF2.setFeedback(m_LPF3.getFeedbackOutput());
		m_LPF1.setFeedback(m_LPF2.getFeedbackOutput());

		// --- form input
		DspFloatType dSigma = m_dSG1*m_LPF1.getFeedbackOutput() + 
						m_dSG2* m_LPF2.getFeedbackOutput() +
						m_dSG3*m_LPF3.getFeedbackOutput() +
						m_dSG4* m_LPF4.getFeedbackOutput();

		// --- for passband gain compensation!
		xn *= 1.0 + m_dAuxControl*m_dK;

		// --- form input
		DspFloatType dU = (xn - m_dK*dSigma)/(1 + m_dK*m_dGamma);

		// ---NLP
		if(m_uNLP == ON)
			dU = fasttanh(m_dSaturation*dU);

		// --- cascade of four filters
		return m_LPF4.doFilter(m_LPF3.doFilter(m_LPF2.doFilter(m_LPF1.doFilter(dU))));
	}
};