#pragma once

#include "CFilter.h"
//#include "fxobjects.h"

class CVAOnePoleFilter : public CFilter
{
public:
	CVAOnePoleFilter(void);
	~CVAOnePoleFilter(void);

	// Trapezoidal Integrator Components
	DspFloatType m_dAlpha;			// Feed Forward coeff

	// -- ADDED for Korg35 and Moog Ladder Filter ---- //
	DspFloatType m_dBeta;

	// -- ADDED for Diode Ladder Filter  ---- //
	DspFloatType m_dGamma;		// Pre-Gain
	DspFloatType m_dDelta;		// FB_IN Coeff
	DspFloatType m_dEpsilon;		// FB_OUT scalar
	DspFloatType m_da0;			// input gain

	// note: this is NOT being used as a z-1 storage register!
	DspFloatType m_dFeedback;		// our own feedback coeff from S 

	// provide access to set feedback input
	void setFeedback(DspFloatType fb){m_dFeedback = fb;}

	// provide access to our feedback output
	// m_dFeedback & m_dDelta = 0 for non-Diode filters
	DspFloatType getFeedbackOutput(){return m_dBeta*(m_dZ1 + m_dFeedback*m_dDelta);}

	// original function
	// DspFloatType getFeedbackOutput(){return m_dZ1*m_dBeta;}
	// ----------------------------------------------- //

	// -- CFilter Overrides ---
	virtual void reset(){m_dZ1 = 0; m_dFeedback = 0;}

	// recalc the coeff
	virtual void update();
	
	// do the filter
	virtual DspFloatType doFilter(DspFloatType xn);

protected:
	DspFloatType m_dZ1;		// our z-1 storage location
};
