#pragma once
#include "CFilter.h"
#include "CVAOnePoleFilter.h"

class CKThreeFiveFilter : public CFilter
{
public:
	CKThreeFiveFilter(void);
	~CKThreeFiveFilter(void);

	// our member filters
	// LPF: LPF1+LPF2+HPF1
	// HPF: HPF1+LPF1+HPF2
	CVAOnePoleFilter m_LPF1;
	CVAOnePoleFilter m_LPF2;
	CVAOnePoleFilter m_HPF1;
	CVAOnePoleFilter m_HPF2;

	// -- CFilter Overrides --
	virtual void reset();
	virtual void setQControl(DspFloatType dQControl);
	virtual void update();
	virtual DspFloatType doFilter(DspFloatType xn);

	// variables
	DspFloatType m_dAlpha0;   // our u scalar value
	DspFloatType m_dK;		// K, set with Q
};