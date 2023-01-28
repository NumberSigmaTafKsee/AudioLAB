#pragma once
#include "CFilter.h"

class CSEMFilter : public CFilter
{
public:
	CSEMFilter(void);
	~CSEMFilter(void);

	// --- Trapezoidal Integrator Components
	DspFloatType m_dAlpha0;		// input scalar
	DspFloatType m_dAlpha;		// alpha is same as VA One Pole
	DspFloatType m_dRho;			// feedback

	// -- CFilter Overrides --
	virtual void reset(){m_dZ11 = 0; m_dZ12 = 0.0;}
	virtual void setQControl(DspFloatType dQControl);
	virtual void update();
	virtual DspFloatType doFilter(DspFloatType xn);

protected:
	DspFloatType m_dZ11;		// our z-1 storage location
	DspFloatType m_dZ12;		// our z-1 storage location # 2
};

