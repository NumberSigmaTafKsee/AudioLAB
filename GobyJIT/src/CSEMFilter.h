#pragma once
#include "CFilter.h"

class CSEMFilter : public CFilter
{
public:
	CSEMFilter(void);
	~CSEMFilter(void);

	// --- Trapezoidal Integrator Components
	double m_dAlpha0;		// input scalar
	double m_dAlpha;		// alpha is same as VA One Pole
	double m_dRho;			// feedback

	// -- CFilter Overrides --
	virtual void reset(){m_dZ11 = 0; m_dZ12 = 0.0;}
	virtual void setQControl(double dQControl);
	virtual void update();
	virtual double doFilter(double xn);

protected:
	double m_dZ11;		// our z-1 storage location
	double m_dZ12;		// our z-1 storage location # 2
};

