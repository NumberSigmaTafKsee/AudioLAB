#pragma once
#include "canalog_filter.hpp"

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

CSEMFilter::CSEMFilter(void)
{
    // --- init
    m_dAlpha0 = 1.0;
    m_dAlpha = 1.0;
    m_dRho = 1.0;
    m_dAuxControl = 0.5; // for BSF this centers it

    //--- our default filter type
    m_uFilterType = LPF2;

    // --- flush registers
    reset();
}

CSEMFilter::~CSEMFilter(void)
{
}

// decode the Q value; Q on UI is 1->10
void CSEMFilter::setQControl(DspFloatType dQControl)
{
    // this maps dQControl = 1->10 to Q = 0.5->25
    m_dQ = (25.0 - 0.5)*(dQControl - 1.0)/(10.0 - 1.0) + 0.5;
}

// recalc the coeffs
void CSEMFilter::update()
{	
    // base class does modulation
    CFilter::update();

    // prewarp the cutoff- these are bilinear-transform filters
    DspFloatType wd = 2*M_PI*m_dFc;          
    DspFloatType T  = 1/m_dSampleRate;             
    DspFloatType wa = (2/T)*tan(wd*T/2); 
    DspFloatType g  = wa*T/2;            

    // note R is the traditional analog damping factor
    DspFloatType R = 1.0/(2.0*m_dQ);

    // set the coeffs
    m_dAlpha0 = 1.0/(1.0 + 2.0*R*g + g*g);
    m_dAlpha = g;
    m_dRho = 2.0*R + g;
}

// do the filter
DspFloatType CSEMFilter::doFilter(DspFloatType xn)
{
    // return xn if filter not supported
    if(m_uFilterType != LPF2 && m_uFilterType != HPF2 && 
    m_uFilterType != BPF2 && m_uFilterType != BSF2)
        return xn;

    // form the HP output first 
    DspFloatType hpf = m_dAlpha0*(xn - m_dRho*m_dZ11 - m_dZ12);

    // BPF Out
    DspFloatType bpf = m_dAlpha*hpf + m_dZ11;

    // for nonlinear proc
    if(m_uNLP == ON)
        bpf = fasttanh(m_dSaturation*bpf);

    // LPF Out
    DspFloatType lpf = m_dAlpha*bpf + m_dZ12;

    // note R is the traditional analog damping factor
    DspFloatType R = 1.0/(2.0*m_dQ);

    // BSF Out
    DspFloatType bsf = xn - 2.0*R*bpf;

    // SEM BPF Output
    // using m_dAuxControl for this one-off control
    DspFloatType semBSF = m_dAuxControl*hpf + (1.0 - m_dAuxControl)*lpf;

    // update memory
    m_dZ11 = m_dAlpha*hpf + bpf;
    m_dZ12 = m_dAlpha*bpf + lpf;

    // return our selected type
    if(m_uFilterType == LPF2)
        return lpf;
    else if(m_uFilterType == HPF2)
        return hpf;
    else if(m_uFilterType == BPF2)
        return bpf; 
    else if(m_uFilterType == BSF2)
        return semBSF;

    // return input if filter not supported
    return xn;
}
