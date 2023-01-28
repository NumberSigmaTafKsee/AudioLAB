#pragma once

#include "CFilter.hpp"
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

CVAOnePoleFilter::CVAOnePoleFilter(void)
{
    // --- init defaults to simple 
    //	   LPF/HPF structure
    m_dAlpha = 1.0;
    m_dBeta = 0.0;
    m_dZ1 = 0.0;
    m_dGamma = 1.0;
    m_dDelta = 0.0;
    m_dEpsilon = 0.0;
    m_da0 = 1.0;
    m_dFeedback = 0.0;

    // --- always set the default!
    m_uFilterType = LPF1;

    // --- flush storage
    reset();
}

CVAOnePoleFilter::~CVAOnePoleFilter(void)
{
}

// recalc coeffs
void CVAOnePoleFilter::update()
{
    // base class does modulation, changes m_fFc
    CFilter::update();

    DspFloatType wd = 2*M_PI*m_dFc;          
    DspFloatType T  = 1/m_dSampleRate;             
    DspFloatType wa = (2/T)*tan(wd*T/2); 
    DspFloatType g  = wa*T/2;            

    m_dAlpha = g/(1.0 + g);
}

// do the filter
DspFloatType CVAOnePoleFilter::doFilter(DspFloatType xn)
{
    // return xn if filter not supported
    if(m_uFilterType != LPF1 && m_uFilterType != HPF1)
        return xn;

    // for diode filter support
    xn = xn*m_dGamma + m_dFeedback + m_dEpsilon*getFeedbackOutput();
    
    // calculate v(n)
    DspFloatType vn = (m_da0*xn - m_dZ1)*m_dAlpha;

    // form LP output
    DspFloatType lpf = vn + m_dZ1;

    // update memory
    m_dZ1 = vn + lpf;

    // do the HPF
    DspFloatType hpf = xn - lpf;

    if(m_uFilterType == LPF1)
        return lpf;
    else if(m_uFilterType == HPF1)
        return hpf;

    return xn; // should never get here
}
