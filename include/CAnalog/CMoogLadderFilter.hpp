#pragma once
#include "CFilter.hpp"
#include "CVAOnePoleFilter.hpp"

class CMoogLadderFilter : public CFilter
{
public:
    CMoogLadderFilter(void);
    ~CMoogLadderFilter(void);

    CVAOnePoleFilter m_LPF1;
    CVAOnePoleFilter m_LPF2;
    CVAOnePoleFilter m_LPF3;
    CVAOnePoleFilter m_LPF4;

    // -- CFilter Overrides --
    virtual void reset();
    virtual void setQControl(DspFloatType dQControl);
    //virtual void update();
    //virtual DspFloatType doFilter(DspFloatType xn);

    // variables
    DspFloatType m_dK;		// K, set with Q
    DspFloatType m_dGamma;	// see block diagram
    DspFloatType m_dAlpha_0;	// see block diagram

    // Oberheim Xpander variations
    DspFloatType m_dA;
    DspFloatType m_dB;
    DspFloatType m_dC;
    DspFloatType m_dD;
    DspFloatType m_dE;

    inline virtual void update()
    {
        // do any modulation first
        CFilter::update();

        // prewarp for BZT
        DspFloatType wd = 2*M_PI*m_dFc;
        DspFloatType T  = 1/m_dSampleRate;
        DspFloatType wa = (2/T)*tan(wd*T/2);
        DspFloatType g  = wa*T/2;

        // G - the feedforward coeff in the VA One Pole
        //     same for LPF, HPF
        DspFloatType G = g/(1.0 + g);

        // set alphas
        m_LPF1.m_dAlpha = G;
        m_LPF2.m_dAlpha = G;
        m_LPF3.m_dAlpha = G;
        m_LPF4.m_dAlpha = G;

        // set betas
        m_LPF1.m_dBeta = G*G*G/(1.0 + g);
        m_LPF2.m_dBeta = G*G/(1.0 + g);
        m_LPF3.m_dBeta = G/(1.0 + g);
        m_LPF4.m_dBeta = 1.0/(1.0 + g);

        m_dGamma = G*G*G*G; // G^4

        m_dAlpha_0 = 1.0/(1.0 + m_dK*m_dGamma);

        // Oberhein variation
        switch(m_uFilterType)
        {
            case LPF4:
                m_dA = 0.0; m_dB = 0.0; m_dC = 0.0; m_dD = 0.0; m_dE = 1.0;
                break;

            case LPF2:
                m_dA = 0.0; m_dB = 0.0; m_dC = 1.0; m_dD = 0.0; m_dE = 0.0;
                break;

            case BPF4:
                m_dA = 0.0; m_dB = 0.0; m_dC = 4.0; m_dD = -8.0; m_dE = 4.0;
                break;

            case BPF2:
                m_dA = 0.0; m_dB = 2.0; m_dC = -2.0; m_dD = 0.0; m_dE = 0.0;
                break;

            case HPF4:
                m_dA = 1.0; m_dB = -4.0; m_dC = 6.0; m_dD = -4.0; m_dE = 1.0;
                break;

            case HPF2:
                m_dA = 1.0; m_dB = -2.0; m_dC = 1.0; m_dD = 0.0; m_dE = 0.0;
                break;

            default: // LPF4
                m_dA = 0.0; m_dB = 0.0; m_dC = 0.0; m_dD = 0.0; m_dE = 1.0;
                break;
        }
    }

    inline virtual DspFloatType doFilter(DspFloatType xn)
    {
        // --- return xn if filter not supported
        if(m_uFilterType == BSF2 || m_uFilterType == LPF1 || m_uFilterType == HPF1)
            return xn;

        DspFloatType dSigma = m_LPF1.getFeedbackOutput() +
                        m_LPF2.getFeedbackOutput() +
                        m_LPF3.getFeedbackOutput() +
                        m_LPF4.getFeedbackOutput();

        // --- for passband gain compensation
        //     you can connect this to a GUI control
        //     and let user choose instead
        xn *= 1.0 + m_dAuxControl*m_dK;

        // calculate input to first filter
        DspFloatType dU = (xn - m_dK*dSigma)*m_dAlpha_0;

        if(m_uNLP == ON)
            dU = tanh(m_dSaturation*dU);

        // --- cascade of 4 filters
        DspFloatType dLP1 = m_LPF1.doFilter(dU);
        DspFloatType dLP2 = m_LPF2.doFilter(dLP1);
        DspFloatType dLP3 = m_LPF3.doFilter(dLP2);
        DspFloatType dLP4 = m_LPF4.doFilter(dLP3);

        // --- Oberheim variations
        return m_dA*dU + m_dB*dLP1 + m_dC*dLP2 + m_dD*dLP3 +  m_dE*dLP4;
    }
 	DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {
		return A * doFilter(I);
	}

    void ProcessSIMD(size_t n, DspFloatType * in, DspFloatType * out) {
		#pragma omp simd aligned(in,out)
		for(size_t i = 0; i < n; i++) {
            DspFloatType xn = in[i];
            // --- return xn if filter not supported
            if(m_uFilterType == BSF2 || m_uFilterType == LPF1 || m_uFilterType == HPF1)
                out[i] = xn;
            else {
                DspFloatType dSigma = m_LPF1.getFeedbackOutput() +
                                m_LPF2.getFeedbackOutput() +
                                m_LPF3.getFeedbackOutput() +
                                m_LPF4.getFeedbackOutput();

                // --- for passband gain compensation
                //     you can connect this to a GUI control
                //     and let user choose instead
                xn *= 1.0 + m_dAuxControl*m_dK;

                // calculate input to first filter
                DspFloatType dU = (xn - m_dK*dSigma)*m_dAlpha_0;

                if(m_uNLP == ON)
                    dU = tanh(m_dSaturation*dU);

                // --- cascade of 4 filters
                DspFloatType dLP1 = m_LPF1.doFilter(dU);
                DspFloatType dLP2 = m_LPF2.doFilter(dLP1);
                DspFloatType dLP3 = m_LPF3.doFilter(dLP2);
                DspFloatType dLP4 = m_LPF4.doFilter(dLP3);

                // --- Oberheim variations
                out[i] = m_dA*dU + m_dB*dLP1 + m_dC*dLP2 + m_dD*dLP3 +  m_dE*dLP4;
            }
        }
    }
};

CMoogLadderFilter::CMoogLadderFilter(void)
{
    // --- init
    m_dK = 0.0;
    m_dAlpha_0 = 1.0;

    // --- Oberheim variables
    m_dA = 0.0; m_dB = 0.0; m_dC = 0.0;
    m_dD = 0.0; m_dE = 0.0;
    
    // --- set all as LPF types
    m_LPF1.m_uFilterType = LPF1;
    m_LPF2.m_uFilterType = LPF1;
    m_LPF3.m_uFilterType = LPF1;
    m_LPF4.m_uFilterType = LPF1;

    // --- set default filter type
    m_uFilterType = LPF4;

    // --- flush everything
    reset();
}

CMoogLadderFilter::~CMoogLadderFilter(void)
{
}

void CMoogLadderFilter::reset()
{
    // flush everything
    m_LPF1.reset();
    m_LPF2.reset();
    m_LPF3.reset();
    m_LPF4.reset();
}

// decode the Q value; Q on UI is 1->10
void CMoogLadderFilter::setQControl(DspFloatType dQControl)
{
    // this maps dQControl = 1->10 to K = 0 -> 4
    m_dK = (4.0)*(dQControl - 1.0)/(10.0 - 1.0);
}
