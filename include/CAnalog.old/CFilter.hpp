#pragma once
//#include "synthfunctions.h"
#include "synthlite.h"
#include <cstdint>


// 46.881879936465680 semitones = semitonesBetweenFrequencies(80, 18000.0)/2.0
#define FILTER_FC_MOD_RANGE 46.881879936465680	
#define FILTER_FC_MIN 80		// 80Hz
#define FILTER_FC_MAX 18000		// 18kHz
#define FILTER_FC_DEFAULT 10000	// 10kHz
#define FILTER_Q_DEFAULT 0.707	// Butterworth

// CFilter Abastract Base Class for all filters
struct CFilter : public FilterProcessor
{
public:
    CFilter(void);
    ~CFilter(void);
    
    // --- ATTRIBUTES
    // --- PUBLIC: these variables may be get/set 
    //             you may make get/set functions for them 
    //             if you like, but will add function call layer

    // --- the user's cutoff frequency control position
    DspFloatType m_dFcControl;
    
    // --- the user's cutoff frequency control position
    DspFloatType m_dQControl;			

    // --- for an aux filter specific like SEM BSF 
    //     control or paasband gain comp (Moog)
    DspFloatType m_dAuxControl;

    // --- for NLP - Non Linear Procssing
    uint32_t m_uNLP;
    enum{OFF,ON};  // switch enum

    // --- to add more distortion
    DspFloatType m_dSaturation;

    // --- NOTE: these are shared; even though some filters won't use some of them
    //           need to maintain the indexing
    enum{LPF1,HPF1,LPF2,HPF2,BPF2,BSF2,LPF4,HPF4,BPF4};

    // --- our selected filter type
    uint32_t m_uFilterType;		


    // --- PROTECTED: generally these are either basic calc variables
    //                and modulation stuff
    //
    DspFloatType m_dSampleRate;	// fs

    // --- the actual cutoff used in the calculation
    DspFloatType m_dFc;			
    
    // --- the current value of Q (internal)
    DspFloatType m_dQ;			

    // --- our cutoff frequency modulation input
    DspFloatType m_dFcMod;		

    // --- add more mods here


    // --- FUNCTIONS: all public
    //	
    inline void setFcMod(DspFloatType d){m_dFcMod = d;}

    // --- VIRTUAL FUNCTIONS ----------------------------------------- //
    //
    // --- PURE ABSTRACT: derived class MUST implement
    virtual DspFloatType doFilter(DspFloatType xn) = 0;

    // --- ABSTRACT: derived class overrides if needed
    //
    inline virtual void setSampleRate(DspFloatType d){m_dSampleRate = d;}

    // --- flush buffers, reset filter
    virtual void reset();

    // --- decode the Q value; this can change from filter to filter
    virtual void setQControl(DspFloatType dQControl);
    
    // --- recalculate the Fc (called after modulations)
    inline virtual void update()
    {
        // --- update Q (filter-dependent)
        setQControl(m_dQControl);
        
        // --- do the modulation freq shift
        m_dFc = m_dFcControl*pitchShiftMultiplier(m_dFcMod);

        // --- bound the final frequency
        if(m_dFc > FILTER_FC_MAX)
            m_dFc = FILTER_FC_MAX;
        if(m_dFc < FILTER_FC_MIN)
            m_dFc = FILTER_FC_MIN;
    }

    DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1)
    {
        DspFloatType f = m_dFc;
        DspFloatType q = m_dQ;
        setFcMod(f*abs(X));
        setQControl(q*abs(Y));
        DspFloatType out =  A*doFilter(I);
        setFcMod(f);
        setQControl(q);
        return out;
    }
    void ProcessBlock(size_t n, DspFloatType * in, DspFloatType * out) {
        #pragma omp simd aligned(in,out)
        for(size_t i = 0; i < n; i++) out[i] = Tick(in[i]);
    }
};

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
