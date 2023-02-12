#pragma once
//#include "pluginconstants.h"
#include "synthlite.h"
#include <cstdint>
#include <cstdlib>
#include <ctime>

#define OSC_FO_MOD_RANGE 2			//2 semitone default
#define OSC_HARD_SYNC_RATIO_RANGE 4	//
#define OSC_PITCHBEND_MOD_RANGE 12	//12 semitone default
#define OSC_FO_MIN 20				//20Hz
#define OSC_FO_MAX 20480			//20.480kHz = 10 octaves up from 20Hz
#define OSC_FO_DEFAULT 440.0		//A5
#define OSC_PULSEWIDTH_MIN 2		//2%
#define OSC_PULSEWIDTH_MAX 98		//98%
#define OSC_PULSEWIDTH_DEFAULT 50	//50%

class COscillator 
{
public:
    COscillator(void);
    virtual ~COscillator(void);	

    // --- ATTRIBUTES
    // --- PUBLIC: these variables may be get/set 
    //             you may make get/set functions for them 
    //             if you like, but will add function call layer	
    // --- oscillator run flag
    bool m_bNoteOn;

    // --- user controls or MIDI 
    DspFloatType m_dOscFo;		// oscillator frequency from MIDI note number
    DspFloatType m_dFoRatio;	    // FM Synth Modulator OR Hard Sync ratio 
    DspFloatType m_dAmplitude;	// 0->1 from GUI
    
    // ---  pulse width in % (sqr only) from GUI
    DspFloatType m_dPulseWidthControl;	

    // --- modulo counter and inc for timebase
    DspFloatType m_dModulo;		// modulo counter 0->1
    DspFloatType m_dInc;			// phase inc = fo/fs

    // --- more pitch mods
    int m_nOctave;			// octave tweak
    int m_nSemitones;		// semitones tweak
    int m_nCents;			// cents tweak
    
    // --- for PITCHED Oscillators
    enum {SINE,SAW1,SAW2,SAW3,TRI,SQUARE,NOISE,PNOISE};
    uint32_t m_uWaveform;	// to store type
    
    // --- for LFOs
    enum {sine,usaw,dsaw,tri,square,expo,rsh,qrsh};

    // --- for LFOs - MODE
    enum {sync,shot,free};
    uint32_t m_uLFOMode;	// to store MODE
        
    // --- MIDI note that is being played
    uint32_t m_uMIDINoteNumber;
    
protected:
    // --- PROTECTED: generally these are either basic calc variables
    //                and modulation stuff
    // --- calculation variables
    DspFloatType m_dSampleRate;	// fs
    DspFloatType m_dFo;			// current (actual) frequency of oscillator	
    DspFloatType m_dPulseWidth;	// pulse width in % for calculation
    
    // --- for noise and random sample/hold
    uint32_t   m_uPNRegister;	// for PN Noise sequence
    int    m_nRSHCounter;	// random sample/hold counter
    DspFloatType m_dRSHValue;		// currnet rsh output
    
    // --- for DPW
    DspFloatType m_dDPWSquareModulator;	// square toggle
    DspFloatType m_dDPW_z1; // memory register for differentiator
    
    // --- mondulation inputs
    DspFloatType m_dFoMod;			/* modulation input -1 to +1 */
    DspFloatType m_dPitchBendMod;	    /* modulation input -1 to +1 */
    DspFloatType m_dFoModLin;			/* FM modulation input -1 to +1 (not actually used in Yamaha FM!) */
    DspFloatType m_dPhaseMod;			/* Phase modulation input -1 to +1 (used for DX synth) */
    DspFloatType m_dPWMod;			/* modulation input for PWM -1 to +1 */
    DspFloatType m_dAmpMod;			/* output amplitude modulation for AM 0 to +1 (not dB)*/
        
public:
    // --- FUNCTIONS: all public
    //
    // --- modulo functions for master/slave operation
    // --- increment the modulo counters
    inline void incModulo(){m_dModulo += m_dInc;}

    // --- check and wrap the modulo
    //     returns true if modulo wrapped
    inline bool checkWrapModulo()
    {	
        // --- for positive frequencies
        if(m_dInc > 0 && m_dModulo >= 1.0) 
        {
            m_dModulo -= 1.0; 
            return true;
        }
        // --- for negative frequencies
        if(m_dInc < 0 && m_dModulo <= 0.0) 
        {
            m_dModulo += 1.0; 
            return true;
        }
        return false;
    }
    
    // --- reset the modulo (required for master->slave operations)
    inline void resetModulo(DspFloatType d = 0.0){m_dModulo = d;}

    // --- modulation functions - NOT needed/used if you implement the Modulation Matrix!
    //
    // --- output amplitude modulation (AM, not tremolo (dB); 0->1, NOT dB
    inline void setAmplitudeMod(DspFloatType dAmp){m_dAmpMod = dAmp;}

    // --- modulation, exponential
    inline void setFoModExp(DspFloatType dMod){m_dFoMod = dMod;}
    inline void setPitchBendMod(DspFloatType dMod){m_dPitchBendMod = dMod;}
    
    // --- for FM only (not used in Yamaha or my DX synths!)
    inline void setFoModLin(DspFloatType dMod){m_dFoModLin = dMod;}

    // --- for Yamaha and my DX Synth
    inline void setPhaseMod(DspFloatType dMod){m_dPhaseMod = dMod;}

    // --- PWM for square waves only
    inline void setPWMod(DspFloatType dMod){m_dPWMod = dMod;}

    // --- VIRTUAL FUNCTIONS ----------------------------------------- //
    //
    // --- PURE ABSTRACT: derived class MUST implement
    // --- start/stop control
    virtual void startOscillator() = 0;
    virtual void stopOscillator() = 0;
    
    // --- render a sample
    //		for LFO:	 pAuxOutput = QuadPhaseOutput
    //			Pitched: pAuxOutput = Right channel (return value is left Channel
    virtual DspFloatType doOscillate(DspFloatType* pAuxOutput = NULL) = 0; 

    // --- ABSTRACT: derived class overrides if needed
    virtual void setSampleRate(DspFloatType dFs){m_dSampleRate = dFs;}
    
    // --- reset counters, etc...
    virtual void reset();

    // INLINE FUNCTIONS: these are inlined because they will be 
    //                   called every sample period
    //					 You may want to move them to the .cpp file and
    //                   enable the compiler Optimization setting for 
    //					 Inline Function Expansion: Any Suitable though
    //					 inlining here forces it.
    //

    // --- update the frequency, amp mod and PWM
    inline virtual void update()
    {		
        // --- ignore LFO mode for noise sources
        if(m_uWaveform == rsh || m_uWaveform == qrsh)
        m_uLFOMode = free;

        // --- do the  complete frequency mod
        m_dFo = m_dOscFo*m_dFoRatio*pitchShiftMultiplier(m_dFoMod + 
                                                        m_dPitchBendMod + 
                                                        m_nOctave*12.0 + 
                                                        m_nSemitones + 
                                                        m_nCents/100.0);
        
        // --- apply linear FM (not used in book projects)
        m_dFo += m_dFoModLin;

        // --- bound Fo (can go outside for FM/PM mod)
        //     +/- 20480 for FM/PM
        if(m_dFo > OSC_FO_MAX)
            m_dFo = OSC_FO_MAX;
        if(m_dFo < -OSC_FO_MAX)
            m_dFo = -OSC_FO_MAX;

        // --- calculate increment (a.k.a. phase a.k.a. phaseIncrement, etc...)
        m_dInc = m_dFo/m_dSampleRate;

        // --- Pulse Width Modulation --- //
        // --- limits are 2% and 98%
        m_dPulseWidth = m_dPulseWidthControl + m_dPWMod*(OSC_PULSEWIDTH_MAX - OSC_PULSEWIDTH_MIN)/OSC_PULSEWIDTH_MIN; 

        // --- bound the PWM to the range
        m_dPulseWidth = fmin(m_dPulseWidth, OSC_PULSEWIDTH_MAX);
        m_dPulseWidth = fmax(m_dPulseWidth, OSC_PULSEWIDTH_MIN);
    }

    DspFloatType Tick(DspFloatType I=0, DspFloatType A=1, DspFloatType X=0, DspFloatType Y=0)
    {
        DspFloatType f = m_dFo;
        DspFloatType p = m_dPulseWidth;
        DspFloatType x = m_dPitchBendMod;
        /*
        m_dPitchBendMod = I;
        m_dFo = abs(X);
        m_dPulseWidth = abs(Y);
        */
        DspFloatType out = doOscillate();
        /*
        m_dFo = f;
        m_dPulseWidth = p;
        m_dPitchBendMod = x;
        */
        return A*out;
    }
    void ProcessSIMD(size_t n, DspFloatType * out) {
		#pragma omp simd
		for(size_t i = 0; i < n; i++) {
            out[i] = doOscillate();
        }
    }
};

// --- construction
COscillator::COscillator(void)
{
    // --- initialize variables
    m_dSampleRate = 44100;	
    m_bNoteOn = false;
    m_uMIDINoteNumber = 0;
    m_dModulo = 0.0;		
    m_dInc = 0.0;			
    m_dOscFo = OSC_FO_DEFAULT; // GUI
    m_dAmplitude = 1.0; // default ON
    m_dPulseWidth = OSC_PULSEWIDTH_DEFAULT;	
    m_dPulseWidthControl = OSC_PULSEWIDTH_DEFAULT; // GUI
    m_dFo = OSC_FO_DEFAULT; 			

    // --- seed the random number generator
    srand(time(NULL));
    m_uPNRegister = rand();

    // --- continue inits
    m_nRSHCounter = -1; // flag for reset condition
    m_dRSHValue = 0.0;
    m_dAmpMod = 1.0; // note default to 1 to avoid silent osc
    m_dFoModLin = 0.0;
    m_dPhaseMod = 0.0;
    m_dFoMod = 0.0;
    m_dPitchBendMod = 0.0;
    m_dPWMod = 0.0;
    m_nOctave = 0.0;
    m_nSemitones = 0.0;
    m_nCents = 0.0;
    m_dFoRatio = 1.0;
    m_uLFOMode = 0;

    // --- pitched
    m_uWaveform = SINE;
}

// --- destruction
COscillator::~COscillator(void)
{
}
    
// --- VIRTUAL FUNCTION; base class implementations
void COscillator::reset()
{
    // --- Pitched modulos, wavetables start at 0.0
    m_dModulo = 0.0;	
        
    // --- needed fror triangle algorithm, DPW
    m_dDPWSquareModulator = -1.0;

    // --- flush DPW registers
    m_dDPW_z1 = 0.0;

    // --- for random stuff
    srand(time(NULL));
    m_uPNRegister = rand();
    m_nRSHCounter = -1; // flag for reset condition
    m_dRSHValue = 0.0;

    // --- modulation variables
    m_dAmpMod = 1.0; // note default to 1 to avoid silent osc
    m_dPWMod = 0.0;
    m_dPitchBendMod = 0.0;
    m_dFoMod = 0.0;
    m_dFoModLin = 0.0;
    m_dPhaseMod = 0.0;
}
