
#include "COscillator.h"

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