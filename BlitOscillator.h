#pragma once

const char * urn = "urn:james5:BlitOscillators";

enum {
	FREQUENCY,
	SEMITONE,
	OCTAVE,
	FINE,	
	HARMONICS,
	DUTY,			
	MIX,	
	CVSYNC,
	CVDUTY,
	CV_PM,
	CV_RM,			
	CV_SAW_OUT,
	CV_SQR_OUT,
	CV_IT_OUT,
	CV_INT_OUT,	
	CV_AUX_IN,		
	LAST
};

enum PortGroups
{
    PORT_AUDIO_IN   = 0,    
    PORT_AUDIO_OUT,
    PORT_MIDI_IN,
    PORT_MIDI_OUT,
    PORT_SAW_OUT,
    PORT_SQUARE_OUT,
    PORT_IMPULSE_OUT,
    PORT_AUX_INPUT,        
    PORT_FREQUENCY,        
    PORT_SEMITONE,
    PORT_OCTAVE,
    PORT_FINE,  
    PORT_HARMONICS,  
    PORT_DUTY,    
    PORT_MIX,    
    PORT_CVSYNC,    
    PORT_CVDUTY,
    PORT_CV_PM,        
    PORT_CV_RM,            
    PORT_CV_SAW_OUT,
    PORT_CV_SQUARE_OUT,
    PORT_CV_IMPULSE_OUT,
    PORT_CV_AUX_IN,
    PORT_NR
};
