#pragma once

enum {
	ATTACK,
	SUSTAIN,
	RELEASE,	
	ENV,
	LAST
};

enum {
	PORT_AUDIO_IN   = 0,    
    PORT_AUDIO_OUT,
    PORT_MIDI_IN,
	PORT_MIDI_OUT,
	PORT_ATTACK,
	PORT_SUSTAIN,
	PORT_RELEASE,
	PORT_ENV,
	PORT_NR,
};

const char * urn = "urn:james5:ASR";
