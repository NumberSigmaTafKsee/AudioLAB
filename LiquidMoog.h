#pragma once

enum {
	CUTOFF,
	Q,
	LAST
};

enum PortGroups
{
    PORT_AUDIO_IN   = 0,    
    PORT_AUDIO_OUT,
    PORT_MIDI_IN,
    PORT_CUTOFF,
    PORT_Q,
    PORT_NR
};

const char * urn = "urn:james5:LiquidMoog";
