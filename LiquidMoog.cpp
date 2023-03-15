// LiquidMoog

// must be define before including LV2Plugin
#include "LiquidMoog.h"
#include "LV2Plugin.hpp"

#include "FX/Amplifiers.hpp"
#include "FX/LiquidMoog.hpp"

using namespace Liquid;

struct LV2LiquidMoog : public DSP
{
	LiquidMoog * moog;
	
	LV2LiquidMoog(float sr) : DSP(sr)
	{
		moog = new LiquidMoog(sr,1000.0,0.5);		
	}
	~LV2LiquidMoog() {
		if(moog) delete moog;
	}
	
	void run(LV2Plugin * m, int sample_count)
	{		
		if(!m) return;   
		if(!m->dsp) return;
		if(!m->audio_in_ptr) return;
		if(!m->audio_out_ptr) return; 
		if(!m->midi_in_ptr) return;
		if(!moog) return;
		
		LOOP(i,0,LAST) {
			if(*m->controls[i] != m->control_values[i]) {
				m->control_values[i] = *m->controls[i];
				m->smoothers[i]->setTarget(m->control_values[i]);
			}
		}
				
		moog->setPort(LiquidMoog::PORT_CUTOFF,m->smoothers[CUTOFF]->process());
		moog->setPort(LiquidMoog::PORT_RESONANCE,m->smoothers[Q]->process());	
		
		moog->ProcessBlock(sample_count,m->audio_in_ptr,m->audio_out_ptr);
		
	}	
	
	void connect_port (LV2Plugin * m, uint32_t port, void *data_location)
	{		
		if (!m) return;		
		switch (port)
		{
		case PORT_AUDIO_IN:
			m->audio_in_ptr = (float*) data_location;
			break;

		case PORT_AUDIO_OUT:
			m->audio_out_ptr = (float*) data_location;
			break;    
		case PORT_MIDI_IN:
			m->midi_in_ptr = static_cast<const LV2_Atom_Sequence*> (data_location);
			break;

		case PORT_CUTOFF:
			m->controls[CUTOFF] = (float*)data_location;
			break;
		
		case PORT_Q:
			m->controls[Q] = (float*)data_location;
			break;
			
		default:
			break;
		}
}

};


DSP* createDSP(float sr) {	
	return new LV2LiquidMoog(sr);
}
