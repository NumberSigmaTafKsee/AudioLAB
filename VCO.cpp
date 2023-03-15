// VCO
// Vectorized Polyblep


#include "VCO.h"
#include "LV2Plugin.hpp"

namespace Analog;


struct LV2VCO : public DSP
{
	VCO * vco;
	
	LV2VCO(float sr) : DSP(sr)
	{
		vco = new VCO(sr,VCOPolyBLEP::Waveform::SAWTOOTH);
	}
	~LV2VCO() {
		if(vco) delete vco;
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
				
		bool allZeros=true;
		LOOP(i,0,sample_count) {
			if(m->audio_in_ptr[i] != 0) {
				allZeros = false;
				break;
			}
		}
		if(allZeros) {
			LOOP(i,0,sample_count) m->audio_in_ptr = 1.0;
		}
		vco->ProcessBlock(sample_count,m->audio_in_ptr,m->audio_out_ptr);
		
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


