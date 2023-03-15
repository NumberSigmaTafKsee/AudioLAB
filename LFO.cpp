/// LFO
/// Audio Modulation
/// CV Output
/// Midi Modulation CC

#include "LFO.h"
#include "LV2Plugin.hpp"

#include "FX/LFO.hpp"

struct LV2LFO : public DSP
{
	LFO * lfo;
	
	LV2LFO(float sr) : DSP(sr)
	{
		lfo = new LFO(sr);
	}
	~LV2LFO() {
		if(adsr) delete adsr;
	}
	void update(LV2Plugin * p)
	{
		if(*p->controls[FREQ] != p->control_values[FREQ])		   
		{			
			lfo->setRate(*p->controls[FREQ]);
		}
		if(*p->controls[WAVEFORM] != p->control_values[WAVEFORM])
		{
			lfo->setWaveform((LFO::waveform_t)*p->controls[WAVEFORM]);
		}
	}
	bool zeros(size_t n, float * memory) {
		LOOP(i,1,n)
			if(memory[i] != 0.0f) return false;
		return true;
	}
	void run(LV2Plugin * m, int sample_count)
	{		
		if(!m) return;   
		if(!m->dsp) return;
		if(!m->audio_in_ptr) return;
		if(!m->audio_out_ptr) return; 
		if(!m->midi_in_ptr) return;
		if(!adsr) return;
		
		update(m);
		
		LOOP(i,0,LAST) {
			if(*m->controls[i] != m->control_values[i]) {
				m->control_values[i] = *m->controls[i];
				m->smoothers[i]->setTarget(m->control_values[i]);
			}
		}
				
		lfo->ProcessInplace(sample_count,m->controls[LFO]);		
		for(size_t i = 0; i < sample_count; i++)
			m->audio_out_ptr[i] = m->audio_in_ptr[i] * m->controls[LFO][i];
		
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

		case PORT_MIDI_OUT:
			m->midi_out_ptr = static_cast<const LV2_Atom_Sequence*> (data_location);
			break;

		case PORT_FREQ:
			m->controls[FREQ] = (float*)data_location;
			break;
		
		case PORT_WAVEFORM:
			m->controls[WAVEFORM] = (float*)data_location;
			break;
		
		case PORT_FREQCV:
			m->controls[FREQCV] = (float*)data_location;
			break;
			
		case PORT_PHASECV:
			m->controls[PHASECV] = (float*)data_location;
			break;
			
		default:
			break;
		}
}

};


DSP* createDSP(float sr) {	
	return new LV2ADSR(sr);
}

