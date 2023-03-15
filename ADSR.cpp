/// ADSR
/// Attack Decay Sustain Release Envelope
/// Audio Modulation
/// CV Output
/// Midi CC Output

#include "ADSR.h"
#include "LV2Plugin.hpp"

#include "FX/ADSR.hpp"

using namespace Envelopes;

const int NUM_FILTERS = 12;


struct LV2ADSR : public DSP
{
	ADSR * adsr;
	
	LV2ADSR(float sr) : DSP(sr)
	{
		adsr = new ADSR(sr);		
	}
	~LV2ADSR() {
		if(adsr) delete adsr;
	}
	void update(LV2Plugin * p)
	{
		if(*p->controls[ATTACK] != p->control_values[ATTACK] ||
		   *p->controls[DECAY] != p->control_values[DECAY] ||
		   *p->controls[SUSTAIN] != p->control_values[SUSTAIN] ||
		   *p->controls[RELEASE] != p->control_values[RELEASE] )
		{			
			adsr->setAllTimes(*p->controls[ATTACK],
							  *p->controls[DECAY],
							  *p->controls[SUSTAIN],
							  *p->controls[RELEASE]);
		}
	}
	bool iszeros(size_t n, float * memory) {
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
		
		uint32_t last_frame = 0;
		LV2_ATOM_SEQUENCE_FOREACH (m->midi_in_ptr, ev)
		{
			/* play frames until event */
			const uint32_t frame = ev->time.frames;
						
			adsr->ProcessInplace(frame, m->controls[ENV] + last_frame);	
			
			last_frame = frame;

			if (ev->body.type == m->urids.midi_MidiEvent)
			{
				const uint8_t* const msg = reinterpret_cast<const uint8_t*> (ev + 1);
				const uint8_t typ = lv2_midi_message_type (msg);

				switch (typ)
				{
				case LV2_MIDI_MSG_NOTE_ON:									
					adsr->noteOn();
					break;

				case LV2_MIDI_MSG_NOTE_OFF:
					adsr->noteOff();
					break;

				case LV2_MIDI_MSG_CONTROLLER:                
					break;
				
				default:
					break;
				}
			}
		}
		
		adsr->ProcessInplace(sample_count,m->controls[ENV]);		
		for(size_t i = last_frame; i < sample_count; i++)
			m->audio_out_ptr[i] = m->audio_in_ptr[i] * m->controls[ENV][i];
		
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

		case PORT_ATTACK:
			m->controls[ATTACK] = (float*)data_location;
			break;
		
		case PORT_DECAY:
			m->controls[DECAY] = (float*)data_location;
			break;
		
		case PORT_SUSTAIN:
			m->controls[SUSTAIN] = (float*)data_location;
			break;
			
		case PORT_RELEASE:
			m->controls[RELEASE] = (float*)data_location;
			break;
		
		case PORT_ENV:
			m->controls[ENV] = (float*)data_location;
			break;
			
		default:
			break;
		}
}

};


DSP* createDSP(float sr) {	
	return new LV2ADSR(sr);
}

