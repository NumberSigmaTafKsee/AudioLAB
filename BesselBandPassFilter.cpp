
// for now
// MIDI CC 103 = Cutoff
// MIDI CC 104 = Bandwidth

#include "BesselBandPassFilter.h"
#include "LV2Plugin.hpp"

#include "Filters/DspBesselBandPass.hpp"

namespace Filters;

struct LV2BesselBandpassFilter : public DSP
{
	BesselBandPassFilter * filter;
	
	LV2BesselBandpassFilter(float sr) : DSP(sr)
	{
		filter = new BesselBandPassFilter(2,100.0,1000.0,sr);		
	}
	~LV2BesselBandpassFilter() {
		if(filter) delete filter;
	}
	void update(LV2Plugin * m) 
	{
		if( *m->controls[PORT_ORDER] != m->control_values[ORDER] || 
		    *m->controls[PORT_CUTOFF] != m->control_values[CUTOFF] ||
		    *m->controls[PORT_BANDWIDTH] != m->control_values[BANDWIDTH]) {
			m->control_values[i] = *m->controls[i];
			filter->order = m->control_values[ORDER];
			filter->bw    = m->control_values[BANDWIDTH];
			filter->setCutoff(m->control_values[CUTOFF]);			
		}
	}	
	void run(LV2Plugin * m, int sample_count)
	{		
		if(!m) return;   
		if(!m->dsp) return;
		if(!m->audio_in_ptr) return;
		if(!m->audio_out_ptr) return; 
		if(!m->midi_in_ptr) return;
		if(!moog) return;
		
		update(m);
		
		uint32_t last_frame = 0;
		LV2_ATOM_SEQUENCE_FOREACH (m->midi_in_ptr, ev)
		{
			/* play frames until event */
			const uint32_t frame = ev->time.frames;			
			filter->ProcessBlock(frame, m->audio_in_ptr+last_frame, m->audio_out_ptr+last_frame);				
			last_frame = frame;

			if (ev->body.type == m->urids.midi_MidiEvent)
			{
				const uint8_t* const msg = reinterpret_cast<const uint8_t*> (ev + 1);
				const uint8_t typ = lv2_midi_message_type (msg);

				switch (typ)
				{
				case LV2_MIDI_MSG_NOTE_ON:									
					
					break;

				case LV2_MIDI_MSG_NOTE_OFF:
					
					break;

				case LV2_MIDI_MSG_CONTROLLER:                
					if( (msg[0] & 0xF0) == 0xB0 )
					{
						if(msg[1] == 103) 
							filter->setCutoff(filter->fc + fc*msg[2]/127.0);
						else if(msg[1] == 104)
							filter->setBandWidth(filter->bw + filter->bw*msg[2]/127.0);
					}
					break;
				
				default:
					break;
				}
			}
		}
		
		LOOP(i,0,LAST) {
			if(*m->controls[i] != m->control_values[i]) {
				m->control_values[i] = *m->controls[i];
				m->smoothers[i]->setTarget(m->control_values[i]);
			}
		}
					
		filter->ProcessBlock(frame, m->audio_in_ptr+last_frame, m->audio_out_ptr+last_frame);				
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
		
		case PORT_BANDWIDTH:
			m->controls[BANDWIDTH] = (float*)data_location;
			break;
			
		case PORT_ORDER:
			m->controls[ORDER] = (float*)data_location;
			break;
			
		default:
			break;
		}
}

};


DSP* createDSP(float sr) {	
	return new LV2BesselBandpassFilter(sr);
}

