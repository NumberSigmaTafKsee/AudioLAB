// Zolzer
#include "ZolzerFilters.h"
#include "LV2Plugin.hpp"

#include "Filters/IIRZolzerFilter.hpp"

using namespace Filters::IIR::ZolzerFilters;

struct LV2ZolzerFilter : public DSP
{
	ZolzerBiquadFilter* filter;
	ZolzerBiquadFilter* stereo[2];	
	DspFloatType fs2;
	DspFloatType kc=0,vel=1;
	
	LV2ZolzerFilter(float sr) : DSP(sr)
	{
		fs2 = sr/2.0;
		filter = new ZolzerBiquadFilter(ZolzerBiquadFilter::ZolzerType::LOWPASS,sr);	
		LOOP(i,0,2) {
			stereo[i] = new ZolzerBiquadFilter(ZolzerBiquadFilter::ZolzerType::LOWPASS,sr);	
		}	
	}
	~LV2ZolzerFilter() {		
		if(filter) delete filter;		
		LOOP(i,0,2) if(stereo[i]) delete stereo[i];		
	}
	void setFilterType(ZolzerBiquadFilter::ZolzerType type)
	{
		filter->setType(type);		
		LOOP(i,0,2) stereo[i]->setType(type);		
	}
	void setCutoff(DspFloatType c) {
		if(c < 30 || c >= fs2) return;
		filter->setCutoff(c);
		LOOP(i,0,2) stereo[i]->setCutoff(c);		
	}
	void setQ(DspFloatType q) {
		filter->setQ(q);
		LOOP(i,0,2) stereo[i]->setQ(q);		
	}
	void setGain(DspFloatType g) {				
		filter->setGain(g);
		LOOP(i,0,2) stereo[i]->setGain(g);		
	}				
	bool zeros(size_t from, size_t to, float * memory) {
		LOOP(i,from,to) if(memory[i] != 0.0f) return false;
		return true;
	}
	void Tick(LV2Plugin * m, uint32_t from, uint32_t to)
	{
		bool zeroPan = zeros(from,to,m->controls[PANCV]);
		LOOP(i,from,to)
		{
			float cut = kc*vel+m->smoothers[CUTOFF1]->process();
			float q   = vel*m->smoothers[Q1]->process();
			float gain= vel*m->smoothers[GAIN1]->process();			
			float pan   = vel*m->smoothers[PAN]->process();
						
			setCutoff(m->controls[CUTOFFCV][i] * cut + cut);
			setQ(m->controls[QCV][i] + q);
			setGain(gain);
						
			m->audio_out_ptr[i] = filter->Tick(m->audio_in_ptr[i]);			
			
			if(m->stereo_in_ptr[0] && m->stereo_out_ptr[0] && m->stereo_in_ptr[1] && m->stereo_out_ptr[1])
			{				
				float p     = m->controls[PANCV][i]*pan;		
				if(zeroPan) p = pan;
				float left  = std::sin((1 - (p)) * M_PI / 2);
				float right = std::cos((p) * M_PI / 2);			
				m->stereo_out_ptr[0][i]  = left*stereo[0]->Tick(m->stereo_in_ptr[0][i]);
				m->stereo_out_ptr[1][i]  = right*stereo[1]->Tick(m->stereo_in_ptr[1][i]);
			}			
		}		
	}
	void run(LV2Plugin * m, int sample_count)
	{		
		if(!m) return;   
		if(!m->dsp) return;
		if(!m->audio_in_ptr) return;
		if(!m->audio_out_ptr) return; 		
		if(!filter) return;
		
		if(*m->controls[TYPE1] != m->control_values[TYPE1]) {
			setFilterType((ZolzerBiquadFilter::ZolzerType)*m->controls[TYPE1]);
		}
		LOOP(i,0,LAST-3) {
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
			Tick(m,last_frame,frame);
			
			if (ev->body.type == m->urids.midi_MidiEvent)
			{
				const uint8_t* const msg = reinterpret_cast<const uint8_t*> (ev + 1);
				const uint8_t typ = lv2_midi_message_type (msg);

				switch (typ)
				{
				case LV2_MIDI_MSG_NOTE_ON:				
					kc = MusicFunctions::midi2freq(msg[0]);				
					vel= msg[1]/127.0f;                
					break;

				case LV2_MIDI_MSG_NOTE_OFF:                
					break;

				case LV2_MIDI_MSG_CONTROLLER:                
					break;
				
				default:
					break;
				}
			}		
			last_frame = frame;	
		}
		Tick(m,last_frame,sample_count);
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
		

		case PORT_LEFTIN:
			m->stereo_in_ptr[0] = (float*) data_location;
			break;    
			
		case PORT_RIGHTIN:
			m->stereo_in_ptr[1] = (float*) data_location;
			break;    

		case PORT_LEFTOUT:
			m->stereo_out_ptr[0] = (float*) data_location;
			break;    
			
		case PORT_RIGHTOUT:
			m->stereo_out_ptr[1] = (float*) data_location;
			break;    

		case PORT_PAN:
			m->controls[PAN] = (float*)data_location;
			break;
			
		case PORT_TYPE1:
			m->controls[TYPE1] = (float*)data_location;
			break;

		case PORT_CUTOFF1:
			m->controls[CUTOFF1] = (float*)data_location;
			break;
		
		case PORT_Q1:
			m->controls[Q1] = (float*)data_location;
			break;
			
		case PORT_GAIN1:
			m->controls[GAIN1] = (float*)data_location;
			break;
		
		case PORT_CUTOFFCV:
			m->controls[CUTOFFCV] = (float*)data_location;
			break;
		
		case PORT_QCV:
			m->controls[QCV] = (float*)data_location;
			break;
					
		case PORT_PANCV:
			m->controls[PANCV] = (float*)data_location;
			break;
		
		default:
			break;
		}
}

};


DSP* createDSP(float sr) {	
	return new LV2ZolzerFilter(sr);
}

