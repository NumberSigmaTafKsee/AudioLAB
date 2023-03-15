// Stereo AnalogSVF
#include "StereoAnalogSVF.h"
#include "LV2Plugin.hpp"

#include "Analog/AnalogSVF.hpp"

using namespace Analog;
using namespace Analog::Filters;
using namespace Analog::Filters::AnalogSVF;

struct LV2StereoAnalogSVF : public DSP
{
	AnalogSVF* filter;
	AnalogSVF* stereo[2];
	DspFloatType fs2;
	
	LV2StereoAnalogSVF(float sr) : DSP(sr)
	{
		fs2 = sr/2.0;
		filter = new RBJBiquadFilter(RBJBiquadFilter::RBJType::LOWPASS,sr);	
		LOOP(i,0,2) {
			stereo[i] = new RBJBiquadFilter(RBJBiquadFilter::RBJType::LOWPASS,sr);	
		}	
	}
	~LV2StereoAnalogSVF() {		
		if(filter) delete filter;		
		LOOP(i,0,2) if(stereo[i]) delete stereo[i];		
	}
	void setFilterType(int type)
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
	void Tick(uint32_t from, uint32_t to )
	{
		bool zeroPan = zeros(sample_count,m->controls[PANCV]);
		LOOP(i,0,sample_count)
		{
			float cut = m->smoothers[CUTOFF1]->process();
			float q   = m->smoothers[Q1]->process();
			float gain= m->smoothers[GAIN1]->process();
			float pan   = m->smoothers[PAN]->process();
						
			setCutoff(m->controls[CUTOFFCV][i] * cut + cut);
			setQ(m->controls[QCV][i] + q);
			setGain(gain);
			setSlope(slope);
			
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
	bool zeros(size_t n, float * memory) {
		LOOP(i,1,n) if(memory[i] != 0.0f) return false;
		return true;
	}
	void run(LV2Plugin * m, int sample_count)
	{		
		if(!m) return;   
		if(!m->dsp) return;
		if(!m->audio_in_ptr) return;
		if(!m->audio_out_ptr) return; 		
		if(!filter) return;
		
		if(*m->controls[TYPE1] != m->control_values[TYPE1]) {
			setFilterType((RBJBiquadFilter::RBJType)*m->controls[TYPE1]);
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
			Tick(last_frame,frame);
			
			if (ev->body.type == m->urids.midi_MidiEvent)
			{
				const uint8_t* const msg = reinterpret_cast<const uint8_t*> (ev + 1);
				const uint8_t typ = lv2_midi_message_type (msg);

				switch (typ)
				{
				case LV2_MIDI_MSG_NOTE_ON:				
					m->kc = MusicFunctions::midi2freq(msg[0]);				
					m->vel= msg[1]/127.0f;                
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
		Tick(last_frame,sample_count);		
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
	return new LV2StereoAnalogSVF(sr);
}

