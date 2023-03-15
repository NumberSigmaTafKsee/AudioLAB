
#include "AnalogSVF.h"
#include "LV2Plugin.hpp"
#include <algorithm>

#include "MusicFunctions.hpp"
#include "Analog/VAAnalogSVF.hpp"
#include "FX/LFO.hpp"
#include "FX/ADSR.hpp"


using namespace Analog::Filters::AnalogSVF;
using namespace FX::Filters::Smoothers;
using namespace Envelopes;

struct LV2AnalogSVF : public DSP
{
    float* lp_out,*hp_out,*bp_out,*ubp_out,*shelf_out,*notch_out,*apf_out,*peak_out;
    float* Acv,*Xcv,*Ycv;    
    double rate;
	double kc = 0;
	double vel= 1.0;
	
	std::array<float,1024> temp1,temp2;
	AnalogSVF * filter;
	LFO * lfo;
	ADSR * adsr;
	
    LV2AnalogSVF(const double sampleRate) : DSP(sampleRate),
    rate (sampleRate),
    filter(nullptr),
    lp_out(nullptr),
    hp_out(nullptr),
    bp_out(nullptr),
    ubp_out(nullptr),
    notch_out(nullptr),
    peak_out(nullptr),
    shelf_out(nullptr),
    apf_out(nullptr)
    {		
		filter = new Analog::Filters::AnalogSVF::AnalogSVF(sampleRate,1000.0/sampleRate,0.5);
		adsr   = new ADSR(sampleRate);
		lfo    = new LFO(sampleRate);				
    }
    ~LV2AnalogSVF() {
		if(filter) delete filter;
		if(adsr) delete adsr;
		if(lfo) delete lfo;		
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
		
		case PORT_TYPE:
			m->controls[TYPE] = (float*)data_location;
			break;
			
		case PORT_GAIN:
			m->controls[GAIN] = (float*)data_location;
			break;
		
		case PORT_MINCLIP:
			m->controls[MINCLIP] = (float*)data_location;
			break;
		
		case PORT_MAXCLIP:
			m->controls[MAXCLIP] = (float*)data_location;
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
		
		case PORT_CUTOFF_ENV:
			m->controls[CUTOFF_ENV] = (float*)data_location;
			break;
		
		case PORT_CUTOFF_LFO:
			m->controls[CUTOFF_LFO] = (float*)data_location;
			break;
			
		case PORT_Q_ENV:
			m->controls[Q_ENV] = (float*)data_location;
			break;
		
		case PORT_Q_LFO:
			m->controls[Q_LFO] = (float*)data_location;
			break;
		
		case PORT_LFO_TYPE:
			m->controls[LFO_TYPE] = (float*)data_location;
			break;
			
		case PORT_LFO_FREQ:
			m->controls[LFO_FREQ] = (float*)data_location;
			break;
		
		case PORT_A:
			Acv = (float*)data_location;				
			break;
		case PORT_X:
			Xcv = (float*)data_location;		
			break;
		case PORT_Y:
			Ycv = (float*)data_location;		
			break;
			
		case PORT_LPOUT:
			lp_out = (float*)data_location;		
			break;
		case PORT_HPOUT:
			hp_out = (float*)data_location;		
			break;
		case PORT_BPOUT:
			bp_out = (float*)data_location;		
			break;
		case PORT_UBPOUT:
			ubp_out = (float*)data_location;		
			break;
		case PORT_SHELFOUT:
			shelf_out = (float*)data_location;		
			break;
		case PORT_NOTCHOUT:
			notch_out = (float*)data_location;		
			break;
		case PORT_APFOUT:
			apf_out = (float*)data_location;		
			break;
		case PORT_PEAKOUT:
			peak_out = (float*)data_location;		
			break;
		
		case PORT_NR:
			break;
			
		default:
			break;
		}
	}


	bool zeros(size_t from, size_t to, float * memory) 
	{	
		LOOP(i,from,to)
			if(memory[i] != 0.0f) return false;
		return true;
	}


	/*streams
		LOOP(i,0,sample_count)
		{	
			m->filter->Aa[i] = 1.0;			
			if(zeros(sample_count,m->controls[A])) std::fill(m->controls[A],m->controls[A]+sample_count,1.0f);
			if(m->controls[A]) 	m->filter->Aa[i] *= std::clamp(m->controls[A][i],-1.0f,1.0f);						
			m->filter->Xa[i] = CE*m->temp1[i];
			m->filter->Xa[i] += CL*m->temp2[i];				
			if(!zeros(sample_count,m->controls[X])) m->filter->Xa[i] *= std::clamp(m->controls[X][i],-1.0f,1.0f); 		
			m->filter->Ya[i] = QE*m->temp1[i];
			m->filter->Ya[i] += QL*m->temp2[i];		
			if(!zeros(sample_count,m->controls[Y])) m->filter->Ya[i] *= std::clamp(m->controls[Y][i],-1.0f,1.0f);		
		}
		m->filter->ProcessBlock(sample_count,m->audio_in_ptr,m->audio_out_ptr);	    
		*/
		
	void Tick(LV2Plugin * m, uint32_t from, uint32_t to)
	{
		// the host should fill these but it doesn't
		if(zeros(from,to,Acv)) std::fill(Acv+from,Acv+to,1.0f);
		
		adsr->ProcessInplace(to-from,temp1.data()+from);	
		lfo->setRate(m->smoothers[LFO_FREQ]->process());	
		lfo->ProcessInplace(to-from,temp2.data()+from);
			
		float eps = std::numeric_limits<float>::epsilon();
				
		LOOP(i,from,to)
		{
			
			float CE = m->smoothers[CUTOFF_ENV]->process();			
			float CL = m->smoothers[CUTOFF_LFO]->process();			
			float QE = m->smoothers[Q_ENV]->process();			
			float QL = m->smoothers[Q_LFO]->process();
				
			float A = 1.0;					
			A *= std::clamp(Acv[i],-1.0f,1.0f);						
			
			float xcv = std::clamp(Xcv[i],-1.0f,1.0f);
			float ycv = std::clamp(Ycv[i],-1.0f,1.0f);
		
			float X = CE*temp1[i] + xcv + CL*temp2[i];
			float Y = QE*temp1[i] + ycv + QL*temp2[i];
								
			float fC = m->smoothers[CUTOFF]->process();
			float fT = kc + vel*fC;
							
			filter->setPort(AnalogSVF::PORT_CUTOFF,fT);
			filter->setPort(AnalogSVF::PORT_Q,vel*m->smoothers[Q]->process());	
			filter->setPort(AnalogSVF::PORT_GAIN,vel*m->smoothers[GAIN]->process());
			
			// these are supposed to be CV ports
			filter->setPort(AnalogSVF::PORT_MINC,m->smoothers[MINCLIP]->process());
			filter->setPort(AnalogSVF::PORT_MAXC,m->smoothers[MAXCLIP]->process());
										
			float r = filter->Tick(m->audio_in_ptr[i],A,3*X,Y);
			m->audio_out_ptr[i] = r;
			lp_out[i] = filter->lp;
			hp_out[i] = filter->hp;
			bp_out[i] = filter->bp;
			ubp_out[i] = filter->ubp;
			shelf_out[i] = filter->shelf;
			peak_out[i] = filter->peak;
			apf_out[i] = filter->apf;
			notch_out[i] = filter->notch;
		}
	}
	void run (LV2Plugin * m, uint32_t sample_count)
	{		
		if(!m) return;   
		if(!m->audio_in_ptr) return;
		if(!m->audio_out_ptr) return; 
		
		if(*m->controls[TYPE] != m->control_values[TYPE]) {
				m->control_values[TYPE] = *m->controls[TYPE];
				filter->setPort(AnalogSVF::PORT_TYPE,m->control_values[TYPE]);
		}
		if(*m->controls[LFO_TYPE] != m->control_values[LFO_TYPE]) {
				m->control_values[LFO_TYPE] = *m->controls[LFO_TYPE];
				lfo->setWaveform((LFO::waveform_t)m->control_values[LFO_TYPE]);			
		}
		if(*m->controls[ATTACK] != m->control_values[ATTACK]) {
				m->control_values[ATTACK] = *m->controls[ATTACK];
				adsr->setAttackTime(m->control_values[ATTACK]);			
		}
		if(*m->controls[DECAY] != m->control_values[DECAY]) {
				m->control_values[DECAY] = *m->controls[DECAY];
				adsr->setDecayTime(m->control_values[DECAY]);			
		}
		if(*m->controls[SUSTAIN] != m->control_values[SUSTAIN]) {
				m->control_values[SUSTAIN] = *m->controls[SUSTAIN];
				adsr->setSustainLevel(m->control_values[SUSTAIN]);			
		}
		if(*m->controls[RELEASE] != m->control_values[RELEASE]) {
				m->control_values[RELEASE] = *m->controls[RELEASE];
				adsr->setReleaseTime(m->control_values[RELEASE]);			
		}
			
		LOOP(j,0,LAST) {
				if( *m->controls[j] != m->control_values[j]) {
					m->control_values[j] = *m->controls[j];
					m->smoothers[j]->setTarget(m->control_values[j]);
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
					adsr->noteOn();
					break;

				case LV2_MIDI_MSG_NOTE_OFF:
					kc = 0;
					vel = 1.0;
					adsr->noteOff();
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
};

DSP* createDSP(float sr) {
	return new LV2AnalogSVF(sr);
}
	
