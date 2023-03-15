#include "BlitOscillator.h"
#include "LV2Plugin.hpp"

using namespace Analog::Oscillators::Blit;
using namespace FX::Filters::Smoothers;


struct LV2BlitOscillator : public DSP
{
    double tune=0;
	double fine=0;	
	Blit    	* blit;
	BlitSaw 	* saw;
	BlitSquare 	* square;	
	float kc  = 0;
	float vel = 1;
	
    LV2BlitOscillator(float sampleRate) :
    DSP(sampleRate),
    midi_in_ptr (nullptr),
    midi_out_ptr (nullptr),
    audio_in_ptr (nullptr),    
    audio_out_ptr (nullptr),    
    map (nullptr),
    rate (sampleRate),
    blit(nullptr),
    saw(nullptr),
    square(nullptr)    
    {		
		blit = new Blit(220,sampleRate);
		saw  = new BlitSaw(220,sampleRate);
		square  = new BlitSquare(220,sampleRate);
		
    }
    ~LV2BlitOscillator() {
		if(blit) delete blit;
		if(saw) delete saw;
		if(square) delete square;	
	}
	void setFrequency(float f) {
		blit->setFrequency(f);
		saw->setFrequency(f);
		square->setFrequency(f);
	}
	void setHarmonics(unsigned int h) {
		blit->setHarmonics(h);
		saw->setHarmonics(h);
		square->setHarmonics(h);
	}
	void reset() {
		blit->reset();
		saw->reset();
		square->reset();
	}
	void savePhase() {
		blit->savePhase();
		saw->savePhase();
		square->savePhase();
	}
	void restorePhase() {
		blit->restorePhase();
		saw->restorePhase();
		square->restorePhase();
	}
	void setPhase(float p) {
		blit->setPhase(p);
		saw->setPhase(p);
		square->setPhase(p);
	}
	void setDuty(float p) {
		square->setDuty(p);
	}
	bool zeros(size_t sample_count, float * memory)
	{
		LOOP(i,1,sample_count) if(memory[i] != 0.0f) return false;
		return true;
	}
	void play(uint32_t frame, uint32_t sample_count) {
		LOOP(i,frame,sample_count) {
			float r;							
			savePhase();
			setPhase(controls[PM][i]);
					
			float br,sr,sqr;
			br  = blit->Tick(); 
			sr  = saw->Tick(); 
			sqr = square->Tick(); 
			r = (1.0-control_values[MIX])*sr + (control_values[MIX])*sqr;			
			r = (1.0-control_values[RM])*r + (control_values[RM])*r*audio_in_ptr[i];			
			audio_out_ptr[i] = r;
			saw_out_ptr[i] = sr;
			square_out_ptr[i] = sqr;
			impulse_out_ptr[i] = br;
			restorePhase();
		}			
	}
	void connect_port (LV2_Handle instance, uint32_t port, void *data_location)
	{
		LV2BlitOscillator *m = (LV2BlitOscillator*) instance;
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
		
		case PORT_FREQUENCY:
			m->controls[FREQUENCY] = (float*)data_location;
			break;

		case PORT_SEMITONE:
			m->controls[SEMITONE] = (float*)data_location;
			break;
		
		case PORT_OCTAVE:
			m->controls[OCTAVE] = (float*)data_location;
			break;
		
		case PORT_FINE:
			m->controls[FINE] = (float*)data_location;
			break;
		
		case PORT_HARMONICS:
			m->controls[HARMONICS] = (float*)data_location;
			break;
		
		case PORT_DUTY:
			m->controls[DUTY] = (float*)data_location;
			break;
		case PORT_MIX:
			m->controls[MIX] = (float*)data_location;
			break;
		
		case PORT_CVSYNC:
			m->controls[CVSYNC] = (float*)data_location;
			break;
			
		case PORT_CVDUTY:
			m->controls[CVDUTY] = (float*)data_location;
			break;
		
		case PORT_PM:
			m->controls[PM] = (float*)data_location;
			break;
		
		case PORT_RM:
			m->controls[RM] = (float*)data_location;
			break;
		
		default:
			break;
		}
	}

	void Tick(LV2Plugin * m, uint32_t from, uint32_t to)
	{
		
	}

	void run (LV2_Handle instance, uint32_t sample_count)
	{
		LV2BlitOscillator* m = (LV2BlitOscillator*) instance;
		if(!m) return;       
		if(!m->audio_in_ptr) return;
		if(!m->audio_out_ptr) return; 
		
		bool setFreq = false;
		bool setHarmonics = false;
		if(*m->controls[FREQUENCY] != m->control_values[FREQUENCY]) setFreq = true;
		if(*m->controls[HARMONICS] != m->control_values[HARMONICS]) setHarmonics = true;
		   
		LOOP(i,0,LAST) {
			if(*m->controls[i] != m->control_values[i]) {
				m->control_values[i] = *m->controls[i];
				m->smoothers[i]->setTarget(m->control_values[i]);
			}
		}
		
		float f = m->kc+m->smoothers[FREQUENCY]->process() + m->smoothers[FINE]->process();
		// + MusicFunctions::semitone(m->control_values[SEMITONE]) + MusicFunctions::octave(m->controls[OCTAVE]);	
		float duty = m->smoothers[DUTY]->process();
		unsigned int harmonics = m->smoothers[HARMONICS]->process();
		
		if(setFreq) m->setFrequency(f);
		m->square->setDuty(duty);
		if(setHarmonics) m->setHarmonics(harmonics);
		
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
		Tick(m,last_frame,sample_count);
	}
};


DSP * createDSP(float sr) {
	return LV2BlitOscillator(sr);
}
