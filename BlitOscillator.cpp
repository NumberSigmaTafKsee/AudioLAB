
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <vector>
#include <array>
#include <string>
#include <iostream>

#include "lv2.h"
#include <lv2/atom/atom.h>
#include <lv2/urid/urid.h>
#include <lv2/midi/midi.h>
#include <lv2/core/lv2_util.h>
#include <lv2/atom/util.h>

#include "SoundObject.hpp"
#include "MusicFunctions.hpp"
#include "Analog/VABlitOscillators.hpp"
#include "FX/CSmoothFilters.hpp"

#define LOOP(X,start,end) for(size_t X = start; X < end; X++)

using namespace Analog::Oscillators::Blit;
using namespace FX::Filters::Smoothers;

struct Urids
{
    LV2_URID midi_MidiEvent;
};

std::string bundlepath;
const char * urn = "urn:james5:BlitOscillators";

enum {
	FREQUENCY,
	SEMITONE,
	OCTAVE,
	FINE,	
	HARMONICS,
	DUTY,			
	MIX,
	IT,
	CVSYNC,
	CVDUTY,
	PM,
	RM,
	LAST
};

enum PortGroups
{
    PORT_AUDIO_IN   = 0,    
    PORT_AUDIO_OUT,
    PORT_MIDI_IN,
    PORT_MIDI_OUT,
    PORT_TYPE,
    PORT_FREQUENCY,        
    PORT_SEMITONE,
    PORT_OCTAVE,
    PORT_FINE,  
    PORT_HARMONICS,  
    PORT_DUTY,    
    PORT_MIX,
    PORT_IT,
    PORT_CVSYNC,    
    PORT_CVDUTY,
    PORT_PM,        
    PORT_RM,            
    PORT_NR
};

struct LV2BlitOscillator
{
    float* audio_in_ptr;
    float* audio_out_ptr;
    const LV2_Atom_Sequence* midi_in_ptr;
    const LV2_Atom_Sequence* midi_out_ptr;    
    LV2_URID_Map* map ;
    Urids urids;
    double rate;
	double tune=0;
	double fine=0;
	double prev_int =0;
	double prev_diff = 0;
	
	std::array<float*,8> controls;
	std::array<float,8> control_values;
	std::array<FX::Filters::Smoothers::CSmoothFilter*,8> smoothers;
	Blit    * blit;
	BlitSaw * saw;
	BlitSquare * square;
	int type = 1;
	float kc  = 0;
	float vel = 1;
	
    LV2BlitOscillator(const double sampleRate, const LV2_Feature *const *features) :    
    midi_in_ptr (nullptr),
    midi_out_ptr (nullptr),
    audio_in_ptr (nullptr),    
    audio_out_ptr (nullptr),    
    map (nullptr),
    rate (sampleRate)
    {
		const char* missing = lv2_features_query
		(
			features,
			LV2_URID__map, &map, true,
			NULL
		);

		if (missing) throw std::invalid_argument ("Feature map not provided by the host. Can't instantiate LV2Lua");

		urids.midi_MidiEvent = map->map (map->handle, LV2_MIDI__MidiEvent);
		
		blit = new Blit(220,sampleRate);
		saw  = new BlitSaw(220,sampleRate);
		square  = new BlitSquare(220,sampleRate);
		LOOP(i,0,LAST) {
			smoothers[i] = new CSmoothFilter(sampleRate,1.0/0.001);
		}
    }
    ~LV2BlitOscillator() {
		if(blit) delete blit;
		if(saw) delete saw;
		if(square) delete square;
		LOOP(i,0,LAST) if(smoothers[i]) delete smoothers[i];
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
			if(!zeros(sample_count,controls[CVDUTY])) setDuty(controls[CVDUTY][i]);		
			if(!zeros(sample_count,controls[CVSYNC])) if( fabs(controls[CVSYNC][i]) < 1e-6 ) reset();		
			savePhase();
			setPhase(controls[PM][i]);
					
			float br,sr,sqr;
			br  = blit->Tick(); 
			sr  = saw->Tick(); 
			sqr = square->Tick(); 
			r = (1.0-control_values[MIX])*sr + (control_values[MIX])*sqr;			
			r = (1.0-control_values[RM])*r + (control_values[RM])*r*audio_in_ptr[i];
			r = (1.0-control_values[IT])*br + (control_values[IT])*r;
			audio_out_ptr[i] = r;
			restorePhase();
		}			
	}
};




static LV2_Handle instantiate (const struct LV2_Descriptor *descriptor, double sample_rate, const char *bundle_path, const LV2_Feature *const *features)
{    	
    LV2BlitOscillator *m = new LV2BlitOscillator(sample_rate, features);
    assert(m != nullptr);
    bundlepath = bundle_path;    
    return m;
}


static void connect_port (LV2_Handle instance, uint32_t port, void *data_location)
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

static void activate (LV2_Handle instance)
{
    /* not needed here */
    
}


static void run (LV2_Handle instance, uint32_t sample_count)
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
		m->play(frame,sample_count);
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
        
	m->play(last_frame,sample_count);	
}

static void deactivate (LV2_Handle instance)
{
    /* not needed here */
}

static void cleanup (LV2_Handle instance)
{
    LV2BlitOscillator * l = (LV2BlitOscillator*)instance;
    if(l) delete l;
}

static const void* extension_data (const char *uri)
{
    return NULL;
}

/* descriptor */
static LV2_Descriptor const descriptor =
{
    urn,
    instantiate,
    connect_port,
    activate /* or NULL */,
    run,
    deactivate /* or NULL */,
    cleanup,
    extension_data /* or NULL */
};

/* interface */
LV2_SYMBOL_EXPORT const LV2_Descriptor* lv2_descriptor (uint32_t index)
{
    if (index == 0) return &descriptor;
    else return NULL;
}
