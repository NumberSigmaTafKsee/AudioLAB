// 2PACSV8
// 8 band State Variable Filter Equalizer


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
#include "Analog/VAAnalogSVF.hpp"
#include "FX/CSmoothFilters.hpp"

#define LOOP(X,start,end) for(size_t X = start; X < end; X++)

using namespace Analog::Filters:AnalogSVF;
using namespace FX::Filters::Smoothers;

struct Urids
{
    LV2_URID midi_MidiEvent;
};

std::string bundlepath;
const char * urn = "urn:james5:MoogLadder";

enum {
	PARALLEL,
	BYPASS1,
	GAIN1,
	CUTOFF1,
	Q1,
	TYPE1,
	BYPASS2,
	GAIN2,		
	CUTOFF2,
	Q2,
	BYPASS3,
	GAIN3,
	TYPE3,		
	CUTOFF3,
	Q3,	
	BYPASS4,
	GAIN4,
	CUTOFF4,
	Q4,
	TYPE4,
	BYPASS5,
	GAIN5,		
	CUTOFF5,
	Q5,
	BYPASS6,
	GAIN6,
	TYPE6,		
	CUTOFF6,
	Q6,	
	BYPASS7,
	GAIN7,
	TYPE7,		
	CUTOFF7,
	Q7,	
	BYPASS8,
	GAIN8,
	TYPE8,		
	CUTOFF8,
	Q8,	
	LAST
};

enum PortGroups
{
    PORT_AUDIO_IN   = 0,    
    PORT_AUDIO_OUT,
    PORT_MIDI_IN,
    PORT_PARALLEL,
	PORT_BYPASS1,
	PORT_GAIN1,
	PORT_CUTOFF1,
	PORT_Q1,
	PORT_TYPE1,
	PORT_BYPASS2,
	PORT_GAIN2,		
	PORT_CUTOFF2,
	PORT_Q2,
	PORT_BYPASS3,
	PORT_GAIN3,
	PORT_TYPE3,		
	PORT_CUTOFF3,
	PORT_Q3,	
	PORT_BYPASS4,
	PORT_GAIN4,
	PORT_CUTOFF4,
	PORT_Q4,
	PORT_TYPE4,
	PORT_BYPASS5,
	PORT_GAIN5,		
	PORT_CUTOFF5,
	PORT_Q5,
	PORT_BYPASS6,
	PORT_GAIN6,
	PORT_TYPE6,		
	PORT_CUTOFF6,
	PORT_Q6,	
	PORT_BYPASS7,
	PORT_GAIN7,
	PORT_TYPE7,		
	PORT_CUTOFF7,
	PORT_Q7,	
	PORT_BYPASS8,
	PORT_GAIN8,
	PORT_TYPE8,		
	PORT_CUTOFF8,
	PORT_Q8,	
    PORT_NR
};

struct LV2Plugin
{
    float* audio_in_ptr;
    float* audio_out_ptr;
    const LV2_Atom_Sequence* midi_in_ptr;
    const LV2_Atom_Sequence* midi_out_ptr;    
    LV2_URID_Map* map ;
    Urids urids;
    double rate;
	
	std::array<float*,LAST> controls;
	std::array<float,LAST> control_values;
	std::array<FX::Filters::Smoothers::CSmoothFilter*,LAST> smoothers;
	
	std::array<AnalogSVF*,8> filters;
	
    LV2Plugin(const double sampleRate, const LV2_Feature *const *features) :    
    midi_in_ptr (nullptr),
    midi_out_ptr (nullptr),
    audio_in_ptr (nullptr),    
    audio_out_ptr (nullptr),    
    map (nullptr),
    rate (sampleRate),
    filter(nullptr)    
    {
		const char* missing = lv2_features_query
		(
			features,
			LV2_URID__map, &map, true,
			NULL
		);

		if (missing) throw std::invalid_argument ("Feature map not provided by the host. Can't instantiate LV2Lua");

		urids.midi_MidiEvent = map->map (map->handle, LV2_MIDI__MidiEvent);
		
		LOOP(i,0,8) {
			filters[i] = new AnalogSVF(44100.0,1000,0.5);
		LOOP(i,0,LAST) {
			smoothers[i] = new CSmoothFilter(sampleRate,1.0/0.001);
		}
    }
    ~LV2Plugin() {
		LOOP(i,0,8) if(filters[i]) delete filters[i];
		LOOP(i,0,LAST) if(smoothers[i]) delete smoothers[i];
	}
} ;




static LV2_Handle instantiate (const struct LV2_Descriptor *descriptor, double sample_rate, const char *bundle_path, const LV2_Feature *const *features)
{    	
    LV2MoogLadder *m = new LV2MoogLadder(sample_rate, features);
    assert(m != nullptr);
    bundlepath = bundle_path;    
    return m;
}


static void connect_port (LV2_Handle instance, uint32_t port, void *data_location)
{
    LV2MoogLadder *m = (LV2MoogLadder*) instance;
    if (!m) return;

    switch (port)
    {
    case PORT_AUDIO_IN:
        m->audio_in_ptr = (float*) data_location;
        break;

    case PORT_AUDIO_OUT:
        m->audio_out_ptr = (float*) data_location;
        break;
    /*
    case PORT_MIDI_IN:
        m->midi_in_ptr = static_cast<const LV2_Atom_Sequence*> (data_location);
        break;
	*/
	case PORT_CUTOFF:
		m->controls[CUTOFF] = (float*)data_location;
		break;

	case PORT_Q:
		m->controls[Q] = (float*)data_location;
		break;
	
	case PORT_TYPE:
		m->controls[TYPE] = (float*)data_location;
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
    LV2MoogLadder* m = (LV2MoogLadder*) instance;
    if(!m) return;   
    if(!m->filter) return;
    if(!m->audio_in_ptr) return;
    if(!m->audio_out_ptr) return; 
    
    
    if(*m->controls[TYPE] != m->control_values[TYPE]) {
		m->control_values[0] = *m->controls[0];
		m->filter->setPort(MoogLadderFilter::PORT_TYPE,m->control_values[TYPE]);
	}
    LOOP(i,0,LAST) {
		if(*m->controls[i] != m->control_values[i]) {
			m->control_values[i] = *m->controls[i];
			m->smoothers[i]->setTarget(m->control_values[i]);
		}
	}
	
	
	m->filter->setPort(MoogLadderFilter::PORT_CUTOFF,m->smoothers[CUTOFF]->process());
	m->filter->setPort(MoogLadderFilter::PORT_RESONANCE,m->smoothers[Q]->process());	
	
	// pre-effect
	// tap-input insert
    m->filter->ProcessBlock(sample_count,m->audio_in_ptr, m->audio_out_ptr);    
    // post-effect
    // tap-output insert
    
}

static void deactivate (LV2_Handle instance)
{
    /* not needed here */
}

static void cleanup (LV2_Handle instance)
{
    LV2MoogLadder * l = (LV2MoogLadder*)instance;
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
