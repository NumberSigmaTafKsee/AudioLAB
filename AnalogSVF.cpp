/* include libs */
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

#include "Analog/VAAnalogSVF.hpp"
#include "FX/CSmoothFilters.hpp"

#define LOOP(X,start,end) for(size_t X = start; X < end; X++)

using namespace Analog::Filters::AnalogSVF;
using namespace FX::Filters::Smoothers;

struct Urids
{
    LV2_URID midi_MidiEvent;
};

std::string bundlepath;

struct LV2AnalogSVF
{
    float* audio_in_ptr;
    float* audio_out_ptr;
    const LV2_Atom_Sequence* midi_in_ptr;
    const LV2_Atom_Sequence* midi_out_ptr;    
    LV2_URID_Map* map ;
    Urids urids;
    double rate;
	
	std::array<float*,11> controls;
	std::array<float,11> control_values;
	std::array<FX::Filters::Smoothers::CSmoothFilter*,11> smoothers;
	Analog::Filters::AnalogSVF::AnalogSVF * filter;
	
    LV2AnalogSVF(const double sampleRate, const LV2_Feature *const *features) :    
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
		
		filter = new Analog::Filters::AnalogSVF::AnalogSVF(sampleRate,1000.0/sampleRate,0.5);
		LOOP(i,0,11) {
			smoothers[i] = new CSmoothFilter(sampleRate,1.0/0.001);
		}
    }
    ~LV2AnalogSVF() {
		if(filter) delete filter;
		LOOP(i,0,11) if(smoothers[i]) delete smoothers[i];
	}
} ;




static LV2_Handle instantiate (const struct LV2_Descriptor *descriptor, double sample_rate, const char *bundle_path, const LV2_Feature *const *features)
{    	
    LV2AnalogSVF *m = new LV2AnalogSVF(sample_rate, features);
    assert(m != nullptr);
    bundlepath = bundle_path;
    std::cout << bundle_path << std::endl;
    return m;
}

enum {
	CUTOFF,
	Q,
	TYPE,
	GAIN,
	MINCLIP,
	MAXCLIP,
	A,
	X,
	Y,	
	LAST
};

enum PortGroups
{
    PORT_AUDIO_IN   = 0,    
    PORT_AUDIO_OUT,
    PORT_CUTOFF,    
    PORT_Q,
    PORT_TYPE,
    PORT_GAIN,
    PORT_MINCLIP,
    PORT_MAXCLIP,
    PORT_A,
    PORT_X,
    PORT_Y,
    PORT_NR
};

static void connect_port (LV2_Handle instance, uint32_t port, void *data_location)
{
    LV2AnalogSVF *m = (LV2AnalogSVF*) instance;
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
		
	case PORT_GAIN:
		m->controls[GAIN] = (float*)data_location;
		break;
	
	case PORT_MINCLIP:
		m->controls[MINCLIP] = (float*)data_location;
		break;
	
	case PORT_MAXCLIP:
		m->controls[MAXCLIP] = (float*)data_location;
		break;
	
	case PORT_A:
		m->controls[A] = (float*)data_location;
		break;
	
	case PORT_X:
		m->controls[X] = (float*)data_location;
		break;
	
	case PORT_Y:
		m->controls[Y] = (float*)data_location;
		break;
	
	case PORT_NR:
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
    LV2AnalogSVF* m = (LV2AnalogSVF*) instance;
    if(!m) return;   
    if(!m->filter) return;
    if(!m->audio_in_ptr) return;
    if(!m->audio_out_ptr) return; 
    
    LOOP(i,0,LAST) {
		if(*m->controls[i] != m->control_values[i]) {
			m->control_values[i] = *m->controls[i];
			m->smoothers[i]->setTarget(m->control_values[i]);
		}
	}
	
	m->filter->setPort(AnalogSVF::PORT_TYPE,m->control_values[TYPE]);
	m->filter->setPort(AnalogSVF::PORT_CUTOFF,m->smoothers[CUTOFF]->process());
	m->filter->setPort(AnalogSVF::PORT_Q,m->smoothers[Q]->process());	
	m->filter->setPort(AnalogSVF::PORT_GAIN,m->smoothers[GAIN]->process());
	m->filter->setPort(AnalogSVF::PORT_MINC,m->smoothers[MINCLIP]->process());
	m->filter->setPort(AnalogSVF::PORT_MAXC,m->smoothers[MAXCLIP]->process());
	
	
	// pre-effect
	// tap-input insert
    m->filter->ProcessSIMD(sample_count,m->audio_in_ptr, m->audio_out_ptr);    
    // post-effect
    // tap-output insert
    
}

static void deactivate (LV2_Handle instance)
{
    /* not needed here */
}

static void cleanup (LV2_Handle instance)
{
    LV2AnalogSVF * l = (LV2AnalogSVF*)instance;
    if(l) delete l;
}

static const void* extension_data (const char *uri)
{
    return NULL;
}

/* descriptor */
static LV2_Descriptor const descriptor =
{
    "urn:james5:AnalogSVF",
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
