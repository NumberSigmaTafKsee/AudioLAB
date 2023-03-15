// it does not clip very much
// the capacitor was set like a filter
// the circuit needs to be tweaked

#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <vector>
#include <array>
#include <iostream>

#include "lv2.h"
#include <lv2/atom/atom.h>
#include <lv2/urid/urid.h>
#include <lv2/midi/midi.h>
#include <lv2/core/lv2_util.h>
#include <lv2/atom/util.h>

#include "Analog/amp_clipper_circuit.hpp"
#include "FX/CSmoothFilters.hpp"

using namespace FX::Distortion::ClipperCircuit;
using namespace FX::Filters::Smoothers;

#define LOOP(X,start,end) for(size_t X = start; X < end; X++)

const char * MYCLASS_URD = "urn:james5:AmpClipperCircuit";

struct Urids
{
    LV2_URID midi_MidiEvent;
};

struct Plugin {

	ClipperCircuit<float> clipper;
	
	Plugin(float sampleRate) : clipper(sampleRate) {
		
	}
	void ProcessBlock(size_t n, float * in, float * out) {
		clipper.ProcessSIMD(n,in,out);
	}
};

struct MYCLASS
{
	// generic code
	float* audio_in_ptr;
    float* audio_out_ptr;
    const LV2_Atom_Sequence* midi_in_ptr;
    const LV2_Atom_Sequence* midi_out_ptr;    
    LV2_URID_Map* map ;
    Urids urids;
    double rate;
	
	// plugin code here
	Plugin * ptr = NULL;
	std::array<float*,3> controls;
	std::array<float,3> control_values;
	std::array<FX::Filters::Smoothers::CSmoothFilter*,3> smoothers;
	
	MYCLASS(Plugin * p, const double sampleRate, const LV2_Feature *const *features) :    
	ptr(p),
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
		
		LOOP(i,0,3) {
			smoothers[i] = new CSmoothFilter(sampleRate,1.0/0.001);
		}
    }
    ~MYCLASS() {
		if(ptr) delete ptr;
		LOOP(i,0,3) if(smoothers[i]) delete smoothers[i];
	}
};


/* internal core methods */
static LV2_Handle instantiate (const struct LV2_Descriptor *descriptor, double sample_rate, const char *bundle_path, const LV2_Feature *const *features)
{    	
	Plugin * plugin = new Plugin(sample_rate);
    MYCLASS *m = new MYCLASS(plugin,sample_rate, features);
    return m;
}

enum {
	IDEALITY,
	CAPACITANCE,
	ASYMMETRY,
	LAST
};

enum PortGroups
{
    PORT_AUDIO_IN   = 0,
    PORT_AUDIO_OUT,
    PORT_IDEALITY,
    PORT_CAPACITANCE,
    PORT_ASYMMETRY,
    PORT_NR,
};

static void connect_port (LV2_Handle instance, uint32_t port, void *data_location)
{
    MYCLASS* m = (MYCLASS*) instance;
    if (!m) return;

    switch (port)
    {
    case PORT_AUDIO_IN:
        m->audio_in_ptr = (float*) data_location;
        break;

    case PORT_AUDIO_OUT:
        m->audio_out_ptr = (float*) data_location;
        break;
    
    // Control go here
	case PORT_IDEALITY:
		m->controls[IDEALITY] = (float*)data_location;
		break;

	case PORT_CAPACITANCE:
		m->controls[CAPACITANCE] = (float*)data_location;
		break;
	
	case PORT_ASYMMETRY:
		m->controls[ASYMMETRY] = (float*)data_location;
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
    MYCLASS* m = (MYCLASS*) instance;
    
    if(!m) return;   
    if(!m->ptr) return;
    if(!m->audio_in_ptr) return;
    if(!m->audio_out_ptr) return; 
    
    
    LOOP(i,0,LAST) {
		if(*m->controls[i] != m->control_values[i]) {
			m->control_values[i] = *m->controls[i];
			m->smoothers[i]->setTarget(m->control_values[i]);
		}
	}
	
	m->ptr->clipper.setIdeality(m->smoothers[IDEALITY]->process());
	m->ptr->clipper.setCapacitance(m->smoothers[CAPACITANCE]->process());
	m->ptr->clipper.setAsymmetry(m->smoothers[ASYMMETRY]->process());
	
	
	// pre-effect
	// tap-input insert
    m->ptr->clipper.ProcessBlock(sample_count,m->audio_in_ptr, m->audio_out_ptr);    
    // post-effect
    // tap-output insert
	// run the plugin code here
}

static void deactivate (LV2_Handle instance)
{
    /* not needed here */
}

static void cleanup (LV2_Handle instance)
{
    MYCLASS * l = (MYCLASS*)instance;
    if(l) delete l;
}

static const void* extension_data (const char *uri)
{
    return NULL;
}

/* descriptor */
static LV2_Descriptor const descriptor =
{
    MYCLASS_URD,
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
