
// Generic Template for LV2 Plugin theory

#include "stddef.h"
#include "stdint.h"
#include "stdlib.h"
#include "lv2.h"
#include <lv2/atom/atom.h>
#include <lv2/urid/urid.h>
#include <lv2/midi/midi.h>
#include <lv2/core/lv2_util.h>
#include <lv2/atom/util.h>

#include "Analog/VAAnalogSVF.hpp"
#include "FX/CSmoothFilters.hpp"


const char * MYCLASS_URD = "urn:james5:lv2luajit";

struct Urids
{
    LV2_URID midi_MidiEvent;
};

struct Plugin {

	void  setParameter(int port, float value) = 0;
	float getParameter(int port) = 0;
	
	void ProcessBlock(size_t n, float * in, float *out = 0);
};

struct MYCLASS

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
		
    }
    ~MYCLASS() {
		if(ptr) delete ptr;
	}
};


/* internal core methods */
static LV2_Handle instantiate (const struct LV2_Descriptor *descriptor, double sample_rate, const char *bundle_path, const LV2_Feature *const *features)
{    	
	Plugin * plugin;
    MYCLASS *m = new MYCLASS(plugin,sample_rate, features);
    return m;
}

enum {
	NOCONTROL=-1,
};

enum PortGroups
{
    PORT_AUDIO_IN   = 0,
    PORT_AUDIO_OUT  = 1,
        
    PORT_NR         = 2,
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
    if (!m) return;    

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
