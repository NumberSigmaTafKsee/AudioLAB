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
#include "FX/CSmoothFilters.hpp"

#define LOOP(X,start,end) for(size_t X = start; X < end; X++)

using namespace FX::Filters::Smoothers;

struct Urids
{
    LV2_URID midi_MidiEvent;
};

std::string bundlepath;



struct LV2Plugin;

struct DSP
{
	float sr;
	DSP(float sample_rate) {
		sr = sample_rate;
	}
		
	virtual void activate() {
	
	}
	virtual void deactivate() {
	
	}
	virtual void run(LV2Plugin * plugin, uint32_t sample_count) = 0;
	virtual void connect_port(LV2Plugin * plugin, uint32_t port, void *data_location) =0;
};

DSP * createDSP(float sr);

struct LV2Plugin
{
    float* audio_in_ptr;
    float* audio_out_ptr;
    float* stereo_in_ptr[2];
    float* stereo_out_ptr[2];
    
    const LV2_Atom_Sequence* midi_in_ptr;
    const LV2_Atom_Sequence* midi_out_ptr;    
    LV2_URID_Map* map ;
    Urids urids;
    double rate;
	
	std::array<float*,LAST> controls;
	std::array<float,LAST> control_values;
	std::array<CSmoothFilter*,LAST> smoothers;

	DSP * dsp;
	
    LV2Plugin(const double sampleRate, const LV2_Feature *const *features) :    
    midi_in_ptr (nullptr),
    midi_out_ptr (nullptr),
    audio_in_ptr (nullptr),    
    audio_out_ptr (nullptr),        
    map (nullptr),
    dsp(nullptr),
    rate (sampleRate) 
    {
		const char* missing = lv2_features_query
		(
			features,
			LV2_URID__map, &map, true,
			NULL
		);

		if (missing) throw std::invalid_argument ("Feature map not provided by the host.");

		urids.midi_MidiEvent = map->map (map->handle, LV2_MIDI__MidiEvent);
		stereo_in_ptr[0] = nullptr;
		stereo_in_ptr[1] = nullptr;
		stereo_out_ptr[0] = nullptr;
		stereo_out_ptr[1] = nullptr;
		dsp = createDSP(sampleRate);
		
		LOOP(i,0,LAST) {
			smoothers[i] = new CSmoothFilter(sampleRate,1.0/0.001);
		}	
    }
    ~LV2Plugin() {		
		if(dsp) delete dsp;
		LOOP(i,0,LAST) if(smoothers[i]) delete smoothers[i];
	}
} ;




static LV2_Handle instantiate (const struct LV2_Descriptor *descriptor, double sample_rate, const char *bundle_path, const LV2_Feature *const *features)
{    	
	bundlepath = bundle_path;    
    LV2Plugin * ptr = new LV2Plugin(sample_rate,features);    
    return ptr;
}


static void connect_port (LV2_Handle instance, uint32_t port, void *data_location)
{
    LV2Plugin *m = (LV2Plugin*) instance;
    if (!m) return;
    if (!m->dsp) return;    
	m->dsp->connect_port(m,port,data_location);    
}

static void activate (LV2_Handle instance)
{	
    LV2Plugin * ptr = (LV2Plugin*)instance;
    if(!ptr) return;
    if(!ptr->dsp) return;
    ptr->dsp->activate();
    
}

static void run (LV2_Handle instance, uint32_t sample_count)
{
    LV2Plugin* m = (LV2Plugin*) instance;
    if(!m) return;
    if (!m->dsp) return;
    m->dsp->run(m,sample_count);       
}

static void deactivate (LV2_Handle instance)
{
    LV2Plugin * ptr = (LV2Plugin*)instance;
    if(!ptr) return;
    if (!ptr->dsp) return;
    ptr->dsp->deactivate();
}

static void cleanup (LV2_Handle instance)
{
    LV2Plugin * l = (LV2Plugin*)instance;
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
