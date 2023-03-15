//2PACRBJ12
// 12-band RBJ equalizer

#include "2PACRBJ12.h"
#include "LV2Plugin.hpp"

#include "Filters/IIRRBJFilters.hpp"

namespace Filters::IIR::RBJFilters

const int NUM_FILTERS = 12;

struct LV22PACRBJ12 : public DSP
{
	std::array<RBJFilter*,NUM_FILTERS> filters;
	const float lowshelf_frequency = 100.0;
	const float hishelf_frequency  = 5280.0;
	LV22PACRBJ12(float sr) : DSP(sr)
	{
		filters[0] = new RBJFilter(RBJFilter::LOWSHELF,sr);
		filters[0]->setFrequency(lowshelf_frequency);
		
		LOOP(i,1,NUM_FILTERS-1) {
			double freq = 440.0*(double)i;
			filters[i] = new RBJFilter(RBJFilter::PEAK,sr);		
		}
		filters[NUM_FILTERS-1] = new RBJFilter(RBJFilter::HIGHSHELF,sr);
	}
	~LV22PACRBJ12() {
		LOOP(i,0,NUM_FILTERS) {
			if(filters[i]) delete filters[i];
		}
	}
	
	void run(LV2Plugin * m, int sample_count)
	{		
		if(!m) return;   
		if(!m->dsp) return;
		if(!m->audio_in_ptr) return;
		if(!m->audio_out_ptr) return; 
		if(!m->midi_in_ptr) return;
		if(!moog) return;
		
		LOOP(i,0,LAST) {
			if(*m->controls[i] != m->control_values[i]) {
				m->control_values[i] = *m->controls[i];
				m->smoothers[i]->setTarget(m->control_values[i]);
			}
		}
				
		if(parallel)
		{
			std::array<float,sample_count> temp;
			std::array<float,sample_count> output;
			memset(output.data(),0,sample_count*sizeof(float));
			LOOP(i,0,NUM_FILTERS)
			
		
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
			
		default:
			break;
		}
}

};


DSP* createDSP(float sr) {	
	return new LV2LiquidMoog(sr);
}

