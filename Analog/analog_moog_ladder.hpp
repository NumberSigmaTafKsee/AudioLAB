#pragma once

#include "analog_moog_ladder_filters.hpp"
#include "analog_moog_ladder_ringbuffer.hpp"
#include "analog_moog_ladder_ladderbase.hpp"
#include "analog_moog_ladder_biquad.hpp"
#include "analog_moog_ladder_rbj.hpp"
#include "analog_moog_ladder_houvilainen.hpp"
#include "analog_moog_ladder_improved.hpp"
#include "analog_moog_ladder_krajeski.hpp"
#include "analog_moog_ladder_microtracker.hpp"
#include "analog_moog_ladder_musicdsp.hpp"
#include "analog_moog_ladder_noise.hpp"
#include "analog_moog_ladder_oberheim.hpp"
#include "analog_moog_ladder_others.hpp"
#include "analog_moog_ladder_rungekutta.hpp"
#include "analog_moog_ladder_simplified.hpp"
#include "analog_moog_ladder_stilson.hpp"


namespace Analog::Moog
{
	enum MoogModelType
	{
		FINN_MOOG,
		IMPROVED_MOOG,
		POLISH_MOOG,
		MICROTRACKER_MOOG,
		MUSICDSP_MOOG,
		OBERHEIM_MOOG,
		RK_MOOG,
		SIMPLIFIED_MOOG,
		STILSON_MOOG
	};

	template<typename DSP>
	struct MoogLadderFilter : public FilterProcessor
	{
		LadderFilterBase * moog;
		DspFloatType sr;
		
		MoogLadderFilter(MoogModelType type, DspFloatType sample_rate=44100) : FilterProcessor() {
			sr = sample_rate;
			switch(type) {
			case FINN_MOOG: moog = new HuovilainenMoog(sample_rate); break;
			case IMPROVED_MOOG: moog = new ImprovedMoog(sample_rate); break;
			case POLISH_MOOG: moog = new KrajeskiMoog(sample_rate); break;
			case MICROTRACKER_MOOG: moog = new MicrotrackerMoog(sample_rate); break;
			case MUSICDSP_MOOG: moog = new MusicDSPMoog(sample_rate); break;
			case OBERHEIM_MOOG: moog = new OberheimVariationMoog(sample_rate); break;
			case RK_MOOG: moog = new RKSimulationMoog(sample_rate); break;
			case STILSON_MOOG: moog = new StilsonMoog(sample_rate); break;
			case SIMPLIFIED_MOOG: moog = new SimplifiedMoog(sample_rate); break;
			}		
			assert(moog != nullptr);
		}
		~MoogLadderFilter() {
			if(moog) delete moog;
		}
		void Process(uint32_t n,DspFloatType * samples) {
			moog->Process(n,samples);
		}    

		void SetResonance(DspFloatType r) {		
			moog->SetResonance(r);
		}
		void SetCutoff(DspFloatType c) {		
			moog->SetCutoff(c);
		}
		void setType(MoogModelType type)
		{
			if(moog) delete moog;
			switch(type) {
			case FINN_MOOG: moog = new HuovilainenMoog(sr); break;
			case IMPROVED_MOOG: moog = new ImprovedMoog(sr); break;
			case POLISH_MOOG: moog = new KrajeskiMoog(sr); break;
			case MICROTRACKER_MOOG: moog = new MicrotrackerMoog(sr); break;
			case MUSICDSP_MOOG: moog = new MusicDSPMoog(sr); break;
			case OBERHEIM_MOOG: moog = new OberheimVariationMoog(sr); break;
			case RK_MOOG: moog = new RKSimulationMoog(sr); break;
			case STILSON_MOOG: moog = new StilsonMoog(sr); break;
			case SIMPLIFIED_MOOG: moog = new SimplifiedMoog(sr); break;
			}
			assert(moog != nullptr);
		}
		enum {
			PORT_CUTOFF,
			PORT_RESONANCE,
			PORT_TYPE
		};
		void setPort(int port, DspFloatType v) {
			switch(port) {
				case PORT_CUTOFF: SetCutoff(v); break;
				case PORT_RESONANCE: SetResonance(v); break;
				case PORT_TYPE: setType((MoogModelType)v); break;
			}
		}
		DspFloatType Tick(DspFloatType I,DspFloatType A=1, DspFloatType X=0, DspFloatType Y=0) {
			DspFloatType o = I;
			DspFloatType c = moog->GetCutoff();
			DspFloatType r = moog->GetResonance();
			SetCutoff(clamp(c + c*(X*0.5),0,sr/2));
			SetResonance(clamp(r + Y,0,1));
			Process(1,&o);
			SetCutoff(c);
			SetResonance(r);
			return o;
		}
	};
}
