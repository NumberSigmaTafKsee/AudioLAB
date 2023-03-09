#pragma once

namespace Analog::Moog
{
	///////////////////////////////////////////////////////////////////////////////////////////
	// Moog Ladder
	///////////////////////////////////////////////////////////////////////////////////////////
	template<typename DSP>
	class LadderFilterBase
	{
	public:

		LadderFilterBase(DspFloatType sampleRate) : sampleRate(sampleRate) {}
		virtual ~LadderFilterBase() {}

		virtual void Process(size_t n,DspFloatType * samples) = 0;	
		virtual void SetResonance(DspFloatType r) = 0;
		virtual void SetCutoff(DspFloatType c) = 0;

		DspFloatType GetResonance() { return resonance; }
		DspFloatType GetCutoff() { return cutoff; }

	protected:

		DspFloatType cutoff;
		DspFloatType resonance;
		DspFloatType sampleRate;
	};
}
