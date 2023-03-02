#pragma once

#include "BogAudioDSP.hpp"


namespace DSP::BogAudio
{
	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Equalizer
	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	struct Equalizer : Filter {
		static constexpr DspFloatType gainDb = 12.0f;
		static constexpr DspFloatType cutDb = -36.0f;

		Amplifier _lowAmp;
		Amplifier _midAmp;
		Amplifier _highAmp;
		FourPoleButtworthLowpassFilter _lowFilter;
		TwoPoleButtworthBandpassFilter _midFilter;
		FourPoleButtworthHighpassFilter _highFilter;

		void setParams(
			DspFloatType sampleRate,
			DspFloatType lowDb,
			DspFloatType midDb,
			DspFloatType highDb
		);
		enum {
			PORT_LOWDB,
			PORT_MIDDB,
			PORT_HIGHDB,
		};
		void setPort(int port, DspFloatType v) {
			switch(port) {
				case PORT_LOWDB: setParams(sampleRate,v,_midAmp,_highAmp); break;
				case PORT_MIDDB: setParams(sampleRate,_lowAmp,v,_highAmp); break;
				case PORT_HIGHDB: setParams(sampleRate,_lowAmp,_midAmp,v); break;
			}
		}
		DspFloatType next(DspFloatType sample) override;
		DspFloatType Tick(DspFloatType I=1, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {
			return A*next();
		}
		void ProcessSIMD(size_t n, DspFloatType * in, DspFloatType * out) {
			#pragma omp simd aligned(in,out)
			for(size_t i = 0; i < n; i++) out[i] = Tick(in[i]);
		}
	};

	void Equalizer::setParams(
		DspFloatType sampleRate,
		DspFloatType lowDb,
		DspFloatType midDb,
		DspFloatType highDb
	) {
		assert(lowDb >= cutDb && lowDb <= gainDb);
		assert(midDb >= cutDb && midDb <= gainDb);
		assert(highDb >= cutDb && highDb <= gainDb);

		_lowAmp.setLevel(lowDb);
		_lowFilter.setParams(sampleRate, 100.0f, 0.0f);

		_midAmp.setLevel(midDb);
		_midFilter.setParams(sampleRate, 350.0f, 0.55f, MultimodeFilter::PITCH_BANDWIDTH_MODE);

		_highAmp.setLevel(highDb);
		_highFilter.setParams(sampleRate, 1000.0f, 0.0f);
	}

	DspFloatType Equalizer::next(DspFloatType sample) {
		DspFloatType low = _lowAmp.next(_lowFilter.next(sample));
		DspFloatType mid = _midAmp.next(_midFilter.next(sample));
		DspFloatType high = _highAmp.next(_highFilter.next(sample));
		return low + mid + high;
	}
}
