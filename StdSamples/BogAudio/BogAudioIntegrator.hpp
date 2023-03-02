#pragma once

#include "BogAudioDSP.hpp"

namespace DSP::BogAudio
{


	struct Integrator : public FunctionProcessor{
		DspFloatType _alpha = 0.0f;
		DspFloatType _last = 0.0f;

		Integrator(DspFloatType alpha = 1.0f) : FunctionProcessor() {
			setParams(alpha);
		}

		void setParams(DspFloatType alpha);
		DspFloatType next(DspFloatType sample);

		DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=1, DspFloatType  Y=1) {
			return next(I);
		}		
		void ProcessSIMD(size_t n, DspFloatType * in, DspFloatType * out) {
			#pragma omp simd aligned(in,out)
			for(size_t i = 0; i < n; i++) out[i] = Tick(in[i]);
		}
	};

 	void Integrator::setParams(DspFloatType alpha) {
		assert(alpha >= 0.0f);
		assert(alpha <= 1.0f);
		_alpha = alpha;
	}

	DspFloatType Integrator::next(DspFloatType sample) {
		// "leaky integrator"
		return _last = (1.0f - _alpha)*_last + _alpha*sample;
	}

}
