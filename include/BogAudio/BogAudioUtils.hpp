#pragma once

#include "BogAudioDSP.hpp"


namespace DSP::BogAudio
{

	// Utilities
	/////////////////////////////////////////////////////////////////////////////////////////////////////////

	struct DCBlocker : Filter {
		DspFloatType _lastIn = 0.0f;
		DspFloatType _lastOut = 0.0f;

		DspFloatType next(DspFloatType sample) override;
	};

    
	struct AverageRectifiedValue : RunningAverage {
		AverageRectifiedValue(DspFloatType sampleRate = 1000.0f, DspFloatType sensitivity = 1.0f, DspFloatType maxDelayMS = 300.0f)
		: RunningAverage(sampleRate, sensitivity, maxDelayMS)
		{
		}

		DspFloatType next(DspFloatType sample) override;
	};

	struct FastRootMeanSquare : AverageRectifiedValue {
		DCBlocker _dcBlocker;

		FastRootMeanSquare(DspFloatType sampleRate = 1000.0f, DspFloatType sensitivity = 1.0f, DspFloatType maxDelayMS = 300.0f)
		: AverageRectifiedValue(sampleRate, sensitivity, maxDelayMS)
		{
		}

		DspFloatType next(DspFloatType sample) override;
	};



	struct PureRootMeanSquare : RunningAverage {
		PureRootMeanSquare(DspFloatType sampleRate = 1000.0f, DspFloatType sensitivity = 1.0f, DspFloatType maxDelayMS = 300.0f)
		: RunningAverage(sampleRate, sensitivity, maxDelayMS)
		{
		}

		DspFloatType next(DspFloatType sample) override;
	};

	typedef FastRootMeanSquare RootMeanSquare;

	// Puckette 2007, "Theory and Technique"
	struct PucketteEnvelopeFollower {
		DCBlocker _dcBlocker;
		LowPassFilter _filter;

		void setParams(DspFloatType sampleRate, DspFloatType sensitivity);
		DspFloatType next(DspFloatType sample);
	};

	typedef PucketteEnvelopeFollower EnvelopeFollower;

	DspFloatType DCBlocker::next(DspFloatType sample) {
		const DspFloatType r = 0.999f;
		_lastOut = sample - _lastIn + r * _lastOut;
		_lastIn = sample;
		return _lastOut;
	}


	DspFloatType AverageRectifiedValue::next(DspFloatType sample) {
		return RunningAverage::next(fabsf(sample));
	}


	DspFloatType FastRootMeanSquare::next(DspFloatType sample) {
		return AverageRectifiedValue::next(_dcBlocker.next(sample));
	}


	DspFloatType PureRootMeanSquare::next(DspFloatType sample) {
		DspFloatType a = RunningAverage::next(sample * sample);
		if (_sum <= 0.0) {
			return 0.0f;
		}
		return sqrtf(a);
	}


	void PucketteEnvelopeFollower::setParams(DspFloatType sampleRate, DspFloatType sensitivity) {
		const DspFloatType maxCutoff = 10000.0f;
		const DspFloatType minCutoff = 100.0f;
		assert(sensitivity >= 0.0f && sensitivity <= 1.0f);
		DspFloatType cutoff = minCutoff + sensitivity * (maxCutoff - minCutoff);
		_filter.setParams(sampleRate, cutoff);
	}

	DspFloatType PucketteEnvelopeFollower::next(DspFloatType sample) {
		return _filter.next(fabsf(_dcBlocker.next(sample)));
	}
}
