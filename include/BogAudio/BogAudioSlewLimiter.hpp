#pragma once

#include "BogAudioDSP.hpp"

namespace DSP::BogAudio
{


	struct SlewLimiter {
		DspFloatType _delta;
		DspFloatType _last = 0.0f;

		SlewLimiter(DspFloatType sampleRate = 1000.0f, DspFloatType milliseconds = 1.0f, DspFloatType range = 10.0f) {
			setParams(sampleRate, milliseconds, range);
		}

		void setParams(DspFloatType sampleRate = 1000.0f, DspFloatType milliseconds = 1.0f, DspFloatType range = 10.0f);
		inline void setLast(DspFloatType last) { _last = last; }
		inline DspFloatType next(DspFloatType sample) {
			return _last = next(sample, _last);
		}
		DspFloatType next(DspFloatType sample, DspFloatType last);
	};

	struct ShapedSlewLimiter {
		const DspFloatType range = 10.0f;
		const DspFloatType minShape = 0.1f;
		const DspFloatType maxShape = 5.0f;
		DspFloatType _sampleTime;
		DspFloatType _time;
		DspFloatType _shapeExponent;
		DspFloatType _inverseShapeExponent;
		DspFloatType _last = 0.0;

		ShapedSlewLimiter(DspFloatType sampleRate = 1000.0f, DspFloatType milliseconds = 1.0f, DspFloatType shape = 1.0f) {
			setParams(sampleRate, milliseconds, shape);
		}

		void setParams(DspFloatType sampleRate, DspFloatType milliseconds, DspFloatType shape);
		DspFloatType next(DspFloatType sample);
	};

    	void SlewLimiter::setParams(DspFloatType sampleRate, DspFloatType milliseconds, DspFloatType range) {
		assert(sampleRate > 0.0f);
		assert(milliseconds >= 0.0f);
		assert(range > 0.0f);
		_delta = range / ((milliseconds / 1000.0f) * sampleRate);
	}

	DspFloatType SlewLimiter::next(DspFloatType sample, DspFloatType last) {
		if (sample > last) {
			return std::min(last + _delta, sample);
		}
		return std::max(last - _delta, sample);
	}


	void ShapedSlewLimiter::setParams(DspFloatType sampleRate, DspFloatType milliseconds, DspFloatType shape) {
		assert(sampleRate > 0.0f);
		assert(milliseconds >= 0.0f);
		assert(shape >= minShape);
		assert(shape <= maxShape);
		_sampleTime = 1.0f / sampleRate;
		_time = milliseconds / 1000.0f;
		_shapeExponent = (shape > -0.05f && shape < 0.05f) ? 0.0f : shape;
		_inverseShapeExponent = 1.0f / _shapeExponent;
	}

	DspFloatType ShapedSlewLimiter::next(DspFloatType sample) {
		DspFloatType difference = sample - _last;
		DspFloatType ttg = fabsf(difference) / range;
		if (_time < 0.0001f || ttg < 0.0001f) {
			return _last = sample;
		}
		if (_shapeExponent != 0.0f) {
			ttg = powf(ttg, _shapeExponent);
		}
		ttg *= _time;
		ttg = std::max(0.0f, ttg - _sampleTime);
		ttg /= _time;
		if (_shapeExponent != 0.0f) {
			ttg = powf(ttg, _inverseShapeExponent);
		}
		DspFloatType y = fabsf(difference) - ttg * range;
		if (difference < 0.0f) {
			return _last = std::max(_last - y, sample);
		}
		return _last = std::min(_last + y, sample);
	}

}
