#pragma once

#include "BogAudioDSP.hpp"


namespace DSP::BogAudio
{


    /////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Resample
	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	struct Decimator {
		Decimator() {}
		virtual ~Decimator() {}

		virtual void setParams(DspFloatType sampleRate, int factor) = 0;
		virtual DspFloatType next(const DspFloatType* buf) = 0;
	};

	struct LPFDecimator : Decimator {
		int _factor;
		MultipoleFilter _filter;

		LPFDecimator(DspFloatType sampleRate = 1000.0f, int factor = 8) {
			setParams(sampleRate, factor);
		}
		void setParams(DspFloatType sampleRate, int factor) override;
		DspFloatType next(const DspFloatType* buf) override;
	};

	struct CICDecimator : Decimator {
		typedef int64_t T;
		static constexpr T scale = ((T)1) << 32;
		int _stages;
		T* _integrators;
		T* _combs;
		int _factor = 0;
		DspFloatType _gainCorrection;

		CICDecimator(int stages = 4, int factor = 8);
		virtual ~CICDecimator();

		void setParams(DspFloatType sampleRate, int factor) override;
		DspFloatType next(const DspFloatType* buf) override;
	};


	struct Interpolator {
		Interpolator() {}
		virtual ~Interpolator() {}

		virtual void setParams(DspFloatType sampleRate, int factor) = 0;
		virtual void next(DspFloatType sample, DspFloatType* buf) = 0;
	};

	struct CICInterpolator : Interpolator {
		typedef int64_t T;
		static constexpr T scale = ((T)1) << 32;
		int _stages;
		T* _integrators;
		T* _combs;
		T* _buffer;
		int _factor = 0;
		DspFloatType _gainCorrection;

		CICInterpolator(int stages = 4, int factor = 8);
		virtual ~CICInterpolator();

		void setParams(DspFloatType sampleRate, int factor) override;
		void next(DspFloatType sample, DspFloatType* buf) override;
	};

	void LPFDecimator::setParams(DspFloatType sampleRate, int factor) {
		_factor = factor;
		_filter.setParams(
			MultipoleFilter::LP_TYPE,
			4,
			factor * sampleRate,
			0.45f * sampleRate,
			0
		);
	}

	DspFloatType LPFDecimator::next(const DspFloatType* buf) {
		DspFloatType s = 0.0f;
		for (int i = 0; i < _factor; ++i) {
			s = _filter.next(buf[i]);
		}
		return s;
	}


	// cascaded integrator–comb decimator: https://en.wikipedia.org/wiki/Cascaded_integrator%E2%80%93comb_filter
	CICDecimator::CICDecimator(int stages, int factor) {
		assert(stages > 0);
		_stages = stages;
		_integrators = new T[_stages + 1] {};
		_combs = new T[_stages] {};
		setParams(0.0f, factor);
	}

	CICDecimator::~CICDecimator() {
		delete[] _integrators;
		delete[] _combs;
	}

	void CICDecimator::setParams(DspFloatType _sampleRate, int factor) {
		assert(factor > 0);
		if (_factor != factor) {
			_factor = factor;
			_gainCorrection = 1.0f / (DspFloatType)(pow(_factor, _stages));
		}
	}

	DspFloatType CICDecimator::next(const DspFloatType* buf) {
		for (int i = 0; i < _factor; ++i) {
			_integrators[0] = buf[i] * scale;
			for (int j = 1; j <= _stages; ++j) {
				_integrators[j] += _integrators[j - 1];
			}
		}
		T s = _integrators[_stages];
		for (int i = 0; i < _stages; ++i) {
			T t = s;
			s -= _combs[i];
			_combs[i] = t;
		}
		return _gainCorrection * (s / (DspFloatType)scale);
	}


	CICInterpolator::CICInterpolator(int stages, int factor) {
		assert(stages > 0);
		_stages = stages;
		_integrators = new T[_stages + 1] {};
		_combs = new T[_stages] {};
		_buffer = NULL;
		setParams(0.0f, factor);
	}

	CICInterpolator::~CICInterpolator() {
		delete[] _integrators;
		delete[] _combs;
		delete[] _buffer;
	}

	void CICInterpolator::setParams(DspFloatType _sampleRate, int factor) {
		assert(factor > 0);
		if (_factor != factor) {
			_factor = factor;
			_gainCorrection = 1.0f / 512.0f; // (DspFloatType)(pow(_factor, _stages / 2));
			if (_buffer) {
				delete[] _buffer;
			}
			_buffer = new T[_factor] {};
		}
	}

	void CICInterpolator::next(DspFloatType sample, DspFloatType* buf) {
		T s = sample * scale;
		for (int i = 0; i < _stages; ++i) {
			T t = s;
			s -= _combs[i];
			_combs[i] = t;
		}

		std::fill(_buffer, _buffer + _factor, (T)0);
		_buffer[0] = s;
		for (int i = 0; i < _factor; ++i) {
			_integrators[0] = _buffer[i];
			for (int j = 1; j <= _stages; ++j) {
				_integrators[j] += _integrators[j - 1];
			}
			buf[i] = _gainCorrection * (_integrators[_stages] / (DspFloatType)scale);
		}
	}
}
