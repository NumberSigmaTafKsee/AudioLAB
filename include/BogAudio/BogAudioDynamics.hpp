#pragma once

#include "BogAudioDSP.hpp"

namespace DSP::BogAudio
{
	struct Limiter : public FunctionProcessor {
		DspFloatType _shape;
		DspFloatType _knee;
		DspFloatType _limit;
		DspFloatType _scale;
		DspFloatType _normalization;
		FastTanhf _tanhf;

		Limiter() : FunctionProcessor() {}

		void setParams(DspFloatType shape = 1.0f, DspFloatType knee = 5.0f, DspFloatType limit = 10.0f, DspFloatType scale = 2.0f);
		DspFloatType next(DspFloatType sample);

		enum {
			PORT_SHAPE,
			PORT_KNEE,
			PORT_LIMIT,
			PORT_SCALE,		
		};
		void setPort(int port, DspFloatType v) {
			switch(port) {
				case PORT_SHAPE: setParams(v,_knee,_limit,_scale); break;
				case PORT_KNEE: setParams(_shape,v,_limit,_scale); break;
				case PORT_LIMIT: setParams(_shape,_knee,v,_scale); break;
				case PORT_SCALE: setParams(_shape,_knee,_limit,v); break;
			}
		}
		DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {
			DspFloatType s = shape;
			DspFloatType k = knee;
			setParams(s*fabs(X),k*fabs(Y));
			DspFloatType r = next(I);
			setParams(s,k);
			return A*r;
		}
		void ProcessSIMD(size_t n, DspFloatType * in, DspFloatType * out) {
			#pragma omp simd aligned(in,out)
			for(size_t i = 0; i < n; i++) out[i] = Tick(in[i]);
		}
	};

    
	struct Saturator : public FunctionProcessor {
		static const DspFloatType limit;

		Saturator() : FunctionProcessor() {

		}
		DspFloatType next(DspFloatType sample);

		DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {
			return A*next(X*Y*I);
		}
	};

	struct Compressor {
		static const DspFloatType maxEffectiveRatio;
		DspFloatType _detectorDb,_thresholdDb,_ratio,_softknee;
		Compressor() {
			_detectorDb = _thresholdDb = _ratio = _softknee =0.0;
		}
		DspFloatType compressionDb(DspFloatType detectorDb, DspFloatType thresholdDb, DspFloatType ratio, bool softKnee);
		enum {
			PORT_DETECTORDB,
			PORT_THRESHDB,
			PORT_RATIO,
			PORT_SOFTKNEE,
		};
		void setPort(int port, DspFloatType v) {
			switch(port) {
				case PORT_DETECTORDB: _detectorDb = v; break;
				case PORT_THRESHDB: _thresholdDb = v; break;
				case PORT_RATIO: _ratio = v; break;
				case PORT_SOFTKNEE: _softknee = v; break;
			}
		}
		DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {
			DspFloatType r = _ratio;
			DspFloatType k = _softknee;
			DspFloatType g = compressionDb(_detectorDb,_thresholdDb,r*fabs(X),k*fabs(Y));
			return A*(g*I);
		}
		void ProcessSIMD(size_t n, DspFloatType * in, DspFloatType * out) {
			#pragma omp simd aligned(in,out)
			for(size_t i = 0; i < n; i++) out[i] = Tick(in[i]);
		}
	};

	struct NoiseGate {
		static const DspFloatType maxEffectiveRatio;
		DspFloatType _detectorDb,_thresholdDb,_ratio,_softknee;
		NoiseGate() {
			_detectorDb = _thresholdDb = _ratio = _softknee =0.0;
		}
		DspFloatType compressionDb(DspFloatType detectorDb, DspFloatType thresholdDb, DspFloatType ratio, bool softKnee);

		enum {
			PORT_DETECTORDB,
			PORT_THRESHDB,
			PORT_RATIO,
			PORT_SOFTKNEE,
		};
		void setPort(int port, DspFloatType v) {
			switch(port) {
				case PORT_DETECTORDB: _detectorDb = v; break;
				case PORT_THRESHDB: _thresholdDb = v; break;
				case PORT_RATIO: _ratio = v; break;
				case PORT_SOFTKNEE: _softknee = v; break;
			}
		}
		DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {
			DspFloatType r = _ratio;
			DspFloatType k = _softknee;
			DspFloatType g = compressionDb(_detectorDb,_thresholdDb,r*fabs(X),k*fabs(Y));
			return A*(g*I);
		}
		void ProcessSIMD(size_t n, DspFloatType * in, DspFloatType * out) {
			#pragma omp simd aligned(in,out)
			for(size_t i = 0; i < n; i++) out[i] = Tick(in[i]);
		}
	};

    void Limiter::setParams(DspFloatType shape, DspFloatType knee, DspFloatType limit, DspFloatType scale) {
		assert(shape >= 0.0f);
		assert(knee >= 0.0f);
		assert(limit >= 0.0f);
		assert(scale >= 1.0f);
		_shape = shape;
		_knee = knee;
		_limit = std::max(knee, limit);
		_scale = scale;

		if (_shape >= 0.1f) {
			if (_shape < 1.0f) {
				_normalization = 1.0f / tanhf(_shape * M_PI);
			}
			else {
				_normalization = 1.0f;
			}
		}
	}

	DspFloatType Limiter::next(DspFloatType sample) {
		DspFloatType out = fabsf(sample);
		if (out > _knee) {
			out -= _knee;
			out /= _scale;
			if (_shape >= 0.1f) {
				// out /= _limit - _knee;
				// out = _tanhf.value(out * _shape * M_PI) * _normalization;
				// out *= _limit - _knee;
				DspFloatType x = out / (_limit - _knee);
				x = _tanhf.value(x * _shape * M_PI) * _normalization;
				x = std::min(x, 1.0f);
				x *= _limit - _knee;
				out = std::min(fabsf(sample) - _knee, x);
			}
			else {
				out = std::min(out, _limit - _knee);
			}
			out += _knee;
		}

		if (sample < 0.0f) {
			return -out;
		}
		return out;
	}


	const DspFloatType Saturator::limit = 12.0f;

	// Zavalishin 2018, "The Art of VA Filter Design", http://www.native-instruments.com/fileadmin/ni_media/downloads/pdf/VAFilterDesign_2.0.0a.pdf
	static inline DspFloatType saturation(DspFloatType x) {
		const DspFloatType y1 = 0.98765f; // (2*x - 1)/x**2 where x is 0.9.
		const DspFloatType offset = 0.075f / Saturator::limit; // magic.
		DspFloatType x1 = (x + 1.0f) * 0.5f;
		return Saturator::limit * (offset + x1 - sqrtf(x1 * x1 - y1 * x) * (1.0f / y1));
	}

	DspFloatType Saturator::next(DspFloatType sample) {
		DspFloatType x = sample * (1.0f / limit);
		if (sample < 0.0f) {
			return -saturation(-x);
		}
		return saturation(x);
	}


	const DspFloatType Compressor::maxEffectiveRatio = 1000.0f;

	DspFloatType Compressor::compressionDb(DspFloatType detectorDb, DspFloatType thresholdDb, DspFloatType ratio, bool softKnee) {
		const DspFloatType softKneeDb = 3.0f;

		if (softKnee) {
			DspFloatType sDb = thresholdDb - softKneeDb;
			if (detectorDb <= sDb) {
				return 0.0f;
			}

			DspFloatType ix = softKneeDb * std::min(ratio, maxEffectiveRatio) + thresholdDb;
			DspFloatType iy = softKneeDb + thresholdDb;
			DspFloatType t = (detectorDb - sDb) / (ix - thresholdDb);
			DspFloatType px = t * (ix - thresholdDb) + thresholdDb;
			DspFloatType py = t * (iy - thresholdDb) + thresholdDb;
			DspFloatType s = (py - sDb) / (px - sDb);
			DspFloatType compressionDb = detectorDb - sDb;
			compressionDb -= s * (detectorDb - sDb);
			return compressionDb;
		}

		if (detectorDb <= thresholdDb) {
			return 0.0f;
		}
		DspFloatType compressionDb = detectorDb - thresholdDb;
		compressionDb -= compressionDb / ratio;
		return compressionDb;
	}


	const DspFloatType NoiseGate::maxEffectiveRatio = Compressor::maxEffectiveRatio;

	DspFloatType NoiseGate::compressionDb(DspFloatType detectorDb, DspFloatType thresholdDb, DspFloatType ratio, bool softKnee) {
		const DspFloatType softKneeDb = 6.0f;

		if (softKnee) {
			// FIXME: this achieves nothing.
			DspFloatType range = thresholdDb - Amplifier::minDecibels;
			DspFloatType ix = thresholdDb + softKneeDb;
			DspFloatType iy = 0;
			if (detectorDb >= ix) {
				return 0.0f;
			}
			DspFloatType ox = thresholdDb - range / ratio;
			if (detectorDb <= ox) {
				return -Amplifier::minDecibels;
			}
			const DspFloatType oy = Amplifier::minDecibels;
			DspFloatType t = (detectorDb - ox) / (ix - ox);
			DspFloatType px = t * (ix - thresholdDb) + thresholdDb;
			DspFloatType py = t * (iy - thresholdDb) + thresholdDb;
			DspFloatType s = (py - oy) / (px - ox);
			return -(oy + s * (detectorDb - ox));
		}

		if (detectorDb >= thresholdDb) {
			return 0.0f;
		}
		DspFloatType differenceDb = thresholdDb - detectorDb;
		DspFloatType compressionDb = differenceDb * ratio - differenceDb;
		return std::min(compressionDb, -Amplifier::minDecibels);
	}
}
