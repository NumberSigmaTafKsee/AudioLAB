
#pragma once

#include "BogAudioDSP.hpp"

namespace DSP::BogAudio
{


    struct CrossFader : public Parameter2Processor {
		DspFloatType _mix = 2.0f;
		DspFloatType _curve = 1.0f;
		bool _linear = true;
		DspFloatType _aMix;
		DspFloatType _bMix;
		Amplifier _aAmp;
		Amplifier _bAmp;

		CrossFader() : Parameter2Processor() {
			setParams(0.0f);
		}

		void setParams(
			DspFloatType mix, // -1 to 1, 0 for equal output of both inputs.
			DspFloatType curve = 1.0f, // -1 to 1: at -1, A will cut fully as mix goes to 0; at 0, A cuts over full mix; at 1, A cuts from 0 to 1.  B symmetric.
			bool linear = true// cut is linear in amplitude if true; linear in decibels otherwise.
		);
		DspFloatType next(DspFloatType a, DspFloatType b);
		enum {
			PORT_MIX,
			PORT_CURVE,
			PORT_LINEAR
		};
		void setPort(int port, DspFloatType v) {
			switch(port) {
				case PORT_MIX: setParams(v); break;
				case PORT_CURVE: setParams(_mix,v); break;
				case PORT_LINEAR: setParams(_mix,_curve,(bool)v); break;
			}
		}
		DspFloatType Tick(DspFloatType a, DspFloatType b) {			
			return next(a,b);
		}
		void ProcessSIMD(size_t n, DspFloatType * a, DspFloatType * b, DspFloatType * out)
		{
			#pragma omp simd aligned(a,b,out)
			for(size_t i = 0; i < n; i++) {
				out[i] = _linear? _aMix * a + _bMix * b : _aAmp.next(a) + _bAmp.next(b);
			}
		}
		
	};

	struct Panner : public StereoSplitterProcessor {
		DspFloatType _pan = 2.0f;
		DspFloatType _lLevel = 0.0f;
		DspFloatType _rLevel = 0.0f;
		const Table& _sineTable;

		Panner() : 
		StereoSplitterProcessor(),
		_sineTable(StaticSineTable::table()) {
			setPan(0.0f);
		}

		void setPan(DspFloatType pan); // -1.0 full left, 0.0 even, 1.0f full right.
		void next(DspFloatType sample, DspFloatType& l, DspFloatType& r);

		enum {
			PORT_PAN,			
		};
		void setPort(int port, DspFloatType v) {
			switch(port) {
				case PORT_PAN: setPan(v); break;
			}
		}
		DspFloatType Tick(DspFloatType in, DspFloatType &L, DspFloatType &R) {
			DspFloatType l=0;
			DspFloatType r=0;
			next(in,l,r);
			L = l;
			R = r;
			return 0.5*(l+r);
		}
	};

    void CrossFader::setParams(DspFloatType mix, DspFloatType curve, bool linear) {
		assert(mix >= -1.0f && mix <= 1.0f);
		assert(curve >= -1.0f && curve <= 1.0f);
		if (_mix != mix || _curve != curve || _linear != linear) {
			_mix = mix;
			_curve = curve;
			_linear = linear;

			DspFloatType aMax, aMin;
			DspFloatType bMax, bMin;
			if (_curve < 0.0f) {
				aMax = 0.0f;
				aMin = _curve + 2.0f;
				bMax = 2.0f;
				bMin = 0.0f - _curve;
			}
			else {
				aMax = _curve;
				aMin = 2.0f;
				bMax = 2.0f - _curve;
				bMin = 0.0f;
			}

			DspFloatType m = _mix + 1.0f;
			if (m < aMax) {
				_aMix = 1.0f;
			}
			else if (m > aMin) {
				_aMix = 0.0f;
			}
			else {
				_aMix = 1.0f - ((m - aMax) / (aMin - aMax));
			}

			if (m > bMax) {
				_bMix = 1.0f;
			}
			else if (m < bMin) {
				_bMix = 0.0f;
			}
			else {
				_bMix = (m - bMin) / (bMax - bMin);
			}

			if (!_linear) {
				_aAmp.setLevel((1.0f - _aMix) * Amplifier::minDecibels);
				_bAmp.setLevel((1.0f - _bMix) * Amplifier::minDecibels);
			}
		}
	}

	DspFloatType CrossFader::next(DspFloatType a, DspFloatType b) {
		if (_linear) {
			return _aMix * a + _bMix * b;
		}
		return _aAmp.next(a) + _bAmp.next(b);
	}


	void Panner::setPan(DspFloatType pan) {
		assert(pan >= -1.0f);
		assert(pan <= 1.0f);
		if (_pan != pan) {
			_pan = pan;
			_lLevel = _sineTable.value(((1.0f + _pan) / 8.0f + 0.25f) * _sineTable.length());
			_rLevel = _sineTable.value(((1.0f + _pan) / 8.0f) * _sineTable.length());
		}
	}

	void Panner::next(DspFloatType sample, DspFloatType& l, DspFloatType& r) {
		l = _lLevel * sample;
		r = _rLevel * sample;
	}


}
