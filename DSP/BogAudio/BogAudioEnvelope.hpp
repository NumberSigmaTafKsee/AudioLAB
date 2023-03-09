#pragma once

#include "BogAudioDSP.hpp"


namespace DSP::BogAudio
{
	////////////////////////////////////////////////////////////////////////////////////////
	// Envelope
	////////////////////////////////////////////////////////////////////////////////////////
	struct EnvelopeGenerator : Generator {
		DspFloatType _sampleRate = -1.0f;
		DspFloatType _sampleTime;

		EnvelopeGenerator(DspFloatType sampleRate = 1000.0f) : Generator() {
			setSampleRate(std::max(1.0f, sampleRate));
		}

		void setSampleRate(DspFloatType sampleRate);
		virtual void _sampleRateChanged() {}
	};


    
	struct ADSR : EnvelopeGenerator {
		enum Stage {
			STOPPED_STAGE,
			ATTACK_STAGE,
			DECAY_STAGE,
			SUSTAIN_STAGE,
			RELEASE_STAGE
		};

		Stage _stage = STOPPED_STAGE;
		bool _gated = false;
		DspFloatType _attack = 0.0f;
		DspFloatType _decay = 0.0f;
		DspFloatType _sustain = 1.0f;
		DspFloatType _release = 0.0f;
		DspFloatType _attackShape;
		DspFloatType _decayShape;
		DspFloatType _releaseShape;
		DspFloatType _stageProgress = 0.0f;
		DspFloatType _releaseLevel = 0.0f;
		DspFloatType _envelope = 0.0f;

		ADSR(bool linear = false, DspFloatType sampleRate = 1000.0f) : EnvelopeGenerator(sampleRate) {
			setLinearShape(linear);
		}

		void reset();
		void setGate(bool high);
		void setAttack(DspFloatType seconds);
		void setDecay(DspFloatType seconds);
		void setSustain(DspFloatType level);
		void setRelease(DspFloatType seconds);
		void setLinearShape(bool linear);
		void setShapes(DspFloatType attackShape, DspFloatType decayShape, DspFloatType releaseShape);
		bool isStage(Stage stage) { return _stage == stage; }
		void retrigger();

		enum {
			PORT_RESET,
			PORT_GATE,
			PORT_ATTACK,
			PORT_DECAY,
			PORT_SUSTAIN,
			PORT_RELEASE,
			PORT_LINEARSHAPE,
			PORT_ATTACKSHAPE,
			PORT_DECAYSHAPE,
			PORT_RELEASESHAPE,
			PORT_RETRUGGER,
		};
		void setPort(int port, DspFloatType v) {
			switch(port) {
				case PORT_RESET: reset(); break;
				case PORT_GATE: setGate((bool)v); break;
				case PORT_ATTACK: setAttack(v); break;
				case PORT_DECAY: setDecay(v); break;
				case PORT_SUSTAIN: setSustain(v); break;
				case PORT_RELEASE: setRelease(v); break;
				case PORT_LINEARSHAPE: setLinearShape((bool)v); break;
				case PORT_ATTACKSHAPE: setShapes(v,_decayShape,_releaseShape); break;
				case PORT_DECAYSHAPE: setShapes(_attackShape,v,_releaseShape); break;
				case PORT_RELEASESHAPE: setShapes(_attackShape,_decayShape,v); break;

			}
		}
		DspFloatType _next() override;
		
		DspFloatType Tick(DspFloatType I=1, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {
			return A*_next();
		}
		void ProcessSIMD(size_t n, DspFloatType * in, DspFloatType * out) {
			#pragma omp simd aligned(in,out)
			for(size_t i = 0; i < n; i++) out[i] = Tick(in[i]);
		}
	};


	void EnvelopeGenerator::setSampleRate(DspFloatType sampleRate) {
		assert(sampleRate >= 1.0f);
		if (_sampleRate != sampleRate) {
			_sampleRate = sampleRate;
			_sampleTime = 1.0f / sampleRate;
			_sampleRateChanged();
		}
	}


	void ADSR::reset() {
		_stage = STOPPED_STAGE;
		_gated = false;
		_envelope = 0.0f;
	}

	void ADSR::setGate(bool high) {
		_gated = high;
	}

	void ADSR::setAttack(DspFloatType seconds) {
		assert(_attack >= 0.0f);
		_attack = std::max(seconds, 0.001f);
	}

	void ADSR::setDecay(DspFloatType seconds) {
		assert(_decay >= 0.0f);
		_decay = std::max(seconds, 0.001f);
	}

	void ADSR::setSustain(DspFloatType level) {
		assert(_sustain >= 0.0f);
		assert(_sustain <= 1.0f);
		_sustain = level;
	}

	void ADSR::setRelease(DspFloatType seconds) {
		assert(_release >= 0.0f);
		_release = std::max(seconds, 0.001f);
	}

	void ADSR::setLinearShape(bool linear) {
		if (linear) {
			setShapes(1.0f, 1.0f, 1.0f);
		}
		else {
			setShapes(0.5f, 2.0f, 2.0f);
		}
	}

	void ADSR::setShapes(DspFloatType attackShape, DspFloatType decayShape, DspFloatType releaseShape) {
		assert(attackShape >= 0.1f && attackShape <= 10.0f);
		assert(decayShape >= 0.0f && decayShape <= 10.0f);
		assert(releaseShape >= 0.0f && releaseShape <= 10.0f);
		_attackShape = attackShape;
		_decayShape = decayShape;
		_releaseShape = releaseShape;
	}

	void ADSR::retrigger() {
		switch (_stage) {
			case STOPPED_STAGE: {
				_stage = ATTACK_STAGE;
				_stageProgress = 0.0f;
				break;
			}
			default: {
				_stage = ATTACK_STAGE;
				DspFloatType e = powf(_envelope, 1.0f / _attackShape);
				_stageProgress = e * _attack;
			}
		}
	}

	DspFloatType ADSR::_next() {
		if (_gated) {
			switch (_stage) {
				case STOPPED_STAGE: {
					_stage = ATTACK_STAGE;
					_stageProgress = 0.0f;
					break;
				}
				case ATTACK_STAGE: {
					if (_envelope >= 1.0) {
						_stage = DECAY_STAGE;
						_stageProgress = 0.0f;
					}
					break;
				}
				case DECAY_STAGE: {
					if (_stageProgress >= _decay) {
						_stage = SUSTAIN_STAGE;
						_stageProgress = 0.0f;
					}
					break;
				}
				case SUSTAIN_STAGE: {
					break;
				}
				case RELEASE_STAGE: {
					_stage = ATTACK_STAGE;
					_stageProgress = _attack * powf(_envelope, _releaseShape);
					break;
				}
			}
		}
		else {
			switch (_stage) {
				case STOPPED_STAGE: {
					break;
				}
				case ATTACK_STAGE:
				case DECAY_STAGE:
				case SUSTAIN_STAGE: {
					_stage = RELEASE_STAGE;
					_stageProgress = 0.0f;
					_releaseLevel = _envelope;
					break;
				}
				case RELEASE_STAGE: {
					if (_stageProgress >= _release) {
						_stage = STOPPED_STAGE;
					}
					break;
				}
			}
		}

		switch (_stage) {
			case STOPPED_STAGE: {
				_envelope = 0.0f;
				break;
			}
			case ATTACK_STAGE: {
				_stageProgress += _sampleTime;
				_envelope = std::min(1.0f, _stageProgress / _attack);
				_envelope = powf(_envelope, _attackShape);
				break;
			}
			case DECAY_STAGE: {
				_stageProgress += _sampleTime;
				_envelope = std::min(1.0f, _stageProgress / _decay);
				_envelope = powf(1.0f - _envelope, _decayShape);
				_envelope *= 1.0f - _sustain;
				_envelope += _sustain;
				break;
			}
			case SUSTAIN_STAGE: {
				_envelope = _sustain;
				break;
			}
			case RELEASE_STAGE: {
				_stageProgress += _sampleTime;
				_envelope = std::min(1.0f, _stageProgress / _release);
				_envelope = powf(1.0f - _envelope, _releaseShape);
				_envelope *= _releaseLevel;
				break;
			}
		}

		return _envelope;
	}
}
