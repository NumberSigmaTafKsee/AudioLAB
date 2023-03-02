#pragma once

#include "BogAudioDSP.hpp"

namespace DSP::BogAudio
{
	struct DelayLine : public FunctionProcessor
	{
		const DspFloatType _maxTimeMS;
		DspFloatType _sampleRate = -1.0f;
		int _bufferN;
		DspFloatType* _buffer = NULL;
		DspFloatType _time = -1.0f;
		bool _initialized = false;
		int _delaySamples;
		int _leadI;
		int _trailI;
		DspFloatType last_out = 0;

		DelayLine(DspFloatType sampleRate = 1000.0f, DspFloatType maxTimeMS = 1000.0f, DspFloatType time = 1.0f) : 			
		FunctionProcessor(),
		_maxTimeMS(maxTimeMS) {
			setSampleRate(sampleRate);
			setTime(time);
		}
		~DelayLine() {
			delete[] _buffer;
		}

		void setSampleRate(DspFloatType sampleRate);
		void setTime(DspFloatType time);
		DspFloatType next(DspFloatType sample);
		int delaySamples();

		enum {
			PORT_TIME,
		};
		void setPort(int port, DspFloatType v) {
			if(port == PORT_TIME) setTime(v);
		}
		DspFloatType Tick(DspFloatType I, DspFloatType A =1 , DspFloatType X=1, DspFloatType Y=1) {
			int t = _time;
			setTime(t * fabs(X));
			DspFloatType r = next(fabs(Y)*I + (1.0-fabs(Y))*last_out);
			last_out = r;
			setTime(t);
			return A*r;
		}
		void ProcessSIMD(size_t n, DspFloatType * in, DspFloatType * out)
		{
			#pragma omp simd aligned(in,out)
			for(size_t i = 0; i < n; i++)
			{
				const DspFloatType samples = in[i];
				DspFloatType delayed = _buffer[_trailI];
				++_trailI;
				_trailI %= _bufferN;
				_buffer[_leadI] = sample;
				++_leadI;
				_leadI %= _bufferN;
				out[i] = delayed;
			}
		}
	};

  	void DelayLine::setSampleRate(DspFloatType sampleRate) {
		assert(sampleRate > 0.0f);
		if (_sampleRate != sampleRate) {
			_sampleRate = sampleRate;
			if (_buffer) {
				delete[] _buffer;
			}
			_bufferN = ceil((_maxTimeMS / 1000.0f) * _sampleRate);
			_buffer = new DspFloatType[_bufferN] {};
			if (_initialized) {
				_initialized = false;
				setTime(_time);
			}
		}
	}

	void DelayLine::setTime(DspFloatType time) {
		assert(time >= 0.0f);
		assert(time <= 1.0f);
		if (_initialized) {
			if (_time != time) {
				_time = time;
				int newDelaySamples = delaySamples();
				int i = newDelaySamples;
				while (i > _delaySamples) {
					--_trailI;
					if (_trailI < 0) {
						_trailI = _bufferN - 1;
					}
					--i;
				}
				while (i < _delaySamples) {
					++_trailI;
					_trailI %= _bufferN;
					++i;
				}
				_delaySamples = newDelaySamples;
			}
		}
		else {
			_initialized = true;
			_time = time;
			_delaySamples = delaySamples();
			_leadI = 0;
			_trailI = _bufferN - _delaySamples;
		}
	}

	DspFloatType DelayLine::next(DspFloatType sample) {
		DspFloatType delayed = _buffer[_trailI];
		++_trailI;
		_trailI %= _bufferN;
		_buffer[_leadI] = sample;
		++_leadI;
		_leadI %= _bufferN;
		return delayed;
	}

	int DelayLine::delaySamples() {
		return std::max((_time * _maxTimeMS / 1000.0f) * _sampleRate, 1.0f);
	}


}
