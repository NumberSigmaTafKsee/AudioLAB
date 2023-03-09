#pragma once

#include "BogAudioDSP.hpp"

namespace DSP::BogAudio
{
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Signal
	/////////////////////////////////////////////////////////////////////////////////////////////////////////

	
	struct Amplifier {
		static const DspFloatType minDecibels;
		static const DspFloatType maxDecibels;
		static const DspFloatType decibelsRange;
		struct LevelTable : Table {
			LevelTable(int n) : Table(n) {}
			void _generate() override;
		};
		struct StaticLevelTable : StaticTable<LevelTable, 13> {};

		DspFloatType _db = 0.0f;
		DspFloatType _level;
		const Table& _table;

		Amplifier() : _table(StaticLevelTable::table())	{
			setLevel(minDecibels);
		}

		void setLevel(DspFloatType db);
		DspFloatType next(DspFloatType s);
	};



	struct RunningAverage {
		DspFloatType _maxDelayMS;
		DspFloatType _sampleRate = -1.0f;
		DspFloatType _sensitivity = -1.0f;

		bool _initialized = false;
		DspFloatType* _buffer = NULL;
		int _bufferN = 0;
		int _sumN = 0;
		DspFloatType _invSumN = 0.0f;
		int _leadI = 0;
		int _trailI = 0;
		DspFloatType _sum = 0;

		RunningAverage(DspFloatType sampleRate = 1000.0f, DspFloatType sensitivity = 1.0f, DspFloatType maxDelayMS = 300.0f) : _maxDelayMS(maxDelayMS) {
			setSampleRate(sampleRate);
			setSensitivity(sensitivity);
		}
		virtual ~RunningAverage() {
			if (_buffer) {
				delete[] _buffer;
			}
		}

		void setSampleRate(DspFloatType sampleRate);
		void setSensitivity(DspFloatType sensitivity);
		void reset();
		virtual DspFloatType next(DspFloatType sample);
	};


        
	const DspFloatType Amplifier::minDecibels = -60.0f;
	const DspFloatType Amplifier::maxDecibels = 20.0f;
	const DspFloatType Amplifier::decibelsRange = maxDecibels - minDecibels;

	void Amplifier::LevelTable::_generate() {
		const DspFloatType rdb = 6.0f;
		const DspFloatType tdb = Amplifier::minDecibels + rdb;
		const DspFloatType ta = decibelsToAmplitude(tdb);
		_table[0] = 0.0f;
		for (int i = 1; i < _length; ++i) {
			DspFloatType db = Amplifier::minDecibels + (i / (DspFloatType)_length) * Amplifier::decibelsRange;
			if (db <= tdb) {
				_table[i] = ((db - minDecibels) / rdb) * ta;
			}
			else {
				_table[i] = decibelsToAmplitude(db);
			}
		}
	}

	void Amplifier::setLevel(DspFloatType db) {
		if (_db != db) {
			_db = db;
			if (_db > minDecibels) {
				if (_db < maxDecibels) {
					_level = _table.value(((_db - minDecibels) / decibelsRange) * _table.length());
				}
				else {
					_level = decibelsToAmplitude(_db);
				}
			}
			else {
				_level = 0.0f;
			}
		}
	}

	DspFloatType Amplifier::next(DspFloatType s) {
		return _level * s;
	}


	void RunningAverage::setSampleRate(DspFloatType sampleRate) {
		assert(sampleRate > 0.0f);
		if (_sampleRate != sampleRate) {
			_sampleRate = sampleRate;
			if (_buffer) {
				delete[] _buffer;
			}
			_bufferN = (_maxDelayMS / 1000.0f) * _sampleRate;
			_buffer = new DspFloatType[_bufferN] {};
			if (_initialized) {
				_initialized = false;
				setSensitivity(_sensitivity);
			}
		}
	}

	void RunningAverage::setSensitivity(DspFloatType sensitivity) {
		assert(sensitivity >= 0.0f);
		assert(sensitivity <= 1.0f);
		if (_initialized) {
			if (_sensitivity != sensitivity) {
				_sensitivity = sensitivity;
				int newSumN = std::max(_sensitivity * _bufferN, 1.0f);
				int i = newSumN;
				while (i > _sumN) {
					--_trailI;
					if (_trailI < 0) {
						_trailI = _bufferN - 1;
					}
					_sum += _buffer[_trailI];
					--i;
				}
				while (i < _sumN) {
					_sum -= _buffer[_trailI];
					++_trailI;
					_trailI %= _bufferN;
					++i;
				}
				_sumN = newSumN;
			}
		}
		else {
			_initialized = true;
			_sensitivity = sensitivity;
			_sumN = std::max(_sensitivity * _bufferN, 1.0f);
			_leadI = 0;
			_trailI = _bufferN - _sumN;
			_sum = 0.0;
		}
		_invSumN = 1.0f / (DspFloatType)_sumN;
	}

	void RunningAverage::reset() {
		_sum = 0.0;
		std::fill(_buffer, _buffer + _bufferN, 0.0);
	}

	DspFloatType RunningAverage::next(DspFloatType sample) {
		_sum -= _buffer[_trailI];
		++_trailI;
		_trailI %= _bufferN;
		_sum += _buffer[_leadI] = sample;
		++_leadI;
		_leadI %= _bufferN;
		return (DspFloatType)_sum * _invSumN;
	}
}
