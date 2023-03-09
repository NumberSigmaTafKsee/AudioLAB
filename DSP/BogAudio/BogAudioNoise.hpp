#pragma once

#include "BogAudioDSP.hpp"


namespace DSP::BogAudio
{
	////////////////////////////////////////////////////////////////////////////////////////
	// Noise
	////////////////////////////////////////////////////////////////////////////////////////
	class Seeds {
	private:
		std::mt19937 _generator;
		Seeds();
		unsigned int _next();

	public:
		Seeds(const Seeds&) = delete;
		void operator=(const Seeds&) = delete;
		static Seeds& getInstance();

		static unsigned int next();
	};

	struct NoiseGenerator : Generator {
		std::minstd_rand _generator; // one of the faster options.

		NoiseGenerator() : _generator(Seeds::next()) {}
	};

	struct WhiteNoiseGenerator : NoiseGenerator {
		std::uniform_real_distribution<DspFloatType> _uniform;

		WhiteNoiseGenerator() : _uniform(-1.0, 1.0) {}

		DspFloatType _next() override {
			return _uniform(_generator);
		}
	};


    

	template<typename G>
	struct BasePinkNoiseGenerator : NoiseGenerator {
		static const int _n = 7;
		G _g;
		G _gs[_n];
		uint32_t _count = _g.next();

		DspFloatType _next() override {
			// See: http://www.firstpr.com.au/dsp/pink-noise/
			DspFloatType sum = _g.next();
			for (int i = 0, bit = 1; i < _n; ++i, bit <<= 1) {
				if (_count & bit) {
					sum += _gs[i].next();
				}
				else {
					sum += _gs[i].current();
				}
			}
			++_count;
			return sum / (DspFloatType)(_n + 1);
		}
	};

	struct PinkNoiseGenerator : BasePinkNoiseGenerator<WhiteNoiseGenerator> {};

	struct RedNoiseGenerator : BasePinkNoiseGenerator<PinkNoiseGenerator> {};

	struct BlueNoiseGenerator : NoiseGenerator {
		PinkNoiseGenerator _pink;
		DspFloatType _last = 0.0f;

		DspFloatType _next() override {
			DspFloatType t = _last;
			_last = _pink.next();
			return _last - t;
		}
	};



	struct GaussianNoiseGenerator : NoiseGenerator {
		std::normal_distribution<DspFloatType> _normal;

		GaussianNoiseGenerator(DspFloatType mean = 0.0f, DspFloatType stdDev = 1.0f) : _normal(mean, stdDev) {}

		DspFloatType _next() override {
			return _normal(_generator);
		}
	};


	struct RandomWalk : Generator {
		DspFloatType _min;
		DspFloatType _max;
		DspFloatType _last = 0.0f;
		DspFloatType _lastOut = 0.0f;
		DspFloatType _damp;
		DspFloatType _bias = 0.0f;
		DspFloatType _biasDamp = 1.0f;
		WhiteNoiseGenerator _noise;
		LowPassFilter _filter;

		RandomWalk(
			DspFloatType min = -5.0f,
			DspFloatType max = 5.0f,
			DspFloatType sampleRate = 1000.0f,
			DspFloatType change = 0.5f
		)
		: _min(min)
		, _max(max)
		{
			assert(_min < _max);
			setParams(sampleRate, change);
		}

		void setParams(DspFloatType sampleRate = 1000.0f, DspFloatType change = 0.5f);
		void jump();
		void tell(DspFloatType v);
		DspFloatType _next() override;
		
		DspFloatType Tick(DspFloatType I=1, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {
			return A*next();
		}
		void ProcessSIMD(size_t n, DspFloatType * in, DspFloatType * out) {
			#pragma omp simd aligned(in,out)
			for(size_t i = 0; i < n; i++) out[i] = Tick(in[i]);
		}
	};

    Seeds::Seeds() {
#ifdef ARCH_WIN
  _generator.seed(time(0));
#else
  std::random_device rd;
  _generator.seed(rd());
#endif
}

unsigned int Seeds::_next() {
  return _generator();
}

Seeds& Seeds::getInstance() {
  static Seeds instance;
  return instance;
}

unsigned int Seeds::next() {
  return getInstance()._next();
};


void RandomWalk::setParams(DspFloatType sampleRate, DspFloatType change) {
	assert(sampleRate > 0.0f);
	assert(change >= 0.0f);
	assert(change <= 1.0f);

	_filter.setParams(sampleRate, std::max(2.0f, change * 0.49f * std::min(44100.0f, sampleRate)));

	const DspFloatType maxDamp = 0.98;
	const DspFloatType minDamp = 0.9999;
	_damp = maxDamp + (1.0f - change)*(minDamp - maxDamp);

	_biasDamp = 1.0f - change*(2.0f / sampleRate);
}

void RandomWalk::jump() {
	DspFloatType x = fabsf(_noise.next()) * (_max - _min);
	x += _min;
	tell(x);
}

void RandomWalk::tell(DspFloatType v) {
	assert(v >= _min && v <= _max);
	_last = _bias = v;
	_filter.reset();
}

DspFloatType RandomWalk::_next() {
	DspFloatType delta = _noise.next();
	if ((_lastOut >= _max && delta > 0.0f) || (_lastOut <= _min && delta < 0.0f)) {
		delta = -delta;
	}
	_last = _damp*_last + delta;
	_bias *= _biasDamp;
	return _lastOut = std::min(std::max(_bias + _filter.next(_last), _min), _max);
}
}
