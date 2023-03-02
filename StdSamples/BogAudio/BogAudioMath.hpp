#pragma once

#include "BogAudioDSP.hpp"


namespace DSP::BogAudio
{

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Math
	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	struct FastTanhf : public FunctionProcessor  {
		struct TanhfTable : Table {
			TanhfTable(int n) : Table(n) {}
			void _generate() override;
		};
		struct StaticTanhfTable : StaticTable<TanhfTable, 11> {};
		const Table& _table;

		FastTanhf() : 
		FunctionProcessor(),
		_table(StaticTanhfTable::table())	{
		}

		DspFloatType value(DspFloatType radians);

		DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=0, DspFloatType Y=0) {
			return value(M_PI*I);
		}		
		void ProcessSIMD(size_t n, DspFloatType * in, DspFloatType * out) {
			#pragma omp simd aligned(in,out)
			for(size_t i = 0; i < n; i++) out[i] = Tick(in[i]);
		}
	};

	void FastTanhf::TanhfTable::_generate() {
		_table[0] = -1.0f;
		_table[_length - 1] = 1.0f;
		for (int i = 1, n = _length - 1; i < n; ++i) {
			_table[i] = tanhf((((i / (DspFloatType)_length) * 2.0f) - 1.0f) * M_PI);
		}
	}

	DspFloatType FastTanhf::value(DspFloatType radians) {
		if (radians <= -M_PI) {
			return -1.0f;
		}
		if (radians >= M_PI) {
			return 1.0f;
		}
		return _table.value(((radians + M_PI) / (2.0f * M_PI)) * _table.length());
	}
}
