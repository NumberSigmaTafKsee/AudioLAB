#pragma once
#include "iir_filters.hpp"
#include "iir_butterworth_filters.hpp"

namespace iir_filters
{
	struct ButterworthBandstopFilter
	{
		int order;
		DspFloatType sampleRate,frequency,q;    
		Filters::BiquadFilterCascade filter;

		ButterworthBandstopFilter(int order, DspFloatType sr)
		{
			this->order = order;
			this->sampleRate = sr;        
			setFilter(1000.0,0.5);
		}
		void setFilter(DspFloatType f, DspFloatType Q) {
			auto x = Filters::butterlp2bs(order,Q);
			frequency = f;
			q = Q;
			auto c = AnalogBiquadCascade(x,f,sampleRate);
			filter.setCoefficients(c);
		}
		void setCutoff(DspFloatType f) {
			if(f <= 0) return;
			if(f >= sampleRate/2) return;
			setFilter(f,q);
		}
		void setQ(DspFloatType Q) {
			if(Q < 0.5) Q = 0.5;
			if(Q > 999) Q = 999;
			setFilter(frequency,Q);
		}
		DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {
			return filter.Tick(I,A,X,Y);
		}
		void ProcessBlock(size_t n, float * in, float * out) {		
			for(size_t i = 0; i < n; i++) in[i] = Tick(out[i]);
		}
		std::vector<DspFloatType> impulse_response(size_t n)
		{
			std::vector<DspFloatType> r(n);
			r[0] = Tick(1.0);
			for(size_t i = 1; i < n; i++) r[i] = Tick(0);
			return r;
		}
	};

}
