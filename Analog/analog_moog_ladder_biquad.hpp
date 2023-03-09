#pragma once

namespace Analog::Moog
{
	template<typename DSP>
	class BiQuadBase : public FilterProcessor
	{
	public:
		
		BiQuadBase() : FilterProcessor()
		{
			bCoef = {{0.0f, 0.0f, 0.0f}};
			aCoef = {{0.0f, 0.0f}};
			w = {{0.0f, 0.0f}};
		}	
		~BiQuadBase()
		{

		}	
		// DF-II impl	
		void Process(uint32_t n, DspFloatType * samples)
		{
			Undenormal denormal;
			DspFloatType out = 0;
			for (int s = 0; s < n; ++s)
			{
				out = bCoef[0] * samples[s] + w[0];
				w[0] = bCoef[1] * samples[s] - aCoef[0] * out + w[1];
				w[1] = bCoef[2] * samples[s] - aCoef[1] * out;
				samples[s] = out;
			}
		}
		DspFloatType Tick(DspFloatType s, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1)
		{
			Undenormal denormal;
			DspFloatType out = bCoef[0] * s + w[0];
			w[0] = bCoef[1] * s - aCoef[0] * out + w[1];
			w[1] = bCoef[2] * s - aCoef[1] * out;
			return out;
		}
		void SetBiquadCoefs(std::array<DspFloatType, 3> b, std::array<DspFloatType, 2> a)
		{
			bCoef = b;
			aCoef = a;
		}
		
	protected:
		std::array<DspFloatType, 3> bCoef; // b0, b1, b2
		std::array<DspFloatType, 2> aCoef; // a1, a2
		std::array<DspFloatType, 2> w; // delays
	};
}
