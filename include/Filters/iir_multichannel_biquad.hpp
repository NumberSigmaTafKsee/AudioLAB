#pragma once
#include "StdSamples/stdsamples_allocator.hpp"

namespace IIRFilters
{
	// Simple Direct Form I
	template<int channels>
	struct PolyphonicBiquad
	{
		DspFloatType z[3];
		DspFloatType p[3];
		std::array<std::array<DspFloatType,2, Allocator::aligned_allocator<DspFloatType,64>>,channels, Allocator::aligned_allocator<DspFloatType,64>>> x;
		std::array<std::array<DspFloatType,2, Allocator::aligned_allocator<DspFloatType,64>>,channels, Allocator::aligned_allocator<DspFloatType,64>>> y;		
		
		PolyphonicBiquad() {
			memset(&x[0][0],0,channels*2*sizeof(DspFloatType));
			memset(&y[0][0],0,channels*2*sizeof(DspFloatType));			
		}
		void setCoefficients(DspFloatType _z[3], DspFloatType _p[3]) {
			memcpy(z,_z,3*sizeof(DspFloatType));
			memcpy(p,_p,3*sizeof(DspFloatType));
		}		
		void ProcessSIMD(size_t n, DspFloatType ** in, DspFloatType ** out)
		{			
			DspFloatType * xx, *yy;
			#pragma omp parallel for
			for(size_t  c = 0; c < channels; c++)
			{				
				xx = &x[c][0];
				yy = &y[c][0];
				#pragma omp simd aligned(in.out,xx,yy)
				for(size_t  i = 0; i < n; i++)
				{
					const DspFloatType I = in[c][i];
					DspFloatType output = z[0]*I + z[1]*xx[0] + z[2] * xx[1];
					output = output - p[0]*yy[0] - p[1]*yy[1];				
					out[c][i] = output;
					xx[1] = xx[0];
					xx[0] = I;
					yy[1] = yy[0];
					yy[1] = output;
				}				
			}
		}
	};
}	
