	

#ifndef _PinkNoise_H
#define _PinkNoise_H

// Technique by Larry "RidgeRat" Trammell 3/2006
// http://home.earthlink.net/~ltrammell/tech/pinkalg.htm
// implementation and optimization by David Lowenfels

#include <cstdlib>
#include <ctime>

#define PINK_NOISE_NUM_STAGES 3

class PinkNoise {
public:
  PinkNoise() {
  srand ( time(NULL) ); // initialize random generator
    clear();
  }

  void clear() {
    #pragma omp simd
    for( size_t i=0; i< PINK_NOISE_NUM_STAGES; i++ )
      state[ i ] = 0.0;
    }

  DspFloatType Tick(DspFloatType I=1,DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {
    static const DspFloatType RMI2 = 2.0 / DspFloatType(RAND_MAX); // + 1.0; // change for range [0,1)
    static const DspFloatType offset = A[0] + A[1] + A[2];

  // unrolled loop
    DspFloatType temp = DspFloatType( rand() );
    state[0] = P[0] * (state[0] - temp) + temp;
    temp = DspFloatType( rand() );
    state[1] = P[1] * (state[1] - temp) + temp;
    temp = DspFloatType( rand() );
    state[2] = P[2] * (state[2] - temp) + temp;
    return ( A[0]*state[0] + A[1]*state[1] + A[2]*state[2] )*RMI2 - offset;
  }
	void ProcessSIMD(size_t n, DspFloatType * in, DspFloatType * out)
	{
		static const DspFloatType RMI2 = 2.0 / DspFloatType(RAND_MAX); // + 1.0; // change for range [0,1)
		static const DspFloatType offset = A[0] + A[1] + A[2];

		#pragma omp simd aligned(in,out)
		for(size_t i = 0; i < n; i++)
		{
			DspFloatType temp = DspFloatType( in[i] );
			state[0] = P[0] * (state[0] - temp) + temp;
			temp = DspFloatType( rand() );
			state[1] = P[1] * (state[1] - temp) + temp;
			temp = DspFloatType( rand() );
			state[2] = P[2] * (state[2] - temp) + temp;
			out[i] = ( A[0]*state[0] + A[1]*state[1] + A[2]*state[2] )*RMI2 - offset;
		}
	}
		 
protected:
  DspFloatType state[ PINK_NOISE_NUM_STAGES ];
  static const DspFloatType A[ PINK_NOISE_NUM_STAGES ];
  static const DspFloatType P[ PINK_NOISE_NUM_STAGES ];
};

const DspFloatType PinkNoise::A[] = { 0.02109238, 0.07113478, 0.68873558 }; // rescaled by (1+P)/(1-P)
const DspFloatType PinkNoise::P[] = { 0.3190,  0.7756,  0.9613  };

#endif
