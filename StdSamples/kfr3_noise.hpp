#pragma once
#include <cmath>

/////////////////////////////////////////////////////////////////////////////////////////
// Noise
/////////////////////////////////////////////////////////////////////////////////////////
#define PINK_NOISE_NUM_STAGES 3

namespace kfr3
{
  double WhiteNoise() {
      return noise.rand();
  }

  class PinkNoise {
  public:
    PinkNoise() {
    srand ( time(NULL) ); // initialize random generator
      clear();
    }

    void clear() {
      for( size_t i=0; i< PINK_NOISE_NUM_STAGES; i++ )
        state[ i ] = 0.0;
      }

    double tick() {
      static const double RMI2 = 2.0 / double(RAND_MAX); // + 1.0; // change for range [0,1)
      static const double offset = A[0] + A[1] + A[2];

    // unrolled loop
      double temp = double( noise.rand() );
      state[0] = P[0] * (state[0] - temp) + temp;
      temp = double( noise.rand() );
      state[1] = P[1] * (state[1] - temp) + temp;
      temp = double( noise.rand() );
      state[2] = P[2] * (state[2] - temp) + temp;
      return ( A[0]*state[0] + A[1]*state[1] + A[2]*state[2] )*RMI2 - offset;
    }

  protected:
    double state[ PINK_NOISE_NUM_STAGES ];
    static const double A[ PINK_NOISE_NUM_STAGES ];
    static const double P[ PINK_NOISE_NUM_STAGES ];
  };

  const double PinkNoise::A[] = { 0.02109238, 0.07113478, 0.68873558 }; // rescaled by (1+P)/(1-P)
  const double PinkNoise::P[] = { 0.3190,  0.7756,  0.9613  };

  double Pink() {
      static PinkNoise pink;
      return pink.tick();
  }
  double GaussianWhiteNoise()
  {
      double R1 = noise.rand();
      double R2 = noise.rand();

      return (double) std::sqrt( -2.0f * std::log( R1 )) * std::cos( 2.0f * M_PI * R2 );
  }
  double GaussRand (int m, double s)
  {
    static int pass = 0;
    static double y2;
    double x1, x2, w, y1;

    if (pass)
    {
        y1 = y2;
    } else  {
        do {
          x1 = 2.0f * noise.rand () - 1.0f;
          x2 = 2.0f * noise.rand () - 1.0f;
          w = x1 * x1 + x2 * x2;
        } while (w >= 1.0f);

        w = (double)std::sqrt (-2.0 * std::log (w) / w);
        y1 = x1 * w;
        y2 = x2 * w;
    }
    pass = !pass;

    return ( (y1 * s + (double) m));
  }

  // +/-0.05dB above 9.2Hz @ 44,100Hz
  class PinkingFilter
  {
    double b0, b1, b2, b3, b4, b5, b6;
  public:
    PinkingFilter() : b0(0), b1(0), b2(0), b3(0), b4(0), b5(0), b6(0) {}
    double process(const double s)
    {
      b0 = 0.99886 * b0 + s * 0.0555179;
      b1 = 0.99332 * b1 + s * 0.0750759;
      b2 = 0.96900 * b2 + s * 0.1538520;
      b3 = 0.86650 * b3 + s * 0.3104856;
      b4 = 0.55000 * b4 + s * 0.5329522;
      b5 = -0.7616 * b5 - s * 0.0168980;
      const double pink = (b0 + b1 + b2 + b3 + b4 + b5 + b6 + (s * 0.5362)) * 0.11;
      b6 = s * 0.115926;
      return pink;
    }
  };

  class BrowningFilter
  {
  double l;
  public:
    BrowningFilter() : l(0) {}
    double process(const double s)
    {
      double brown = (l + (0.02f * s)) / 1.02f;
      l = brown;
      return brown * 3.5f; // compensate for gain
    }
  };


  double PinkNoiseGenerator()
  {
      static PinkingFilter pink;
      return pink.process(noise.rand());
  }
  double BrownNoiseGenerator()
  {
      static BrowningFilter brown;
      return brown.process(noise.rand());
  }
}
