#pragma once

#include "Vc/vectorclass.h"

#include <cstdint>


// Numerical Recipes In C p.284.
struct alignas(64) White16 {
  Vec16ui buffer{0};

  White16(uint32_t seed) { setSeed(seed); }

  void setSeed(uint32_t seed)
  {
    for (int idx = 0; idx < 16; ++idx) {
      seed = 1664525L * seed + 1013904223L;
      buffer.insert(idx, seed);
    }
  }

  Vec16f process()
  {
    buffer = 1664525L * buffer + 1013904223L;
    return to_float(buffer) / float(UINT32_MAX);
  }
};

// Paul Kellet's refined method in Allan's analysis.
// http://www.firstpr.com.au/dsp/pink-noise/
template<typename Sample> class Pink {
public:
  Pink(int32_t seed) : rng(seed) {}

  Sample process()
  {
    const Sample gain = 0.125;
    auto white = rng.process();
    b0 = Sample(0.99886) * b0 + white * Sample(0.0555179) * gain;
    b1 = Sample(0.99332) * b1 + white * Sample(0.0750759) * gain;
    b2 = Sample(0.96900) * b2 + white * Sample(0.1538520) * gain;
    b3 = Sample(0.86650) * b3 + white * Sample(0.3104856) * gain;
    b4 = Sample(0.55000) * b4 + white * Sample(0.5329522) * gain;
    b5 = Sample(-0.7616) * b5 - white * Sample(0.0168980) * gain;
    auto pink = b0 + b1 + b2 + b3 + b4 + b5 + b6 + white * Sample(0.5362) * gain;
    b6 = white * Sample(0.115926);
    return pink;
  }

private:
  White<Sample> rng;
  Sample b0 = 0.0;
  Sample b1 = 0.0;
  Sample b2 = 0.0;
  Sample b3 = 0.0;
  Sample b4 = 0.0;
  Sample b5 = 0.0;
  Sample b6 = 0.0;
};