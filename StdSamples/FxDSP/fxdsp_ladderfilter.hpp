#pragma once

#include <cstddef>
#include <cstdlib>
#include <cmath>

#include "fxdsp.hpp"

namespace FXDSP
{
    /** Boltzman's constant */
    static const float BOLTZMANS_CONSTANT = 1.38e-23;

    /** Magnitude of electron charge */
    static const float Q = 1.609e-19;

    /* LadderFilter ********************************************************/
    template<typename T>
    struct LadderFilter
    {
        T y[4];
        T w[4];
        T Vt;           // transistor treshold voltage [V]
        T sample_rate;
        T cutoff;
        T resonance;

        LadderFilter(T sample_rate)
        {
            ClearBuffer(y, 4);
            ClearBuffer(w, 4);
            Vt = 0.026;
            cutoff = 0;
            resonance = 0;
            this->sample_rate = _sample_rate;
        }
        void setCutoff(T c) {
            // work with cv
            cutoff = c;
        }
        void setResonance(T r) {
            resonance = r;
        }
        void ProcessBlock(size_t n, T * inBuffer, T * outBuffer)
        {
            // Pre-calculate Scalars
            T TWO_VT_INV = 1.0 / (2 * Vt);
            T TWO_VT_G = 2 * Vt * (1 - exp(-TWO_PI * cutoff / sample_rate));
            
            // Filter audio 
            #pragma omp simd aligned(inBuffer,outBuffer)
            for (unsigned i = 0; i < n; ++i)
            {

                // Stage 1 output
                y[0] = y[0] + TWO_VT_G * (f_tanh(inBuffer[i] - 4 * \
                            f_tanh(2 * resonance * y[3]) * \
                            TWO_VT_INV) - w[0]);

                w[0] = f_tanh(y[0] * TWO_VT_INV);

                // Stage 2 output
                y[1] = y[1] + TWO_VT_G * (w[0]- w[1]);
                w[1] = f_tanh(y[1] * TWO_VT_INV);

                // Stage 3 output
                y[2] = y[2] + TWO_VT_G * (w[1]- w[2]);
                w[2] = f_tanh(y[2] * TWO_VT_INV);

                // Stage 4 output
                y[3] = y[3] + TWO_VT_G * (w[2]- w[3]);
                w[3] = f_tanh(y[3] * TWO_VT_INV);

                // Write output
                outBuffer[i] = y[3];
            }
        }
        T Tick(T I, T A=1, T X=1, T Y=1)
        {
            T c = cutoff;
            T r = resonance;            
            setCutoff(c + X*c);
            setResonance(r + Y*r);

            T TWO_VT_INV = 1.0 / (2 * Vt);
            T TWO_VT_G = 2 * Vt * (1 - exp(-TWO_PI * cutoff / sample_rate));
            // Stage 1 output            
            y[0] = y[0] + TWO_VT_G * (f_tanh(inBuffer[i] - 4 * \
                        f_tanh(2 * resonance * y[3]) * \
                        TWO_VT_INV) - w[0]);

            w[0] = f_tanh(y[0] * TWO_VT_INV);

            // Stage 2 output
            y[1] = y[1] + TWO_VT_G * (w[0]- w[1]);
            w[1] = f_tanh(y[1] * TWO_VT_INV);

            // Stage 3 output
            y[2] = y[2] + TWO_VT_G * (w[1]- w[2]);
            w[2] = f_tanh(y[2] * TWO_VT_INV);

            // Stage 4 output
            y[3] = y[3] + TWO_VT_G * (w[2]- w[3]);
            w[3] = f_tanh(y[3] * TWO_VT_INV);

            setCutoff(c);
            setResonance(r);

            return A*y[3];
        }
        void setTemperature(T tempC)
        {
            T t = tempC + 273.15;
            Vt = BOLTZMANS_CONSTANT * t / Q;
        }
    };

}
