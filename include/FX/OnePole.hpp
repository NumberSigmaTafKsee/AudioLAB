#pragma once

#include <cmath>

namespace FX::Filters
{
    class OnePole {
    public:
        OnePole() {a0 = 1.0; b1 = 0.0; z1 = 0.0;};
        OnePole(DspFloatType Fc) {z1 = 0.0; setFc(Fc);};
        ~OnePole() = default;

        void setFc(DspFloatType Fc);
        void setHighPass(DspFloatType Fc);
        DspFloatType process(DspFloatType in);

        DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {
            return A*process(I);
        }
        void ProcessSIMD(size_t n, DspFloatType * in, DspFloatType * out) {
            #pragma omp simd
            for(size_t i = 0; i < n; i++)    
            {
                z1 = in[i] * a0 + z1 * b1;
                out[i] = z1;
            }
        }
        void ProcessBlock(size_t n, DspFloatType * in, DspFloatType * out) { 
            ProcessSIMD(n,in,out);
        }
        void ProcessInplace(size_t n, DspFloatType * in) { 
            ProcessSIMD(n,in,in);
        }

    protected:
        DspFloatType a0, b1, z1;
    };

    // low pass
    inline void OnePole::setFc(DspFloatType Fc) {
        b1 = exp(-2.0 * M_PI * Fc);
        a0 = 1.0 - b1;
    }

    inline void OnePole::setHighPass(DspFloatType Fc) {
        b1 = -exp(-2.0 * M_PI * (0.5 - Fc));
        a0 = 1.0 + b1;
    }

    inline DspFloatType OnePole::process(DspFloatType in) {
        return z1 = in * a0 + z1 * b1;
    }
    
}