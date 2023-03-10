//
//  OnePole.h
//

#ifndef OnePole_h
#define OnePole_h

#include <math.h>

class OnePole {
public:
    
    OnePole() {a0 = 1.0; b1 = 0.0; z1 = 0.0;};
    OnePole(DspFloatType Fc, DspFloatType sample_rate=44100) { 
        z1 = 0.0; 
        setFc(Fc);
        sr = sample_rate;
    };
    ~OnePole();
    
    void  setLowPassFc(DspFloatType Fc) {
        b1 = exp(-2.0 * M_PI * Fc/sr);
        a0 = 1.0 - b1;
    }
    void  setHighPassFc(DspFloatType Fc) {
        b1 = -exp(-2.0 * M_PI * (0.5 - Fc/sr));
        a0 = 1.0 + b1;
    }
    void setDCBlocker(DspFloatType f = 10.0) {
        setLowPass(f/sampleRate);
    }
    DspFloatType process(DspFloatType in)
    {
        x = in;
        return z1 = in * a0 + z1 * b1;
    }
    DspFloatType hi() {
        return x - z1;
    }
    void ProcessSIMD(size_t n, DspFloatType * in, DspFloatType * out) {
        #pragma omp simd aligned(in,out)
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
    DspFloatType x;
    DspFloatType a0, b1, z1;
    DspFloatType sr;
};


