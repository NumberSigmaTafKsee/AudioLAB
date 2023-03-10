#pragma once

namespace FX::Filters
{
    class CDCFilter : public FilterProcessor {
    public:
        
        CDCFilter(DspFloatType Fs=44100.0f, DspFloatType Fc=10)
        : FilterProcessor()
        {
            a0 = 1.0; b1 = 0.0; z1 = 0.0;
            setCutoff(Fc/Fs);
        };
        
        void setCutoff(DspFloatType Fc)
        {
            b1 = exp(-2.0 * M_PI * Fc);
            a0 = 1.0 - b1;
        }
        DspFloatType process(DspFloatType in)
        {
            z1 = in * a0 + z1 * b1;
            return in - z1;
        }
        DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=0, DspFloatType Y=0)
        {
            return process(I);
        }
        void ProcessSIMD(size_t num, DspFloatType * in, DspFloatType * out)
        {
            #pragma omp simd aligned(in,out)
            for(size_t i = 0; i < num ; i++) {
                z1 = in[i] * a0 + z1 * b1;
                out[i] = in[i] - z1;
            }
        }
        void ProcessBlock(size_t num, DspFloatType * in, DspFloatType * out)
        {
            ProcessSIMD(num,in,out);            
        }
        void ProcessInplace(size_t n, DspFloatType * out) {
            ProcessSIMD(n,out,out);
        }
    protected:    
        DspFloatType a0, b1, z1;
    };


    class CDCRemover : public FilterProcessor {
    public:

        CDCRemover(DspFloatType Fs=44100.0f, DspFloatType Fc=10.0f) : FilterProcessor() {
            a0 = 1.0; b1 = 0.0; z1 = 0.0;
            setCutoff(Fc/Fs);
        };
        
        // this is a high pass
        void setCutoff(DspFloatType Fc)
        {
            b1 = -exp(-2.0 * M_PI * (0.5 - Fc));
            a0 = 1.0 + b1;
        }
        DspFloatType process(DspFloatType in)
        {
            return z1 = in * a0 + z1 * b1;
        }
        DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=0, DspFloatType Y=0)
        {
            return process(I);
        }
        void ProcessSIMD(size_t num, DspFloatType * in, DspFloatType * out)
        {
            #pragma omp simd
            for(size_t i = 0; i < num ; i++) {
                z1 = in[i] * a0 + z1 * b1;
                out[i] = z1;
            }
        }
        void ProcessBlock(size_t num, DspFloatType * in, DspFloatType * out)
        {
            ProcessSIMD(num,in,out);            
        }
        void ProcessInplace(size_t n, DspFloatType * out) {
            ProcessSIMD(n,out,out);
        }
    protected:    
        DspFloatType a0, b1, z1;
    };
}
