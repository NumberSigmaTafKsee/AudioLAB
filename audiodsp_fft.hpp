#pragma once

namespace AudioDSP
{    
    struct FFTPlan
    {
        virtual void set_input(const DspFloatType * ptr) {
            assert(1==0);
        }
        virtual void get_output(DspFloatType * ptr) {
            assert(1==0);
        }
        virtual void set_complex_input(const std::complex<DspFloatType> * ptr) {
            assert(1==0);
        }
        virtual void get_complex_output(std::complex<DspFloatType> * ptr) {
            assert(1==0);
        }

        virtual void forward() = 0;
        virtual void backward() = 0;
    };

    void fft(FFTPlan& plan, const std::complex<DspFloatType> * in, std::complex<DspFloatType>* out)
    {
        plan.set_complex_input(in);
        plan.forward();
        plan.get_complex_output(out);
    }

    void ifft(FFTPlan& plan, const  std::complex<DspFloatType> * in, std::complex<DspFloatType>* out)
    {
        plan.set_complex_input(in);
        plan.backward();
        plan.get_complex_output(out);
    }
    void fft(FFTPlan& plan, const  DspFloatType * in, std::complex<DspFloatType>* out)
    {
        plan.set_input(in);
        plan.forward();
        plan.get_complex_output(out);
    }

    void ifft(FFTPlan& plan, std::complex<DspFloatType> * in, DspFloatType* out)
    {
        plan.set_complex_input(in);
        plan.backward();
        plan.get_output(out);
    }
}