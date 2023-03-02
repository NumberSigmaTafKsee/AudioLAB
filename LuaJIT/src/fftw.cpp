#include <fftw.h>
#include <complex>
#include "StdVector.h"


struct FFTWReal1D
{
    rfftw_plan plan_forward;
    rfftw_plan plan_backward;
    size_t num_bins;
    
    FFTWReal1D(size_t n) {
        num_bins = n;
        plan_forward = rfftw_create_plan(n,FFTW_REAL_TO_COMPLEX,FFTW_ESTIMATE);
        plan_backward = rfftw_create_plan(n,FFTW_COMPLEX_TO_REAL,FFTW_ESTIMATE);
    }
    ~FFTWReal1D() {
        rfftw_plan_destroy(plan_forward);
        rfftw_plan_destroy(plan_backward);
    }

    void Forward( Std::ValVector<float> & input, Std::ValVector<std::complex<float>> & output) {
        fftw_complex * A = (fftw_complex*)(&output[0]);
        rfftw_one_real_to_complex(plan_forward,input.data(),A);
    }
    void Backward( Std::ValVector<std::complex<float>> & input, Std::ValVector<float> & output) {
        fftw_complex * A = (fftw_complex*)(&input[0]);
        rfftw_one_complex_to_real(plan_backward,)A,output.data());
    }
};

int main()
{
    Std::ValVector<float> input(1024), r(1024);
    Std::ValVector<std::complex<float>> output(1024);
    FFTWReal1D fft(1024);
    input.random();    
    fft.Forward(input,output);
    Std::println(input);
    fft.Backward(output,r);
    Std::println(r);
}