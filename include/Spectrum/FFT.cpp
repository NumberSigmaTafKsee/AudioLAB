#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include "samples/sample_dsp.hpp"
#include "fftwpp/fftw++.h"
#include "Spectrum/vcvfft.hpp"


void print_complex(std::vector<std::complex<float>> r)
{
    for(size_t i = 0; i < r.size(); i++) {
        std::cout << r[i] << ",";            
    }
    std::cout << std::endl;
}

void print_complex2(DSP::complex_vector<float> r)
{
    for(size_t i = 0; i < r.size(); i++) {
        std::cout << r[i] << ",";            
    }
    std::cout << std::endl;
}

void print_real(std::vector<std::complex<float>> r)
{
    for(size_t i = 0; i < r.size(); i++) {
        std::cout << abs(r[i]) << ",";            
    }
    std::cout << std::endl;
}
void print_vector(std::vector<float> v)
{
    for(size_t i = 0; i < v.size(); i++) {
        std::cout << v[i] << ",";            
    }
    std::cout << std::endl;
}
void print_vector(std::vector<double> v)
{
    for(size_t i = 0; i < v.size(); i++) {
        std::cout << v[i] << ",";            
    }
    std::cout << std::endl;
}
void print_vector2(sample_vector<float> v)
{
    for(size_t i = 0; i < v.size(); i++) {
        std::cout << v[i] << ",";            
    }
    std::cout << std::endl;
}
void tests()
{
    sample_vector<float> test2 = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};
    sample_vector<float> temp1;
    DSP::complex_vector<float> c;

    std::cout << std::endl;
    std::cout << "audiofft" << std::endl;
    c = DSP::audiofft_forward(test2);
    for(size_t i = 0; i < 16; i++)
        std::cout << c[i] << ",";        
    std::cout << std::endl;
    temp1 = DSP::audiofft_inverse(c);
    for(size_t i = 0; i < 16; i++)
        std::cout << temp1[i] << ",";
    std::cout << std::endl;

    std::cout << std::endl;
    std::cout << "vcvfft" << std::endl;
    
    if (32 < pffft_min_fft_size(PFFFT_REAL))
    {
        fprintf(stderr, "Error: minimum FFT transformation length is %d\n", pffft_min_fft_size(PFFFT_REAL));
        return;
    }


    VCVFFT::RealFFT real_fft(32);    
    float * x1 = (float*)pffft_aligned_malloc(32 * sizeof(float));
    float * x2 = (float*)pffft_aligned_malloc(64 * sizeof(float));
    memcpy(x1,test2.data(),16*sizeof(float));
    memcpy(x1+16,test2.data(),16*sizeof(float));
    real_fft.rfftUnordered(x1,x2);    
    for(size_t i = 0; i < 32; i++)
        x2[i] /= 32;    
    real_fft.irfftUnordered(x2,x1);
    for(size_t i = 0; i < 32; i++)
        std::cout << x1[i] << ",";
    std::cout << std::endl;

    std::cout << std::endl;
    std::cout << "fftwpp" << std::endl;

    fftwpp::fft1d forward(16,-1);
    fftwpp::fft1d backward(16,1);
    Complex * f = utils::ComplexAlign(16);

    for(size_t i = 0; i < 16; i++)
    {
        f[i] = Complex(test2[i],0);        
    }
    forward.fft(f);
    for(size_t i = 0; i < 16; i++)
    {
        f[i] /= 16;
    }
    for(size_t i = 0; i < 16; i++)
        std::cout << f[i] << ",";
    std::cout << std::endl;
    backward.fft(f);
    for(size_t i = 0; i < 16; i++)
        std::cout << abs(f[i]) << ",";
    std::cout << std::endl;

    std::cout << std::endl;
    std::cout << "samples/dsp" << std::endl;

    DSP::R2CF rfft(16);
    rfft.set_input(test2);
    rfft.Execute();
    c = rfft.get_output();
    for(size_t i = 0; i < c.size(); i++)
        c[i] /= test2.size();
    print_complex2(c);
    DSP::C2RF ifft(16);
    ifft.set_input(c);
    ifft.Execute();
    test2 = ifft.get_output();
    print_vector2(test2);
}




int main()
{
    std::vector<float>  test2 = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};    
    //tests();
    VCVFFT::RealFFT real_fft(32);    
    //float * x1 = (float*)pffft_aligned_malloc(32 * sizeof(float));
    //float * x2 = (float*)pffft_aligned_malloc(64 * sizeof(float));
    sample_vector<float> x1(32),x2(64);
    memcpy(x1.data(),test2.data(),16*sizeof(float));
    memcpy(x1.data()+16,test2.data(),16*sizeof(float));
    real_fft.rfftUnordered(x1.data(),x2.data());    
    for(size_t i = 0; i < 32; i++)
        x2[i] /= 32;    
    real_fft.irfftUnordered(x2.data(),x1.data());
    for(size_t i = 0; i < 32; i++)
        std::cout << x1[i] << ",";
    std::cout << std::endl;
   
}