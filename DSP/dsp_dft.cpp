#include <complex>
#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>

template<typename T>
struct DFT 
{        
    void Forward(std::vector<T> & input, std::vector<std::complex<T>> & bins)
    {        
        #pragma omp parallel for
        for(size_t k = 0; k < bins.size(); k++)
        {
            bins[k].real(0);
            bins[k].imag(0);
            
            T a = 0.0;
            T b = 0.0;
            for(size_t j = 0; j < input.size(); j++) 
            {
                T t = ((T)j*(T)k*2.0f*M_PI)/(T)input.size();                                
                T value = input[j];
                a += std::cos(t)*value;
                b += std::sin(t)*value;
            } 
            bins[k] = std::complex<T>(a,-b);
        }
    }

    void Backward(std::vector<std::complex<T>> & bins, std::vector<T> & output)
    {                
        #pragma omp parallel for
        for(size_t i = 0; i < output.size(); i++)
        {
            std::complex<T> result;
            std::complex<T> ci(0,1);
            for(size_t j = 0; j < bins.size(); j++)
            {     
                // seperate real and imaginary       
                //auto phase = (2 * M_PI * k * n) / N;
                //a += cos(phase) * output[k].real() - sin(phase) * output[k].imag();                

                T  t = (2*M_PI*(T)i*(T)j)/(T)bins.size();
                result += std::exp(t*ci) * bins[j];
            }
            result /= bins.size();            
            output[i] = std::abs(result);
        }
    }
};

void TestDFT() {
    DFT<float> dft;    
    std::vector<float> x(1024);
    std::vector<std::complex<float>> o(1024);
    for(size_t i = 0; i < 1024; i++) {
        x[i] = (float)rand()/(float)RAND_MAX;
    }
    dft.Forward(x,o);
    for(size_t i = 0; i < 1024; i++)
        std::cout << o[i] << std::endl;    
    std::vector<float> r(1024);
    dft.Backward(o,r);
    for(size_t i = 0; i < 1024; i++)
    {
        std::cout << "x=" << x[i] << ", r=" << r[i] << std::endl;
    }
}
#include "SndFile.hpp"
int main() 
{
    SoundWave::SndFileReader wav("test.wav");
    size_t size = wav.size();
    std::cout << size << std::endl;
    std::vector<float> samples(size);
    wav >> samples;
    std::vector<std::complex<float>> bins(size);
    DFT<float> dft;
    dft.Forward(samples,bins);
    std::vector<float> mag(bins.size()/2);
    for(size_t i = 0; i < bins.size()/2; i++)
        mag[i] = std::abs(bins[i])/22050.0f;
    typename std::vector<float>::iterator x = std::max_element(mag.begin(),mag.end());
    std::cout << "mag[" << x - mag.begin() << "] =" << *x << std::endl;
    
    for(; x != mag.end(); x++) {
        bins[*x].real(0);
        bins[*x].imag(0);
    }
    dft.Backward(bins,samples);
    SoundWave::SndFileWriter out("out.wav",wav.format(),wav.channels(),wav.samplerate());
    out << samples;

}
