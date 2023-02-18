%module kissfftr
%{ 
#include <kiss_fft.h>
#include <kiss_fftr.h>
#include <cassert>
#include <complex>
using namespace std;
%}

%include "std_math.i"
%include "stdint.i"
%include "std_vector.i"

%import <kiss_fft.h>
%include <kiss_fftr.h>

%template(float_vector) std::vector<float>;
%template(complex_vector) std::vector<std::complex<float>>;

%inline 
%{
struct KissFFTR
{
    kiss_fftr_cfg  config;    
    kiss_fftr_cfg  iconfig;    

    KissFFTR(size_t nfft) 
    {
        config = kiss_fftr_alloc(nfft, 0,NULL,NULL);
        assert(config != NULL);
        iconfig = kiss_fftr_alloc(nfft, 1,NULL,NULL);
        assert(iconfig != NULL);
    }
    ~KissFFTR()
    {
        if(config) free(config);
        if(iconfig) free(iconfig);
    }        
    void forward( float* in, std::complex<float> * out)
    {                        
        kiss_fftr(config, in, (kiss_fft_cpx*)out);
    }
    void inverse(std::complex<float>* in, float* out)
    {        
        kiss_fftri(iconfig, (const kiss_fft_cpx*) in, out);
    }

};
%}