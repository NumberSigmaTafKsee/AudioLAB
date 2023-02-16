#pragma once 

#include <complex>
#include <vector>
#include <fftw3.h>

struct R2CD
{
    double       * in;    
    fftw_complex * out;
    size_t size;
    fftw_plan p;

    R2CD(size_t n) {
        in = fftw_alloc_real(n);
        out= fftw_alloc_complex(n);        
        size = n;
        p = fftw_plan_dft_r2c_1d(n, in, out, FFTW_ESTIMATE);
    }
    ~R2CD() {
        fftw_destroy_plan(p);
        fftw_free(in);
        fftw_free(out);    
    }
    void set_input(std::vector<double> & input) {
        memcpy(in,input.data(),size*sizeof(double));
    }
    std::vector<std::complex<double>> get_output() {
        std::vector<std::complex<double>> r(size);
        for(size_t i = 0; i < size; i++) {
            r[i].real(out[i][0]);
            r[i].imag(out[i][1]);
        }
        return r;
    }
    void Execute() {
        fftw_execute(p);
    }
};

struct R2CF
{
    float * in;    
    fftwf_complex * out;
    size_t size;
    fftwf_plan p;

    R2CF(size_t n) {
        in = fftwf_alloc_real(n);
        out= fftwf_alloc_complex(n);        
        size = n;
        p = p = fftwf_plan_dft_r2c_1d(n, in, out, FFTW_ESTIMATE);
    }
    ~R2CF() {
        fftwf_destroy_plan(p);
        fftwf_free(in);
        fftwf_free(out);    
    }
    void set_input(std::vector<float> & input) {
        memcpy(in,input.data(),size*sizeof(float));
    }
    std::vector<std::complex<float>> get_output() {
        std::vector<std::complex<float>> r(size);
        for(size_t i = 0; i < size; i++) {
            r[i].real(out[i][0]);
            r[i].imag(out[i][1]);
        }
            
        return r;
    }    
    void Execute() {
        fftwf_execute(p);
    }
};

