#pragma once 

#include <complex>
#include <vector>
#include <fftw3.h>
#include <cstring>

struct C2RD
{
    fftw_complex * in;    
    double * out;
    size_t size;
    fftw_plan p;

    C2RD(size_t n) {
        in = fftw_alloc_complex(n);
        out= fftw_alloc_real(n);        
        size = n;
        p = fftw_plan_dft_c2r_1d(n, in, out, FFTW_ESTIMATE);
    }
    ~C2RD() {
        fftw_destroy_plan(p);
        fftw_free(in);
        fftw_free(out);    
    }
    void set_input(std::vector<std::complex<double>> & input) {
        for(size_t i = 0; i < i; i++) {
            in[i][0] = input[i].real();
            in[i][1] = input[i].imag();
        }
    }
    std::vector<double> get_output() {
        std::vector<double> r(size);
        memcpy(r.data(),out, size * sizeof(double));
        return r;
    }
    void Execute() {
        fftw_execute(p);
    }
};

struct C2RF
{
    fftwf_complex * in;    
    float * out;
    size_t size;
    fftwf_plan p;

    C2RF(size_t n) {
        in = fftwf_alloc_complex(n);
        out= fftwf_alloc_real(n);        
        size = n;
        p = fftwf_plan_dft_c2r_1d(n, in, out, FFTW_ESTIMATE);
    }
    ~C2RF() {
        fftwf_destroy_plan(p);
        fftwf_free(in);
        fftwf_free(out);    
    }
    void set_input(std::vector<std::complex<float>> & input) {
        for(size_t i = 0; i < i; i++) {
            in[i][0] = input[i].real();
            in[i][1] = input[i].imag();
        }
    }
    std::vector<float> get_output() {
        std::vector<float> r(size);
        memcpy(r.data(),out, size*sizeof(float));
        return r;
    }
    void Execute() {
        fftwf_execute(p);
    }
};

