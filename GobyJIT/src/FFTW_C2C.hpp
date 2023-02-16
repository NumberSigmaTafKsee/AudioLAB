#pragma once 

#include <complex>
#include <vector>
#include <fftw3.h>

struct C2CD
{
    fftw_complex * in;    
    fftw_complex * out;
    size_t size;
    fftw_plan p;

    enum Direction {
        BACKWARD= FFTW_BACKWARD,
        FORWARD = FFTW_FORWARD,
    };

    C2CD(size_t n, Direction dir = FORWARD) {
        in = fftw_alloc_complex(n);
        out= fftw_alloc_complex(n);        
        size = n;
        p = fftw_plan_dft_1d(n, in, out, dir, FFTW_ESTIMATE);
    }
    ~C2CD() {
        fftw_destroy_plan(p);
        fftw_free(in);
        fftw_free(out);    
    }
    void set_input(std::vector<std::complex<double>> & input) {
        for(size_t i = 0; i < size; i++) {
            in[i][0] = input[i].real();
            in[i][1] = input[i].imag();
        }
    }
    std::vector<std::complex<double>> get_output() {
        std::vector<std::complex<double>> r(size);
        for(size_t i = 0; i < size; i++ )
        {
            r[i].real(out[i][0]);
            r[i].imag(out[i][1]);
        }
        return r;
    }
    void Execute() {
        fftw_execute(p);
    }
};

struct C2CF
{
    fftwf_complex * in;    
    fftwf_complex * out;
    size_t size;
    fftwf_plan p;

    enum Direction {
        BACKWARD=FFTW_BACKWARD,
        FORWARD=FFTW_FORWARD,
    };

    C2CF(size_t n, Direction dir = FORWARD) {
        in = fftwf_alloc_complex(n);
        out= fftwf_alloc_complex(n);        
        size = n;
        p = fftwf_plan_dft_1d(n, in, out, dir, FFTW_ESTIMATE);
    }
    ~C2CF() {
        fftwf_destroy_plan(p);
        fftwf_free(in);
        fftwf_free(out);    
    }
    void set_input(std::vector<std::complex<float>> & input) {
        for(size_t i = 0; i < size; i++) {
            in[i][0] = input[i].real();
            in[i][1] = input[i].imag();
        }
    }
    std::vector<std::complex<float>> get_output() {
        std::vector<std::complex<float>> r(size);
        for(size_t i = 0; i < size; i++ )
        {
            r[i].real(out[i][0]);
            r[i].imag(out[i][1]);
        }
        return r;
    }
    void Execute() {
        fftwf_execute(p);
    }
};

