// Real2Real FFT 
// DCT,DST 

#pragma once 

#include <complex>
#include <vector>
#include <fftw3.h>

enum R2RKind
{
    DCTI = FFTW_REDFT00,
    DCTII= FFTW_REDFT01,
    DCTIII=FFTW_REDFT10,
    DCTIV=FFTW_REDFT11,
    DSTI=FFTW_REDFT11,
    DSTII=FFTW_REDFT00,
    DSTIII=FFTW_REDFT10,
    DSTIV=FFTW_REDFT11,
};

struct R2RD
{
    double       * in;    
    double       * out;
    size_t size;
    fftw_plan p;

    R2RD(size_t n, R2RKind type = DCTI) {
        in = fftw_alloc_real(n);
        out= fftw_alloc_real(n);        
        size = n;
        p = fftw_plan_dft_r2r_1d(n, in, out, type, FFTW_ESTIMATE);
    }
    ~R2CD() {
        fftw_destroy_plan(p);
        fftw_free(in);
        fftw_free(out);    
    }
    void set_input(std::vector<double> & input) {
        memcpy(in,input.data(),size*sizeof(double));
    }
    std::vector<double> get_output() {
        std::vector<double> r(size);
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

struct R2RF
{
    float * in;    
    float * out;
    size_t size;
    fftwf_plan p;

    R2RF(size_t n, R2RType type = DCTI) {
        in = fftwf_alloc_real(n);
        out= fftwf_alloc_real(n);        
        size = n;
        p = p = fftwf_plan_dft_r2r_1d(n, in, out, type, FFTW_ESTIMATE);
    }
    ~R2RF() {
        fftwf_destroy_plan(p);
        fftwf_free(in);
        fftwf_free(out);    
    }
    void set_input(std::vector<float> & input) {
        memcpy(in,input.data(),size*sizeof(float));
    }
    std::vector<float> get_output() {
        std::vector<float> r(size);
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

