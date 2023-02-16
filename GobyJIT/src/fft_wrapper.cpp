#include <algorithm>
#include <fftw3.h>

#include "fft_wrapper.hpp"

FFTwrapper::FFTwrapper(int n){
    n_r_ = n;
    n_c_ = n/2+1;
    real_ = (double*)fftw_malloc(sizeof(double)*n_r_);
    complex_ = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*n_c_);
    r2c_ = fftw_plan_dft_r2c_1d(n_r_, real_, complex_, FFTW_ESTIMATE);
    c2r_ = fftw_plan_dft_c2r_1d(n_r_, complex_, real_, FFTW_ESTIMATE);
}
FFTwrapper::FFTwrapper(int n, bool is_complex_state): FFTwrapper(n){
    is_complex_state_ = is_complex_state;
}
FFTwrapper::FFTwrapper(FFTwrapper &&fft){
    n_r_ = fft.n_r_;
    n_c_ = fft.n_c_;
    real_ = fft.real_;
    complex_ = fft.complex_;
    r2c_ = fft.r2c_;
    c2r_ = fft.c2r_;

    fft.real_ = nullptr;
    fft.complex_ = nullptr;
    fft.r2c_ = nullptr;
    fft.c2r_ = nullptr;
}
FFTwrapper &FFTwrapper::operator=(FFTwrapper &&fft){
    n_r_ = fft.n_r_;
    n_c_ = fft.n_c_;
    std::swap(real_, fft.real_);
    std::swap(complex_, fft.complex_);
    std::swap(r2c_, fft.r2c_);
    std::swap(c2r_, fft.c2r_);
    return *this;
}

FFTwrapper::~FFTwrapper(){
    //fftw_plan is just a typedefed pointer
    if(r2c_)
        fftw_destroy_plan(r2c_);
    if(c2r_)
        fftw_destroy_plan(c2r_);
    if(real_)
        fftw_free(real_);
    if(complex_)
        fftw_free(complex_);
}
double *FFTwrapper::get_real_ptr(){
    return real_;
}
size_t FFTwrapper::get_real_c() const{
    return n_r_;
}
fftw_complex *FFTwrapper::get_complex_ptr(){
    return complex_;
}
size_t FFTwrapper::get_complex_c() const{
    return n_c_;
}
void FFTwrapper::r2c_transform(){
    if(!is_complex_state_){
        fftw_execute(r2c_);
        is_complex_state_ = true;
    }
}
void FFTwrapper::c2r_transform(){
    if(is_complex_state_){
        fftw_execute(c2r_);
        is_complex_state_ = false;
        //Rescaling
        for(size_t i=0; i < (size_t)n_r_; ++i)
            real_[i] /= n_r_;
    }
}