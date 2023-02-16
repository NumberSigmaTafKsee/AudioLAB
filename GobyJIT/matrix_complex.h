#pragma once
#include <cuda/std/ccomplex>

template<typename T> using complex = cuda::std::complex<T>;

template<typename T>
void cmatrix_r_realf(complex<T> * x, T * y, int m,int n);
template<typename T>
void cmatrix_r_imagf(complex<T> * x, T * y, int m,int n);

template<typename T> void cmatrix_r_argf(complex<T> * devPtr, T * x, int m,int n);
template<typename T> void cmatrix_r_normf(complex<T> * devPtr, complex<T> * x, int m,int n);
template<typename T> void cmatrix_r_conjf(complex<T> * devPtr, complex<T> * x, int m,int n);
template<typename T> void cmatrix_r_projf(complex<T> * devPtr, complex<T> * x,int m,int n);

template<typename T>
void cmatrix_2d_r_addf(complex<T> * x, complex<T> * y, complex<T>* z,int m,int n);
template<typename T>
void cmatrix_2d_r_subf(complex<T> * x, complex<T> * y, complex<T>* z, int m,int n);
template<typename T>
void cmatrix_2d_r_mulf(complex<T> * x, complex<T> * y, complex<T>* z, int m,int n);
template<typename T>
void cmatrix_2d_r_divf(complex<T> * x, complex<T> * y, complex<T>* z, int m,int n);

template<typename T>
void cmatrix_powf_scalar(complex<T> * x, complex<T> *y,  complex<T>* z, int m,int n);

template<typename T>
void cmatrix_r_acosf(complex<T> * devPtr, complex<T> * outputs, int m,int n);
template<typename T>
void cmatrix_r_asinf(complex<T> * devPtr, complex<T> * outputs, int m,int n);
template<typename T>
void cmatrix_r_atanf(complex<T> * devPtr, complex<T> * outputs, int m,int n);
template<typename T>
void cmatrix_r_acoshf(complex<T> * devPtr, complex<T> * outputs, int m,int n);
template<typename T>
void cmatrix_r_asinhf(complex<T> * devPtr, complex<T> * outputs, int m,int n);
template<typename T>
void cmatrix_r_atanhf(complex<T> * devPtr, complex<T> * outputs, int m,int n);
template<typename T>
void cmatrix_r_cosf(complex<T> * devPtr, complex<T> * outputs, int m,int n);
template<typename T>
void cmatrix_r_sinf(complex<T> * devPtr, complex<T> * outputs, int m,int n);
template<typename T>
void cmatrix_r_tanf(complex<T> * devPtr, complex<T> * outputs, int m,int n);
template<typename T>
void cmatrix_r_coshf(complex<T> * devPtr, complex<T> * outputs, int m,int n);
template<typename T>
void cmatrix_r_sinhf(complex<T> * devPtr, complex<T> * outputs, int m,int n);
template<typename T>
void cmatrix_r_tanhf(complex<T> * devPtr, complex<T> * outputs, int m,int n);
template<typename T>
void cmatrix_r_sinf(complex<T> * x, complex<T> *output, int m,int n);
template<typename T>
void cmatrix_r_sinhf(complex<T> * x, complex<T> *output, int m,int n);
template<typename T>
void cmatrix_r_tanf(complex<T> * x, complex<T> *output, int m,int n);
template<typename T>
void cmatrix_r_tanhf(complex<T> * x, complex<T> *output, int m,int n);

template<typename T>
void cmatrix_r_expf(complex<T> * devPtr, complex<T> * output, int m,int n);
template<typename T>
void cmatrix_r_logf(complex<T> * x, complex<T> *output, int m,int n);
template<typename T>
void cmatrix_r_log10f(complex<T> * x, complex<T> *output, int m,int n);
template<typename T>
void cmatrix_r_powf(complex<T> * x, complex<T> * y, complex<T> *output, int m,int n);
template<typename T>
void cmatrix_r_sqrtf(complex<T> * x, complex<T> *output, int m,int n);

template<typename T>
void cmatrix_r_sigmoidf(complex<T> * x, complex<T> *output, int m,int n);
template<typename T>
void cmatrix_r_sigmoid_gradf(complex<T> * x, complex<T> *output, int m,int n);
template<typename T>
void cmatrix_r_tanh_gradf(complex<T> * x, complex<T> *output, int m,int n);

template<typename T>
void cmatrix_r_addf_const(complex<T> * x, complex<T>  y, complex<T> *output, int m,int n);
template<typename T>
void cmatrix_r_subf_const(complex<T> * x, complex<T>  y, complex<T> *output,int m,int n);
template<typename T>
void cmatrix_r_mulf_const(complex<T> * x, complex<T>  y, complex<T> *output,int m,int n);
template<typename T>
void cmatrix_r_divf_const(complex<T> * x, complex<T>  y, complex<T> *output,int m,int n);
template<typename T>
void cmatrix_r_modf_const(complex<T> * x, complex<T>  y, complex<T> *output,int m,int n);
template<typename T>
void cmatrix_r_powf_const(complex<T> * x, complex<T> y, complex<T> *output,int m,int n);

template<typename T>
void cmatrix_r_addf_scalar(complex<T> * x, complex<T> * y, complex<T> *output,int m,int n);
template<typename T>
void cmatrix_r_subf_scalar(complex<T> * x, complex<T> * y, complex<T> *output,int m,int n);
template<typename T>
void cmatrix_r_mulf_scalar(complex<T> * x, complex<T> * y, complex<T> *output,int m,int n);
template<typename T>
void cmatrix_r_divf_scalar(complex<T> * x, complex<T> * y, complex<T> *output,int m,int n);


template<typename T>
void cmatrix_r_powf_scalar(complex<T> * x, complex<T> *y, complex<T> *output,int m,int n);

template<typename T> 
void cmatrix_r_hadamardf( complex<T> * a, complex<T> * b, complex<T> *output, int M, int N, int K);

template<typename T>
void cmatrix_r_multiplyf(complex<T> * a, complex<T> * b, complex<T> * output, int M, int N, int K);
template<typename T>
void cmatrix_r_transposef( complex<T> * input, complex<T> *output, int M, int N);
