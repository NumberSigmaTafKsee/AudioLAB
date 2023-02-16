
#pragma once
#include <cuda/std/ccomplex>

template<typename T> using complex = cuda::std::complex<T>;

template<typename T>
void cmatrix_2d_r_mulf(T * a, T * b, T *output, int M, int N);
template<typename T>
void cmatrix_2d_r_addf(T * a, T * b, T *output,int M, int N);
template<typename T>
void cmatrix_2d_r_subf(T * a, T * b, T *output,int M, int N);
template<typename T>
void cmatrix_2d_r_divf(T * a, T * b, T *output,int M, int N);
template<typename T>
void cmatrix_2d_r_modf(T * a, T * b, T *output,int M, int N);


template<typename T> complex<T>* cvector_addf(complex<T> * x, complex<T> * y, int n);
template<typename T> complex<T>* cvector_subf(complex<T> * x, complex<T> * y, int n);
template<typename T> complex<T>* cvector_mulf(complex<T> * x, complex<T> * y, int n);
template<typename T> complex<T>* cvector_divf(complex<T> * x, complex<T> * y, int n);

template<typename T> void cvector_realf(complex<T> * x, T * y, int n);
template<typename T> void cvector_imagf(complex<T> * x, T * y, int n);
template<typename T>  void cvector_argf(complex<T> * devPtr, T * y, int n);
template<typename T>  void cvector_normf(complex<T> * devPtr, complex<T> * y, int n);
template<typename T>  void cvector_conjf(complex<T> * devPtr, complex<T> * y, int n);
template<typename T>  void cvector_projf(complex<T> * devPtr, complex<T> * y, int n);
//template<typename T>  void cvector_polar(T * r, T *theta,  complex<T> * y, int n);

template<typename T> complex<T>* cvector_expf(complex<T> * devPtr, int n);
template<typename T> complex<T>* cvector_logf(complex<T> * x, int n);
template<typename T> complex<T>* cvector_log10f(complex<T> * x, int n);

template<typename T> complex<T>* cvector_acosf(complex<T> * devPtr, int n);
template<typename T> complex<T>* cvector_acoshf(complex<T> * devPtr, int n);
template<typename T> complex<T>* cvector_asinf(complex<T> * devPtr, int n);
template<typename T> complex<T>* cvector_asinhf(complex<T> * devPtr, int n);
template<typename T> complex<T>* cvector_atanf(complex<T> * devPtr, int n);
template<typename T> complex<T>* cvector_atanhf(complex<T> * devPtr, int n);
template<typename T> complex<T>* cvector_cosf(complex<T> * devPtr, int n);
template<typename T> complex<T>* cvector_coshf(complex<T> * devPtr, int n);
template<typename T> complex<T>* cvector_sinf(complex<T> * x, int n);
template<typename T> complex<T>* cvector_sinhf(complex<T> * x, int n);
template<typename T> complex<T>* cvector_tanf(complex<T> * x, int n);
template<typename T> complex<T>* cvector_tanhf(complex<T> * x, int n);
template<typename T> complex<T>* cvector_powf(complex<T> * x, complex<T> * y, int n);
template<typename T> complex<T>* cvector_sqrtf(complex<T> * x, int n);

template<typename T> complex<T>* cvector_sigmoidf(complex<T> * devPtr, int n);
template<typename T> complex<T>* cvector_sigmoid_gradf(complex<T> * devPtr, int n);
template<typename T> complex<T>* cvector_tanh_gradf(complex<T> * devPtr, int n);

template<typename T> complex<T>* cvector_addf_const(complex<T> * x, complex<T>  y, int n);
template<typename T> complex<T>* cvector_subf_const(complex<T> * x, complex<T>  y, int n);
template<typename T> complex<T>* cvector_mulf_const(complex<T> * x, complex<T>  y, int n);
template<typename T> complex<T>* cvector_divf_const(complex<T> * x, complex<T>  y, int n);

template<typename T> complex<T>* cvector_powf_const(complex<T> * x, complex<T> y, int n);

template<typename T> complex<T>* cvector_addf_scalar(complex<T> * x, complex<T> * y, int n);
template<typename T> complex<T>* cvector_subf_scalar(complex<T> * x, complex<T> * y, int n);
template<typename T> complex<T>* cvector_mulf_scalar(complex<T> * x, complex<T> * y, int n);
template<typename T> complex<T>* cvector_divf_scalar(complex<T> * x, complex<T> * y, int n);

template<typename T> complex<T>* cvector_powf_scalar(complex<T> * x, complex<T> *y, int n);


template<typename T> void cvector_addf_row(complex<T> * x, int row, complex<T> * y, int src,size_t n);
template<typename T> void cvector_subf_row(complex<T> * x, int row, complex<T> * y, int src,size_t n);
template<typename T> void cvector_mulf_row(complex<T> * x, int row, complex<T> * y, int src,size_t n);
template<typename T> void cvector_divf_row(complex<T> * x, int row, complex<T> * y, int src,size_t n);

template<typename T> void cvector_setrowf(complex<T> * dst, int row, complex<T> * sc, int row_src, size_t n);
template<typename T> void cvector_copyrowf(complex<T> * dst, int row, complex<T> * src, int row_src, int n);

template<typename T> void cvector_r_addf(complex<T> * x, complex<T> * y, complex<T> * output, int n);
template<typename T> void cvector_r_subf(complex<T> * x, complex<T> * y, complex<T> * output, int n);
template<typename T> void cvector_r_mulf(complex<T> * x, complex<T> * y, complex<T> * output, int n);
template<typename T> void cvector_r_divf(complex<T> * x, complex<T> * y, complex<T> * output, int n);


template<typename T> void cvector_r_acosf(complex<T> * devPtr, complex<T> * outputs, int n);
template<typename T> void cvector_r_asinf(complex<T> * devPtr, complex<T> * outputs, int n);
template<typename T> void cvector_r_atanf(complex<T> * devPtr, complex<T> * outputs, int n);
template<typename T> void cvector_r_acoshf(complex<T> * devPtr, complex<T> * outputs, int n);
template<typename T> void cvector_r_asinhf(complex<T> * devPtr, complex<T> * outputs, int n);
template<typename T> void cvector_r_atanhf(complex<T> * devPtr, complex<T> * outputs, int n);
template<typename T> void cvector_r_cosf(complex<T> * devPtr, complex<T> * outputs, int n);
template<typename T> void cvector_r_sinf(complex<T> * devPtr, complex<T> * outputs, int n);
template<typename T> void cvector_r_tanf(complex<T> * devPtr, complex<T> * outputs, int n);
template<typename T> void cvector_r_coshf(complex<T> * devPtr, complex<T> * outputs, int n);
template<typename T>  void cvector_r_sinhf(complex<T> * devPtr, complex<T> * outputs, int n);
template<typename T>  void cvector_r_tanhf(complex<T> * devPtr, complex<T> * outputs, int n);
template<typename T>  void cvector_r_sinf(complex<T> * x, complex<T> *output, int n);
template<typename T>  void cvector_r_sinhf(complex<T> * x, complex<T> *output, int n);
template<typename T>  void cvector_r_tanf(complex<T> * x, complex<T> *output, int n);
template<typename T>  void cvector_r_tanhf(complex<T> * x, complex<T> *output, int n);

template<typename T>  void cvector_r_expf(complex<T> * devPtr, complex<T> * output, int n);
template<typename T>  void cvector_r_log10f(complex<T> * x, complex<T> *output, int n);
template<typename T>  void cvector_r_powf(complex<T> * x, complex<T> * y, complex<T> *output, int n);
template<typename T>  void cvector_r_sqrtf(complex<T> * x, complex<T> *output, int n);

template<typename T>  void cvector_r_sigmoidf(complex<T> * x, complex<T> *output, int n);
template<typename T>  void cvector_r_sigmoid_gradf(complex<T> * x, complex<T> *output, int n);
template<typename T>  void cvector_r_tanh_gradf(complex<T> * x, complex<T> *output, int n);

template<typename T>  void cvector_r_addf_const(complex<T> * x, complex<T>  y, complex<T> *output, int n);
template<typename T>  void cvector_r_subf_const(complex<T> * x, complex<T>  y, complex<T> *output,int n);
template<typename T>  void cvector_r_mulf_const(complex<T> * x, complex<T>  y, complex<T> *output,int n);
template<typename T>  void cvector_r_divf_const(complex<T> * x, complex<T>  y, complex<T> *output,int n);
template<typename T>  void cvector_r_modf_const(complex<T> * x, complex<T>  y, complex<T> *output,int n);

template<typename T>  void cvector_r_powf_const(complex<T> * x, complex<T> y, complex<T> *output,int n);
template<typename T>  void cvector_r_addf_scalar(complex<T> * x, complex<T> * y, complex<T> *output,int n);
template<typename T>  void cvector_r_subf_scalar(complex<T> * x, complex<T> * y, complex<T> *output,int n);
template<typename T>  void cvector_r_mulf_scalar(complex<T> * x, complex<T> * y, complex<T> *output,int n);
template<typename T>  void cvector_r_divf_scalar(complex<T> * x, complex<T> * y, complex<T> *output,int n);
template<typename T>  void cvector_r_powf_scalar(complex<T> * x, complex<T> *y, complex<T> *output,int n);


template<typename T>  void cvector_r_hadamardf(complex<T> * a, complex<T> * b, complex<T> *output, int N);
