#pragma once
#include <cuda/std/ccomplex>
//#include "cuda_complex.hpp"

template<typename T>
using complex = cuda::std::complex<T>;


template<typename T>
T vector_sumf(T * in, int size);
template<typename T>
T vector_prodf(T * in, int size);

template<typename T> void vector_addf_row(T * x, int row, T * y, int src,size_t n);
template<typename T> void vector_subf_row(T * x, int row, T * y, int src,size_t n);
template<typename T> void vector_mulf_row(T * x, int row, T * y, int src,size_t n);
template<typename T> void vector_divf_row(T * x, int row, T * y, int src,size_t n);
template<typename T> void vector_modf_row(T * x, int row, T * y, int src,size_t n);


template<typename T> void vector_setrowf(T * dst, int row, T * sc, int row_src, size_t n);
template<typename T> void vector_copyrowf(T * dst, int row, T * src, int row_src, int n);


// r registers are used now
template<typename T> void vector_r_truncf(T * x, T *output, int n);
template<typename T> void vector_r_copysignf(T * X, T *Y, T *output, int n);

template<typename T> void vector_r_addf(T * x, T * y, T * output, int n);
template<typename T> void vector_r_subf(T * x, T * y, T * output, int n);
template<typename T> void vector_r_mulf(T * x, T * y, T * output, int n);
template<typename T> void vector_r_divf(T * x, T * y, T * output, int n);
template<typename T> void vector_r_modf(T * x, T * y, T * output, int n);
template<typename T> void vector_r_acosf(T * devPtr, T * outputs, int n);
template<typename T> void vector_r_asinf(T * devPtr, T * outputs, int n);
template<typename T> void vector_r_atanf(T * devPtr, T * outputs, int n);
template<typename T> void vector_r_atan2f(T * a, T * b, T * output, int n);
template<typename T> void vector_r_acoshf(T * devPtr, T * outputs, int n);
template<typename T> void vector_r_asinhf(T * devPtr, T * outputs, int n);
template<typename T> void vector_r_atanhf(T * devPtr, T * outputs, int n);
template<typename T> void vector_r_cosf(T * devPtr, T * outputs, int n);
template<typename T> void vector_r_sinf(T * devPtr, T * outputs, int n);
template<typename T> void vector_r_tanf(T * devPtr, T * outputs, int n);
template<typename T> void vector_r_coshf(T * devPtr, T * outputs, int n);
template<typename T> void vector_r_sinhf(T * devPtr, T * outputs, int n);
template<typename T> void vector_r_tanhf(T * devPtr, T * outputs, int n);
template<typename T> void vector_r_ceilf(T * devPtr, T * output, int n);
template<typename T> void vector_r_exp10f(T * devPtr, T * outputs, int n);
template<typename T> void vector_r_exp2f(T * devPtr, T * output, int n);
template<typename T> void vector_r_expf(T * devPtr, T * output, int n);
template<typename T> void vector_r_expm1f(T * devPtr, T * output, int n);
template<typename T> void vector_r_fabsf(T * devPtr, T * output, int n);
template<typename T> void vector_r_floorf(T * devPtr, T * output, int n);
template<typename T> void vector_r_fmaxf(T * x, T * y, T * output, int n);
template<typename T> void vector_r_fminf(T * x, T * y, T * output, int n);
template<typename T> void vector_r_fmodf(T * x, T * y, T * output, int n);
template<typename T> void vector_r_log10f(T * x, T *output, int n);
template<typename T> void vector_r_log1pf(T * x, T *output, int n);
template<typename T> void vector_r_log2f(T * x, T * output, int n);
template<typename T> void vector_r_logbf(T * x, T *output, int n);
template<typename T> void vector_r_powf(T * x, T * y, T *output, int n);
template<typename T> void vector_r_rsqrtf(T * x, T * output, int n);
template<typename T> void vector_r_sinf(T * x, T *output, int n);
template<typename T> void vector_r_sinhf(T * x, T *output, int n);
template<typename T> void vector_r_sqrtf(T * x, T *output, int n);
template<typename T> void vector_r_tanf(T * x, T *output, int n);
template<typename T> void vector_r_tanhf(T * x, T *output, int n);
template<typename T> void vector_r_softmaxf(T * x, T *output, int n);
template<typename T> void vector_r_sigmoidf(T * x, T *output, int n);
template<typename T> void vector_r_sigmoid_gradf(T * x, T *output, int n);
template<typename T> void vector_r_tanh_gradf(T * x, T *output, int n);
template<typename T> void vector_r_reluf(T * x, T *output, int n);
template<typename T> void vector_r_relu_gradf(T * x, T *output, int n);
template<typename T> void vector_r_cbrtf(T * devPtr, T * output, int n);
template<typename T> void vector_r_cospif(T * devPtr, T * output, int n);
template<typename T> void vector_r_cyl_bessel_i0f(T * devPtr, T * output, int n);
template<typename T> void vector_r_cyl_bessel_i1f(T * devPtr, T * output, int n);
template<typename T> void vector_r_erfcf(T * devPtr, T * output, int n);
template<typename T> void vector_r_erfcinvf(T * devPtr, T * output, int n);
template<typename T> void vector_r_erfcxf(T * devPtr, T * output, int n);
template<typename T> void vector_r_erff(T * devPtr, T * output, int n);
template<typename T> void vector_r_erfinvf(T * devPtr, T * output, int n);
template<typename T> void vector_r_fdimf(T * a, T * b, T * output, int n);
template<typename T> void vector_r_fdividef(T * a, T * b, T * output, int n);
template<typename T> void vector_r_fmaf(T * x, T * y, T * z, T *output, int n);
template<typename T> void vector_r_hypotf(T * x, T * y, T * output, int n);
template<typename T> void vector_r_ilogbf(T * x, T *output, int n);
template<typename T> void vector_r_j0f(T * x, T *output, int n);
template<typename T> void vector_r_j1f(T * x, T *output, int n);
template<typename T> void vector_r_jnf(T * x, T * output, int M, int n);
template<typename T> void vector_r_ldexpf(T * x, T * output, int exp, int n);
template<typename T> void vector_r_lgammaf(T * x, T *output, int n);
template<typename T> void vector_r_nearbyintf(T * x, T *output, int n);
template<typename T> void vector_r_norm3df(T * x, T * y, T * z, T * output, int n);
template<typename T> void vector_r_norm4df(T * x, T * y, T * z, T * q, T * output, int n);
template<typename T> void vector_r_normcdff(T * x, T * output, int n);
template<typename T> void vector_r_normcdfinvf(T * x, T *output, int n);
template<typename T> void vector_r_normf(int dim, T * x, T * output, int n);
template<typename T> void vector_r_rcbrtf(T * x, T *output, int n);
template<typename T> void vector_r_remainderf(T * x, T * y, T *output, int n);
template<typename T> void vector_r_rhypotf(T * x, T * y, T *output, int n);
template<typename T> void vector_r_rnorm3df(T * x, T * y, T * z, T * output, int n);
template<typename T> void vector_r_rnorm4df(T * x, T * y, T * z, T * q, T *output, int n);
template<typename T> void vector_r_rnormf(int dim, T * x, T *output, int n);
template<typename T> void vector_r_scalblnf(T * x, long int M, T * output, int n);
template<typename T> void vector_r_tgammaf(T * x, T * output, int n);
template<typename T> void vector_r_truncf(T * x, T *output, int n);
template<typename T> void vector_r_y0f(T * x, T *output, int n);
template<typename T> void vector_r_y1f(T * x, T * output, int n);
template<typename T> void vector_r_ynf(int M, T * x, T *output, int n);
template<typename T> void vector_r_sinpif(T * x, T *output, int n);


template<typename T> void vector_r_addf_const(T * x, T  y, T *output, int n);
template<typename T> void vector_r_subf_const(T * x, T  y, T *output,int n);
template<typename T> void vector_r_mulf_const(T * x, T  y, T *output,int n);
template<typename T> void vector_r_divf_const(T * x, T  y, T *output,int n);
template<typename T> void vector_r_modf_const(T * x, T  y, T *output,int n);
template<typename T> void vector_r_atan2f_const(T * a, T b, T *output,int n);
template<typename T> void vector_r_fmaxf_const(T * x, T  y, T *output,int n);
template<typename T> void vector_r_fminf_const(T * x, T  y, T *output,int n);
template<typename T> void vector_r_fmodf_const(T * x, T  y, T *output,int n);
template<typename T> void vector_r_powf_const(T * x, T y, T *output,int n);


template<typename T> void vector_r_addf_scalar(T * x, T * y, T *output,int n);
template<typename T> void vector_r_subf_scalar(T * x, T * y, T *output,int n);
template<typename T> void vector_r_mulf_scalar(T * x, T * y, T *output,int n);
template<typename T> void vector_r_divf_scalar(T * x, T * y, T *output,int n);
template<typename T> void vector_r_modf_scalar(T * x, T * y, T *output,int n);
template<typename T> void vector_r_atan2f_scalar(T * a, T *b, T *output,int n);
template<typename T> void vector_r_fmaxf_scalar(T * x, T  *y, T *output,int n);
template<typename T> void vector_r_fminf_scalar(T * x, T  *y, T *output,int n);
template<typename T> void vector_r_fmodf_scalar(T * x, T  *y, T *output,int n);
template<typename T> void vector_r_powf_scalar(T * x, T *y, T *output,int n);

template<typename T> void vector_r_fdimf_const(T * a, T  b, T *output,int n);
template<typename T> void vector_r_fdividef_const(T * a, T  b, T *output,int n);
template<typename T> void vector_r_hypotf_const(T * x, T  y, T *output,int n);
template<typename T> void vector_r_remainderf_const(T * x, T y, T *output,int n);
template<typename T> void vector_r_rhypotf_const(T * x, T y, T *output,int n);

template<typename T> void vector_r_fdimf_scalar(T * a, T  *b, T *output,int n);
template<typename T> void vector_r_fdividef_scalar(T * a, T *b, T *output,int n);
template<typename T> void vector_r_hypotf_scalar(T * x, T  *y, T *output,int n);
template<typename T> void vector_r_remainderf_scalar(T * x, T *y, T *output,int n);
template<typename T> void vector_r_rhypotf_scalar(T * x, T *y, T *output,int n);

/*
// i do not t
template<typename T>
T vector_realf(complex<T> * x, T * y, int n);
template<typename T>
T vector_imagf(complex<T> * x, T * y, int n);

template<typename T>
complex<T>* vector_argf(complex<T> * devPtr, T * x, int n);
template<typename T>
complex<T>* vector_normf(complex<T> * devPtr, complex<T> * x, int n);
template<typename T> 
void vector_r_normf(complex<T> * x, complex<T> * output, int n);
template<typename T>
complex<T>* vector_conjf(complex<T> * devPtr, complex<T> * x,int n);
template<typename T>
complex<T>* vector_projf(complex<T> * devPtr, complex<T> * x, int n);
template<typename T>
complex<T>* vector_polarf(T * r, T *theta, complex<T> * out, int n);
*/
template<typename T>
void    return_memory(int length, T *fp);
template<typename T>
T*  find_memory(int length);
template<typename T>
void    add_memory(int length, T * ptr);
void    clear_cache();
void    calcSize(int N,int * gridSize, int * blockSize);

template<typename T>
void cuda_zero(T * dst, int n);
template<typename T>
void cuda_memcpy(T * dst, T * src, int n);
