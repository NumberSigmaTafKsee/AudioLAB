#ifndef __MATRIXFLOAT_H
#define __MATRIXFLOAT_H

#ifdef __cplusplus 
extern "C" {
#endif 


float* _2d_mulf(float * a, float * b, int M, int N);
float* _2d_addf(float * a, float * b, int M, int N);
float* _2d_subf(float * a, float * b, int M, int N);
float* _2d_divf(float * a, float * b, int M, int N);
float* _2d_modf(float * a, float * b, int M, int N);

void _2d_r_mulf(float * a, float * b, float *output, int M, int N);
void _2d_r_addf(float * a, float * b, float *output,int M, int N);
void _2d_r_subf(float * a, float * b, float *output,int M, int N);
void _2d_r_divf(float * a, float * b, float *output,int M, int N);
void _2d_r_modf(float * a, float * b, float *output,int M, int N);


// goes in matrix_float.h
float* matrix_transposef(float * input, int M, int N);
float* matrix_hadamardf(float * a, float * b, int M, int N, int K);
float* matrix_multiplyf(float * a, float * b, int M, int N, int K);
float* matrix_addf(float * a, float * b, int M, int N, int K);
float* matrix_subf(float * a, float * b, int M, int N, int K);


float* matrix_acosf(float * input, int M, int N);
float* matrix_asinf(float * input, int M, int N);
float* matrix_atanf(float * input, int M, int N);
float* matrix_atan2f(float * x, float * y, int M, int N);

float* matrix_cosf(float * input, int M, int N);
float* matrix_sinf(float * input, int M, int N);
float* matrix_tanf(float * input, int M, int N);

float* matrix_acoshf(float * input, int M, int N);
float* matrix_asinhf(float * input, int M, int N);
float* matrix_atanhf(float * input, int M, int N);

float* matrix_coshf(float * input, int M, int N);
float* matrix_sinhf(float * input, int M, int N);
float* matrix_tanhf(float * input, int M, int N);

float* matrix_ceilf(float * devPtr, int M, int N);

float* matrix_exp10f(float * devPtr, int M, int N);
float* matrix_exp2f(float * devPtr, int M, int N);
float* matrix_expf(float * devPtr, int M, int N);
float* matrix_expm1f(float * devPtr,int M, int N);

float* matrix_fabsf(float * devPtr,int M, int N);
float* matrix_floorf(float * devPtr,int M, int N);
float* matrix_fmaxf(float * x, float * y, int M, int N);
float* matrix_fminf(float * x, float * y,int M, int N);
float* matrix_fmodf(float * x, float * y,int M, int N);
float* matrix_log10f(float * x, int M, int N);
float* matrix_log1pf(float * x, int M, int N);
float* matrix_log2f(float * x, int M, int N);
float* matrix_logbf(float * x, int M, int N);
float* matrix_powf(float * x, float * y, int M, int N);
float* matrix_rsqrtf(float * x, int M, int N);
float* matrix_sqrtf(float * x, int M, int N);
float* matrix_sigmoidf(float * devPtr,int M, int N);
float* matrix_sigmoid_gradf(float * devPtr, int M, int N);
float* matrix_tanh_gradf(float * devPtr, int M, int N);
float* matrix_reluf(float * devPtr, int M, int N);
float* matrix_relu_gradf(float * devPtr, int M, int N);
float* matrix_softmaxf(float * x, int M, int N);



float* matrix_cbrtf(float * devPtr, int M, int N);
float* matrix_cospif(float * devPtr, int M, int N);
float* matrix_cyl_bessel_i0f(float * devPtr, int M, int N);
float* matrix_cyl_bessel_i1f(float * devPtr, int M, int N);
float* matrix_erfcf(float * devPtr,int M, int N);
float* matrix_erfcinvf(float * devPtr,int M, int N);
float* matrix_erfcxf(float * devPtr, int M, int N);
float* matrix_erff(float * devPtr, int M, int N);
float* matrix_erfinvf(float * devPtr, int M, int N);
float* matrix_fdimf(float * a, float * b, int M, int N);
float* matrix_fdividef(float * a, float * b, int M, int N);
float* matrix_fmaf(float * x, float * y, float * z, int M, int N);
float* matrix_hypotf(float * x, float * y,int M, int N);
float* matrix_ilogbf(float * x,int M, int N);
float* matrix_j0f(float * x,int M, int N);
float* matrix_j1f(float * x,int M, int N);
float* matrix_jnf(float * x, int n, int M, int N);
float* matrix_ldexpf(float * x, int exp, int M, int N);
float* matrix_lgammaf(float * x, int M, int N);
float* matrix_nearbyintf(float * x, int m, int n);
float* matrix_norm3df(float * x, float * y, float * z, int M, int N);
float* matrix_norm4df(float * x, float * y, float * z, float * q, int M, int N);
float* matrix_normcdff(float * x,int M, int N);
float* matrix_normcdfinvf(float * x,int M, int N);
float* matrix_normf(int dim, float * x, int M, int N);
float* matrix_remainderf(float * x, float * y,int M, int N);
float* matrix_rcbrtf(float * x,int M, int N);
float* matrix_rhypotf(float * x, float * y,int M, int N);
float* matrix_rnorm3df(float * x, float * y, float * z,int M, int N);
float* matrix_rnorm4df(float * x, float * y, float * z, float * q, int M, int N);
float* matrix_rnormf(int dim, float * x, int M, int N);
float* matrix_scalblnf(float * x, long int n, int M, int N);
float* matrix_sinpif(float * x, int M, int N);
float* matrix_tgammaf(float * x, int M, int N);
float* matrix_y0f(float * x, int M, int N);
float* matrix_y1f(float * x, int M, int N);
float* matrix_ynf(int n, float * x, int M, int N);


void matrix_r_addf(float * a, float * b, float *output, int M, int N, int K);
void matrix_r_subf(float * a, float * b, float *output, int M, int N, int K);
void matrix_r_hadamardf(float * a, float * b, float *output, int M, int N, int K);
void matrix_r_multiplyf(float * a, float * b, float * output, int M, int N, int K);
void matrix_r_transposef(float * input, float *output, int M, int N);
void matrix_r_acosf(float * input, float * output,  int M, int N);
void matrix_r_asinf(float * input, float * output,  int M, int N);
void matrix_r_atanf(float * input, float * output,  int M, int N);
void matrix_r_atan2f(float * x, float * y, float * output, int M, int N);
void matrix_r_acoshf(float * input, float * output,  int M, int N);
void matrix_r_asinhf(float * input, float * output,  int M, int N);
void matrix_r_atanhf(float * input, float * output,  int M, int N);
void matrix_r_cosf(float * input, float * output,  int M, int N);
void matrix_r_sinf(float * input, float * output,  int M, int N);
void matrix_r_tanf(float * input, float * output, int M, int N);
void matrix_r_coshf(float * input, float * output,  int M, int N);
void matrix_r_sinhf(float * input, float * output,  int M, int N);
void matrix_r_tanhf(float * input, float * output, int M, int N);
void matrix_r_atan2f_const(float * input, float b, float * output, int M, int N);
void matrix_r_ceilf(float * input, float *output, int M, int N);
void matrix_r_exp10f(float * input, float *output, int M, int N);
void matrix_r_exp2f(float * input, float *output, int M, int N);
void matrix_r_expf(float * input, float *output, int M, int N);
void matrix_r_expm1f(float * input, float *output, int M, int N);
void matrix_r_fabsf(float * input, float *output, int M, int N);
void matrix_r_floorf(float * input, float *output, int M, int N);
void matrix_r_fmaxf(float * x, float *y, float *output, int M, int N);
void matrix_r_fminf(float * x, float *y, float *output, int M, int N);
void matrix_r_fmodf(float * x, float *y, float *output, int M, int N);
void matrix_r_log10f(float * input, float *output, int M, int N);
void matrix_r_log1pf(float * input, float *output, int M, int N);
void matrix_r_log2f(float * input, float *output, int M, int N);
void matrix_r_logbf(float * input, float *output, int M, int N);
void matrix_r_powf(float * x, float *y, float *output, int M, int N);
void matrix_r_rsqrtf(float * input, float *output, int M, int N);
void matrix_r_sqrtf(float * input, float *output, int M, int N);
void matrix_r_cbrtf(float * input, float *output, int M, int N);
void matrix_r_cospif(float * input, float *output, int M, int N);
void matrix_r_cyl_bessel_i0f(float * input, float *output, int M, int N);
void matrix_r_cyl_bessel_i1f(float * input, float *output, int M, int N);
void matrix_r_erfcf(float * input, float *output, int M, int N);
void matrix_r_erfcinvf(float * input, float *output, int M, int N);
void matrix_r_erfcxf(float * input, float *output, int M, int N);
void matrix_r_erff(float * input, float *output, int M, int N);
void matrix_r_erfinvf(float * input, float * outputs, int M, int N);
void matrix_r_fdimf(float * x, float * y, float *output, int M, int N);
void matrix_r_fdividef(float * x, float *y, float *output, int M, int N);
void matrix_r_fmaf(float * x, float *y, float *z, float *output, int M, int N);
void matrix_r_hypotf(float * x, float *y, float *output, int M, int N);
void matrix_r_ilogbf(float * input, float *output, int M, int N);
void matrix_r_j0f(float * input, float *output, int M, int N);
void matrix_r_j1f(float * input, float *output, int M, int N);
void matrix_r_jnf(float * input, float *output, int n, int M, int N);
void matrix_r_ldexpf(float * input, float *output, int exp, int M, int N);
void matrix_r_lgammaf(float * input, float *output, int M, int N);
void matrix_r_nearbyintf(float * input, float *output, int M, int N);
void matrix_r_norm3df(float * x, float *y, float *z, float *output, int M, int N);
void matrix_r_norm4df(float * x, float *y, float *z, float * w, float *output,int M, int N);
void matrix_r_normcdff(float * input, float *output, int M, int N);
void matrix_r_normcdfinvf(float * input, float *output, int M, int N);
void matrix_r_normf(int dim, float *input, float *output, int M, int N);
void matrix_r_remainderf(float * x, float *y, float *output, int M, int N);
void matrix_r_rcbrtf(float * input, float * output, int M, int N);
void matrix_r_rhypotf(float * x, float *y, float * output, int M, int N);
void matrix_r_rnorm3df(float * x, float *y, float *z, float * output, int M, int N);
void matrix_r_rnorm4df(float * x, float *y, float *z, float *w, float *output, int M, int N);
void matrix_r_rnormf(int dim, float *input, float *output, int M, int N);
void matrix_r_scalblnf(float * input, long int n,  float *output, int M, int N);
void matrix_r_sinpif(float * input, float *output, int M, int N);
void matrix_r_tgammaf(float * input, float *output, int M, int N);
void matrix_r_y0f(float * input, float *output, int M, int N);
void matrix_r_y1f(float * input, float *output, int M, int N);
void matrix_r_ynf(int n, float * input, float *output, int M, int N);

void matrix_r_softmaxf(float * input, float *output, int M, int N);
void matrix_r_sigmoidf(float * x, float *output, int M, int N);
void matrix_r_sigmoid_gradf(float * x, float *output,int M, int N);
void matrix_r_tanh_gradf(float * x, float *output, int M, int N);
void matrix_r_reluf(float * x, float *output, int M, int N);
void matrix_r_relu_gradf(float * x, float *output, int M, int N);


float* matrix_addf_const(float * x, float  y, int M, int N);
float* matrix_subf_const(float * x, float  y, int M, int N);
float* matrix_mulf_const(float * x, float  y, int M, int N);
float* matrix_divf_const(float * x, float  y, int M, int N);
float* matrix_modf_const(float * x, float  y, int M, int N);
float* matrix_atan2f_const(float * a, float b, int M, int N);
float* matrix_fmaxf_const(float * x, float  y, int M, int N);
float* matrix_fminf_const(float * x, float  y, int M, int N);
float* matrix_fmodf_const(float * x, float  y, int M, int N);
float* matrix_powf_const(float * x, float y, int M, int N);
float* matrix_fdimf_const(float * a, float  b, int M, int N);
float* matrix_fdividef_const(float * a, float  b, int M, int N);
float* matrix_hypotf_const(float * x, float  y, int M, int N);
float* matrix_remainderf_const(float * x, float y, int M, int N);
float* matrix_rhypotf_const(float * x, float y,int M, int N);

float* matrix_addf_scalar(float * x, float * y, int M, int N);
float* matrix_subf_scalar(float * x, float * y, int M, int N);
float* matrix_mulf_scalar(float * x, float * y, int M, int N);
float* matrix_divf_scalar(float * x, float * y, int M, int N);
float* matrix_modf_scalar(float * x, float * y, int M, int N);
float* matrix_atan2f_scalar(float * a, float *b, int M, int N);
float* matrix_fmaxf_scalar(float * x, float *y, int M, int N);
float* matrix_fminf_scalar(float * x, float *y, int M, int N);
float* matrix_fmodf_scalar(float * x, float *y, int M, int N);
float* matrix_powf_scalar(float * x, float *y, int M, int N);
float *matrix_fdimf_scalar(float * x, float *y, int M, int N);
float *matrix_fdividef_scalar(float * x, float *y, int M, int N);
float *matrix_hypotf_scalar(float * x, float *y,  int M, int N);
float* matrix_remainderf_scalar(float * x, float *y, int M, int N);
float* matrix_rhypotf_scalar(float * x, float *y, int M, int N);

void matrix_r_addf_const(float * x, float y, float *output, int M, int N);
void matrix_r_subf_const(float * x, float y, float *output, int M, int N);
void matrix_r_mulf_const(float * x, float y, float *output,int M, int N);
void matrix_r_divf_const(float * x, float y, float *output, int M, int N);
void matrix_r_modf_const(float * x, float y, float *output, int M, int N);
void matrix_r_fmaxf_const(float * x, float y, float *output, int M, int N);
void matrix_r_fminf_const(float * x, float y, float *output, int M, int N);
void matrix_r_fmodf_const(float * x, float y, float *output, int M, int N);
void matrix_r_powf_const(float * x, float y, float *output, int M, int N);
void matrix_r_fdimf_const(float * x, float y, float *output, int M, int N);
void matrix_r_fdividef_const(float * x, float y, float *output, int M, int N);
void matrix_r_hypotf_const(float * x, float y, float *output, int M, int N);
void matrix_r_remainderf_const(float * x, float y, float * output, int M, int N);
void matrix_r_rhypotf_const(float * x, float y, float *output, int M, int N);

void matrix_r_addf_scalar(float * x, float *y, float *output, int M, int N);
void matrix_r_subf_scalar(float * x, float *y, float *output, int M, int N);
void matrix_r_mulf_scalar(float * x, float *y, float *output, int M, int N);
void matrix_r_divf_scalar(float * x, float *y, float *output, int M, int N);
void matrix_r_modf_scalar(float * x, float *y, float *output, int M, int N);
void matrix_r_atan2f_scalar(float * a, float *b, float *output,int M, int N);
void matrix_r_fmaxf_scalar(float * x, float *y, float *output,int M, int N);
void matrix_r_fminf_scalar(float * x, float *y, float *output,int M, int N);
void matrix_r_fmodf_scalar(float * x, float *y, float *output,int M, int N);
void matrix_r_powf_scalar(float * x, float *y, float *output,int M, int N);
void matrix_r_fdimf_scalar(float * x, float *y, float *output, int M, int N);
void matrix_r_fdividef_scalar(float * x, float *y, float *output, int M, int N);
void matrix_r_hypotf_scalar(float * x, float *y, float *output, int M, int N);
void matrix_r_remainderf_scalar(float * x, float *y, float * output, int M, int N);
void matrix_r_rhypotf_scalar(float * x, float *y, float *output, int M, int N);

float* matrix_copysignf(float * X, float *Y, int M, int N);
void matrix_r_copysignf(float * X, float *Y, float *output, int M, int N);

float* matrix_truncf(float * x, int M, int N);
void matrix_r_truncf(float * x, float *output, int M, int N);


#ifdef __cplusplus
}
#endif 

#endif