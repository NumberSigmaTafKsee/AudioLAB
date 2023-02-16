#pragma once


template<typename T>
void _2d_r_mulf(T * a, T * b, T *output, int M, int N);
template<typename T>
void _2d_r_addf(T * a, T * b, T *output,int M, int N);
template<typename T>
void _2d_r_subf(T * a, T * b, T *output,int M, int N);
template<typename T>
void _2d_r_divf(T * a, T * b, T *output,int M, int N);
template<typename T>
void _2d_r_modf(T * a, T * b, T *output,int M, int N);



template<typename T>
void matrix_r_addf(T * a, T * b, T *output, int M, int N, int K);
template<typename T>
void matrix_r_subf(T * a, T * b, T *output, int M, int N, int K);
template<typename T>
void matrix_r_hadamardf(T * a, T * b, T *output, int M, int N, int K);
template<typename T>
void matrix_r_multiplyf(T * a, T * b, T * output, int M, int N, int K);
template<typename T>
void matrix_r_transposef(T * input, T *output, int M, int N);
template<typename T>
void matrix_r_acosf(T * input, T * output,  int M, int N);
template<typename T>
void matrix_r_asinf(T * input, T * output,  int M, int N);
template<typename T>
void matrix_r_atanf(T * input, T * output,  int M, int N);
template<typename T>
void matrix_r_atan2f(T * x, T * y, T * output, int M, int N);
template<typename T>
void matrix_r_acoshf(T * input, T * output,  int M, int N);
template<typename T>
void matrix_r_asinhf(T * input, T * output,  int M, int N);
template<typename T>
void matrix_r_atanhf(T * input, T * output,  int M, int N);
template<typename T>
void matrix_r_cosf(T * input, T * output,  int M, int N);
template<typename T>
void matrix_r_sinf(T * input, T * output,  int M, int N);
template<typename T>
void matrix_r_tanf(T * input, T * output, int M, int N);
template<typename T>
void matrix_r_coshf(T * input, T * output,  int M, int N);
template<typename T>
void matrix_r_sinhf(T * input, T * output,  int M, int N);
template<typename T>
void matrix_r_tanhf(T * input, T * output, int M, int N);
template<typename T>
void matrix_r_atan2f_const(T * input, T b, T * output, int M, int N);
template<typename T>
void matrix_r_ceilf(T * input, T *output, int M, int N);
template<typename T>
void matrix_r_exp10f(T * input, T *output, int M, int N);
template<typename T>
void matrix_r_exp2f(T * input, T *output, int M, int N);
template<typename T>
void matrix_r_expf(T * input, T *output, int M, int N);
template<typename T>
void matrix_r_expm1f(T * input, T *output, int M, int N);
template<typename T>
void matrix_r_fabsf(T * input, T *output, int M, int N);
template<typename T>
void matrix_r_floorf(T * input, T *output, int M, int N);
template<typename T>
void matrix_r_fmaxf(T * x, T *y, T *output, int M, int N);
template<typename T>
void matrix_r_fminf(T * x, T *y, T *output, int M, int N);
template<typename T>
void matrix_r_fmodf(T * x, T *y, T *output, int M, int N);
template<typename T>
void matrix_r_log10f(T * input, T *output, int M, int N);
template<typename T>
void matrix_r_log1pf(T * input, T *output, int M, int N);
template<typename T>
void matrix_r_log2f(T * input, T *output, int M, int N);
template<typename T>
void matrix_r_logbf(T * input, T *output, int M, int N);
template<typename T>
void matrix_r_powf(T * x, T *y, T *output, int M, int N);
template<typename T>
void matrix_r_rsqrtf(T * input, T *output, int M, int N);
template<typename T>
void matrix_r_sqrtf(T * input, T *output, int M, int N);
template<typename T>
void matrix_r_cbrtf(T * input, T *output, int M, int N);
template<typename T>
void matrix_r_cospif(T * input, T *output, int M, int N);
template<typename T>
void matrix_r_cyl_bessel_i0f(T * input, T *output, int M, int N);
template<typename T>
void matrix_r_cyl_bessel_i1f(T * input, T *output, int M, int N);
template<typename T>
void matrix_r_erfcf(T * input, T *output, int M, int N);
template<typename T>
void matrix_r_erfcinvf(T * input, T *output, int M, int N);
template<typename T>
void matrix_r_erfcxf(T * input, T *output, int M, int N);
template<typename T>
void matrix_r_erff(T * input, T *output, int M, int N);
template<typename T>
void matrix_r_erfinvf(T * input, T * outputs, int M, int N);
template<typename T>
void matrix_r_fdimf(T * x, T * y, T *output, int M, int N);
template<typename T>
void matrix_r_fdividef(T * x, T *y, T *output, int M, int N);
template<typename T>
void matrix_r_fmaf(T * x, T *y, T *z, T *output, int M, int N);
template<typename T>
void matrix_r_hypotf(T * x, T *y, T *output, int M, int N);
template<typename T>
void matrix_r_ilogbf(T * input, T *output, int M, int N);
template<typename T>
void matrix_r_j0f(T * input, T *output, int M, int N);
template<typename T>
void matrix_r_j1f(T * input, T *output, int M, int N);
template<typename T>
void matrix_r_jnf(T * input, T *output, int n, int M, int N);
template<typename T>
void matrix_r_ldexpf(T * input, T *output, int exp, int M, int N);
template<typename T>
void matrix_r_lgammaf(T * input, T *output, int M, int N);
template<typename T>
void matrix_r_nearbyintf(T * input, T *output, int M, int N);
template<typename T>
void matrix_r_norm3df(T * x, T *y, T *z, T *output, int M, int N);
template<typename T>
void matrix_r_norm4df(T * x, T *y, T *z, T * w, T *output,int M, int N);
template<typename T>
void matrix_r_normcdff(T * input, T *output, int M, int N);
template<typename T>
void matrix_r_normcdfinvf(T * input, T *output, int M, int N);
template<typename T>
void matrix_r_normf(int dim, T *input, T *output, int M, int N);
template<typename T>
void matrix_r_remainderf(T * x, T *y, T *output, int M, int N);
template<typename T>
void matrix_r_rcbrtf(T * input, T * output, int M, int N);
template<typename T>
void matrix_r_rhypotf(T * x, T *y, T * output, int M, int N);
template<typename T>
void matrix_r_rnorm3df(T * x, T *y, T *z, T * output, int M, int N);
template<typename T>
void matrix_r_rnorm4df(T * x, T *y, T *z, T *w, T *output, int M, int N);
template<typename T>
void matrix_r_rnormf(int dim, T *input, T *output, int M, int N);
template<typename T>
void matrix_r_scalblnf(T * input, long int n,  T *output, int M, int N);
template<typename T>
void matrix_r_sinpif(T * input, T *output, int M, int N);
template<typename T>
void matrix_r_tgammaf(T * input, T *output, int M, int N);
template<typename T>
void matrix_r_y0f(T * input, T *output, int M, int N);
template<typename T>
void matrix_r_y1f(T * input, T *output, int M, int N);
template<typename T>
void matrix_r_ynf(int n, T * input, T *output, int M, int N);

template<typename T>
void matrix_r_softmaxf(T * input, T *output, int M, int N);
template<typename T>
void matrix_r_sigmoidf(T * x, T *output, int M, int N);
template<typename T>
void matrix_r_sigmoid_gradf(T * x, T *output,int M, int N);
template<typename T>
void matrix_r_tanh_gradf(T * x, T *output, int M, int N);
template<typename T>
void matrix_r_reluf(T * x, T *output, int M, int N);
template<typename T>
void matrix_r_relu_gradf(T * x, T *output, int M, int N);



template<typename T>
void matrix_r_addf_const(T * x, T y, T *output, int M, int N);
template<typename T>
void matrix_r_subf_const(T * x, T y, T *output, int M, int N);
template<typename T>
void matrix_r_mulf_const(T * x, T y, T *output,int M, int N);
template<typename T>
void matrix_r_divf_const(T * x, T y, T *output, int M, int N);
template<typename T>
void matrix_r_divf_matrix_const(T * x, T y, T *output, int M, int N);
template<typename T>
void matrix_r_modf_const(T * x, T y, T *output, int M, int N);
template<typename T>
void matrix_r_fmaxf_const(T * x, T y, T *output, int M, int N);
template<typename T>
void matrix_r_fminf_const(T * x, T y, T *output, int M, int N);
template<typename T>
void matrix_r_fmodf_const(T * x, T y, T *output, int M, int N);
template<typename T>
void matrix_r_powf_const(T * x, T y, T *output, int M, int N);
template<typename T>
void matrix_r_fdimf_const(T * x, T y, T *output, int M, int N);
template<typename T>
void matrix_r_fdividef_const(T * x, T y, T *output, int M, int N);
template<typename T>
void matrix_r_hypotf_const(T * x, T y, T *output, int M, int N);
template<typename T>
void matrix_r_remainderf_const(T * x, T y, T * output, int M, int N);
template<typename T>
void matrix_r_rhypotf_const(T * x, T y, T *output, int M, int N);

template<typename T>
void matrix_r_addf_scalar(T * x, T *y, T *output, int M, int N);
template<typename T>
void matrix_r_subf_scalar(T * x, T *y, T *output, int M, int N);
template<typename T>
void matrix_r_mulf_scalar(T * x, T *y, T *output, int M, int N);
template<typename T>
void matrix_r_divf_scalar(T * x, T *y, T *output, int M, int N);
template<typename T>
void matrix_r_modf_scalar(T * x, T *y, T *output, int M, int N);
template<typename T>
void matrix_r_atan2f_scalar(T * a, T *b, T *output,int M, int N);
template<typename T>
void matrix_r_fmaxf_scalar(T * x, T *y, T *output,int M, int N);
template<typename T>
void matrix_r_fminf_scalar(T * x, T *y, T *output,int M, int N);
template<typename T>
void matrix_r_fmodf_scalar(T * x, T *y, T *output,int M, int N);
template<typename T>
void matrix_r_powf_scalar(T * x, T *y, T *output,int M, int N);
template<typename T>
void matrix_r_fdimf_scalar(T * x, T *y, T *output, int M, int N);
template<typename T>
void matrix_r_fdividef_scalar(T * x, T *y, T *output, int M, int N);
template<typename T>
void matrix_r_hypotf_scalar(T * x, T *y, T *output, int M, int N);
template<typename T>
void matrix_r_remainderf_scalar(T * x, T *y, T * output, int M, int N);
template<typename T>
void matrix_r_rhypotf_scalar(T * x, T *y, T *output, int M, int N);


template<typename T>
void matrix_r_copysignf(T * X, T *Y, T *output, int M, int N);

template<typename T>
void matrix_r_truncf(T * x, T *output, int M, int N);

