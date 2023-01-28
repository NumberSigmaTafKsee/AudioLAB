#pragma once
// these can be C callable.
#ifdef __cplusplus
extern "C" {
#endif 


double* vector_addf(double * x, double * y, int n);
double* vector_subf(double * x, double * y, int n);
double* vector_mulf(double * x, double * y, int n);
double* vector_divf(double * x, double * y, int n);
double* vector_modf(double * x, double * y, int n);

double* vector_acosf(double * devPtr, int n);
double* vector_acoshf(double * devPtr, int n);
double* vector_asinf(double * devPtr, int n);
double* vector_asinhf(double * devPtr, int n);
double* vector_atan2f(double * a, double * b, int n);
double* vector_atanf(double * devPtr, int n);
double* vector_atanhf(double * devPtr, int n);
double* vector_cbrtf(double * devPtr, int n);
double* vector_ceilf(double * devPtr, int n);
double* vector_cosf(double * devPtr, int n);
double* vector_coshf(double * devPtr, int n);
double* vector_exp10f(double * devPtr, int n);
double* vector_exp2f(double * devPtr, int n);
double* vector_expf(double * devPtr, int n);
double* vector_expm1f(double * devPtr, int n);
double* vector_fabsf(double * devPtr, int n);
double* vector_floorf(double * devPtr, int n);
double* vector_fmaxf(double * x, double * y, int n);
double* vector_fminf(double * x, double * y, int n);
double* vector_fmodf(double * x, double * y, int n);
double* vector_hypotf(double * x, double * y, int n);
double* vector_log10f(double * x, int n);
double* vector_log1pf(double * x, int n);
double* vector_log2f(double * x, int n);
double* vector_logbf(double * x, int n);
double* vector_powf(double * x, double * y, int n);
double* vector_rsqrtf(double * x, int n);
double* vector_sinf(double * x, int n);
double* vector_sinhf(double * x, int n);
double* vector_sinpif(double * x, int n);
double* vector_sqrtf(double * x, int n);
double* vector_tanf(double * x, int n);
double* vector_tanhf(double * x, int n);

double* vector_sigmoidf(double * devPtr, int n);
double* vector_sigmoid_gradf(double * devPtr, int n);
double* vector_tanh_gradf(double * devPtr, int n);
double* vector_reluf(double * devPtr, int n);
double* vector_relu_gradf(double * devPtr, int n);
double* vector_softmaxf(double * x, int n);

double  vector_sumf(double * devPtr, int n);
double  vector_prodf(double * devPtr, int n);

double* vector_addf_const(double * x, double  y, int n);
double* vector_subf_const(double * x, double  y, int n);
double* vector_mulf_const(double * x, double  y, int n);
double* vector_divf_const(double * x, double  y, int n);
double* vector_modf_const(double * x, double  y, int n);
double* vector_atan2f_const(double * a, double b, int n);
double* vector_fmaxf_const(double * x, double  y, int n);
double* vector_fminf_const(double * x, double  y, int n);
double* vector_fmodf_const(double * x, double  y, int n);
double* vector_powf_const(double * x, double y, int n);
double* vector_fdimf_const(double * a, double  b, int n);
double* vector_fdividef_const(double * a, double  b, int n);
double* vector_remainderf_const(double * x, double y, int n);
double* vector_hypotf_const(double * x, double  y, int n);
double* vector_rhypotf_const(double * x, double y, int n);


double* vector_addf_scalar(double * x, double * y, int n);
double* vector_subf_scalar(double * x, double * y, int n);
double* vector_mulf_scalar(double * x, double * y, int n);
double* vector_divf_scalar(double * x, double * y, int n);
double* vector_modf_scalar(double * x, double * y, int n);
double* vector_atan2f_scalar(double * a, double *b, int n);
double* vector_fmaxf_scalar(double * x, double  *y, int n);
double* vector_fminf_scalar(double * x, double  *y, int n);
double* vector_fmodf_scalar(double * x, double  *y, int n);
double* vector_powf_scalar(double * x, double *y, int n);
double* vector_fdimf_scalar(double * a, double  *b, int n);
double* vector_fdividef_scalar(double * a, double *b, int n);
double* vector_hypotf_scalar(double * x, double  *y, int n);
double* vector_remainderf_scalar(double * x, double *y, int n);
double* vector_rhypotf_scalar(double * x, double *y, int n);

double* vector_copysignf(double * X, double *Y, int n);
double* vector_cospif(double * devPtr, int n);
double* vector_cyl_bessel_i0f(double * devPtr, int n);
double* vector_cyl_bessel_i1f(double * devPtr, int n);
double* vector_erfcf(double * devPtr, int n);
double* vector_erfcinvf(double * devPtr, int n);
double* vector_erfcxf(double * devPtr, int n);
double* vector_erff(double * devPtr, int n);
double* vector_erfinvf(double * devPtr, int n);
double* vector_fdimf(double * a, double * b, int n);
double* vector_fdividef(double * a, double * b, int n);
double* vector_fmaf(double * x, double * y, double * z, int n);
double* vector_ilogbf(double * x, int n);
double* vector_j0f(double * x, int n);
double* vector_j1f(double * x, int n);
double* vector_jnf(double * x, int N, int n);
double* vector_ldexpf(double * x, int exp, int n);
double* vector_lgammaf(double * x, int n);
long long* vector_llrintf(double * x, int n);
long long* vector_llroundf(double * x, int n);
long* vector_lrintf(double * x, int n);
long* vector_lroundf(double * x, int n);
double* vector_nearbyintf(double * x, int n);
double* vector_norm3df(double * x, double * y, double * z, int n);
double* vector_norm4df(double * x, double * y, double * z, double * q, int n);
double* vector_normcdff(double * x, int n);
double* vector_normcdfinvf(double * x, int n);
double* vector_normf(int dim, double * x, int n);
double* vector_rcbrtf(double * x, int n);
double* vector_remainderf(double * x, double * y, int n);
double* vector_rhypotf(double * x, double * y, int n);
double* vector_rnorm3df(double * x, double * y, double * z, int n);
double* vector_rnorm4df(double * x, double * y, double * z, double * q, int n);
double* vector_rnormf(int dim, double * x, int n);
double* vector_tgammaf(double * x, int n);
double* vector_y0f(double * x, int n);
double* vector_y1f(double * x, int n);
double* vector_ynf(int N, double * x, int n);
double* vector_scalblnf(double * x, long int M, int n);

double* vector_truncf(double * x, int n);
void vector_r_truncf(double * x, double *output, int n);
void vector_r_copysignf(double * X, double *Y, double *output, int n);



void vector_addf_row(double * x, int row, double * y, int src,size_t n);
void vector_subf_row(double * x, int row, double * y, int src,size_t n);
void vector_mulf_row(double * x, int row, double * y, int src,size_t n);
void vector_divf_row(double * x, int row, double * y, int src,size_t n);
void vector_modf_row(double * x, int row, double * y, int src,size_t n);


void vector_setrowf(double * dst, int row, double * sc, int row_src, size_t n);
void vector_copyrowf(double * dst, int row, double * src, int row_src, int n);
    
void vector_r_addf(double * x, double * y, double * output, int n);
void vector_r_subf(double * x, double * y, double * output, int n);
void vector_r_mulf(double * x, double * y, double * output, int n);
void vector_r_divf(double * x, double * y, double * output, int n);
void vector_r_modf(double * x, double * y, double * output, int n);
void vector_r_acosf(double * devPtr, double * outputs, int n);
void vector_r_asinf(double * devPtr, double * outputs, int n);
void vector_r_atanf(double * devPtr, double * outputs, int n);
void vector_r_atan2f(double * a, double * b, double * output, int n);
void vector_r_acoshf(double * devPtr, double * outputs, int n);
void vector_r_asinhf(double * devPtr, double * outputs, int n);
void vector_r_atanhf(double * devPtr, double * outputs, int n);
void vector_r_cosf(double * devPtr, double * outputs, int n);
void vector_r_sinf(double * devPtr, double * outputs, int n);
void vector_r_tanf(double * devPtr, double * outputs, int n);
void vector_r_coshf(double * devPtr, double * outputs, int n);
void vector_r_sinhf(double * devPtr, double * outputs, int n);
void vector_r_tanhf(double * devPtr, double * outputs, int n);
void vector_r_ceilf(double * devPtr, double * output, int n);
void vector_r_exp10f(double * devPtr, double * outputs, int n);
void vector_r_exp2f(double * devPtr, double * output, int n);
void vector_r_expf(double * devPtr, double * output, int n);
void vector_r_expm1f(double * devPtr, double * output, int n);
void vector_r_fabsf(double * devPtr, double * output, int n);
void vector_r_floorf(double * devPtr, double * output, int n);
void vector_r_fmaxf(double * x, double * y, double * output, int n);
void vector_r_fminf(double * x, double * y, double * output, int n);
void vector_r_fmodf(double * x, double * y, double * output, int n);
void vector_r_log10f(double * x, double *output, int n);
void vector_r_log1pf(double * x, double *output, int n);
void vector_r_log2f(double * x, double * output, int n);
void vector_r_logbf(double * x, double *output, int n);
void vector_r_powf(double * x, double * y, double *output, int n);
void vector_r_rsqrtf(double * x, double * output, int n);
void vector_r_sinf(double * x, double *output, int n);
void vector_r_sinhf(double * x, double *output, int n);
void vector_r_sqrtf(double * x, double *output, int n);
void vector_r_tanf(double * x, double *output, int n);
void vector_r_tanhf(double * x, double *output, int n);
void vector_r_softmaxf(double * x, double *output, int n);
void vector_r_sigmoidf(double * x, double *output, int n);
void vector_r_sigmoid_gradf(double * x, double *output, int n);
void vector_r_tanh_gradf(double * x, double *output, int n);
void vector_r_reluf(double * x, double *output, int n);
void vector_r_relu_gradf(double * x, double *output, int n);
void vector_r_cbrtf(double * devPtr, double * output, int n);
void vector_r_cospif(double * devPtr, double * output, int n);
void vector_r_cyl_bessel_i0f(double * devPtr, double * output, int n);
void vector_r_cyl_bessel_i1f(double * devPtr, double * output, int n);
void vector_r_erfcf(double * devPtr, double * output, int n);
void vector_r_erfcinvf(double * devPtr, double * output, int n);
void vector_r_erfcxf(double * devPtr, double * output, int n);
void vector_r_erff(double * devPtr, double * output, int n);
void vector_r_erfinvf(double * devPtr, double * output, int n);
void vector_r_fdimf(double * a, double * b, double * output, int n);
void vector_r_fdividef(double * a, double * b, double * output, int n);
void vector_r_fmaf(double * x, double * y, double * z, double *output, int n);
void vector_r_hypotf(double * x, double * y, double * output, int n);
void vector_r_ilogbf(double * x, double *output, int n);
void vector_r_j0f(double * x, double *output, int n);
void vector_r_j1f(double * x, double *output, int n);
void vector_r_jnf(double * x, double * output, int M, int n);
void vector_r_ldexpf(double * x, double * output, int exp, int n);
void vector_r_lgammaf(double * x, double *output, int n);
void vector_r_nearbyintf(double * x, double *output, int n);
void vector_r_norm3df(double * x, double * y, double * z, double * output, int n);
void vector_r_norm4df(double * x, double * y, double * z, double * q, double * output, int n);
void vector_r_normcdff(double * x, double * output, int n);
void vector_r_normcdfinvf(double * x, double *output, int n);
void vector_r_normf(int dim, double * x, double * output, int n);
void vector_r_rcbrtf(double * x, double *output, int n);
void vector_r_remainderf(double * x, double * y, double *output, int n);
void vector_r_rhypotf(double * x, double * y, double *output, int n);
void vector_r_rnorm3df(double * x, double * y, double * z, double * output, int n);
void vector_r_rnorm4df(double * x, double * y, double * z, double * q, double *output, int n);
void vector_r_rnormf(int dim, double * x, double *output, int n);
void vector_r_scalblnf(double * x, long int M, double * output, int n);
void vector_r_tgammaf(double * x, double * output, int n);
void vector_r_truncf(double * x, double *output, int n);
void vector_r_y0f(double * x, double *output, int n);
void vector_r_y1f(double * x, double * output, int n);
void vector_r_ynf(int M, double * x, double *output, int n);
void vector_r_sinpif(double * x, double *output, int n);


void vector_r_addf_const(double * x, double  y, double *output, int n);
void vector_r_subf_const(double * x, double  y, double *output,int n);
void vector_r_mulf_const(double * x, double  y, double *output,int n);
void vector_r_divf_const(double * x, double  y, double *output,int n);
void vector_r_modf_const(double * x, double  y, double *output,int n);
void vector_r_atan2f_const(double * a, double b, double *output,int n);
void vector_r_fmaxf_const(double * x, double  y, double *output,int n);
void vector_r_fminf_const(double * x, double  y, double *output,int n);
void vector_r_fmodf_const(double * x, double  y, double *output,int n);
void vector_r_powf_const(double * x, double y, double *output,int n);


void vector_r_addf_scalar(double * x, double * y, double *output,int n);
void vector_r_subf_scalar(double * x, double * y, double *output,int n);
void vector_r_mulf_scalar(double * x, double * y, double *output,int n);
void vector_r_divf_scalar(double * x, double * y, double *output,int n);
void vector_r_modf_scalar(double * x, double * y, double *output,int n);
void vector_r_atan2f_scalar(double * a, double *b, double *output,int n);
void vector_r_fmaxf_scalar(double * x, double  *y, double *output,int n);
void vector_r_fminf_scalar(double * x, double  *y, double *output,int n);
void vector_r_fmodf_scalar(double * x, double  *y, double *output,int n);
void vector_r_powf_scalar(double * x, double *y, double *output,int n);

void vector_r_fdimf_const(double * a, double  b, double *output,int n);
void vector_r_fdividef_const(double * a, double  b, double *output,int n);
void vector_r_hypotf_const(double * x, double  y, double *output,int n);
void vector_r_remainderf_const(double * x, double y, double *output,int n);
void vector_r_rhypotf_const(double * x, double y, double *output,int n);

void vector_r_fdimf_scalar(double * a, double  *b, double *output,int n);
void vector_r_fdividef_scalar(double * a, double *b, double *output,int n);
void vector_r_hypotf_scalar(double * x, double  *y, double *output,int n);
void vector_r_remainderf_scalar(double * x, double *y, double *output,int n);
void vector_r_rhypotf_scalar(double * x, double *y, double *output,int n);


void    return_memory(int length, double *fp);
double*  find_memory(int length);
void    add_memory(int length, double * ptr);
void    clear_cache();

void    calcSize(int N,int * gridSize, int * blockSize);

void cuda_zero(double * dst, int n);
void cuda_memcpy(double * dst, double * src, int n);


/*
// the kernel must wrap the call to the GPU kernel.
typedef double* (*vector_kernel1)(double * input, int n);
typedef double* (*vector_kernel2)(double * x, double * y, int n);
typedef double* (*vector_kernel3)(double * x, double * y, double * z, int n);
typedef double* (*vector_kernel4)(double * x, double * y, double * z, double * w, int n);

void    register_vector_kernel1(const char* name, vector_kernel1 kernel);
void    register_vector_kernel2(const char* name, vector_kernel2 kernel);
void    register_vector_kernel3(const char* name, vector_kernel3 kernel);
void    register_vector_kernel4(const char* name, vector_kernel4 kernel);

double*  execute_vector_kernel1(const char* name, double * input, int n);
double*  execute_vector_kernel2(const char* name, double * i, double * j, int n);
double*  execute_vector_kernel3(const char* name, double * i, double * j, double *k, int n);
double*  execute_vector_kernel4(const char* name, double * i, double * j, double *k, double *w, int n);
*/


#ifdef __cplusplus 
}
#endif
