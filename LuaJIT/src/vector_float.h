#ifndef __VECTORFLOAT_H 
#define __VECTORFLOAT_H 

// these can be C callable.
#ifdef __cplusplus
extern "C" {
#endif 


float* vector_addf(float * x, float * y, int n);
float* vector_subf(float * x, float * y, int n);
float* vector_mulf(float * x, float * y, int n);
float* vector_divf(float * x, float * y, int n);
float* vector_modf(float * x, float * y, int n);

float* vector_acosf(float * devPtr, int n);
float* vector_acoshf(float * devPtr, int n);
float* vector_asinf(float * devPtr, int n);
float* vector_asinhf(float * devPtr, int n);
float* vector_atan2f(float * a, float * b, int n);
float* vector_atanf(float * devPtr, int n);
float* vector_atanhf(float * devPtr, int n);
float* vector_cbrtf(float * devPtr, int n);
float* vector_ceilf(float * devPtr, int n);
float* vector_cosf(float * devPtr, int n);
float* vector_coshf(float * devPtr, int n);
float* vector_exp10f(float * devPtr, int n);
float* vector_exp2f(float * devPtr, int n);
float* vector_expf(float * devPtr, int n);
float* vector_expm1f(float * devPtr, int n);
float* vector_fabsf(float * devPtr, int n);
float* vector_floorf(float * devPtr, int n);
float* vector_fmaxf(float * x, float * y, int n);
float* vector_fminf(float * x, float * y, int n);
float* vector_fmodf(float * x, float * y, int n);
float* vector_hypotf(float * x, float * y, int n);
float* vector_log10f(float * x, int n);
float* vector_log1pf(float * x, int n);
float* vector_log2f(float * x, int n);
float* vector_logbf(float * x, int n);
float* vector_powf(float * x, float * y, int n);
float* vector_rsqrtf(float * x, int n);
float* vector_sinf(float * x, int n);
float* vector_sinhf(float * x, int n);
float* vector_sinpif(float * x, int n);
float* vector_sqrtf(float * x, int n);
float* vector_tanf(float * x, int n);
float* vector_tanhf(float * x, int n);

float* vector_sigmoidf(float * devPtr, int n);
float* vector_sigmoid_gradf(float * devPtr, int n);
float* vector_tanh_gradf(float * devPtr, int n);
float* vector_reluf(float * devPtr, int n);
float* vector_relu_gradf(float * devPtr, int n);
float* vector_softmaxf(float * x, int n);

float  vector_sumf(float * devPtr, int n);
float  vector_prodf(float * devPtr, int n);

float* vector_addf_const(float * x, float  y, int n);
float* vector_subf_const(float * x, float  y, int n);
float* vector_mulf_const(float * x, float  y, int n);
float* vector_divf_const(float * x, float  y, int n);
float* vector_modf_const(float * x, float  y, int n);
float* vector_atan2f_const(float * a, float b, int n);
float* vector_fmaxf_const(float * x, float  y, int n);
float* vector_fminf_const(float * x, float  y, int n);
float* vector_fmodf_const(float * x, float  y, int n);
float* vector_powf_const(float * x, float y, int n);
float* vector_fdimf_const(float * a, float  b, int n);
float* vector_fdividef_const(float * a, float  b, int n);
float* vector_remainderf_const(float * x, float y, int n);
float* vector_hypotf_const(float * x, float  y, int n);
float* vector_rhypotf_const(float * x, float y, int n);


float* vector_addf_scalar(float * x, float * y, int n);
float* vector_subf_scalar(float * x, float * y, int n);
float* vector_mulf_scalar(float * x, float * y, int n);
float* vector_divf_scalar(float * x, float * y, int n);
float* vector_modf_scalar(float * x, float * y, int n);
float* vector_atan2f_scalar(float * a, float *b, int n);
float* vector_fmaxf_scalar(float * x, float  *y, int n);
float* vector_fminf_scalar(float * x, float  *y, int n);
float* vector_fmodf_scalar(float * x, float  *y, int n);
float* vector_powf_scalar(float * x, float *y, int n);
float* vector_fdimf_scalar(float * a, float  *b, int n);
float* vector_fdividef_scalar(float * a, float *b, int n);
float* vector_hypotf_scalar(float * x, float  *y, int n);
float* vector_remainderf_scalar(float * x, float *y, int n);
float* vector_rhypotf_scalar(float * x, float *y, int n);

float* vector_copysignf(float * X, float *Y, int n);
float* vector_cospif(float * devPtr, int n);
float* vector_cyl_bessel_i0f(float * devPtr, int n);
float* vector_cyl_bessel_i1f(float * devPtr, int n);
float* vector_erfcf(float * devPtr, int n);
float* vector_erfcinvf(float * devPtr, int n);
float* vector_erfcxf(float * devPtr, int n);
float* vector_erff(float * devPtr, int n);
float* vector_erfinvf(float * devPtr, int n);
float* vector_fdimf(float * a, float * b, int n);
float* vector_fdividef(float * a, float * b, int n);
float* vector_fmaf(float * x, float * y, float * z, int n);
float* vector_ilogbf(float * x, int n);
float* vector_j0f(float * x, int n);
float* vector_j1f(float * x, int n);
float* vector_jnf(float * x, int N, int n);
float* vector_ldexpf(float * x, int exp, int n);
float* vector_lgammaf(float * x, int n);
long long* vector_llrintf(float * x, int n);
long long* vector_llroundf(float * x, int n);
long* vector_lrintf(float * x, int n);
long* vector_lroundf(float * x, int n);
float* vector_nearbyintf(float * x, int n);
float* vector_norm3df(float * x, float * y, float * z, int n);
float* vector_norm4df(float * x, float * y, float * z, float * q, int n);
float* vector_normcdff(float * x, int n);
float* vector_normcdfinvf(float * x, int n);
float* vector_normf(int dim, float * x, int n);
float* vector_rcbrtf(float * x, int n);
float* vector_remainderf(float * x, float * y, int n);
float* vector_rhypotf(float * x, float * y, int n);
float* vector_rnorm3df(float * x, float * y, float * z, int n);
float* vector_rnorm4df(float * x, float * y, float * z, float * q, int n);
float* vector_rnormf(int dim, float * x, int n);
float* vector_tgammaf(float * x, int n);
float* vector_y0f(float * x, int n);
float* vector_y1f(float * x, int n);
float* vector_ynf(int N, float * x, int n);
float* vector_scalblnf(float * x, long int M, int n);

float* vector_truncf(float * x, int n);
void vector_r_truncf(float * x, float *output, int n);
void vector_r_copysignf(float * X, float *Y, float *output, int n);



void vector_addf_row(float * x, int row, float * y, int src,size_t n);
void vector_subf_row(float * x, int row, float * y, int src,size_t n);
void vector_mulf_row(float * x, int row, float * y, int src,size_t n);
void vector_divf_row(float * x, int row, float * y, int src,size_t n);
void vector_modf_row(float * x, int row, float * y, int src,size_t n);


void vector_setrowf(float * dst, int row, float * sc, int row_src, size_t n);
void vector_copyrowf(float * dst, int row, float * src, int row_src, int n);
    
void vector_r_addf(float * x, float * y, float * output, int n);
void vector_r_subf(float * x, float * y, float * output, int n);
void vector_r_mulf(float * x, float * y, float * output, int n);
void vector_r_divf(float * x, float * y, float * output, int n);
void vector_r_modf(float * x, float * y, float * output, int n);
void vector_r_acosf(float * devPtr, float * outputs, int n);
void vector_r_asinf(float * devPtr, float * outputs, int n);
void vector_r_atanf(float * devPtr, float * outputs, int n);
void vector_r_atan2f(float * a, float * b, float * output, int n);
void vector_r_acoshf(float * devPtr, float * outputs, int n);
void vector_r_asinhf(float * devPtr, float * outputs, int n);
void vector_r_atanhf(float * devPtr, float * outputs, int n);
void vector_r_cosf(float * devPtr, float * outputs, int n);
void vector_r_sinf(float * devPtr, float * outputs, int n);
void vector_r_tanf(float * devPtr, float * outputs, int n);
void vector_r_coshf(float * devPtr, float * outputs, int n);
void vector_r_sinhf(float * devPtr, float * outputs, int n);
void vector_r_tanhf(float * devPtr, float * outputs, int n);
void vector_r_ceilf(float * devPtr, float * output, int n);
void vector_r_exp10f(float * devPtr, float * outputs, int n);
void vector_r_exp2f(float * devPtr, float * output, int n);
void vector_r_expf(float * devPtr, float * output, int n);
void vector_r_expm1f(float * devPtr, float * output, int n);
void vector_r_fabsf(float * devPtr, float * output, int n);
void vector_r_floorf(float * devPtr, float * output, int n);
void vector_r_fmaxf(float * x, float * y, float * output, int n);
void vector_r_fminf(float * x, float * y, float * output, int n);
void vector_r_fmodf(float * x, float * y, float * output, int n);
void vector_r_log10f(float * x, float *output, int n);
void vector_r_log1pf(float * x, float *output, int n);
void vector_r_log2f(float * x, float * output, int n);
void vector_r_logbf(float * x, float *output, int n);
void vector_r_powf(float * x, float * y, float *output, int n);
void vector_r_rsqrtf(float * x, float * output, int n);
void vector_r_sinf(float * x, float *output, int n);
void vector_r_sinhf(float * x, float *output, int n);
void vector_r_sqrtf(float * x, float *output, int n);
void vector_r_tanf(float * x, float *output, int n);
void vector_r_tanhf(float * x, float *output, int n);
void vector_r_softmaxf(float * x, float *output, int n);
void vector_r_sigmoidf(float * x, float *output, int n);
void vector_r_sigmoid_gradf(float * x, float *output, int n);
void vector_r_tanh_gradf(float * x, float *output, int n);
void vector_r_reluf(float * x, float *output, int n);
void vector_r_relu_gradf(float * x, float *output, int n);
void vector_r_cbrtf(float * devPtr, float * output, int n);
void vector_r_cospif(float * devPtr, float * output, int n);
void vector_r_cyl_bessel_i0f(float * devPtr, float * output, int n);
void vector_r_cyl_bessel_i1f(float * devPtr, float * output, int n);
void vector_r_erfcf(float * devPtr, float * output, int n);
void vector_r_erfcinvf(float * devPtr, float * output, int n);
void vector_r_erfcxf(float * devPtr, float * output, int n);
void vector_r_erff(float * devPtr, float * output, int n);
void vector_r_erfinvf(float * devPtr, float * output, int n);
void vector_r_fdimf(float * a, float * b, float * output, int n);
void vector_r_fdividef(float * a, float * b, float * output, int n);
void vector_r_fmaf(float * x, float * y, float * z, float *output, int n);
void vector_r_hypotf(float * x, float * y, float * output, int n);
void vector_r_ilogbf(float * x, float *output, int n);
void vector_r_j0f(float * x, float *output, int n);
void vector_r_j1f(float * x, float *output, int n);
void vector_r_jnf(float * x, float * output, int M, int n);
void vector_r_ldexpf(float * x, float * output, int exp, int n);
void vector_r_lgammaf(float * x, float *output, int n);
void vector_r_nearbyintf(float * x, float *output, int n);
void vector_r_norm3df(float * x, float * y, float * z, float * output, int n);
void vector_r_norm4df(float * x, float * y, float * z, float * q, float * output, int n);
void vector_r_normcdff(float * x, float * output, int n);
void vector_r_normcdfinvf(float * x, float *output, int n);
void vector_r_normf(int dim, float * x, float * output, int n);
void vector_r_rcbrtf(float * x, float *output, int n);
void vector_r_remainderf(float * x, float * y, float *output, int n);
void vector_r_rhypotf(float * x, float * y, float *output, int n);
void vector_r_rnorm3df(float * x, float * y, float * z, float * output, int n);
void vector_r_rnorm4df(float * x, float * y, float * z, float * q, float *output, int n);
void vector_r_rnormf(int dim, float * x, float *output, int n);
void vector_r_scalblnf(float * x, long int M, float * output, int n);
void vector_r_tgammaf(float * x, float * output, int n);
void vector_r_truncf(float * x, float *output, int n);
void vector_r_y0f(float * x, float *output, int n);
void vector_r_y1f(float * x, float * output, int n);
void vector_r_ynf(int M, float * x, float *output, int n);
void vector_r_sinpif(float * x, float *output, int n);


void vector_r_addf_const(float * x, float  y, float *output, int n);
void vector_r_subf_const(float * x, float  y, float *output,int n);
void vector_r_mulf_const(float * x, float  y, float *output,int n);
void vector_r_divf_const(float * x, float  y, float *output,int n);
void vector_r_modf_const(float * x, float  y, float *output,int n);
void vector_r_atan2f_const(float * a, float b, float *output,int n);
void vector_r_fmaxf_const(float * x, float  y, float *output,int n);
void vector_r_fminf_const(float * x, float  y, float *output,int n);
void vector_r_fmodf_const(float * x, float  y, float *output,int n);
void vector_r_powf_const(float * x, float y, float *output,int n);


void vector_r_addf_scalar(float * x, float * y, float *output,int n);
void vector_r_subf_scalar(float * x, float * y, float *output,int n);
void vector_r_mulf_scalar(float * x, float * y, float *output,int n);
void vector_r_divf_scalar(float * x, float * y, float *output,int n);
void vector_r_modf_scalar(float * x, float * y, float *output,int n);
void vector_r_atan2f_scalar(float * a, float *b, float *output,int n);
void vector_r_fmaxf_scalar(float * x, float  *y, float *output,int n);
void vector_r_fminf_scalar(float * x, float  *y, float *output,int n);
void vector_r_fmodf_scalar(float * x, float  *y, float *output,int n);
void vector_r_powf_scalar(float * x, float *y, float *output,int n);

void vector_r_fdimf_const(float * a, float  b, float *output,int n);
void vector_r_fdividef_const(float * a, float  b, float *output,int n);
void vector_r_hypotf_const(float * x, float  y, float *output,int n);
void vector_r_remainderf_const(float * x, float y, float *output,int n);
void vector_r_rhypotf_const(float * x, float y, float *output,int n);

void vector_r_fdimf_scalar(float * a, float  *b, float *output,int n);
void vector_r_fdividef_scalar(float * a, float *b, float *output,int n);
void vector_r_hypotf_scalar(float * x, float  *y, float *output,int n);
void vector_r_remainderf_scalar(float * x, float *y, float *output,int n);
void vector_r_rhypotf_scalar(float * x, float *y, float *output,int n);


void    return_memory(int length, float *fp);
float*  find_memory(int length);
void    add_memory(int length, float * ptr);
void    clear_cache();

void    calcSize(int N,int * gridSize, int * blockSize);

void cuda_zero(float * dst, int n);
void cuda_memcpy(float * dst, float * src, int n);


/*
// the kernel must wrap the call to the GPU kernel.
typedef float* (*vector_kernel1)(float * input, int n);
typedef float* (*vector_kernel2)(float * x, float * y, int n);
typedef float* (*vector_kernel3)(float * x, float * y, float * z, int n);
typedef float* (*vector_kernel4)(float * x, float * y, float * z, float * w, int n);

void    register_vector_kernel1(const char* name, vector_kernel1 kernel);
void    register_vector_kernel2(const char* name, vector_kernel2 kernel);
void    register_vector_kernel3(const char* name, vector_kernel3 kernel);
void    register_vector_kernel4(const char* name, vector_kernel4 kernel);

float*  execute_vector_kernel1(const char* name, float * input, int n);
float*  execute_vector_kernel2(const char* name, float * i, float * j, int n);
float*  execute_vector_kernel3(const char* name, float * i, float * j, float *k, int n);
float*  execute_vector_kernel4(const char* name, float * i, float * j, float *k, float *w, int n);
*/


#ifdef __cplusplus 
}
#endif

#endif
