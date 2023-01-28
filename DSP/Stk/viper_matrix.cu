//////////////////////////////////////////////////////////////////////////
// matrix
//////////////////////////////////////////////////////////////////////////
#include <cassert>
#include "cuda_runtime.h"
#include "math_constants.h"
#include "viper_matrix.hpp"
#include "viper_vector.hpp"
//////////////////////////////////////////////////////////////////////////
#include "viper_ops.h"

static const int BLOCK_SIZE =16;
cudaStream_t get_cuda_stream();

template void _2d_r_mulf(float * a, float * b, float *output, int M, int N);
template void _2d_r_addf(float * a, float * b, float *output,int M, int N);
template void _2d_r_subf(float * a, float * b, float *output,int M, int N);
template void _2d_r_divf(float * a, float * b, float *output,int M, int N);
template void _2d_r_modf(float * a, float * b, float *output,int M, int N);

template void matrix_r_addf(float * a, float * b, float *output, int M, int N, int K);
template void matrix_r_subf(float * a, float * b, float *output, int M, int N, int K);
template void matrix_r_hadamardf(float * a, float * b, float *output, int M, int N, int K);
template void matrix_r_multiplyf(float * a, float * b, float * output, int M, int N, int K);
template void matrix_r_transposef(float * input, float *output, int M, int N);
template void matrix_r_acosf(float * input, float * output,  int M, int N);
template void matrix_r_asinf(float * input, float * output,  int M, int N);
template void matrix_r_atanf(float * input, float * output,  int M, int N);
template void matrix_r_atan2f(float * x, float * y, float * output, int M, int N);
template void matrix_r_acoshf(float * input, float * output,  int M, int N);
template void matrix_r_asinhf(float * input, float * output,  int M, int N);
template void matrix_r_atanhf(float * input, float * output,  int M, int N);
template void matrix_r_cosf(float * input, float * output,  int M, int N);
template void matrix_r_sinf(float * input, float * output,  int M, int N);
template void matrix_r_tanf(float * input, float * output, int M, int N);
template void matrix_r_coshf(float * input, float * output,  int M, int N);
template void matrix_r_sinhf(float * input, float * output,  int M, int N);
template void matrix_r_tanhf(float * input, float * output, int M, int N);
template void matrix_r_atan2f_const(float * input, float b, float * output, int M, int N);
template void matrix_r_ceilf(float * input, float *output, int M, int N);
template void matrix_r_exp10f(float * input, float *output, int M, int N);
template void matrix_r_exp2f(float * input, float *output, int M, int N);
template void matrix_r_expf(float * input, float *output, int M, int N);
template void matrix_r_expm1f(float * input, float *output, int M, int N);
template void matrix_r_fabsf(float * input, float *output, int M, int N);
template void matrix_r_floorf(float * input, float *output, int M, int N);
template void matrix_r_fmaxf(float * x, float *y, float *output, int M, int N);
template void matrix_r_fminf(float * x, float *y, float *output, int M, int N);
template void matrix_r_fmodf(float * x, float *y, float *output, int M, int N);
template void matrix_r_log10f(float * input, float *output, int M, int N);
template void matrix_r_log1pf(float * input, float *output, int M, int N);
template void matrix_r_log2f(float * input, float *output, int M, int N);
template void matrix_r_logbf(float * input, float *output, int M, int N);
template void matrix_r_powf(float * x, float *y, float *output, int M, int N);
template void matrix_r_rsqrtf(float * input, float *output, int M, int N);
template void matrix_r_sqrtf(float * input, float *output, int M, int N);
template void matrix_r_cbrtf(float * input, float *output, int M, int N);
template void matrix_r_cospif(float * input, float *output, int M, int N);
template void matrix_r_cyl_bessel_i0f(float * input, float *output, int M, int N);
template void matrix_r_cyl_bessel_i1f(float * input, float *output, int M, int N);
template void matrix_r_erfcf(float * input, float *output, int M, int N);
template void matrix_r_erfcinvf(float * input, float *output, int M, int N);
template void matrix_r_erfcxf(float * input, float *output, int M, int N);
template void matrix_r_erff(float * input, float *output, int M, int N);
template void matrix_r_erfinvf(float * input, float * outputs, int M, int N);
template void matrix_r_fdimf(float * x, float * y, float *output, int M, int N);
template void matrix_r_fdividef(float * x, float *y, float *output, int M, int N);
template void matrix_r_fmaf(float * x, float *y, float *z, float *output, int M, int N);
template void matrix_r_hypotf(float * x, float *y, float *output, int M, int N);
template void matrix_r_ilogbf(float * input, float *output, int M, int N);
template void matrix_r_j0f(float * input, float *output, int M, int N);
template void matrix_r_j1f(float * input, float *output, int M, int N);
template void matrix_r_jnf(float * input, float *output, int n, int M, int N);
template void matrix_r_ldexpf(float * input, float *output, int exp, int M, int N);
template void matrix_r_lgammaf(float * input, float *output, int M, int N);
template void matrix_r_nearbyintf(float * input, float *output, int M, int N);
template void matrix_r_norm3df(float * x, float *y, float *z, float *output, int M, int N);
template void matrix_r_norm4df(float * x, float *y, float *z, float * w, float *output,int M, int N);
template void matrix_r_normcdff(float * input, float *output, int M, int N);
template void matrix_r_normcdfinvf(float * input, float *output, int M, int N);
template void matrix_r_normf(int dim, float *input, float *output, int M, int N);
template void matrix_r_remainderf(float * x, float *y, float *output, int M, int N);
template void matrix_r_rcbrtf(float * input, float * output, int M, int N);
template void matrix_r_rhypotf(float * x, float *y, float * output, int M, int N);
template void matrix_r_rnorm3df(float * x, float *y, float *z, float * output, int M, int N);
template void matrix_r_rnorm4df(float * x, float *y, float *z, float *w, float *output, int M, int N);
template void matrix_r_rnormf(int dim, float *input, float *output, int M, int N);
template void matrix_r_scalblnf(float * input, long int n,  float *output, int M, int N);
template void matrix_r_sinpif(float * input, float *output, int M, int N);
template void matrix_r_tgammaf(float * input, float *output, int M, int N);
template void matrix_r_y0f(float * input, float *output, int M, int N);
template void matrix_r_y1f(float * input, float *output, int M, int N);
template void matrix_r_ynf(int n, float * input, float *output, int M, int N);
template void matrix_r_softmaxf(float * input, float *output, int M, int N);
template void matrix_r_sigmoidf(float * x, float *output, int M, int N);
template void matrix_r_sigmoid_gradf(float * x, float *output,int M, int N);
template void matrix_r_tanh_gradf(float * x, float *output, int M, int N);
template void matrix_r_reluf(float * x, float *output, int M, int N);
template void matrix_r_relu_gradf(float * x, float *output, int M, int N);
template void matrix_r_addf_const(float * x, float y, float *output, int M, int N);
template void matrix_r_subf_const(float * x, float y, float *output, int M, int N);
template void matrix_r_mulf_const(float * x, float y, float *output,int M, int N);
template void matrix_r_divf_const(float * x, float y, float *output, int M, int N);
template void matrix_r_divf_matrix_const(float * x, float y, float *output, int M, int N);
template void matrix_r_modf_const(float * x, float y, float *output, int M, int N);
template void matrix_r_fmaxf_const(float * x, float y, float *output, int M, int N);
template void matrix_r_fminf_const(float * x, float y, float *output, int M, int N);
template void matrix_r_fmodf_const(float * x, float y, float *output, int M, int N);
template void matrix_r_powf_const(float * x, float y, float *output, int M, int N);
template void matrix_r_fdimf_const(float * x, float y, float *output, int M, int N);
template void matrix_r_fdividef_const(float * x, float y, float *output, int M, int N);
template void matrix_r_hypotf_const(float * x, float y, float *output, int M, int N);
template void matrix_r_remainderf_const(float * x, float y, float * output, int M, int N);
template void matrix_r_rhypotf_const(float * x, float y, float *output, int M, int N);
template void matrix_r_addf_scalar(float * x, float *y, float *output, int M, int N);
template void matrix_r_subf_scalar(float * x, float *y, float *output, int M, int N);
template void matrix_r_mulf_scalar(float * x, float *y, float *output, int M, int N);
template void matrix_r_divf_scalar(float * x, float *y, float *output, int M, int N);
template void matrix_r_modf_scalar(float * x, float *y, float *output, int M, int N);
template void matrix_r_atan2f_scalar(float * a, float *b, float *output,int M, int N);
template void matrix_r_fmaxf_scalar(float * x, float *y, float *output,int M, int N);
template void matrix_r_fminf_scalar(float * x, float *y, float *output,int M, int N);
template void matrix_r_fmodf_scalar(float * x, float *y, float *output,int M, int N);
template void matrix_r_powf_scalar(float * x, float *y, float *output,int M, int N);
template void matrix_r_fdimf_scalar(float * x, float *y, float *output, int M, int N);
template void matrix_r_fdividef_scalar(float * x, float *y, float *output, int M, int N);
template void matrix_r_hypotf_scalar(float * x, float *y, float *output, int M, int N);
template void matrix_r_remainderf_scalar(float * x, float *y, float * output, int M, int N);
template void matrix_r_rhypotf_scalar(float * x, float *y, float *output, int M, int N);
template void matrix_r_copysignf(float * X, float *Y, float *output, int M, int N);
template void matrix_r_truncf(float * x, float *output, int M, int N);

template void _2d_r_mulf(double * a, double * b, double *output, int M, int N);
template void _2d_r_addf(double * a, double * b, double *output,int M, int N);
template void _2d_r_subf(double * a, double * b, double *output,int M, int N);
template void _2d_r_divf(double * a, double * b, double *output,int M, int N);
template void _2d_r_modf(double * a, double * b, double *output,int M, int N);

template void matrix_r_addf(double * a, double * b, double *output, int M, int N, int K);
template void matrix_r_subf(double * a, double * b, double *output, int M, int N, int K);
template void matrix_r_hadamardf(double * a, double * b, double *output, int M, int N, int K);
template void matrix_r_multiplyf(double * a, double * b, double * output, int M, int N, int K);
template void matrix_r_transposef(double * input, double *output, int M, int N);
template void matrix_r_acosf(double * input, double * output,  int M, int N);
template void matrix_r_asinf(double * input, double * output,  int M, int N);
template void matrix_r_atanf(double * input, double * output,  int M, int N);
template void matrix_r_atan2f(double * x, double * y, double * output, int M, int N);
template void matrix_r_acoshf(double * input, double * output,  int M, int N);
template void matrix_r_asinhf(double * input, double * output,  int M, int N);
template void matrix_r_atanhf(double * input, double * output,  int M, int N);
template void matrix_r_cosf(double * input, double * output,  int M, int N);
template void matrix_r_sinf(double * input, double * output,  int M, int N);
template void matrix_r_tanf(double * input, double * output, int M, int N);
template void matrix_r_coshf(double * input, double * output,  int M, int N);
template void matrix_r_sinhf(double * input, double * output,  int M, int N);
template void matrix_r_tanhf(double * input, double * output, int M, int N);
template void matrix_r_atan2f_const(double * input, double b, double * output, int M, int N);
template void matrix_r_ceilf(double * input, double *output, int M, int N);
template void matrix_r_exp10f(double * input, double *output, int M, int N);
template void matrix_r_exp2f(double * input, double *output, int M, int N);
template void matrix_r_expf(double * input, double *output, int M, int N);
template void matrix_r_expm1f(double * input, double *output, int M, int N);
template void matrix_r_fabsf(double * input, double *output, int M, int N);
template void matrix_r_floorf(double * input, double *output, int M, int N);
template void matrix_r_fmaxf(double * x, double *y, double *output, int M, int N);
template void matrix_r_fminf(double * x, double *y, double *output, int M, int N);
template void matrix_r_fmodf(double * x, double *y, double *output, int M, int N);
template void matrix_r_log10f(double * input, double *output, int M, int N);
template void matrix_r_log1pf(double * input, double *output, int M, int N);
template void matrix_r_log2f(double * input, double *output, int M, int N);
template void matrix_r_logbf(double * input, double *output, int M, int N);
template void matrix_r_powf(double * x, double *y, double *output, int M, int N);
template void matrix_r_rsqrtf(double * input, double *output, int M, int N);
template void matrix_r_sqrtf(double * input, double *output, int M, int N);
template void matrix_r_cbrtf(double * input, double *output, int M, int N);
template void matrix_r_cospif(double * input, double *output, int M, int N);
template void matrix_r_cyl_bessel_i0f(double * input, double *output, int M, int N);
template void matrix_r_cyl_bessel_i1f(double * input, double *output, int M, int N);
template void matrix_r_erfcf(double * input, double *output, int M, int N);
template void matrix_r_erfcinvf(double * input, double *output, int M, int N);
template void matrix_r_erfcxf(double * input, double *output, int M, int N);
template void matrix_r_erff(double * input, double *output, int M, int N);
template void matrix_r_erfinvf(double * input, double * outputs, int M, int N);
template void matrix_r_fdimf(double * x, double * y, double *output, int M, int N);
template void matrix_r_fdividef(double * x, double *y, double *output, int M, int N);
template void matrix_r_fmaf(double * x, double *y, double *z, double *output, int M, int N);
template void matrix_r_hypotf(double * x, double *y, double *output, int M, int N);
template void matrix_r_ilogbf(double * input, double *output, int M, int N);
template void matrix_r_j0f(double * input, double *output, int M, int N);
template void matrix_r_j1f(double * input, double *output, int M, int N);
template void matrix_r_jnf(double * input, double *output, int n, int M, int N);
template void matrix_r_ldexpf(double * input, double *output, int exp, int M, int N);
template void matrix_r_lgammaf(double * input, double *output, int M, int N);
template void matrix_r_nearbyintf(double * input, double *output, int M, int N);
template void matrix_r_norm3df(double * x, double *y, double *z, double *output, int M, int N);
template void matrix_r_norm4df(double * x, double *y, double *z, double * w, double *output,int M, int N);
template void matrix_r_normcdff(double * input, double *output, int M, int N);
template void matrix_r_normcdfinvf(double * input, double *output, int M, int N);
template void matrix_r_normf(int dim, double *input, double *output, int M, int N);
template void matrix_r_remainderf(double * x, double *y, double *output, int M, int N);
template void matrix_r_rcbrtf(double * input, double * output, int M, int N);
template void matrix_r_rhypotf(double * x, double *y, double * output, int M, int N);
template void matrix_r_rnorm3df(double * x, double *y, double *z, double * output, int M, int N);
template void matrix_r_rnorm4df(double * x, double *y, double *z, double *w, double *output, int M, int N);
template void matrix_r_rnormf(int dim, double *input, double *output, int M, int N);
template void matrix_r_scalblnf(double * input, long int n,  double *output, int M, int N);
template void matrix_r_sinpif(double * input, double *output, int M, int N);
template void matrix_r_tgammaf(double * input, double *output, int M, int N);
template void matrix_r_y0f(double * input, double *output, int M, int N);
template void matrix_r_y1f(double * input, double *output, int M, int N);
template void matrix_r_ynf(int n, double * input, double *output, int M, int N);
template void matrix_r_softmaxf(double * input, double *output, int M, int N);
template void matrix_r_sigmoidf(double * x, double *output, int M, int N);
template void matrix_r_sigmoid_gradf(double * x, double *output,int M, int N);
template void matrix_r_tanh_gradf(double * x, double *output, int M, int N);
template void matrix_r_reluf(double * x, double *output, int M, int N);
template void matrix_r_relu_gradf(double * x, double *output, int M, int N);
template void matrix_r_addf_const(double * x, double y, double *output, int M, int N);
template void matrix_r_subf_const(double * x, double y, double *output, int M, int N);
template void matrix_r_mulf_const(double * x, double y, double *output,int M, int N);
template void matrix_r_divf_const(double * x, double y, double *output, int M, int N);
template void matrix_r_divf_matrix_const(double * x, double y, double *output, int M, int N);
template void matrix_r_modf_const(double * x, double y, double *output, int M, int N);
template void matrix_r_fmaxf_const(double * x, double y, double *output, int M, int N);
template void matrix_r_fminf_const(double * x, double y, double *output, int M, int N);
template void matrix_r_fmodf_const(double * x, double y, double *output, int M, int N);
template void matrix_r_powf_const(double * x, double y, double *output, int M, int N);
template void matrix_r_fdimf_const(double * x, double y, double *output, int M, int N);
template void matrix_r_fdividef_const(double * x, double y, double *output, int M, int N);
template void matrix_r_hypotf_const(double * x, double y, double *output, int M, int N);
template void matrix_r_remainderf_const(double * x, double y, double * output, int M, int N);
template void matrix_r_rhypotf_const(double * x, double y, double *output, int M, int N);
template void matrix_r_addf_scalar(double * x, double *y, double *output, int M, int N);
template void matrix_r_subf_scalar(double * x, double *y, double *output, int M, int N);
template void matrix_r_mulf_scalar(double * x, double *y, double *output, int M, int N);
template void matrix_r_divf_scalar(double * x, double *y, double *output, int M, int N);
template void matrix_r_modf_scalar(double * x, double *y, double *output, int M, int N);
template void matrix_r_atan2f_scalar(double * a, double *b, double *output,int M, int N);
template void matrix_r_fmaxf_scalar(double * x, double *y, double *output,int M, int N);
template void matrix_r_fminf_scalar(double * x, double *y, double *output,int M, int N);
template void matrix_r_fmodf_scalar(double * x, double *y, double *output,int M, int N);
template void matrix_r_powf_scalar(double * x, double *y, double *output,int M, int N);
template void matrix_r_fdimf_scalar(double * x, double *y, double *output, int M, int N);
template void matrix_r_fdividef_scalar(double * x, double *y, double *output, int M, int N);
template void matrix_r_hypotf_scalar(double * x, double *y, double *output, int M, int N);
template void matrix_r_remainderf_scalar(double * x, double *y, double * output, int M, int N);
template void matrix_r_rhypotf_scalar(double * x, double *y, double *output, int M, int N);
template void matrix_r_copysignf(double * X, double *Y, double *output, int M, int N);
template void matrix_r_truncf(double * x, double *output, int M, int N);


template<typename T> __global__ void gpu_2d_addf(T *a,T *b, T *c, int m, int n)
{ 
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;    
    if( col < n && row < m) 
    {        
        c[row * n + col] = a[row * n + col] + b[row*n + col];
    }
} 

template<typename T> void _2d_r_addf(T * a, T * b, T *output, int M, int N)
{
       
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);
    gpu_2d_addf<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(a,b,output,M,N);
}

template<typename T> __global__ void gpu_2d_mulf(T *a,T *b, T *c, int m, int n)
{ 
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;    
    if( col < n && row < m) 
    {        
        c[row * n + col] = a[row * n + col] * b[row*n + col];
    }
} 


template<typename T> void _2d_r_mulf(T * a, T * b, T *output, int M, int N)
{
    
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);
    gpu_2d_mulf<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(a,b,output,M,N);
}


template<typename T> __global__ void gpu_2d_subf(T *a,T *b, T *c, int m, int n)
{ 
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;    
    if( col < n && row < m) 
    {        
        c[row * n + col] = a[row * n + col] - b[row*n + col];
    }
} 


template<typename T>
void _2d_r_subf(T * a, T * b, T *output, int M, int N)
{   
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);
    gpu_2d_subf<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(a,b,output,M,N);
}

template<typename T> __global__ void gpu_2d_divf(T *a,T *b, T *c, int m, int n)
{ 
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;    
    if( col < n && row < m) 
    {        
        c[row * n + col] = a[row * n + col] / b[row*n + col];
    }
} 


template<typename T>
void _2d_r_divf(T * a, T * b, T *output, int M, int N)
{
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);
    gpu_2d_divf<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(a,b,output,M,N);
}

template<typename T> __global__ void gpu_2d_modf(T *a,T *b, T *c, int m, int n)
{ 
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;    
    if( col < n && row < m) 
    {        
        c[row * n + col] = fmodf(a[row * n + col],b[row*n + col]);
    }
} 


template<typename T> void _2d_r_modf(T * a, T * b, T *output, int M, int N)
{
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);
    gpu_2d_modf<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(a,b,output,M,N);
}



template<typename T> __global__ void gpu_matrix_addf(T *a,T *b, T *c, int m, int n, int k)
{ 
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;    
    if( col < k && row < m) 
    {        
        c[row * n + col] = a[row * n + col] + b[row*n + col];
    }
} 

template<typename T> void matrix_r_addf(T * a, T * b, T *output, int M, int N, int K)
{
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (K + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);
    gpu_matrix_addf<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(a,b,output,M,N,K);
}

template<typename T> __global__ void gpu_matrix_subf(T *a,T *b, T *c, int m, int n, int k)
{ 
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;    
    if( col < k && row < m) 
    {                
        c[row * n + col] = a[row * n + col] - b[row*n + col];
    }
} 


template<typename T> void matrix_r_subf(T * a, T * b, T *output, int M, int N, int K)
{  
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (K + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);
    gpu_matrix_subf<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(a,b,output,M,N,K); 
}



template<typename T> __global__ void gpu_matrix_hadamardf(T *a,T *b, T *c, int m, int n, int k)
{ 
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;    
    if( col < k && row < m) 
    {                
        c[row * n + col] = a[row * n + col] * b[row * k + col];
    }
} 


template<typename T> void matrix_r_hadamardf(T * a, T * b, T *output, int M, int N, int K)
{   
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (K + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);
    gpu_matrix_hadamardf<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(a,b,output,M,N,K);
}


template<typename T> __global__ void gpu_matrix_multiplyf(T *a,T *b, T *c, int m, int n, int k)
{ 
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    T sum = 0;    
    
    if( col < k && row < m) 
    {     
        for(int i = 0; i < n; i++) 
        {            
            sum += a[row*n+i]*b[i*k+col];
        }        
        c[row*n + col] = sum;
    }

} 

template<typename T> void matrix_r_multiplyf(T * a, T * b, T * output, int M, int N, int K)
{    
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (K + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_rows,grid_cols);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_matrix_multiplyf<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(a,b,output,M,N,K);
}

template<typename T> __global__ void gpu_matrix_transposef(T* mat_in, T* mat_out, unsigned int rows, unsigned int cols) 
{
    unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int idy = blockIdx.y * blockDim.y + threadIdx.y;

    if (idx < cols && idy < rows) 
    {
        
        unsigned int pos        =  idx * rows + idy;
        unsigned int trans_pos  =  idy * cols + idx;
        mat_out[trans_pos] = mat_in[pos];    
    }
}



template<typename T> void matrix_r_transposef(T * input, T *output, int M, int N)
{
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);
    gpu_matrix_transposef<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(input,output,M,N);
}

template<typename T> __global__ void gpu_matrix_softmaxf(T* x, T* c, T sum, int m, int n) 
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (col < n && row < m) 
    {        
        c[row*n + col] = Exp(x[row*n + col])/sum;        
    }
}


template<typename T> void matrix_r_softmaxf(T * input, T *output, int M, int N)
{   
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);
    T sum = vector_sumf(input,M*N);
    gpu_matrix_softmaxf<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(input,output,sum,M,N);
}


template<typename T> __global__ void gpu_matrix_acosf(T* x, T* c, int m, int n) 
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (col < n && row < m) 
    {        
        c[row*n + col] = ACos(x[row*n + col]);
    }
}

template<typename T> void matrix_r_acosf(T * input, T *output, int M, int N)
{
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_matrix_acosf<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(input,output,M,N);
}


template<typename T> __global__ void gpu_matrix_asinf(T* x, T* c, int m, int n) 
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (col < n && row < m) 
    {        
        c[row*n + col] = ASin(x[row*n + col]);
    }
}


template<typename T> void matrix_r_asinf(T * input, T * output,  int M, int N)
{    
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_matrix_asinf<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(input,output,M,N);
}

template<typename T> __global__ void gpu_matrix_atanf(T* x, T* c, int m, int n) 
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (col < n && row < m) 
    {        
        c[row*n + col] = ATan(x[row*n + col]);
    }
}


template<typename T> void matrix_r_atanf(T * input, T * output, int M, int N)
{    
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_matrix_atanf<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(input,output,M,N);    
}

template<typename T> __global__ void gpu_matrix_atan2f(T* x, T *y, T* c, int m, int n) 
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (col < n && row < m) 
    {        
        c[row*n + col] = ATan2(x[row*n + col],y[row*n+col]);
    }
}


template<typename T> void matrix_r_atan2f(T * x, T * y, T * output, int M, int N)
{    
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_matrix_atan2f<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(x,y,output,M,N); 
}


template<typename T> __global__ void gpu_matrix_cosf(T* x, T* c, int m, int n) 
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (col < n && row < m) 
    {        
        c[row*n + col] = ACos(x[row*n + col]);
    }
}


template<typename T> void matrix_r_cosf(T * input, T * output, int M, int N)
{   
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_matrix_cosf<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(input,output,M,N); 
}

template<typename T> __global__ void gpu_matrix_sinf(T* x, T* c, int m, int n) 
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (col < n && row < m) 
    {        
        c[row*n + col] = ASin(x[row*n + col]);
    }
}


template<typename T> void matrix_r_sinf(T * input, T * output, int M, int N)
{   
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_matrix_sinf<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(input,output,M,N); 
}

template<typename T> __global__ void gpu_matrix_tanf(T* x, T* c, int m, int n) 
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (col < n && row < m) 
    {        
        c[row*n + col] = ATan(x[row*n + col]);
    }
}

template<typename T> void matrix_r_tanf(T * input, T * output, int M, int N)
{
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_matrix_tanf<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(input,output,M,N);
}

template<typename T> __global__ void gpu_matrix_coshf(T* x, T* c, int m, int n) 
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (col < n && row < m) 
    {        
        c[row*n + col] = ACosh(x[row*n + col]);
    }
}

template<typename T> void matrix_r_coshf(T * input, T * output, int M, int N)
{    
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_matrix_coshf<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(input,output,M,N); 
}

template<typename T> __global__ void gpu_matrix_sinhf(T* x, T* c, int m, int n) 
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (col < n && row < m) 
    {        
        c[row*n + col] = ASinh(x[row*n + col]);
    }
}

template<typename T> void matrix_r_sinhf(T * input, T *output, int M, int N)
{
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_matrix_sinhf<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(input,output,M,N);
}

template<typename T> __global__ void gpu_matrix_tanhf(T* x, T* c, int m, int n) 
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (col < n && row < m) 
    {        
        c[row*n + col] = Tanh(x[row*n + col]);
    }
}

template<typename T> void matrix_r_tanhf(T * input, T * output, int M, int N)
{
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_matrix_tanhf<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(input,output,M,N);
}

template<typename T> __global__ void gpu_matrix_acoshf(T* x, T* c, int m, int n) 
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (col < n && row < m) 
    {        
        c[row*n + col] = ACosh(x[row*n + col]);
    }
}


template<typename T> void matrix_r_acoshf(T * input, T * output, int M, int N)
{    
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_matrix_acoshf<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(input,output,M,N); 
}

template<typename T> __global__ void gpu_matrix_asinhf(T* x, T* c, int m, int n) 
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (col < n && row < m) 
    {        
        c[row*n + col] = ASinh(x[row*n + col]);
    }
}

template<typename T> void matrix_r_asinhf(T * input, T * output, int M, int N)
{   
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_matrix_asinhf<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(input,output,M,N); 
}

template<typename T> __global__ void gpu_matrix_atanhf(T* x, T* c, int m, int n) 
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (col < n && row < m) 
    {        
        c[row*n + col] = ATanh(x[row*n + col]);
    }
}

template<typename T> void matrix_r_atanhf(T * input, T *output, int M, int N)
{  
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_matrix_atanhf<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(input,output,M,N); 
}


template<typename T> __global__ void gpu_matrix_ceilf(T* x, T* c, int m, int n) 
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (col < n && row < m) 
    {        
        c[row*n + col] = Ceil(x[row*n + col]);
    }
}

template<typename T>
void matrix_r_ceilf(T * input, T *output, int M, int N)
{  
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_matrix_ceilf<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(input,output,M,N); 
}



template<typename T> __global__ void gpu_matrix_exp10f(T* x, T* c, int m, int n) 
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (col < n && row < m) 
    {        
        c[row*n + col] = Exp10(x[row*n + col]);
    }
}


template<typename T>
void matrix_r_exp10f(T * input, T *output, int M, int N)
{
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_matrix_exp10f<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(input,output,M,N);
    
}

template<typename T> __global__ void gpu_matrix_exp2f(T* x, T* c, int m, int n) 
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (col < n && row < m) 
    {        
        c[row*n + col] = Exp2(x[row*n + col]);
    }
}

template<typename T>
void matrix_r_exp2f(T * input, T *output, int M, int N)
{  
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_matrix_exp2f<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(input,output,M,N); 
}

template<typename T> __global__ void gpu_matrix_expf(T* x, T* c, int m, int n) 
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    int idx = row*n + col;
    if (col < n && row < m)     {        
        c[idx] = expf(x[idx]);
    }    
}

template<typename T>
void matrix_r_expf(T * input, T *output, int M, int N)
{
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_matrix_expf<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(input,output,M,N);
}


template<typename T> __global__ void gpu_matrix_expm1f(T* x, T* c, int m, int n) 
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (col < n && row < m) 
    {        
        c[row*n + col] = Expm1(x[row*n + col]);
    }
}

template<typename T>
void matrix_r_expm1f(T * input, T *output, int M, int N)
{
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_matrix_expm1f<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(input,output,M,N);
}

template<typename T> __global__ void gpu_matrix_fabsf(T* x, T* c, int m, int n) 
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (col < n && row < m) 
    {        
        c[row*n + col] = Fabs(x[row*n + col]);
    }
}

template<typename T>
void matrix_r_fabsf(T * input, T *output, int M, int N)
{   
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_matrix_fabsf<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(input,output,M,N); 
}


template<typename T> __global__ void gpu_matrix_floorf(T* x, T* c, int m, int n) 
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (col < n && row < m) 
    {        
        c[row*n + col] = Floor(x[row*n + col]);
    }
}

template<typename T>
void matrix_r_floorf(T * input, T *output, int M, int N)
{
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_matrix_floorf<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(input,output,M,N);
}


template<typename T> __global__ void gpu_matrix_fmaxf(T* x, T * y, T* c, int m, int n) 
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (col < n && row < m) 
    {        
        c[row*n + col] = Fmax(x[row*n + col],y[row*n + col]);
    }
}

template<typename T>
void matrix_r_fmaxf(T * x, T *y, T *output, int M, int N)
{
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_matrix_fmaxf<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(x,y,output,M,N);
}

template<typename T> __global__ void gpu_matrix_fminf(T* x, T * y, T* c, int m, int n) 
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (col < n && row < m) 
    {        
        c[row*n + col] = Fmin(x[row*n + col],y[row*n + col]);
    }
}

template<typename T>
void matrix_r_fminf(T * x, T *y, T *output, int M, int N)
{
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_matrix_fminf<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(x,y,output,M,N);    
}


template<typename T> __global__ void gpu_matrix_fmodf(T* x, T * y, T* c, int m, int n) 
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (col < n && row < m) 
    {        
        c[row*n + col] = Fmod(x[row*n + col],y[row*n + col]);
    }
}

template<typename T>
void matrix_r_fmodf(T * x, T *y, T *output, int M, int N)
{
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_matrix_fmodf<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(x,y,output,M,N);
}



template<typename T> __global__ void gpu_matrix_log10f(T* x, T* c, int m, int n) 
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (col < n && row < m) 
    {        
        c[row*n + col] = Log10(x[row*n + col]);
    }
}

template<typename T>
void matrix_r_log10f(T * input, T *output, int M, int N)
{
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_matrix_log10f<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(input,output,M,N);
}


template<typename T> __global__ void gpu_matrix_log1pf(T* x, T* c, int m, int n) 
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (col < n && row < m) 
    {        
        c[row*n + col] = Log1p(x[row*n + col]);
    }
}

template<typename T> 
void matrix_r_log1pf(T * input, T *output, int M, int N)
{
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_matrix_log1pf<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(input,output,M,N);
}

template<typename T> __global__ void gpu_matrix_log2f(T* x, T* c, int m, int n) 
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (col < n && row < m) 
    {        
        c[row*n + col] = Log2(x[row*n + col]);
    }
}

template<typename T> 
void matrix_r_log2f(T * input, T *output, int M, int N)
{
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_matrix_log2f<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(input,output,M,N);
}


template<typename T> __global__ void gpu_matrix_logbf(T* x, T* c, int m, int n) 
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (col < n && row < m) 
    {        
        c[row*n + col] = Logb(x[row*n + col]);
    }
}

template<typename T> 
void matrix_r_logbf(T * input, T *output, int M, int N)
{
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_matrix_logbf<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(input,output,M,N);
}



template<typename T> __global__ void gpu_matrix_powf(T* x, T * y, T* c, int m, int n) 
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    int idx = row*n + col;
    if (col < n && row < m) 
    {        
        if(x[idx] == 0) c[idx]=0;
        else if(x[idx]==1) c[idx]=1;
        else if(x[idx] < 0) c[idx]=-Pow(Fabs(x[idx]),y[idx]);
        else c[idx] = Pow(x[idx],y[idx]);
    }
}

template<typename T> 
void matrix_r_powf(T * x, T *y, T *output, int M, int N)
{
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_matrix_powf<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(x,y,output,M,N);
}


template<typename T> __global__ void gpu_matrix_rsqrtf(T* x, T* c, int m, int n) 
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (col < n && row < m) 
    {        
        if(x[row*n + col] <= 0) c[row*n + col] = 0;
        else c[row*n + col] = RSqrt(x[row*n + col]);
    }
}

template<typename T> 
void matrix_r_rsqrtf(T * input, T *output, int M, int N)
{
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_matrix_rsqrtf<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(input,output,M,N);
}

template<typename T> __global__ void gpu_matrix_sqrtf(T* x, T* c, int m, int n) 
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (col < n && row < m) 
    {        
        if(x[row*n + col] <= 0) c[row*n + col] = 0;
        else c[row*n + col] = Sqrt(x[row*n + col]);
    }
}

template<typename T> 
void matrix_r_sqrtf(T * input, T *output, int M, int N)
{
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_matrix_sqrtf<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(input,output,M,N);
}






template<typename T> __global__ void matrix_sigmoid_device(T * x, T * out, int m, int n)
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    int idx = row*n + col; 
    if (col < n && row < m) 
    {        
        T r = 1.0f / (1.0f + Exp(-x[idx]));
        out[idx] = r;
    }
}

template<typename T> 
void matrix_r_sigmoidf(T * x, T *output, int M, int N)
{
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);       
    matrix_sigmoid_device<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(x,output,M,N);        
}

template<typename T> __global__ void matrix_sigmoid_grad_device(T * x, T * out, int m, int n)
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    int idx = row*n + col; 
    if (col < n && row < m) 
    {        
        T r = x[idx] * (1.0f - x[idx]);
        out[idx] = r;
    }    
}


template<typename T> 
void matrix_r_sigmoid_gradf(T * x, T *output,int M, int N)
{
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    matrix_sigmoid_grad_device<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(x,output,M,N);        
}


template<typename T> __global__ void matrix_tanh_grad_device(T * x, T * out, int m, int n)
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    int idx = row*n + col; 
    if (col < n && row < m) 
    {        
        out[idx] = 1.0 - (x[idx] *x[idx]);
    }    
}


template<typename T>
void matrix_r_tanh_gradf(T * x, T *output, int M, int N)
{
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    matrix_tanh_grad_device<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(x,output,M,N);        
}

template<typename T> __global__ void matrix_relu_device(T * x, T * out, int m, int n)
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    int idx = row*n + col; 
    if (col < n && row < m) 
    {   
        if(x[idx] < 0) out[idx] = 0.0f;        
        else out[idx] = x[idx];                
    }
}    

template<typename T> 
void matrix_r_reluf(T * x, T *output, int M, int N)
{    
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    matrix_relu_device<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(x,output,M,N);         
}

template<typename T> __global__ void matrix_relu_grad_device(T * x, T * out, int m, int n)
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    int idx = row*n + col; 
    if (col < n && row < m) 
    {        
        if(x[idx] > 0) out[idx] = 1.0;
        else out[idx] = 0.0f;
    }
}

template<typename T> 
void matrix_r_relu_gradf(T * x, T *output, int M, int N)
{
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    matrix_relu_grad_device<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(x,output,M,N);        
}

template<typename T> __global__ void matrix_add_const_device(T * x, T y, T * out, int m, int n)
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    int idx = row*n + col; 
    if (col < n && row < m) 
    {        
        out[idx] = x[idx] + y;
    }
}


template<typename T> 
void matrix_r_addf_const(T * x, T y, T *output, int M, int N)
{    
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    matrix_add_const_device<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(x,y,output,M,N);     
}

template<typename T> __global__ void matrix_sub_const_device(T * x, T y, T * out, int m, int n)
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    int idx = row*n + col; 
    if (col < n && row < m) 
    {        
        out[idx] = x[idx] - y;
    }
}


template<typename T> 
void matrix_r_subf_const(T * x, T y, T *output, int M, int N)
{
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    matrix_sub_const_device<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(x,y,output,M,N);    
}

template<typename T> __global__ void matrix_mul_const_device(T * x, T y, T * out, int m, int n)
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    int idx = row*n + col; 
    if (col < n && row < m) 
    {        
        out[idx] = x[idx] * y;
    }
}

template<typename T> 
void matrix_r_mulf_const(T * x, T y, T *output,int M, int N)
{
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    matrix_mul_const_device<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(x,y,output,M,N);    
}


template<typename T> __global__ void matrix_div_const_device(T * x, T y, T * out, int m, int n)
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    int idx = row*n + col; 
    if (col < n && row < m && y != 0.0f) 
    {        
        out[idx] = x[idx] / y;
    }
}

template<typename T> 
void matrix_r_divf_const(T * x, T y, T *output, int M, int N)
{
    assert(y != 0.0f);
    
    // if(N < 1024) BLOCK_SIZE=N;    
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    matrix_div_const_device<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(x,y,output,M,N);    
}

template<typename T> __global__ void matrix_div_matrix_const_device(T * x, T y, T * out, int m, int n)
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    int idx = row*n + col; 
    if (col < n && row < m && y != 0.0f) 
    {        
        out[idx] = y / x[idx];
    }
}

template<typename T> 
void matrix_r_divf_matrix_const(T * x, T y, T *output, int M, int N)
{
    assert(y != 0.0f);
    
    // if(N < 1024) BLOCK_SIZE=N;    
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    matrix_div_matrix_const_device<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(x,y,output,M,N);    
}

template<typename T> __global__ void matrix_mod_const_device(T * x, T y, T * out, int m, int n)
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    int idx = row*n + col; 
    if (col < n && row < m) 
    {        
        out[idx] = Fmod(x[idx],y);
    }
}

template<typename T> 
void matrix_r_modf_const(T * x, T y, T *output, int M, int N)
{
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    matrix_mod_const_device<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(x,y,output,M,N);    
}


template<typename T> __global__ void matrix_atan2f_const_device(T * x, T y, T * out, int m, int n)
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    int idx = row*n + col; 
    if (col < n && row < m) 
    {        
        out[idx] = ATan2(x[idx],y);
    }
}


template<typename T> 
void matrix_r_atan2f_const(T * x, T y, T * output, int M, int N)
{    
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    matrix_atan2f_const_device<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(x,y,output,M,N);     
}


template<typename T> __global__ void matrix_fmaxf_const_device(T * x, T y, T * out, int m, int n)
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    int idx = row*n + col; 
    if (col < n && row < m) 
    {        
        out[idx] = Fmax(x[idx],y);
    }
}


template<typename T> 
void matrix_r_fmaxf_const(T * x, T y, T *output, int M, int N)
{
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    matrix_fmaxf_const_device<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(x,y,output,M,N);    
}


template<typename T> __global__ void matrix_fminf_const_device(T * x, T y, T * out, int m, int n)
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    int idx = row*n + col; 
    if (col < n && row < m) 
    {        
        out[idx] = Fmin(x[idx],y);
    }
}

template<typename T> 
void matrix_r_fminf_const(T * x, T y, T *output, int M, int N)
{
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    matrix_fminf_const_device<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(x,y,output,M,N);    
}

template<typename T> __global__ void matrix_fmodf_const_device(T * x, T y, T * out, int m, int n)
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    int idx = row*n + col; 
    if (col < n && row < m) 
    {        
        out[idx] = Fmod(x[idx],y);
    }
}

template<typename T> 
void matrix_r_fmodf_const(T * x, T y, T *output, int M, int N)
{
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    matrix_fmodf_const_device<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(x,y,output,M,N);    
}

template<typename T> __global__ void matrix_powf_const_device(T * x, T y, T * out, int m, int n)
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    int idx = row*n + col; 
    if (col < n && row < m) 
    {        
        if(x[idx] == 0) out[idx]=0;
        else if(x[idx]==1) out[idx]=1;
        else if(x[idx] < 0) out[idx]=-Pow(Fabs(x[idx]),y);
        else out[idx] = Pow(x[idx],y);        
    }
}

template<typename T> 
void matrix_r_powf_const(T * x, T y, T *output, int M, int N)
{       
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    matrix_powf_const_device<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(x,y,output,M,N);    
}



template<typename T> __global__ void matrix_add_scalar_device(T * x, T *y, T * out, int m, int n)
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    int idx = row*n + col; 
    if (col < n && row < m) 
    {        
        out[idx] = x[idx] +y[0];
    }
}

template<typename T> 
void matrix_r_addf_scalar(T * x, T *y, T *output, int M, int N)
{
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    matrix_add_scalar_device<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(x,y,output,M,N);    
}

template<typename T> __global__ void matrix_subf_scalar_device(T * x, T *y, T * out, int m, int n)
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    int idx = row*n + col; 
    if (col < n && row < m) 
    {        
        out[idx] = x[idx] - y[0];
    }
}

template<typename T> 
void matrix_r_subf_scalar(T * x, T *y, T *output, int M, int N)
{
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    matrix_add_scalar_device<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(x,y,output,M,N);    
}

template<typename T> __global__ void matrix_mul_scalar_device(T * x, T *y, T * out, int m, int n)
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    int idx = row*n + col; 
    if (col < n && row < m) 
    {        
        out[idx] = x[idx] * y[0];
    }
}

template<typename T> 
void matrix_r_mulf_scalar(T * x, T *y, T *output, int M, int N)
{
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    matrix_mul_scalar_device<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(x,y,output,M,N);    
}

template<typename T> __global__ void matrix_div_scalar_device(T * x, T *y, T * out, int m, int n)
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    int idx = row*n + col; 
    if (col < n && row < m) 
    {        
        if(y[0] != 0.0)
            out[idx] = x[idx] / y[0];
        // if there is some reason not to do this, it can be changed.
        else 
            out[idx] = CUDART_NAN_F;
    }
}

template<typename T> void matrix_r_divf_scalar(T * x, T *y, T *output, int M, int N)
{
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    matrix_div_scalar_device<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(x,y,output,M,N);    
}


template<typename T> __global__ void matrix_mod_scalar_device(T * x, T *y, T * out, int m, int n)
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    int idx = row*n + col; 
    if (col < n && row < m) 
    {            
        out[idx] = Fmod(x[idx],y[0]);        
    }
}

template<typename T> void matrix_r_modf_scalar(T * x, T *y, T *output, int M, int N)
{
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    matrix_mod_scalar_device<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(x,y,output,M,N);    
}


template<typename T> __global__ void matrix_fmax_scalar_device(T * x, T *y, T * out, int m, int n)
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    int idx = row*n + col; 
    if (col < n && row < m) 
    {            
        out[idx] = Fmax(x[idx],y[0]);        
    }
}

template<typename T> void matrix_r_fmaxf_scalar(T * x, T *y, T *output, int M, int N)
{
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    matrix_fmax_scalar_device<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(x,y,output,M,N);    
}

template<typename T> __global__ void matrix_fmin_scalar_device(T * x, T *y, T * out, int m, int n)
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    int idx = row*n + col; 
    if (col < n && row < m) 
    {            
        out[idx] = Fmin(x[idx],y[0]);        
    }
}

template<typename T> void matrix_r_fminf_scalar(T * x, T *y, T *output, int M, int N)
{
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    matrix_fmin_scalar_device<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(x,y,output,M,N);    
}


template<typename T> __global__ void matrix_pow_scalar_device(T * x, T *y, T * out, int m, int n)
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    int idx = row*n + col; 
    if (col < n && row < m) 
    {        
        if(x[idx] == 0) out[idx]=0;
        else if(x[idx]==1) out[idx]=1;
        else if(x[idx] < 0) out[idx]=-Pow(Fabs(x[idx]),y[0]);    
        else out[idx] = Pow(x[idx],y[0]);        
    }
}

template<typename T> void matrix_r_powf_scalar(T * x, T *y, T *output, int M, int N)
{
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    matrix_pow_scalar_device<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(x,y,output,M,N);    
}

template<typename T> __global__ void matrix_hypot_scalar_device(T * x, T *y, T * out, int m, int n)
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    int idx = row*n + col; 
    if (col < n && row < m) 
    {            
        out[idx] = Hypot(x[idx],y[0]);        
    }
}

template<typename T> void matrix_r_hypotf_scalar(T * x, T *y, T *output, int M, int N)
{
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    matrix_hypot_scalar_device<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(x,y,output,M,N);    
}


template<typename T> __global__ void matrix_rhypot_scalar_device(T * x, T *y, T * out, int m, int n)
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    int idx = row*n + col; 
    if (col < n && row < m) 
    {            
        out[idx] = RHypot(x[idx],y[0]);        
    }
}

template<typename T> void matrix_r_rhypotf_scalar(T * x, T *y, T *output, int M, int N)
{
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    matrix_rhypot_scalar_device<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(x,y,output,M,N);    
}

template<typename T> __global__ void matrix_fdividef_scalar_device(T * x, T *y, T * out, int m, int n)
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    int idx = row*n + col; 
    if (col < n && row < m) 
    {            
        out[idx] = FDivide(x[idx],y[0]);        
    }
}

template<typename T> void matrix_r_fdividef_scalar(T * x, T *y, T *output, int M, int N)
{
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    matrix_fdividef_scalar_device<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(x,y,output,M,N);    
}


template<typename T> __global__ void matrix_fmodf_scalar_device(T * x, T *y, T * out, int m, int n)
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    int idx = row*n + col; 
    if (col < n && row < m) 
    {            
        out[idx] = Fmod(x[idx],y[0]);        
    }
}

template<typename T> void matrix_r_fmodf_scalar(T * x, T *y, T *output, int M, int N)
{
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    matrix_fmodf_scalar_device<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(x,y,output,M,N);    
}


template<typename T> __global__ void matrix_remainderf_scalar_device(T * x, T *y, T * out, int m, int n)
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    int idx = row*n + col; 
    if (col < n && row < m) 
    {            
        out[idx] = Remainder(x[idx],y[0]);        
    }
}

template<typename T> void matrix_r_remainderf_scalar(T * x, T *y, T *output, int M, int N)
{
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    matrix_remainderf_scalar_device<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(x,y,output,M,N);    
}


template<typename T> __global__ void matrix_fdimf_scalar_device(T * x, T *y, T * out, int m, int n)
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    int idx = row*n + col; 
    if (col < n && row < m) 
    {            
        out[idx] = FDim(x[idx],y[0]);        
    }
}

template<typename T> void matrix_r_fdimf_scalar(T * x, T *y, T *output, int M, int N)
{
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    matrix_fdimf_scalar_device<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(x,y,output,M,N);    
}

template<typename T> __global__ void matrix_atan2f_scalar_device(T * x, T *y, T * out, int m, int n)
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    int idx = row*n + col; 
    if (col < n && row < m) 
    {            
        out[idx] = ATan2(x[idx],y[0]);        
    }
}

template<typename T> void matrix_r_atan2f_scalar(T * x, T *y, T *output, int M, int N)
{
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    matrix_atan2f_scalar_device<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(x,y,output,M,N);    
}

////////////////////////////////////////////////////////////////////////////////
// matrix

template<typename T> __global__ void gpu_matrix_cbrtf(T* x, T* c, int m, int n) 
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (col < n && row < m) 
    {        
        c[row*n + col] = Cbrt(x[row*n + col]);
    }
}


template<typename T> void matrix_r_cbrtf(T * input, T *output, int M, int N)
{
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_matrix_cbrtf<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(input,output,M,N);
}

template<typename T> __global__ void gpu_matrix_cospif(T* x, T* c, int m, int n) 
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (col < n && row < m) 
    {        
        c[row*n + col] = Cospi(x[row*n + col]);
    }
}


template<typename T> void matrix_r_cospif(T * input, T *output, int M, int N)
{
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_matrix_cospif<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(input,output,M,N);
}


template<typename T> __global__ void gpu_matrix_cyl_bessel_i0f(T* x, T* c, int m, int n) 
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (col < n && row < m) 
    {        
        c[row*n + col] = Besseli0(x[row*n + col]);
    }
}


template<typename T> 
void matrix_r_cyl_bessel_i0f(T * input, T *output, int M, int N)
{
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_matrix_cyl_bessel_i0f<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(input,output,M,N);
}

template<typename T> __global__ void gpu_matrix_cyl_bessel_i1f(T* x, T* c, int m, int n) 
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (col < n && row < m) 
    {        
        c[row*n + col] = Besseli1(x[row*n + col]);
    }
}


template<typename T> void matrix_r_cyl_bessel_i1f(T * input, T *output, int M, int N)
{
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_matrix_cyl_bessel_i1f<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(input,output,M,N);
}


template<typename T> __global__ void gpu_matrix_erfcf(T* x, T* c, int m, int n) 
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (col < n && row < m) 
    {        
        c[row*n + col] = Erfc(x[row*n + col]);
    }
}

template<typename T> void matrix_r_erfcf(T * input, T *output, int M, int N)
{
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_matrix_erfcf<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(input,output,M,N);
}

template<typename T> __global__ void gpu_matrix_erfcinvf(T* x, T* c, int m, int n) 
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (col < n && row < m) 
    {        
        c[row*n + col] = Erfcinv(x[row*n + col]);
    }
}


template<typename T> void matrix_r_erfcinvf(T * input, T *output, int M, int N)
{
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_matrix_erfcinvf<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(input,output,M,N);
}

template<typename T> __global__ void gpu_matrix_erfcxf(T* x, T* c, int m, int n) 
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (col < n && row < m) 
    {        
        c[row*n + col] = Erfcx(x[row*n + col]);
    }
}



template<typename T> void matrix_r_erfcxf(T * input, T *output, int M, int N)
{
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_matrix_erfcxf<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(input,output,M,N);
}


template<typename T> __global__ void gpu_matrix_erff(T* x, T* c, int m, int n) 
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (col < n && row < m) 
    {        
        c[row*n + col] = Erf(x[row*n + col]);
    }
}


template<typename T> 
void matrix_r_erff(T * input, T *output, int M, int N)
{
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_matrix_erff<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(input,output,M,N);
}


template<typename T> __global__ void gpu_matrix_erfinvf(T* x, T* c, int m, int n) 
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (col < n && row < m) 
    {        
        c[row*n + col] = Erfinv(x[row*n + col]);
    }
}

template<typename T> 
void matrix_r_erfinvf(T * input, T * output, int M, int N)
{
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_matrix_erfinvf<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(input,output,M,N);
}

template<typename T> __global__ void gpu_matrix_fdimf(T* x, T * y, T* c, int m, int n) 
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (col < n && row < m) 
    {        
        c[row*n + col] = FDim(x[row*n + col],y[row*n + col]);
    }
}


template<typename T> 
void matrix_r_fdimf(T * x, T * y, T *output, int M, int N)
{
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_matrix_fdimf<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(x,y,output,M,N);
}

template<typename T> __global__ void gpu_matrix_fdividef(T* x, T * y, T* c, int m, int n) 
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (col < n && row < m) 
    {        
        c[row*n + col] = FDivide(x[row*n + col],y[row*n + col]);
    }
}

template<typename T> 
void matrix_r_fdividef(T * x, T *y, T *output, int M, int N)
{
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_matrix_fdividef<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(x,y,output,M,N);
}

template<typename T> __global__ void gpu_matrix_fmaf(T* x, T * y, T *z, T* c, int m, int n) 
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (col < n && row < m) 
    {        
        c[row*n + col] = FMa(x[row*n + col],y[row*n + col],z[row*n + col]);
    }
}


template<typename T>
void matrix_r_fmaf(T * x, T *y, T *z, T *output, int M, int N)
{
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_matrix_fmaf<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(x,y,z,output,M,N);
}

template<typename T> __global__ void gpu_matrix_hypotf(T* x, T * y, T* c, int m, int n) 
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (col < n && row < m) 
    {        
        c[row*n + col] = Hypot(x[row*n + col],y[row*n + col]);
    }
}

template<typename T> void matrix_r_hypotf(T * x, T *y, T *output, int M, int N)
{
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_matrix_hypotf<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(x,y,output,M,N);
}

template<typename T> __global__ void gpu_matrix_ilogbf(T* x, T* c, int m, int n) 
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (col < n && row < m) 
    {        
        c[row*n + col] = ILogb(x[row*n + col]);
    }
}


template<typename T> void matrix_r_ilogbf(T * input, T *output, int M, int N)
{
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_matrix_ilogbf<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(input,output,M,N);
}

template<typename T> __global__ void gpu_matrix_j0f(T* x, T* c, int m, int n) 
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (col < n && row < m) 
    {        
        c[row*n + col] = J0(x[row*n + col]);
    }
}

template<typename T> void matrix_r_j0f(T * input, T *output, int M, int N)
{
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_matrix_j0f<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(input,output,M,N);
}

template<typename T> __global__ void gpu_matrix_j1f(T* x, T* c, int m, int n) 
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (col < n && row < m) 
    {        
        c[row*n + col] = J1(x[row*n + col]);
    }
}


template<typename T> 
void matrix_r_j1f(T * input, T *output, int M, int N)
{
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_matrix_j1f<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(input,output,M,N);
}


template<typename T> __global__ void gpu_matrix_jnf(T* x, int N, T* c, int m, int n) 
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (col < n && row < m) 
    {        
        c[row*n + col] = Jn(x[row*n + col],N);
    }
}


template<typename T> 
void matrix_r_jnf(T * input, T *output, int n, int M, int N)
{
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_matrix_jnf<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(input,n,output,M,N);
}


template<typename T> __global__ void gpu_matrix_ldexpf(T* x, int exp, T* c, int m, int n) 
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (col < n && row < m) 
    {        
        c[row*n + col] = Ldexp(x[row*n + col],exp);
    }
}



template<typename T> 
void matrix_r_ldexpf(T * input, T *output, int exp, int M, int N)
{
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_matrix_ldexpf<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(input,exp,output,M,N);
}

template<typename T> __global__ void gpu_matrix_lgammaf(T* x, T* c, int m, int n) 
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (col < n && row < m) 
    {        
        c[row*n + col] = LGamma(x[row*n + col]);
    }
}


template<typename T> 
void matrix_r_lgammaf(T * input, T *output, int M, int N)
{
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_matrix_lgammaf<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(input,output,M,N);
}

template<typename T> __global__ void gpu_matrix_copysign(T* x, T *y, T* c, int m, int n) 
{    
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    int idx = row*n + col;
    if (col < n && row < m) 
    {        
        c[idx] = Copysign(x[idx],y[idx]);
    }
}


template<typename T> 
void matrix_r_copysignf(T * X, T *Y, T *output, int M, int N)
{    
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_matrix_copysign<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(X,Y,output,M,N);
}


template<typename T> __global__ void gpu_matrix_nearbyintf(T* x, T* c, int m, int n) 
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (col < n && row < m) 
    {        
        c[row*n + col] = Nearbyint(x[row*n + col]);
    }
}

template<typename T> 
void matrix_r_nearbyintf(T * input, T *output, int M, int N)
{
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_matrix_nearbyintf<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(input,output,M,N);
}

template<typename T> __global__ void gpu_matrix_norm3df(T* x, T *y, T *z, T* c, int m, int n) 
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (col < n && row < m) 
    {        
        c[row*n + col] = Norm3d(x[row*n + col],y[row*n + col],z[row*n + col]);
    }
}


template<typename T> 
void matrix_r_norm3df(T * x, T *y, T *z, T *output, int M, int N)
{   
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_matrix_norm3df<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(x,y,z,output,M,N); 
}


template<typename T> __global__ void gpu_matrix_norm4df(T* x, T *y, T *z, T * w, T* c, int m, int n) 
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (col < n && row < m) 
    {        
        c[row*n + col] = Norm4d(x[row*n + col],y[row*n + col],z[row*n + col],w[row*n + col]);
    }
}


template<typename T> 
void matrix_r_norm4df(T * x, T *y, T *z, T * w, T *output,int M, int N)
{
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_matrix_norm4df<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(x,y,z,w,output,M,N);
}

template<typename T> __global__ void gpu_matrix_normcdff(T* x, T* c, int m, int n) 
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (col < n && row < m) 
    {        
        c[row*n + col] = Normcdf(x[row*n + col]);
    }
}



template<typename T> 
void matrix_r_normcdff(T * input, T *output, int M, int N)
{
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_matrix_normcdff<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(input,output,M,N);
}

template<typename T> __global__ void gpu_matrix_normcdfinvf(T* x, T* c, int m, int n) 
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (col < n && row < m) 
    {        
        c[row*n + col] = Normcdfinv(x[row*n + col]);
    }
}


template<typename T> 
void matrix_r_normcdfinvf(T * input, T *output, int M, int N)
{
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_matrix_normcdfinvf<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(input,output,M,N);
}


template<typename T> __global__ void gpu_matrix_normf(int dim, T* x, T* c, int m, int n) 
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (col < n && row < m) 
    {        
        c[row*n + col] = Norm(dim, x);
    }
}


template<typename T> 
void matrix_r_normf(int dim, T *input, T *output, int M, int N)
{
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_matrix_normf<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(dim,input,output,M,N);
}

template<typename T> __global__ void gpu_matrix_remainderf(T* x, T * y, T* c, int m, int n) 
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (col < n && row < m) 
    {        
        c[row*n + col] = Remainder(x[row*n + col],y[row*n + col]);
    }
}



template<typename T> 
void matrix_r_remainderf(T * x, T *y, T *output, int M, int N)
{
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_matrix_remainderf<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(x,y,output,M,N);
}




template<typename T> __global__ void gpu_matrix_rcbrtf(T* x, T* c, int m, int n) 
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (col < n && row < m) 
    {        
        c[row*n + col] = Rcbrt(x[row*n + col]);
    }
}

template<typename T> 
void matrix_r_rcbrtf(T * input, T * output, int M, int N)
{
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_matrix_rcbrtf<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(input,output,M,N);
}



template<typename T> __global__ void gpu_matrix_rhypotf(T* x, T *y, T* c, int m, int n) 
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (col < n && row < m) 
    {        
        c[row*n + col] = RHypot(x[row*n + col],y[row*n + col]);
    }
}


template<typename T> 
void matrix_r_rhypotf(T * x, T *y, T * output, int M, int N)
{
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_matrix_rhypotf<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(x,y,output,M,N);
}

template<typename T> __global__ void gpu_matrix_rnorm3df(T* x, T *y, T *z, T* c, int m, int n) 
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (col < n && row < m) 
    {        
        c[row*n + col] = RNorm3d(x[row*n + col],y[row*n + col],z[row*n + col]);
    }
}




template<typename T> 
void matrix_r_rnorm3df(T * x, T *y, T *z, T * output, int M, int N)
{    
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_matrix_rnorm3df<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(x,y,z,output,M,N); 
}


template<typename T> __global__ void gpu_matrix_rnorm4df(T* x, T *y, T *z, T *w, T* c, int m, int n) 
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (col < n && row < m) 
    {        
        c[row*n + col] = RNorm4d(x[row*n + col],y[row*n + col],z[row*n + col],w[row*n + col]);
    }
}



template<typename T> 
void matrix_r_rnorm4df(T * x, T *y, T *z, T *w, T *output, int M, int N)
{
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_matrix_rnorm4df<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(x,y,z,w,output,M,N);
}


// i dont know if this works this way yet
template<typename T> __global__ void gpu_matrix_rnormf(int dim, T* x, T* c, int m, int n) 
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (col < n && row < m) 
    {        
        c[row*n + col] = RNorm(dim, x);
    }
}

template<typename T> 
void matrix_r_rnormf(int dim, T *input, T *output, int M, int N)
{
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_matrix_rnormf<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(dim,input,output,M,N);    
}


template<typename T> __global__ void gpu_matrix_scalblnf(T* x, long int N, T* c, int m, int n) 
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (col < n && row < m) 
    {        
        c[row*n + col] = Scalbln(x[row*n + col],N);
    }
}       


template<typename T> 
void matrix_r_scalblnf(T * input, long int n, T *output, int M, int N)
{
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_matrix_scalblnf<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(input,n,output,M,N);
}

template<typename T> __global__ void gpu_matrix_sinpif(T* x, T* c, int m, int n) 
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (col < n && row < m) 
    {        
        c[row*n + col] = Sinpi(x[row*n + col]);
    }
}



template<typename T> 
void matrix_r_sinpif(T * input, T *output, int M, int N)
{
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_matrix_sinpif<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(input,output,M,N);
}

template<typename T> __global__ void gpu_matrix_tgammaf(T* x, T* c, int m, int n) 
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (col < n && row < m) 
    {        
        c[row*n + col] = TGamma(x[row*n + col]);
    }
}

template<typename T> 
void matrix_r_tgammaf(T * input, T *output, int M, int N)
{
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_matrix_tgammaf<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(input,output,M,N);
}


template<typename T> __global__ void gpu_matrix_y0f(T* x, T* c, int m, int n) 
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (col < n && row < m) 
    {        
        c[row*n + col] = y0(x[row*n + col]);
    }
}


template<typename T> 
void matrix_r_y0f(T * input, T *output, int M, int N)
{
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_matrix_y0f<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(input,output,M,N);
}


template<typename T> __global__ void gpu_matrix_y1f(T* x, T* c, int m, int n) 
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (col < n && row < m) 
    {        
        c[row*n + col] = y1(x[row*n + col]);
    }
}



template<typename T> 
void matrix_r_y1f(T * input, T *output, int M, int N)
{
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_matrix_y1f<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(input,output,M,N);
}

template<typename T> __global__ void gpu_matrix_ynf(int N,T* x, T* c, int m, int n) 
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (col < n && row < m) 
    {        
        c[row*n + col] = yn(N,x[row*n + col]);
    }
}


template<typename T> 
void matrix_r_ynf(int n, T * input, T *output, int M, int N)
{
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_matrix_ynf<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(n,input,output,M,N);
}

template<typename T> __global__ void matrix_fdimf_const_device(T * x, T y, T * out, int m, int n)
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    int idx = row*n + col; 
    if (col < n && row < m) 
    {        
        out[idx] = FDim(x[idx],y);
    }
}



template<typename T> 
void matrix_r_fdimf_const(T * x, T y, T *output, int M, int N)
{
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    matrix_fdimf_const_device<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(x,y,output,M,N);    
}


template<typename T> __global__ void matrix_fdividef_const_device(T * x, T y, T * out, int m, int n)
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    int idx = row*n + col; 
    if (col < n && row < m) 
    {        
        out[idx] = FDivide(x[idx],y);
    }
}



template<typename T> 
void matrix_r_fdividef_const(T * x, T y, T *output, int M, int N)
{
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    matrix_fdividef_const_device<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(x,y,output,M,N);    
}

template<typename T> __global__ void matrix_hypotf_const_device(T * x, T y, T * out, int m, int n)
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    int idx = row*n + col; 
    if (col < n && row < m) 
    {        
        out[idx] = Hypot(x[idx],y);
    }
}

template<typename T> 
void matrix_r_hypotf_const(T * x, T y, T *output, int M, int N)
{
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    matrix_hypotf_const_device<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(x,y,output,M,N);    
}


template<typename T> __global__ void matrix_remainderf_const_device(T * x, T y, T * out, int m, int n)
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    int idx = row*n + col; 
    if (col < n && row < m) 
    {        
        out[idx] = Remainder(x[idx],y);
    }
}

template<typename T> 
void matrix_r_remainderf_const(T * x, T y, T * output, int M, int N)
{
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    matrix_remainderf_const_device<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(x,y,output,M,N);    
}

template<typename T> __global__ void matrix_rhypotf_const_device(T * x, T y, T * out, int m, int n)
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    int idx = row*n + col; 
    if (col < n && row < m) 
    {        
        out[idx] = Fmax(x[idx],y);
    }
}

template<typename T> 
void matrix_r_rhypotf_const(T * x, T y, T *output, int M, int N)
{
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    matrix_rhypotf_const_device<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(x,y,output,M,N);    
}

template<typename T> __global__ void matrix_truncf_device(T * a, T * out, int m, int n)
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    int idx = row*n + col; 
    if (col < n && row < m) 
    {        
        out[idx] = Trunc(a[idx]);
    }
}

template<typename T> 
void matrix_r_truncf(T * x, T *output, int M, int N)
{
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    matrix_truncf_device<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(x,output,M,N);
}
