////////////////////////////////////////////////////////////////////////
// vector
////////////////////////////////////////////////////////////////////////


#include <cassert>
#include <map>
#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <complex>

#include "cuda_runtime.h"
#include "math_constants.h"

#include "viper_vector.hpp"
#include "viper_ops.h"

template void vector_addf_row(float * x, int row, float * y, int src,size_t n);
template void vector_subf_row(float * x, int row, float * y, int src,size_t n);
template void vector_mulf_row(float * x, int row, float * y, int src,size_t n);
template void vector_divf_row(float * x, int row, float * y, int src,size_t n);
template void vector_modf_row(float * x, int row, float * y, int src,size_t n);
template void vector_r_truncf(float * x, float *output, int n);
template void vector_r_copysignf(float * X, float *Y, float *output, int n);
template void vector_setrowf(float * dst, int row, float * sc, int row_src, size_t n);
template void vector_copyrowf(float * dst, int row, float * src, int row_src, int n);
template void vector_r_addf(float * x, float * y, float * output, int n);
template void vector_r_subf(float * x, float * y, float * output, int n);
template void vector_r_mulf(float * x, float * y, float * output, int n);
template void vector_r_divf(float * x, float * y, float * output, int n);
template void vector_r_modf(float * x, float * y, float * output, int n);
template void vector_r_acosf(float * devPtr, float * outputs, int n);
template void vector_r_asinf(float * devPtr, float * outputs, int n);
template void vector_r_atanf(float * devPtr, float * outputs, int n);
template void vector_r_atan2f(float * a, float * b, float * output, int n);
template void vector_r_acoshf(float * devPtr, float * outputs, int n);
template void vector_r_asinhf(float * devPtr, float * outputs, int n);
template void vector_r_atanhf(float * devPtr, float * outputs, int n);
template void vector_r_cosf(float * devPtr, float * outputs, int n);
template void vector_r_sinf(float * devPtr, float * outputs, int n);
template void vector_r_tanf(float * devPtr, float * outputs, int n);
template void vector_r_coshf(float * devPtr, float * outputs, int n);
template void vector_r_sinhf(float * devPtr, float * outputs, int n);
template void vector_r_tanhf(float * devPtr, float * outputs, int n);
template void vector_r_ceilf(float * devPtr, float * output, int n);
template void vector_r_exp10f(float * devPtr, float * outputs, int n);
template void vector_r_exp2f(float * devPtr, float * output, int n);
template void vector_r_expf(float * devPtr, float * output, int n);
template void vector_r_expm1f(float * devPtr, float * output, int n);
template void vector_r_fabsf(float * devPtr, float * output, int n);
template void vector_r_floorf(float * devPtr, float * output, int n);
template void vector_r_fmaxf(float * x, float * y, float * output, int n);
template void vector_r_fminf(float * x, float * y, float * output, int n);
template void vector_r_fmodf(float * x, float * y, float * output, int n);
template void vector_r_log10f(float * x, float *output, int n);
template void vector_r_log1pf(float * x, float *output, int n);
template void vector_r_log2f(float * x, float * output, int n);
template void vector_r_logbf(float * x, float *output, int n);
template void vector_r_powf(float * x, float * y, float *output, int n);
template void vector_r_rsqrtf(float * x, float * output, int n);
template void vector_r_sqrtf(float * x, float *output, int n);
template void vector_r_softmaxf(float * x, float *output, int n);
template void vector_r_sigmoidf(float * x, float *output, int n);
template void vector_r_sigmoid_gradf(float * x, float *output, int n);
template void vector_r_tanh_gradf(float * x, float *output, int n);
template void vector_r_reluf(float * x, float *output, int n);
template void vector_r_relu_gradf(float * x, float *output, int n);
template void vector_r_cbrtf(float * devPtr, float * output, int n);
template void vector_r_cospif(float * devPtr, float * output, int n);
template void vector_r_cyl_bessel_i0f(float * devPtr, float * output, int n);
template void vector_r_cyl_bessel_i1f(float * devPtr, float * output, int n);
template void vector_r_erfcf(float * devPtr, float * output, int n);
template void vector_r_erfcinvf(float * devPtr, float * output, int n);
template void vector_r_erfcxf(float * devPtr, float * output, int n);
template void vector_r_erff(float * devPtr, float * output, int n);
template void vector_r_erfinvf(float * devPtr, float * output, int n);
template void vector_r_fdimf(float * a, float * b, float * output, int n);
template void vector_r_fdividef(float * a, float * b, float * output, int n);
template void vector_r_fmaf(float * x, float * y, float * z, float *output, int n);
template void vector_r_hypotf(float * x, float * y, float * output, int n);
template void vector_r_ilogbf(float * x, float *output, int n);
template void vector_r_j0f(float * x, float *output, int n);
template void vector_r_j1f(float * x, float *output, int n);
template void vector_r_jnf(float * x, float * output, int M, int n);
template void vector_r_ldexpf(float * x, float * output, int exp, int n);
template void vector_r_lgammaf(float * x, float *output, int n);
template void vector_r_nearbyintf(float * x, float *output, int n);
template void vector_r_norm3df(float * x, float * y, float * z, float * output, int n);
template void vector_r_norm4df(float * x, float * y, float * z, float * q, float * output, int n);
template void vector_r_normcdff(float * x, float * output, int n);
template void vector_r_normcdfinvf(float * x, float *output, int n);
template void vector_r_normf(int dim, float * x, float * output, int n);
template void vector_r_rcbrtf(float * x, float *output, int n);
template void vector_r_remainderf(float * x, float * y, float *output, int n);
template void vector_r_rhypotf(float * x, float * y, float *output, int n);
template void vector_r_rnorm3df(float * x, float * y, float * z, float * output, int n);
template void vector_r_rnorm4df(float * x, float * y, float * z, float * q, float *output, int n);
template void vector_r_rnormf(int dim, float * x, float *output, int n);
template void vector_r_scalblnf(float * x, long int M, float * output, int n);
template void vector_r_tgammaf(float * x, float * output, int n);
template void vector_r_y0f(float * x, float *output, int n);
template void vector_r_y1f(float * x, float * output, int n);
template void vector_r_ynf(int M, float * x, float *output, int n);
template void vector_r_sinpif(float * x, float *output, int n);
template void vector_r_addf_const(float * x, float  y, float *output, int n);
template void vector_r_subf_const(float * x, float  y, float *output,int n);
template void vector_r_mulf_const(float * x, float  y, float *output,int n);
template void vector_r_divf_const(float * x, float  y, float *output,int n);
template void vector_r_modf_const(float * x, float  y, float *output,int n);
template void vector_r_atan2f_const(float * a, float b, float *output,int n);
template void vector_r_fmaxf_const(float * x, float  y, float *output,int n);
template void vector_r_fminf_const(float * x, float  y, float *output,int n);
template void vector_r_fmodf_const(float * x, float  y, float *output,int n);
template void vector_r_powf_const(float * x, float y, float *output,int n);
template void vector_r_addf_scalar(float * x, float * y, float *output,int n);
template void vector_r_subf_scalar(float * x, float * y, float *output,int n);
template void vector_r_mulf_scalar(float * x, float * y, float *output,int n);
template void vector_r_divf_scalar(float * x, float * y, float *output,int n);
template void vector_r_modf_scalar(float * x, float * y, float *output,int n);
template void vector_r_atan2f_scalar(float * a, float *b, float *output,int n);
template void vector_r_fmaxf_scalar(float * x, float  *y, float *output,int n);
template void vector_r_fminf_scalar(float * x, float  *y, float *output,int n);
template void vector_r_fmodf_scalar(float * x, float  *y, float *output,int n);
template void vector_r_powf_scalar(float * x, float *y, float *output,int n);
template void vector_r_fdimf_const(float * a, float  b, float *output,int n);
template void vector_r_fdividef_const(float * a, float  b, float *output,int n);
template void vector_r_hypotf_const(float * x, float  y, float *output,int n);
template void vector_r_remainderf_const(float * x, float y, float *output,int n);
template void vector_r_rhypotf_const(float * x, float y, float *output,int n);
template void vector_r_fdimf_scalar(float * a, float  *b, float *output,int n);
template void vector_r_fdividef_scalar(float * a, float *b, float *output,int n);
template void vector_r_hypotf_scalar(float * x, float  *y, float *output,int n);
template void vector_r_remainderf_scalar(float * x, float *y, float *output,int n);
template void vector_r_rhypotf_scalar(float * x, float *y, float *output,int n);
template float vector_sumf(float * in, int size);
template float vector_prodf(float * in, int size);

template void vector_addf_row(double * x, int row, double * y, int src,size_t n);
template void vector_subf_row(double * x, int row, double * y, int src,size_t n);
template void vector_mulf_row(double * x, int row, double * y, int src,size_t n);
template void vector_divf_row(double * x, int row, double * y, int src,size_t n);
template void vector_modf_row(double * x, int row, double * y, int src,size_t n);
template void vector_r_truncf(double * x, double *output, int n);
template void vector_r_copysignf(double * X, double *Y, double *output, int n);
template void vector_setrowf(double * dst, int row, double * sc, int row_src, size_t n);
template void vector_copyrowf(double * dst, int row, double * src, int row_src, int n);
template void vector_r_addf(double * x, double * y, double * output, int n);
template void vector_r_subf(double * x, double * y, double * output, int n);
template void vector_r_mulf(double * x, double * y, double * output, int n);
template void vector_r_divf(double * x, double * y, double * output, int n);
template void vector_r_modf(double * x, double * y, double * output, int n);
template void vector_r_acosf(double * devPtr, double * outputs, int n);
template void vector_r_asinf(double * devPtr, double * outputs, int n);
template void vector_r_atanf(double * devPtr, double * outputs, int n);
template void vector_r_atan2f(double * a, double * b, double * output, int n);
template void vector_r_acoshf(double * devPtr, double * outputs, int n);
template void vector_r_asinhf(double * devPtr, double * outputs, int n);
template void vector_r_atanhf(double * devPtr, double * outputs, int n);
template void vector_r_cosf(double * devPtr, double * outputs, int n);
template void vector_r_sinf(double * devPtr, double * outputs, int n);
template void vector_r_tanf(double * devPtr, double * outputs, int n);
template void vector_r_coshf(double * devPtr, double * outputs, int n);
template void vector_r_sinhf(double * devPtr, double * outputs, int n);
template void vector_r_tanhf(double * devPtr, double * outputs, int n);
template void vector_r_ceilf(double * devPtr, double * output, int n);
template void vector_r_exp10f(double * devPtr, double * outputs, int n);
template void vector_r_exp2f(double * devPtr, double * output, int n);
template void vector_r_expf(double * devPtr, double * output, int n);
template void vector_r_expm1f(double * devPtr, double * output, int n);
template void vector_r_fabsf(double * devPtr, double * output, int n);
template void vector_r_floorf(double * devPtr, double * output, int n);
template void vector_r_fmaxf(double * x, double * y, double * output, int n);
template void vector_r_fminf(double * x, double * y, double * output, int n);
template void vector_r_fmodf(double * x, double * y, double * output, int n);
template void vector_r_log10f(double * x, double *output, int n);
template void vector_r_log1pf(double * x, double *output, int n);
template void vector_r_log2f(double * x, double * output, int n);
template void vector_r_logbf(double * x, double *output, int n);
template void vector_r_powf(double * x, double * y, double *output, int n);
template void vector_r_rsqrtf(double * x, double * output, int n);
template void vector_r_sqrtf(double * x, double *output, int n);
template void vector_r_softmaxf(double * x, double *output, int n);
template void vector_r_sigmoidf(double * x, double *output, int n);
template void vector_r_sigmoid_gradf(double * x, double *output, int n);
template void vector_r_tanh_gradf(double * x, double *output, int n);
template void vector_r_reluf(double * x, double *output, int n);
template void vector_r_relu_gradf(double * x, double *output, int n);
template void vector_r_cbrtf(double * devPtr, double * output, int n);
template void vector_r_cospif(double * devPtr, double * output, int n);
template void vector_r_cyl_bessel_i0f(double * devPtr, double * output, int n);
template void vector_r_cyl_bessel_i1f(double * devPtr, double * output, int n);
template void vector_r_erfcf(double * devPtr, double * output, int n);
template void vector_r_erfcinvf(double * devPtr, double * output, int n);
template void vector_r_erfcxf(double * devPtr, double * output, int n);
template void vector_r_erff(double * devPtr, double * output, int n);
template void vector_r_erfinvf(double * devPtr, double * output, int n);
template void vector_r_fdimf(double * a, double * b, double * output, int n);
template void vector_r_fdividef(double * a, double * b, double * output, int n);
template void vector_r_fmaf(double * x, double * y, double * z, double *output, int n);
template void vector_r_hypotf(double * x, double * y, double * output, int n);
template void vector_r_ilogbf(double * x, double *output, int n);
template void vector_r_j0f(double * x, double *output, int n);
template void vector_r_j1f(double * x, double *output, int n);
template void vector_r_jnf(double * x, double * output, int M, int n);
template void vector_r_ldexpf(double * x, double * output, int exp, int n);
template void vector_r_lgammaf(double * x, double *output, int n);
template void vector_r_nearbyintf(double * x, double *output, int n);
template void vector_r_norm3df(double * x, double * y, double * z, double * output, int n);
template void vector_r_norm4df(double * x, double * y, double * z, double * q, double * output, int n);
template void vector_r_normcdff(double * x, double * output, int n);
template void vector_r_normcdfinvf(double * x, double *output, int n);
template void vector_r_normf(int dim, double * x, double * output, int n);
template void vector_r_rcbrtf(double * x, double *output, int n);
template void vector_r_remainderf(double * x, double * y, double *output, int n);
template void vector_r_rhypotf(double * x, double * y, double *output, int n);
template void vector_r_rnorm3df(double * x, double * y, double * z, double * output, int n);
template void vector_r_rnorm4df(double * x, double * y, double * z, double * q, double *output, int n);
template void vector_r_rnormf(int dim, double * x, double *output, int n);
template void vector_r_scalblnf(double * x, long int M, double * output, int n);
template void vector_r_tgammaf(double * x, double * output, int n);
template void vector_r_y0f(double * x, double *output, int n);
template void vector_r_y1f(double * x, double * output, int n);
template void vector_r_ynf(int M, double * x, double *output, int n);
template void vector_r_sinpif(double * x, double *output, int n);
template void vector_r_addf_const(double * x, double  y, double *output, int n);
template void vector_r_subf_const(double * x, double  y, double *output,int n);
template void vector_r_mulf_const(double * x, double  y, double *output,int n);
template void vector_r_divf_const(double * x, double  y, double *output,int n);
template void vector_r_modf_const(double * x, double  y, double *output,int n);
template void vector_r_atan2f_const(double * a, double b, double *output,int n);
template void vector_r_fmaxf_const(double * x, double  y, double *output,int n);
template void vector_r_fminf_const(double * x, double  y, double *output,int n);
template void vector_r_fmodf_const(double * x, double  y, double *output,int n);
template void vector_r_powf_const(double * x, double y, double *output,int n);
template void vector_r_addf_scalar(double * x, double * y, double *output,int n);
template void vector_r_subf_scalar(double * x, double * y, double *output,int n);
template void vector_r_mulf_scalar(double * x, double * y, double *output,int n);
template void vector_r_divf_scalar(double * x, double * y, double *output,int n);
template void vector_r_modf_scalar(double * x, double * y, double *output,int n);
template void vector_r_atan2f_scalar(double * a, double *b, double *output,int n);
template void vector_r_fmaxf_scalar(double * x, double  *y, double *output,int n);
template void vector_r_fminf_scalar(double * x, double  *y, double *output,int n);
template void vector_r_fmodf_scalar(double * x, double  *y, double *output,int n);
template void vector_r_powf_scalar(double * x, double *y, double *output,int n);
template void vector_r_fdimf_const(double * a, double  b, double *output,int n);
template void vector_r_fdividef_const(double * a, double  b, double *output,int n);
template void vector_r_hypotf_const(double * x, double  y, double *output,int n);
template void vector_r_remainderf_const(double * x, double y, double *output,int n);
template void vector_r_rhypotf_const(double * x, double y, double *output,int n);
template void vector_r_fdimf_scalar(double * a, double  *b, double *output,int n);
template void vector_r_fdividef_scalar(double * a, double *b, double *output,int n);
template void vector_r_hypotf_scalar(double * x, double  *y, double *output,int n);
template void vector_r_remainderf_scalar(double * x, double *y, double *output,int n);
template void vector_r_rhypotf_scalar(double * x, double *y, double *output,int n);
template double vector_sumf(double * in, int size);
template double vector_prodf(double * in, int size);

#define BLOCK_SIZE 16

#include <pthread.h>

extern pthread_mutex_t mutex;

cudaStream_t get_cuda_stream();
cudaStream_t random_stream();

template void    return_memory(int length, float *fp);
template float*  find_memory(int length);
template void    add_memory(int length, float * ptr);
template void cuda_zero(float * dst, int n);
template void cuda_memcpy(float * dst, float * src, int n);

template void    return_memory(int length, double *fp);
template double*  find_memory(int length);
template void    add_memory(int length, double * ptr);
template void cuda_zero(double * dst, int n);
template void cuda_memcpy(double * dst, double * src, int n);

template void    return_memory(int length, complex<float> *fp);
template void    return_memory(int length, complex<double> *fp);
template complex<float>*  find_memory(int length);
template complex<double>*  find_memory(int length);
template void    add_memory(int length, complex<float> * ptr);
template void    add_memory(int length, complex<double> * ptr);

template void    return_memory(int length, std::complex<float> *fp);
template void    return_memory(int length, std::complex<double> *fp);
template std::complex<float>*  find_memory(int length);
template std::complex<double>*  find_memory(int length);
template void    add_memory(int length, std::complex<float> * ptr);
template void    add_memory(int length, std::complex<double> * ptr);

// memory cache
std::multimap<int,void*> cuda_memory;

template<typename T>
void add_memory(int length, T * f) {
    pthread_mutex_lock(&mutex);
    cuda_memory.insert(std::pair<int,T*>(length*sizeof(T),f));
    pthread_mutex_unlock(&mutex);
}

template<typename T>
void return_memory(int length, T *fp) {        
    pthread_mutex_lock(&mutex);
    cuda_memory.insert(std::pair<int,T*>(length*sizeof(T),fp));    
    pthread_mutex_unlock(&mutex);
}

template<typename T>
T* find_memory(int length) {   
    pthread_mutex_lock(&mutex);                
    typename std::multimap<int,void*>::iterator i = cuda_memory.find(length*sizeof(T));
    if(i == cuda_memory.end()) 
    {
        pthread_mutex_unlock(&mutex);
        return nullptr;                
    }
    //cuda_zero((T*)i->second,length);
    cuda_memory.erase(i);
    pthread_mutex_unlock(&mutex);
    return static_cast<T*>(i->second);
}

void clear_cache() {
    typename std::multimap<int,void*>::iterator i = cuda_memory.begin();
    while(i != cuda_memory.end()) {
        cudaFree(i->second);       
        i++; 
    }
    cuda_memory.clear();    
}




// obviously this is crap but reduction on GPU is a giant pain
template<typename T>
T vector_sumf(T * in, int size)
{
    T r = 0;
    for(size_t i = 0; i < size; i++) r += in[i];
    return r;
}



template<typename T>
T vector_prodf(T * in, int size)
{
    T r = 0;
    for(size_t i = 0; i < size; i++) r *= in[i];
    return r;
}


template<typename T> __global__ void vector_dummy(T * x, T * out, int N) {
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) out[idx] =x[idx];
}


void calcSize(int N,int * gridSize, int * blockSize) {
    
    //int minGridSize = 0;    
    //cudaOccupancyMaxPotentialBlockSize(&minGridSize, blockSize, vector_dummy, 0, N); 
    //*gridSize = (N + *blockSize - 1) / *blockSize; 
    
    
    *blockSize=BLOCK_SIZE;    
    *gridSize=(N+*blockSize)/ *blockSize;    
}


template<typename T> __global__ void cuda_memcpy_device(T * dst, T * src, size_t N) {
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) dst[idx] = src[idx];
}

template<typename T>
void cuda_memcpy(T * dst, T * src, int n)
{    
    //int gridSize,blockSize;
    //calcSize(n,&gridSize,&blockSize);        
    //cuda_memcpy_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(dst,src,n);        
    //cudaDeviceSynchronize();
    cudaMemcpyAsync(dst,src,n*sizeof(T),cudaMemcpyDeviceToDevice, get_cuda_stream());
}



template<typename T> __global__ void cuda_zero_device(T * dst, size_t N) {
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) dst[idx] = 0;
}


template<typename T>
void cuda_zero(T * dst, int n)
{    
    /*
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    cuda_zero_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(dst,n);        
    */    
    cudaMemsetAsync(dst,0,n*sizeof(T), get_cuda_stream());    
}




template<typename T> __global__ void vector_addf_device(T * x, T * y, T * out, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) out[idx] = x[idx] + y[idx];
}

template<typename T> void vector_r_addf(T * x, T * y, T * output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);    
    vector_addf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);    
}




template<typename T> __global__ void vector_subf_device(T * x, T * y, T * out, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) out[idx] = x[idx] - y[idx];
}


template<typename T> void vector_r_subf(T * x, T * y, T * output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);    
    vector_subf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);    
}




template<typename T> __global__ void vector_mulf_device(T * x, T * y, T * out, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) out[idx] = x[idx] * y[idx];
}


template<typename T> void vector_r_mulf(T * x, T * y, T * output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);    
    vector_mulf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);    
}


template<typename T> __global__ void vector_divf_device(T * x, T * y, T * out, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) out[idx] = x[idx] / y[idx];
}

template<typename T> void vector_r_divf(T * x, T * y, T * output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);    
    vector_divf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);    
}


template<typename T> __global__ void vector_modf_device(T * x, T * y, T * out, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) out[idx] = fmodf(x[idx],y[idx]);
}

template<typename T> void vector_r_modf(T * x, T * y, T * output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_modf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);    
}

template<typename T> __global__ void vector_acosf_device(T * in, T * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n) out[idx] = ACos(in[idx]);
}


template<typename T> void vector_r_acosf(T * devPtr, T * output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);    
    vector_acosf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n);        
}


template<typename T> __global__ void vector_acoshf_device(T * in, T * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n) out[idx] = ACosh(in[idx]);
}


template<typename T> void vector_r_acoshf(T * devPtr, T * output, int n)
{    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);    
    vector_acoshf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n); 
}




template<typename T> __global__ void vector_asinhf_device(T * in, T * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = ASinh(in[idx]);
}


template<typename T> void vector_r_asinhf(T * devPtr, T * output, int n)
{    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);    
    vector_asinhf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n); 
}



template<typename T> __global__ void vector_asinf_device(T * in, T * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = asinf(in[idx]);
}

template<typename T> void vector_r_asinf(T * devPtr, T * output, int n)
{    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);    
    vector_asinf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n);
}



template<typename T> __global__ void vector_atan2f_device(T * a, T * b, T * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = ATan2(a[idx],b[idx]);
}

template<typename T> void vector_r_atan2f(T * a, T * b, T * output, int n)
{    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);    
    vector_atan2f_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(a,b,output,n); 
}


template<typename T> __global__ void vector_atanf_device(T * in, T * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = atanf(in[idx]);
}

template<typename T> void vector_r_atanf(T * devPtr, T * output, int n)
{    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);    
    vector_atanf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n);    
}



template<typename T> __global__ void vector_atanhf_device(T * in, T * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = ATanh(in[idx]);
}


template<typename T> void vector_r_atanhf(T * devPtr, T * output, int n)
{    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);    
    vector_atanhf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n); 
}



template<typename T> __global__ void vector_ceilf_device(T * in, T * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = Ceil(in[idx]);
}

template<typename T> void vector_r_ceilf(T * devPtr, T * output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);    
    vector_ceilf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n);
}



template<typename T> __global__ void vector_cosf_device(T * in, T * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = Cos(in[idx]);
}


template<typename T> void vector_r_cosf(T * devPtr, T * output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_cosf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n);    
}



template<typename T> __global__ void vector_coshf_device(T * in, T * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = Cosh(in[idx]);
}

template<typename T> void vector_r_coshf(T * devPtr, T * output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_coshf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n);
}



template<typename T> __global__ void vector_exp10f_device(T * in, T * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = Exp10(in[idx]);
}

template<typename T> void vector_r_exp10f(T * devPtr, T * output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_exp10f_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n);
}



template<typename T> __global__ void vector_exp2f_device(T * in, T * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = Exp2(in[idx]);
}


template<typename T> void vector_r_exp2f(T * devPtr, T * output, int n)
{   
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_exp2f_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n); 
}




template<typename T> __global__ void vector_expf_device(T * in, T * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = Exp(in[idx]);
}


template<typename T> void vector_r_expf(T * devPtr, T * output, int n)
{   
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_expf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n); 
}


template<typename T> __global__ void vector_expm1f_device(T * in, T * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = Expm1(in[idx]);
}

template<typename T> void vector_r_expm1f(T * devPtr, T * output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_expm1f_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n);
}

template<typename T> __global__ void vector_fabsf_device(T * in, T * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = Fabs(in[idx]);
}

template<typename T> void vector_r_fabsf(T * devPtr, T * output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_fabsf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n);
}

template<typename T> __global__ void vector_floorf_device(T * a,T * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = Floor(a[idx]);
}

template<typename T> void vector_r_floorf(T * devPtr, T * output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_floorf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n);    
}

template<typename T> __global__ void vector_fmaxf_device(T * a,T * b, T * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = Fmax(a[idx],b[idx]);
}

template<typename T> void vector_r_fmaxf(T * x, T * y, T * output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_fmaxf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);
}
template<typename T> __global__ void vector_fminf_device(T * a,T * b, T * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = Fmin(a[idx],b[idx]);
}


template<typename T> void vector_r_fminf(T * x, T * y, T * output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_fminf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);
}




template<typename T> __global__ void vector_fmodf_device(T * a,T * b, T * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = Fmod(a[idx],b[idx]);
}



template<typename T> void vector_r_fmodf(T * x, T * y, T * output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_fmodf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);
}




template<typename T> __global__ void vector_log10f_device(T * a, T * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = Log10(a[idx]);
}

template<typename T> void vector_r_log10f(T * x, T *output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_log10f_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);
}

template<typename T> __global__ void vector_log1pf_device(T * a, T * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = Log1p(a[idx]);
}


template<typename T> void vector_r_log1pf(T * x, T *output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_log1pf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);
}



template<typename T> __global__ void vector_log2f_device(T * a, T * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = Log2(a[idx]);
}



template<typename T> void vector_r_log2f(T * x, T * output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_log2f_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);
}



template<typename T> __global__ void vector_logbf_device(T * a, T * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = Logb(a[idx]);
}


template<typename T> void vector_r_logbf(T * x, T *output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_logbf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);
}



template<typename T> __global__ void vector_powf_device(T * a, T * b, T * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
    {
        if(a[idx] == 0) out[idx]=0;
        else if(a[idx]==1) out[idx]=1;
        else if(a[idx] < 0) out[idx]=-Pow(Fabs(a[idx]),b[idx]);
        else out[idx] = Pow(a[idx],b[idx]);
    }
}

template<typename T> void vector_r_powf(T * x, T * y, T *output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_powf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);
}



template<typename T> __global__ void vector_rsqrtf_device(T * a, T * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        if(a[idx] <= 0) out[idx] = 0;
        else out[idx] = RSqrt(a[idx]);
}

template<typename T> void vector_r_rsqrtf(T * x, T * output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_rsqrtf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);
}




template<typename T> __global__ void vector_sinf_device(T * a, T * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = Sin(a[idx]);
}


template<typename T> void vector_r_sinf(T * x, T *output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_sinf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);
}




template<typename T> __global__ void vector_sinhf_device(T * a, T * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = Sinh(a[idx]);
}


template<typename T> void vector_r_sinhf(T * x, T *output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_sinhf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);    
}
template<typename T> __global__ void vector_sqrtf_device(T * a, T * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
    {
        if(a[idx] <= 0) out[idx] = 0;
        else out[idx] = Sqrt(a[idx]);
    }
}

template<typename T> void vector_r_sqrtf(T * x, T *output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_sqrtf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);
}



template<typename T> __global__ void vector_tanf_device(T * a, T * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = Tan(a[idx]);
}


template<typename T> void vector_r_tanf(T * x, T *output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_tanf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);
}



template<typename T> __global__ void vector_tanhf_device(T * a, T * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = Tanh(a[idx]);
}

template<typename T> void vector_r_tanhf(T * x, T *output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_tanhf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);
}

template<typename T> __global__ void vector_softmax_device(T * x,T *out, T sum, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) out[idx] = Exp(x[idx]) / sum;
}

template<typename T> void vector_r_softmaxf(T * x, T *output, int n)
{
    T sum = vector_sumf(x,n);
    assert(sum != 0.0);
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);
    vector_softmax_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,sum,n);            
}



template<typename T> __global__ void vector_sigmoid_device(T * x, T * out, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) out[idx] = 1.0 / (1.0 + Exp(-x[idx]));
}

template<typename T> void vector_r_sigmoidf(T * x, T *output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);            
    vector_sigmoid_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);        
}

template<typename T> __global__ void vector_sigmoid_grad_device(T * x, T * out, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) out[idx] = x[idx] * (1.0 - x[idx]);
}

template<typename T> void vector_r_sigmoid_gradf(T * x, T *output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_sigmoid_grad_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);        
}




template<typename T> __global__ void vector_tanh_grad_device(T * x, T * out, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) out[idx] = 1.0 - (x[idx]*x[idx]);
}


template<typename T> void vector_r_tanh_gradf(T * x, T *output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_tanh_grad_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);        
}



template<typename T> __global__ void vector_relu_device(T * x, T * out, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) 
    {
        if(x[idx] < 0) out[idx] = 0.0f;        
        else out[idx] = x[idx]; 
    }
}

template<typename T> void vector_r_reluf(T * x, T *output, int n)
{    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_relu_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);        
}



template<typename T> __global__ void vector_relu_grad_device(T * x, T * out, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) 
    {
        if(x[idx] > 0) out[idx] = 1.0;        
        else out[idx] = 0.0f;
    }
}

template<typename T> void vector_r_relu_gradf(T * x, T *output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_relu_grad_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);        
}

template<typename T> __global__ void vector_add_const_device(T * x, T y, T * out, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) out[idx] = x[idx] + y;
}

template<typename T> void vector_r_addf_const(T * x, T y, T * output, int n)
{    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_add_const_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);     
}



template<typename T> __global__ void vector_sub_const_device(T * x, T y, T * out, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) out[idx] = x[idx] - y;
}


template<typename T> void vector_r_subf_const(T * x, T y, T *output, int n)
{  
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_sub_const_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);          
}



template<typename T> __global__ void vector_mul_const_device(T * x, T y, T * out, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) out[idx] = x[idx] * y;
}

template<typename T> void vector_r_mulf_const(T * x, T y, T * output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_mul_const_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);        
}

template<typename T> __global__ void vector_div_const_device(T * x, T y, T * out, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) out[idx] = x[idx] / y;
}

template<typename T> void vector_r_divf_const(T * x, T y, T *output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_div_const_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);        
}


template<typename T> __global__ void vector_mod_const_device(T * x, T y, T * out, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) out[idx] = fmodf(x[idx],y);
}

template<typename T> void vector_r_modf_const(T * x, T y, T *output, int n)
{    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_mod_const_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);        
}



template<typename T> __global__ void vector_atan2f_const_device(T * a, T b, T * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = ATan2(a[idx],b);
}

template<typename T> void vector_r_atan2f_const(T * a, T  b, T * output, int n)
{   
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_atan2f_const_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(a,b,output,n);    
}





template<typename T> __global__ void vector_fmaxf_const_device(T * a,T b, T * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = Fmax(a[idx],b);
}

template<typename T> void vector_r_fmaxf_const(T * x, T y, T *p, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_fmaxf_const_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,p,n);        
}


template<typename T> __global__ void vector_fminf_const_device(T * a,T b, T * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = Fmin(a[idx],b);
}

template<typename T> void vector_r_fminf_const(T * x, T y, T *p, int n)
{    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_fminf_const_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,p,n);        
}

template<typename T> __global__ void vector_fmodf_const_device(T * a,T b, T * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = Fmod(a[idx],b);
}

template<typename T> void vector_r_fmodf_const(T * x, T y, T * p, int n)
{    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_fmodf_const_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,p,n);       
}



template<typename T> __global__ void vector_powf_const_device(T * a, T b, T * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
    {
        if(a[idx] == 0) out[idx]=0;
        else if(a[idx]==1) out[idx]=1;
        else if(a[idx] < 0) out[idx]=-Pow(Fabs(a[idx]),b);
        else out[idx] = Pow(a[idx],b);
    }
}

template<typename T> void vector_r_powf_const(T * x, T y, T *p, int n)
{    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_powf_const_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,p,n);    
}






/////////////////////////////////
// const/scalar
/////////////////////////////////
template<typename T> __global__ void vector_add_scalar_device(T * x, T * y, T * out, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) out[idx] = x[idx] + y[0];
}
template<typename T> void vector_r_addf_scalar(T * x, T * y, T *output, int n)
{      
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_add_scalar_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);            
}

template<typename T> __global__ void vector_sub_scalar_device(T * x, T * y, T * out, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) out[idx] = x[idx] - y[0];
}

template<typename T> void vector_r_subf_scalar(T * x, T * y, T *output, int n)
{    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_sub_scalar_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);        
}


template<typename T> __global__ void vector_mul_scalar_device(T * x, T * y, T * out, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) out[idx] = x[idx] * y[0];
}
template<typename T> void vector_r_mulf_scalar(T * x, T * y, T * output, int n)
{        
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_mul_scalar_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);            
}


template<typename T> __global__ void vector_div_scalar_device(T * x, T * y, T * out, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N && y[0] != 0.0f) out[idx] = x[idx] / y[0];
}

template<typename T> void vector_r_divf_scalar(T * x, T * y, T *output, int n)
{    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_div_scalar_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);        
}



template<typename T> __global__ void vector_mod_scalar_device(T * x, T * y, T * out, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) out[idx] = Fmod(x[idx],y[0]);
}

template<typename T> void vector_r_modf_scalar(T * x, T * y, T *output, int n)
{    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_mod_scalar_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);        
}
template<typename T> void vector_r_fmodf_scalar(T * x, T * y, T *output, int n)
{    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_mod_scalar_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);        
}



template<typename T> __global__ void vector_fmaxf_scalar_device(T * x, T * y, T * out, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) out[idx] = Fmax(x[idx],y[0]);
}

template<typename T> void vector_r_fmaxf_scalar(T * x, T * y, T *output, int n)
{    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_fmaxf_scalar_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);        
}


template<typename T> __global__ void vector_fminf_scalar_device(T * x, T * y, T * out, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) out[idx] = Fmin(x[idx],y[0]);
}

template<typename T> void vector_r_fminf_scalar(T * x, T * y, T *output, int n)
{    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_fminf_scalar_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);        
}

template<typename T> __global__ void vector_powf_scalar_device(T * x, T * y, T * out, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) 
    {
        if(x[idx] == 0) out[idx]=0;
        else if(x[idx]==1) out[idx]=1;
        else if(x[idx] < 0) out[idx]=-Pow(Fabs(x[idx]),y[0]);
        out[idx] = Pow(x[idx],y[0]);
    }
}

template<typename T> void vector_r_powf_scalar(T * x, T * y, T *output, int n)
{    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_powf_scalar_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);        
}

template<typename T> __global__ void vector_atan2f_scalar_device(T * x, T * y, T * out, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) out[idx] = ATan2(x[idx],y[0]);
}

template<typename T> void vector_r_atan2f_scalar(T * x, T * y, T *output, int n)
{    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_atan2f_scalar_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);        
}

template<typename T> __global__ void vector_fdimf_scalar_device(T * x, T * y, T * out, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) out[idx] = FDim(x[idx],y[0]);
}

template<typename T> void vector_r_fdimf_scalar(T * x, T * y, T *output, int n)
{    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_fdimf_scalar_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);        
}

template<typename T> __global__ void vector_fdividef_scalar_device(T * x, T * y, T * out, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N && y[0] != 0.0) out[idx] = fdividef(x[idx],y[0]);
}

template<typename T> void vector_r_fdividef_scalar(T * x, T * y, T *output, int n)
{    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_fdividef_scalar_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);        
}


template<typename T> __global__ void vector_remainderf_scalar_device(T * x, T * y, T * out, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) out[idx] = Remainder(x[idx],y[0]);
}

template<typename T> void vector_r_remainderf_scalar(T * x, T * y, T *output, int n)
{    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_remainderf_scalar_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);        
}


template<typename T> __global__ void vector_hypot_scalar_device(T * x, T * y, T * out, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) out[idx] = Hypot(x[idx],y[0]);
}

template<typename T> void vector_r_hypotf_scalar(T * x, T * y, T *output, int n)
{    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_hypot_scalar_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);        
}

template<typename T> __global__ void vector_rhypot_scalar_device(T * x, T * y, T * out, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) out[idx] = RHypot(x[idx],y[0]);
}

template<typename T> void vector_r_rhypotf_scalar(T * x, T * y, T *output, int n)
{    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_rhypot_scalar_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);        
}





template<typename T> __global__ void vector_setrowf_device(T * dst, T * src, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) dst[idx] = src[idx];
}

template<typename T> void vector_setrowf(T * dst, int dst_row, T * src, int row_src, size_t n)
{    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);    
    vector_setrowf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(dst+dst_row*n,src+row_src*n,n);        
}

template<typename T> void vector_copyrowf(T * dst, int dst_row, T * src, int row_src, int n) {
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);    
    vector_setrowf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(dst+dst_row*n,src+row_src*n,n);        
}

template<typename T> __global__ void vector_add_rowf_device(T * x,T * y, T * out, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) out[idx] = x[idx] + y[idx];
}

template<typename T> void vector_addf_row(T * x, int row, T * y, int row_src, size_t n)
{    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);    
    vector_add_rowf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x+row,y+row_src*n,x+row,n);        
}

template<typename T> __global__ void vector_sub_rowf_device(T * x, T * y, T * out, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) out[idx] = x[idx] - y[idx];
}

template<typename T> void vector_subf_row(T * x, int row, T * y, int row_src, size_t n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);    
    vector_sub_rowf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x+row,y+row_src,x+row,n);        
}

template<typename T> __global__ void vector_mul_rowf_device(T * x, T * y, T * out, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) out[idx] = x[idx] * y[idx];    
}

template<typename T> void vector_mulf_row(T * x, int row, T * y, int row_src, size_t n)
{    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);    
    vector_mul_rowf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x+row,y+row_src,x+row,n);            

}


template<typename T> __global__ void vector_div_rowf_device(T * x, T * y, T * out, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) out[idx] = x[idx] / y[idx];
}

template<typename T> void vector_divf_row(T * x,int row, T * y, int row_src, size_t n)
{   
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);    
    vector_div_rowf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x+row,y+row_src,x+row,n);            
}

template<typename T> __global__ void vector_mod_rowf_device(T * x, T * y, T * out, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) out[idx] = Fmod(x[idx],y[idx]);
}

template<typename T> void vector_modf_row(T * x, int row, T * y, int row_src, size_t n)
{    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);    
    vector_mod_rowf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x+row,y+row_src,x+row,n);            
}
// vector.cu
template<typename T> __global__ void vector_cbrtf_device(T * in, T * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = Cbrt(in[idx]);
}

template<typename T> void vector_r_cbrtf(T * devPtr, T * output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);    
    vector_cbrtf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n);
}



template<typename T> __global__ void vector_copysignf_device(T * x, T * y, T * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = Copysign(x[idx],y[idx]);
}

template<typename T> void vector_r_copysignf(T * X, T *Y, T *output, int n)
{    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_copysignf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(X,Y,output,n); 
}



template<typename T> __global__ void vector_cospif_device(T * in, T * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = Cospi(in[idx]);
}

template<typename T> void vector_r_cospif(T * devPtr, T * output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_cospif_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n);
}



template<typename T> __global__ void vector_cyl_bessel_i0f_device(T * in, T * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = Besseli0(in[idx]);
}

template<typename T> void vector_r_cyl_bessel_i0f(T * devPtr, T * output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_cyl_bessel_i0f_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n);
}


template<typename T> __global__ void vector_cyl_bessel_i1f_device(T * in, T * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = Besseli1(in[idx]);
}

template<typename T> void vector_r_cyl_bessel_i1f(T * devPtr, T * output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_cyl_bessel_i1f_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n);
}

template<typename T> __global__ void vector_erfcf_device(T * in, T * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = Erfc(in[idx]);
}

template<typename T> void vector_r_erfcf(T * devPtr, T * output, int n)
{ 
    int gridSize,blockSize;    
    calcSize(n,&gridSize,&blockSize);        
    vector_erfcf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n);
}

template<typename T> __global__ void vector_erfcinvf_device(T * in, T * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = Erfcinv(in[idx]);
}


template<typename T> void vector_r_erfcinvf(T * devPtr, T * output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_erfcinvf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n);
}

template<typename T> __global__ void vector_erfcxf_device(T * in, T * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = Erfcx(in[idx]);
}

template<typename T> void vector_r_erfcxf(T * devPtr, T * output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_erfcxf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n);
}
template<typename T> __global__ void vector_erff_device(T * in, T * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = Erf(in[idx]);
}

template<typename T> void vector_r_erff(T * devPtr, T * output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_erff_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n);
}

template<typename T> __global__ void vector_erfinvf_device(T * in, T * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = Erfinv(in[idx]);
}


template<typename T> void vector_r_erfinvf(T * devPtr, T * output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_erfinvf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n);
}


template<typename T> __global__ void vector_fdimf_device(T * a, T * b, T * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = FDim(a[idx],b[idx]);
}

template<typename T> void vector_r_fdimf(T * a, T * b, T * output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_fdimf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(a,b,output,n);
}

template<typename T> __global__ void vector_fdividef_device(T * a, T * b, T * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = FDivide(a[idx],b[idx]);
}

template<typename T> void vector_r_fdividef(T * a, T * b, T * output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_fdividef_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(a,b,output,n);
}

template<typename T> __global__ void vector_fmaf_device(T * a, T * b, T * c, T * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = FMa(a[idx],b[idx],c[idx]);
}

template<typename T> void vector_r_fmaf(T * x, T * y, T * z, T *output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_fmaf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,z,output,n);
}


template<typename T> __global__ void vector_hypotf_device(T * a,T * b, T * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = Hypot(a[idx],b[idx]);
}

template<typename T> void vector_r_hypotf(T * x, T * y, T * output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_hypotf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);
}
template<typename T> __global__ void vector_ilogbf_device(T * a, T * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = ILogb(a[idx]);
}

template<typename T> void vector_r_ilogbf(T * x, T *output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_ilogbf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);
}

template<typename T> __global__ void vector_j0f_device(T * a, T * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = J0(a[idx]);
}

template<typename T> void vector_r_j0f(T * x, T *output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_j0f_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);
}


template<typename T> __global__ void vector_j1f_device(T * a, T * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = j1f(a[idx]);
}

template<typename T> void vector_r_j1f(T * x, T *output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_j1f_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);
}


template<typename T> __global__ void vector_jnf_device(T * a, int N, T * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = jnf(a[idx],N);
}

template<typename T> void vector_r_jnf(T * x, T * output, int M, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_jnf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,M,output,n);
}

template<typename T> __global__ void vector_ldexpf_device(T * a, int exp, T * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = Ldexp(a[idx],exp);
}

template<typename T> void vector_r_ldexpf(T * x, T * output, int exp, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_ldexpf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,exp,output,n);
}

template<typename T> __global__ void vector_lgammaf_device(T * a, T * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = LGamma(a[idx]);
}

template<typename T> void vector_r_lgammaf(T * x, T *output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_lgammaf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);
}


/*
template<typename T> __global__ void vector_llrintf_device(T * a, long long * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = llrintf(a[idx]);
}


template<typename T> 
long long* vector_llrintf(T * x, int n)
{
    long long * p;
    cudaMalloc((void**)&p,sizeof(long long)*n);
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_llrintf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,p,n);
    return p;
}

template<typename T> __global__ void vector_llroundf_device(T * a, long long * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = llroundf(a[idx]);
}

template<typename T> 
long long* vector_llroundf(T * x, int n)
{
    long long * p;
    cudaMalloc((void**)&p,sizeof(long long)*n);
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_llroundf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,p,n);
    return p;
}

template<typename T> __global__ void vector_lrintf_device(T * a, long * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = lrintf(a[idx]);
}

template<typename T> 
long * vector_lrintf(T * x, int n)
{
    long * p;
    cudaMalloc((void**)&p,sizeof(long )*n);
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_lrintf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,p,n);
    return p;
}

template<typename T> __global__ void vector_lroundf_device(T * a, long * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = lroundf(a[idx]);
}

template<typename T> 
long * vector_lroundf(T * x, int n)
{
    long * p;
    cudaMalloc((void**)&p,sizeof(long )*n);
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_lroundf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,p,n);
    return p;
}
*/

template<typename T> __global__ void vector_nearbyintf_device(T * a, T * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = Nearbyint(a[idx]);
}

template<typename T> void vector_r_nearbyintf(T * x, T *output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_nearbyintf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);
}
template<typename T> __global__ void vector_norm3df_device(T * a, T * b, T * c, T* out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = Norm3d(a[idx],b[idx],c[idx]);
}

template<typename T> void vector_r_norm3df(T * x, T * y, T * z, T * output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_norm3df_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,z,output,n);
}

template<typename T> __global__ void vector_norm4df_device(T * a, T * b, T * c, T * d, T* out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = Norm4d(a[idx],b[idx],c[idx],d[idx]);
}

template<typename T> void vector_r_norm4df(T * x, T * y, T * z, T * q, T * output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_norm4df_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,z,q,output,n);
}


template<typename T> __global__ void vector_normcdff_device(T * a, T * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = Normcdf(a[idx]);
}

template<typename T> void vector_r_normcdff(T * x, T * output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_normcdff_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);    
}



template<typename T> __global__ void vector_normcdfinvf_device(T * a, T * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = Normcdfinv(a[idx]);
}

template<typename T> void vector_r_normcdfinvf(T * x, T *output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_normcdfinvf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);
}


template<typename T> __global__ void vector_normf_device(T * a, complex<T> * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = Norm(a[idx]);
}


template<typename T> void vector_r_normf(complex<T> * x, complex<T> * output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_normf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);
}

template<typename T> __global__ void vector_normf_device(int dim, T * a, T * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = Norm(dim,a);
}

template<typename T> void vector_r_normf(int dim, T * x, T * output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_normf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(dim,x,output,n);
}

template<typename T> __global__ void vector_rcbrtf_device(T * a, T * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = Rcbrt(a[idx]);
}

template<typename T> void vector_r_rcbrtf(T * x, T *output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_rcbrtf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);
}



template<typename T> __global__ void vector_remainderf_device(T * a, T * b, T * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = Remainder(a[idx],b[idx]);
}

template<typename T> void vector_r_remainderf(T * x, T * y, T *output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_remainderf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);
}


template<typename T> __global__ void vector_rhypotf_device(T * a, T * b, T * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = RHypot(a[idx],b[idx]);
}

template<typename T> void vector_r_rhypotf(T * x, T * y, T *output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_rhypotf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n); 
}





template<typename T> __global__ void vector_rnorm3df_device(T * a, T * b, T * c, T* out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = RNorm3d(a[idx],b[idx],c[idx]);
}

template<typename T> void vector_r_rnorm3df(T * x, T * y, T * z, T * output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_rnorm3df_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,z,output,n);
}


template<typename T> __global__ void vector_rnorm4df_device(T * a, T * b, T * c, T * d, T* out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = RNorm4d(a[idx],b[idx],c[idx],d[idx]);
}


template<typename T> void vector_r_rnorm4df(T * x, T * y, T * z, T * q, T *output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_rnorm4df_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,z,q,output,n);
}


template<typename T> __global__ void vector_rnormf_device(int dim, T * a, T * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = RNorm(dim,a);
}

template<typename T> void vector_r_rnormf(int dim, T * x, T *output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_rnormf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(dim,x,output,n);
}

template<typename T> __global__ void vector_scalblnf_device(T * a, long int N,T * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = Scalbln(a[idx],N);
}


template<typename T> void vector_r_scalblnf(T * x, long int M, T * output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_scalblnf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,M,output,n);
}

template<typename T> __global__ void vector_sinpif_device(T * a, T * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = Sinpi(a[idx]);
}

template<typename T> void vector_r_sinpif(T * x, T *output, int n)
{    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_sinpif_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n); 
}

template<typename T> __global__ void vector_tgammaf_device(T * a, T * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = TGamma(a[idx]);
}

template<typename T> void vector_r_tgammaf(T * x, T * output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_tgammaf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);
}

template<typename T> __global__ void vector_truncf_device(T * a, T * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = Trunc(a[idx]);
}


template<typename T> void vector_r_truncf(T * x, T *output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_truncf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);
}

template<typename T> __global__ void vector_y0f_device(T * a, T * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = Y0(a[idx]);
}

template<typename T> void vector_r_y0f(T * x, T *output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_y0f_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);
}



template<typename T> __global__ void vector_y1f_device(T * a, T * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = Y1(a[idx]);
}

template<typename T> void vector_r_y1f(T * x, T * output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_y1f_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);
}


template<typename T> __global__ void vector_ynf_device(int N, T * a, T * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = Yn(n,a[idx]);
}

template<typename T> void vector_r_ynf(int M, T * x, T *output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_ynf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(M,x,output,n);
}


template<typename T> __global__ void vector_fdimf_const_device(T * a, T b, T * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = FDim(a[idx],b);
}

template<typename T> void vector_r_fdimf_const(T * a, T  b, T * output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_fdimf_const_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(a,b,output,n);        
}


template<typename T> __global__ void vector_fdividef_const_device(T * a, T b, T * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = FDivide(a[idx],b);
}

template<typename T> void vector_r_fdividef_const(T * a, T b, T *p, int n)
{    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_fdividef_const_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(a,b,p,n);     
}


template<typename T> __global__ void vector_hypotf_const_device(T * a,T b, T * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = Hypot(a[idx],b);
}

template<typename T> void vector_r_hypotf_const(T * x, T y, T *p, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_hypotf_const_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,p,n);    
}



template<typename T> __global__ void vector_remainderf_const_device(T * a, T b, T * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = Remainder(a[idx],b);
}


template<typename T> void vector_r_remainderf_const(T * x, T y, T *p, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_remainderf_const_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,p,n);        
}

template<typename T> __global__ void vector_rhypotf_const_device(T * a, T b, T * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = RHypot(a[idx],b);
}

template<typename T> void vector_r_rhypotf_const(T * x, T y, T *p, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_rhypotf_const_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,p,n);    
}

