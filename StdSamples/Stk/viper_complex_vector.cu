/////////////////////////////////////////////////////////////////////////////////////////
// vector
/////////////////////////////////////////////////////////////////////////////////////////


#include <cassert>
#include "cuda_runtime.h"
#include "math_constants.h"
#include <map>
#include <iostream>
#include <vector>
#include <map>
#include <string>
#include "viper_vector.hpp"
#include "vector_complex.h"


template  complex<float>* cvector_addf(complex<float> * x, complex<float> * y, int n);
template  complex<float>* cvector_subf(complex<float> * x, complex<float> * y, int n);
template  complex<float>* cvector_mulf(complex<float> * x, complex<float> * y, int n);
template  complex<float>* cvector_divf(complex<float> * x, complex<float> * y, int n);


template  void cvector_realf(complex<float> * x, float * y, int n);
template  void cvector_imagf(complex<float> * x, float * y, int n);

template  void cvector_argf(complex<float> * devPtr, float * y, int n);
template  void cvector_normf(complex<float> * devPtr, complex<float> * y, int n);
template  void cvector_conjf(complex<float> * devPtr, complex<float> * y, int n);
template  void cvector_projf(complex<float> * devPtr, complex<float> * y, int n);
//template  void cvector_polar(float * r, float *theta,  complex<float> * y, int n);

template  complex<float>* cvector_expf(complex<float> * devPtr, int n);
template  complex<float>* cvector_logf(complex<float> * x, int n);
template  complex<float>* cvector_log10f(complex<float> * x, int n);

template  complex<float>* cvector_acosf(complex<float> * devPtr, int n);
template  complex<float>* cvector_acoshf(complex<float> * devPtr, int n);
template  complex<float>* cvector_asinf(complex<float> * devPtr, int n);
template  complex<float>* cvector_asinhf(complex<float> * devPtr, int n);
template  complex<float>* cvector_atanf(complex<float> * devPtr, int n);
template  complex<float>* cvector_atanhf(complex<float> * devPtr, int n);
template  complex<float>* cvector_cosf(complex<float> * devPtr, int n);
template  complex<float>* cvector_coshf(complex<float> * devPtr, int n);
template  complex<float>* cvector_sinf(complex<float> * x, int n);
template  complex<float>* cvector_sinhf(complex<float> * x, int n);
template  complex<float>* cvector_tanf(complex<float> * x, int n);
template  complex<float>* cvector_tanhf(complex<float> * x, int n);
template  complex<float>* cvector_powf(complex<float> * x, complex<float> * y, int n);
template  complex<float>* cvector_sqrtf(complex<float> * x, int n);

template  complex<float>* cvector_sigmoidf(complex<float> * devPtr, int n);
template  complex<float>* cvector_sigmoid_gradf(complex<float> * devPtr, int n);
template  complex<float>* cvector_tanh_gradf(complex<float> * devPtr, int n);

template  complex<float>* cvector_addf_const(complex<float> * x, complex<float>  y, int n);
template  complex<float>* cvector_subf_const(complex<float> * x, complex<float>  y, int n);
template  complex<float>* cvector_mulf_const(complex<float> * x, complex<float>  y, int n);
template  complex<float>* cvector_divf_const(complex<float> * x, complex<float>  y, int n);

template  complex<float>* cvector_powf_const(complex<float> * x, complex<float> y, int n);

template  complex<float>* cvector_addf_scalar(complex<float> * x, complex<float> * y, int n);
template  complex<float>* cvector_subf_scalar(complex<float> * x, complex<float> * y, int n);
template  complex<float>* cvector_mulf_scalar(complex<float> * x, complex<float> * y, int n);
template  complex<float>* cvector_divf_scalar(complex<float> * x, complex<float> * y, int n);

template  complex<float>* cvector_powf_scalar(complex<float> * x, complex<float> *y, int n);


template  void cvector_addf_row(complex<float> * x, int row, complex<float> * y, int src,size_t n);
template  void cvector_subf_row(complex<float> * x, int row, complex<float> * y, int src,size_t n);
template  void cvector_mulf_row(complex<float> * x, int row, complex<float> * y, int src,size_t n);
template  void cvector_divf_row(complex<float> * x, int row, complex<float> * y, int src,size_t n);

template  void cvector_setrowf(complex<float> * dst, int row, complex<float> * sc, int row_src, size_t n);
template  void cvector_copyrowf(complex<float> * dst, int row, complex<float> * src, int row_src, int n);

template  void cvector_r_addf(complex<float> * x, complex<float> * y, complex<float> * output, int n);
template  void cvector_r_subf(complex<float> * x, complex<float> * y, complex<float> * output, int n);
template  void cvector_r_mulf(complex<float> * x, complex<float> * y, complex<float> * output, int n);
template  void cvector_r_divf(complex<float> * x, complex<float> * y, complex<float> * output, int n);


template  void cvector_r_acosf(complex<float> * devPtr, complex<float> * outputs, int n);
template  void cvector_r_asinf(complex<float> * devPtr, complex<float> * outputs, int n);
template  void cvector_r_atanf(complex<float> * devPtr, complex<float> * outputs, int n);
template  void cvector_r_acoshf(complex<float> * devPtr, complex<float> * outputs, int n);
template  void cvector_r_asinhf(complex<float> * devPtr, complex<float> * outputs, int n);
template  void cvector_r_atanhf(complex<float> * devPtr, complex<float> * outputs, int n);
template  void cvector_r_cosf(complex<float> * devPtr, complex<float> * outputs, int n);
template  void cvector_r_sinf(complex<float> * devPtr, complex<float> * outputs, int n);
template  void cvector_r_tanf(complex<float> * devPtr, complex<float> * outputs, int n);
template  void cvector_r_coshf(complex<float> * devPtr, complex<float> * outputs, int n);
template   void cvector_r_sinhf(complex<float> * devPtr, complex<float> * outputs, int n);
template   void cvector_r_tanhf(complex<float> * devPtr, complex<float> * outputs, int n);

template   void cvector_r_expf(complex<float> * devPtr, complex<float> * output, int n);
template   void cvector_r_log10f(complex<float> * x, complex<float> *output, int n);
template   void cvector_r_powf(complex<float> * x, complex<float> * y, complex<float> *output, int n);
template   void cvector_r_sqrtf(complex<float> * x, complex<float> *output, int n);

template   void cvector_r_sigmoidf(complex<float> * x, complex<float> *output, int n);
template   void cvector_r_sigmoid_gradf(complex<float> * x, complex<float> *output, int n);
template   void cvector_r_tanh_gradf(complex<float> * x, complex<float> *output, int n);

template   void cvector_r_addf_const(complex<float> * x, complex<float>  y, complex<float> *output, int n);
template   void cvector_r_subf_const(complex<float> * x, complex<float>  y, complex<float> *output,int n);
template   void cvector_r_mulf_const(complex<float> * x, complex<float>  y, complex<float> *output,int n);
template   void cvector_r_divf_const(complex<float> * x, complex<float>  y, complex<float> *output,int n);


template   void cvector_r_powf_const(complex<float> * x, complex<float> y, complex<float> *output,int n);
template   void cvector_r_addf_scalar(complex<float> * x, complex<float> * y, complex<float> *output,int n);
template   void cvector_r_subf_scalar(complex<float> * x, complex<float> * y, complex<float> *output,int n);
template   void cvector_r_mulf_scalar(complex<float> * x, complex<float> * y, complex<float> *output,int n);
template   void cvector_r_divf_scalar(complex<float> * x, complex<float> * y, complex<float> *output,int n);
template   void cvector_r_powf_scalar(complex<float> * x, complex<float> *y, complex<float> *output,int n);

template   void cvector_r_hadamardf(complex<float> * a, complex<float> * b, complex<float> *output, int N);

template  complex<double>* cvector_addf(complex<double> * x, complex<double> * y, int n);
template  complex<double>* cvector_subf(complex<double> * x, complex<double> * y, int n);
template  complex<double>* cvector_mulf(complex<double> * x, complex<double> * y, int n);
template  complex<double>* cvector_divf(complex<double> * x, complex<double> * y, int n);

template  void cvector_realf(complex<double> * x, double * y, int n);
template  void cvector_imagf(complex<double> * x, double * y, int n);

template  void cvector_argf(complex<double> * devPtr, double * y, int n);
template  void cvector_normf(complex<double> * devPtr, complex<double> * y, int n);
template  void cvector_conjf(complex<double> * devPtr, complex<double> * y, int n);
template  void cvector_projf(complex<double> * devPtr, complex<double> * y, int n);
//template  void cvector_polar(double * r, double *theta,  complex<double> * y, int n);

template  complex<double>* cvector_expf(complex<double> * devPtr, int n);
template  complex<double>* cvector_logf(complex<double> * x, int n);
template  complex<double>* cvector_log10f(complex<double> * x, int n);

template  complex<double>* cvector_acosf(complex<double> * devPtr, int n);
template  complex<double>* cvector_acoshf(complex<double> * devPtr, int n);
template  complex<double>* cvector_asinf(complex<double> * devPtr, int n);
template  complex<double>* cvector_asinhf(complex<double> * devPtr, int n);
template  complex<double>* cvector_atanf(complex<double> * devPtr, int n);
template  complex<double>* cvector_atanhf(complex<double> * devPtr, int n);
template  complex<double>* cvector_cosf(complex<double> * devPtr, int n);
template  complex<double>* cvector_coshf(complex<double> * devPtr, int n);
template  complex<double>* cvector_sinf(complex<double> * x, int n);
template  complex<double>* cvector_sinhf(complex<double> * x, int n);
template  complex<double>* cvector_tanf(complex<double> * x, int n);
template  complex<double>* cvector_tanhf(complex<double> * x, int n);
template  complex<double>* cvector_powf(complex<double> * x, complex<double> * y, int n);
template  complex<double>* cvector_sqrtf(complex<double> * x, int n);

template  complex<double>* cvector_sigmoidf(complex<double> * devPtr, int n);
template  complex<double>* cvector_sigmoid_gradf(complex<double> * devPtr, int n);
template  complex<double>* cvector_tanh_gradf(complex<double> * devPtr, int n);

template  complex<double>* cvector_addf_const(complex<double> * x, complex<double>  y, int n);
template  complex<double>* cvector_subf_const(complex<double> * x, complex<double>  y, int n);
template  complex<double>* cvector_mulf_const(complex<double> * x, complex<double>  y, int n);
template  complex<double>* cvector_divf_const(complex<double> * x, complex<double>  y, int n);

template  complex<double>* cvector_powf_const(complex<double> * x, complex<double> y, int n);

template  complex<double>* cvector_addf_scalar(complex<double> * x, complex<double> * y, int n);
template  complex<double>* cvector_subf_scalar(complex<double> * x, complex<double> * y, int n);
template  complex<double>* cvector_mulf_scalar(complex<double> * x, complex<double> * y, int n);
template  complex<double>* cvector_divf_scalar(complex<double> * x, complex<double> * y, int n);

template  complex<double>* cvector_powf_scalar(complex<double> * x, complex<double> *y, int n);


template  void cvector_addf_row(complex<double> * x, int row, complex<double> * y, int src,size_t n);
template  void cvector_subf_row(complex<double> * x, int row, complex<double> * y, int src,size_t n);
template  void cvector_mulf_row(complex<double> * x, int row, complex<double> * y, int src,size_t n);
template  void cvector_divf_row(complex<double> * x, int row, complex<double> * y, int src,size_t n);

template  void cvector_setrowf(complex<double> * dst, int row, complex<double> * sc, int row_src, size_t n);
template  void cvector_copyrowf(complex<double> * dst, int row, complex<double> * src, int row_src, int n);

template  void cvector_r_addf(complex<double> * x, complex<double> * y, complex<double> * output, int n);
template  void cvector_r_subf(complex<double> * x, complex<double> * y, complex<double> * output, int n);
template  void cvector_r_mulf(complex<double> * x, complex<double> * y, complex<double> * output, int n);
template  void cvector_r_divf(complex<double> * x, complex<double> * y, complex<double> * output, int n);


template  void cvector_r_acosf(complex<double> * devPtr, complex<double> * outputs, int n);
template  void cvector_r_asinf(complex<double> * devPtr, complex<double> * outputs, int n);
template  void cvector_r_atanf(complex<double> * devPtr, complex<double> * outputs, int n);
template  void cvector_r_acoshf(complex<double> * devPtr, complex<double> * outputs, int n);
template  void cvector_r_asinhf(complex<double> * devPtr, complex<double> * outputs, int n);
template  void cvector_r_atanhf(complex<double> * devPtr, complex<double> * outputs, int n);
template  void cvector_r_cosf(complex<double> * devPtr, complex<double> * outputs, int n);
template  void cvector_r_sinf(complex<double> * devPtr, complex<double> * outputs, int n);
template  void cvector_r_tanf(complex<double> * devPtr, complex<double> * outputs, int n);
template  void cvector_r_coshf(complex<double> * devPtr, complex<double> * outputs, int n);
template   void cvector_r_sinhf(complex<double> * devPtr, complex<double> * outputs, int n);
template   void cvector_r_tanhf(complex<double> * devPtr, complex<double> * outputs, int n);

template   void cvector_r_expf(complex<double> * devPtr, complex<double> * output, int n);
template   void cvector_r_log10f(complex<double> * x, complex<double> *output, int n);
template   void cvector_r_powf(complex<double> * x, complex<double> * y, complex<double> *output, int n);
template   void cvector_r_sqrtf(complex<double> * x, complex<double> *output, int n);

template   void cvector_r_sigmoidf(complex<double> * x, complex<double> *output, int n);
template   void cvector_r_sigmoid_gradf(complex<double> * x, complex<double> *output, int n);
template   void cvector_r_tanh_gradf(complex<double> * x, complex<double> *output, int n);

template   void cvector_r_addf_const(complex<double> * x, complex<double>  y, complex<double> *output, int n);
template   void cvector_r_subf_const(complex<double> * x, complex<double>  y, complex<double> *output,int n);
template   void cvector_r_mulf_const(complex<double> * x, complex<double>  y, complex<double> *output,int n);
template   void cvector_r_divf_const(complex<double> * x, complex<double>  y, complex<double> *output,int n);


template   void cvector_r_powf_const(complex<double> * x, complex<double> y, complex<double> *output,int n);
template   void cvector_r_addf_scalar(complex<double> * x, complex<double> * y, complex<double> *output,int n);
template   void cvector_r_subf_scalar(complex<double> * x, complex<double> * y, complex<double> *output,int n);
template   void cvector_r_mulf_scalar(complex<double> * x, complex<double> * y, complex<double> *output,int n);
template   void cvector_r_divf_scalar(complex<double> * x, complex<double> * y, complex<double> *output,int n);
template   void cvector_r_powf_scalar(complex<double> * x, complex<double> *y, complex<double> *output,int n);

template   void cvector_r_hadamardf(complex<double> * a, complex<double> * b, complex<double> *output, int N);

cudaStream_t get_cuda_stream();
cudaStream_t random_stream();
void    calcSize(int N,int * gridSize, int * blockSize);

template<typename T>
__global__ void cvector_addf_device(complex<T> * x, complex<T> * y, complex<T> * out, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) out[idx] = x[idx] + y[idx];
}

template<typename T>
complex<T>* cvector_addf(complex<T> * x, complex<T> * y, int n)
{
    complex<T> * output = find_memory<complex<T>>(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(complex<T>));
    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);    
    cvector_addf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);    
    return output;
}


template<typename T>
void cvector_r_addf(complex<T> * x, complex<T> * y, complex<T> * output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);    
    cvector_addf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);    
}



template<typename T>
__global__ void cvector_subf_device(complex<T> * x, complex<T> * y, complex<T> * out, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) out[idx] = x[idx] - y[idx];
}

template<typename T>
complex<T>* cvector_subf(complex<T> * x, complex<T> * y, int n)
{
    complex<T> * output = find_memory<complex<T>>(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(complex<T>));
    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);    
    cvector_subf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);    
    return output;
}


template<typename T>
void cvector_r_subf(complex<T> * x, complex<T> * y, complex<T> * output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);    
    cvector_subf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);    
}


template<typename T>
__global__ void cvector_mulf_device(complex<T> * x, complex<T> * y, complex<T> * out, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) out[idx] = x[idx] * y[idx];
}

template<typename T>
complex<T>* cvector_mulf(complex<T> * x, complex<T> * y, int n)
{
    complex<T> * output = find_memory<complex<T>>(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(complex<T>));
    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);    
    cvector_mulf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);    
    return output;
}

template<typename T>
void cvector_r_mulf(complex<T> * x, complex<T> * y, complex<T> * output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);    
    cvector_mulf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);    
}

template<typename T>
void cvector_r_hadamardf(complex<T> * x, complex<T> * y, complex<T> * output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);    
    cvector_mulf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);    
}

template<typename T>
__global__ void cvector_divf_device(complex<T> * x, complex<T> * y, complex<T> * out, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) out[idx] = x[idx] / y[idx];
}


template<typename T>
complex<T>* cvector_divf(complex<T> * x, complex<T> * y, int n)
{
    complex<T> * output = find_memory<complex<T>>(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(complex<T>));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);    
    cvector_divf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);    
    return output;
}

template<typename T>
void cvector_r_divf(complex<T> * x, complex<T> * y, complex<T> * output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);    
    cvector_divf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);    
}


template<typename T>
__global__ void cvector_acosf_device(complex<T> * in, complex<T> * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n) out[idx] = acos(in[idx]);
}

template<typename T>
complex<T>* cvector_acosf(complex<T> * devPtr, int n)
{
    complex<T> * output = find_memory<complex<T>>(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(complex<T>));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);    
    cvector_acosf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n);    
    return output;
}

template<typename T>
void cvector_r_acosf(complex<T> * devPtr, complex<T> * output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);    
    cvector_acosf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n);        
}

template<typename T>
__global__ void cvector_acoshf_device(complex<T> * in, complex<T> * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n) out[idx] = acosh(in[idx]);
}

template<typename T>
complex<T>* cvector_acoshf(complex<T> * devPtr, int n)
{
    complex<T> * output = find_memory<complex<T>>(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(complex<T>));
    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);    
    cvector_acoshf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n);
    return output;
}

template<typename T>
void cvector_r_acoshf(complex<T> * devPtr, complex<T> * output, int n)
{    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);    
    cvector_acoshf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n); 
}

template<typename T>
__global__ void cvector_asinhf_device(complex<T> * in, complex<T> * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = asinh(in[idx]);
}


template<typename T>
complex<T>* cvector_asinhf(complex<T> * devPtr, int n)
{
    complex<T> * output = find_memory<complex<T>>(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(complex<T>));
    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);    
    cvector_asinhf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n);
    return output;
}

template<typename T>
void cvector_r_asinhf(complex<T> * devPtr, complex<T> * output, int n)
{    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);    
    cvector_asinhf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n); 
}


template<typename T>
__global__ void cvector_asinf_device(complex<T> * in, complex<T> * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = asin(in[idx]);
}


template<typename T>
complex<T>* cvector_asinf(complex<T> * devPtr, int n)
{
    complex<T> * output = find_memory<complex<T>>(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(complex<T>));
    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);    
    cvector_asinf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n);
    return output;
}

template<typename T>
void cvector_r_asinf(complex<T> * devPtr, complex<T> * output, int n)
{    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);    
    cvector_asinf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n);
}

template<typename T>
__global__ void cvector_atanf_device(complex<T> * in, complex<T> * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = atan(in[idx]);
}


template<typename T>
complex<T>* cvector_atanf(complex<T> * devPtr, int n)
{
    complex<T> * output = find_memory<complex<T>>(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(complex<T>));
    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);    
    cvector_atanf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n);
    return output;
}

template<typename T>
void cvector_r_atanf(complex<T> * devPtr, complex<T> * output, int n)
{    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);    
    cvector_atanf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n);    
}


template<typename T>
__global__ void cvector_atanhf_device(complex<T> * in, complex<T> * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = atanh(in[idx]);
}

template<typename T>
complex<T>* cvector_atanhf(complex<T> * devPtr, int n)
{
    complex<T> * output = find_memory<complex<T>>(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(complex<T>));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);    
    cvector_atanhf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n);
    return output;
}


template<typename T>
void cvector_r_atanhf(complex<T> * devPtr, complex<T> * output, int n)
{    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);    
    cvector_atanhf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n); 
}


template<typename T>
__global__ void cvector_cosf_device(complex<T> * in, complex<T> * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = cos(in[idx]);
}

template<typename T>
complex<T>* cvector_cosf(complex<T> * devPtr, int n)
{
    complex<T> * output = find_memory<complex<T>>(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(complex<T>));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    cvector_cosf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n);
    return output;
}

template<typename T>
void cvector_r_cosf(complex<T> * devPtr, complex<T> * output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    cvector_cosf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n);    
}


template<typename T>
__global__ void cvector_coshf_device(complex<T> * in, complex<T> * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = cosh(in[idx]);
}

template<typename T>
complex<T>* cvector_coshf(complex<T> * devPtr, int n)
{
    complex<T> * output = find_memory<complex<T>>(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(complex<T>));
    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    cvector_coshf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n);
    return output;
}

template<typename T>
void cvector_r_coshf(complex<T> * devPtr, complex<T> * output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    cvector_coshf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n);
}


template<typename T>
__global__ void cvector_expf_device(complex<T> * in, complex<T> * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = exp(in[idx]);
}

template<typename T>
complex<T>* cvector_expf(complex<T> * devPtr, int n)
{
    complex<T> * output = find_memory<complex<T>>(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(complex<T>));
    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    cvector_expf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n);
    return output;
}

template<typename T>
void cvector_r_expf(complex<T> * devPtr, complex<T> * output, int n)
{   
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    cvector_expf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n); 
}


template<typename T>
__global__ void cvector_log10f_device(complex<T> * a, complex<T> * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = log10(a[idx]);
}


template<typename T>
complex<T>* cvector_log10f(complex<T> * x, int n)
{
    complex<T> * output = find_memory<complex<T>>(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(complex<T>));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    cvector_log10f_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);
    return output;
}

template<typename T>
void cvector_r_log10f(complex<T> * x, complex<T> *output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    cvector_log10f_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);
}

template<typename T>
__global__ void cvector_logf_device(complex<T> * a, complex<T> * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = log(a[idx]);
}


template<typename T>
complex<T>* cvector_logf(complex<T> * x, int n)
{
    complex<T> * output = find_memory<complex<T>>(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(complex<T>));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    cvector_logf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);
    return output;
}

template<typename T>
void cvector_r_logf(complex<T> * x, complex<T> *output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    cvector_logf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);
}

template<typename T> 
__global__ void cvector_powf_device(complex<T> * a, complex<T> * b, complex<T> * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = pow(a[idx],b[idx]);
}


template<typename T>
complex<T>* cvector_powf(complex<T> * x, complex<T> * y, int n)
{
    complex<T> * output = find_memory<complex<T>>(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(complex<T>));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    cvector_powf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);
    return output;
}

template<typename T>
void cvector_r_powf(complex<T> * x, complex<T> * y, complex<T> *output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    cvector_powf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);
}

template<typename T> 
__global__ void cvector_sinf_device(complex<T> * a, complex<T> * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = sin(a[idx]);
}

template<typename T>
complex<T>* cvector_sinf(complex<T> * x, int n)
{
    complex<T> * output = find_memory<complex<T>>(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(complex<T>));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    cvector_sinf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);
    return output;
}

template<typename T>
void cvector_r_sinf(complex<T> * x, complex<T> *output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    cvector_sinf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);
}


template<typename T>
__global__ void cvector_sinhf_device(complex<T> * a, complex<T> * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = sinh(a[idx]);
}

template<typename T>
complex<T>* cvector_sinhf(complex<T> * x, int n)
{
    complex<T> * output = find_memory<complex<T>>(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(complex<T>));
    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    cvector_sinhf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);
    return output;
}

template<typename T>
void cvector_r_sinhf(complex<T> * x, complex<T> *output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    cvector_sinhf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);    
}

template<typename T>
__global__ void cvector_sqrtf_device(complex<T> * a, complex<T> * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = sqrt(a[idx]);
}

template<typename T>
complex<T>* cvector_sqrtf(complex<T> * x, int n)
{
    complex<T> * output = find_memory<complex<T>>(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(complex<T>));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    cvector_sqrtf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);
    return output;
}

template<typename T>
void cvector_r_sqrtf(complex<T> * x, complex<T> *output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    cvector_sqrtf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);
}


template<typename T>
__global__ void cvector_tanf_device(complex<T> * a, complex<T> * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = tan(a[idx]);
}


template<typename T>
complex<T>* cvector_tanf(complex<T> * x, int n)
{
    complex<T> * output = find_memory<complex<T>>(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(complex<T>));
    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    cvector_tanf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);
    return output;
}

template<typename T>
void cvector_r_tanf(complex<T> * x, complex<T> *output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    cvector_tanf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);
}

template<typename T>
__global__ void cvector_tanhf_device(complex<T> * a, complex<T> * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = tanh(a[idx]);
}

template<typename T>
complex<T>* cvector_tanhf(complex<T> * x, int n)
{
    complex<T> * output = find_memory<complex<T>>(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(complex<T>));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    cvector_tanhf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);
    return output;
}

template<typename T>
void cvector_r_tanhf(complex<T> * x, complex<T> *output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    cvector_tanhf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);
}

template<typename T>
__global__ void cvector_sigmoid_device(complex<T> * x, complex<T> * out, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;         
    complex<T> one(1.0,0.0);
    if(idx < N) out[idx] = one / (one + exp(-x[idx]));
}

template<typename T>
complex<T>* cvector_sigmoidf(complex<T> * x, int n)
{
    complex<T> * output = find_memory<complex<T>>(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(complex<T>));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);            
    cvector_sigmoid_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);        
    return output;
}

template<typename T>
void cvector_r_sigmoidf(complex<T> * x, complex<T> *output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);            
    cvector_sigmoid_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);        
}

template<typename T>
__global__ void cvector_sigmoid_grad_device(complex<T> * x, complex<T> * out, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    complex<T> one(1.0,0.0);
    if(idx < N) out[idx] = x[idx] * (one - x[idx]);
}

template<typename T>
complex<T>* cvector_sigmoid_gradf(complex<T> * x, int n)
{
    complex<T> * output = find_memory<complex<T>>(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(complex<T>));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    cvector_sigmoid_grad_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);        
    return output;
}

template<typename T>
void cvector_r_sigmoid_gradf(complex<T> * x, complex<T> *output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    cvector_sigmoid_grad_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);        
}

template<typename T>
__global__ void cvector_tanh_grad_device(complex<T> * x, complex<T> * out, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    complex<T> one(1.0,0.0);
    if(idx < N) out[idx] = one - (x[idx]*x[idx]);
}


template<typename T>
complex<T>* cvector_tanh_gradf(complex<T> * x, int n)
{
    complex<T> * output = find_memory<complex<T>>(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(complex<T>));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    cvector_tanh_grad_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);        
    return output;
}


template<typename T>
void cvector_r_tanh_gradf(complex<T> * x, complex<T> *output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    cvector_tanh_grad_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);        
}


template<typename T>
__global__ void cvector_add_const_device(complex<T> * x, complex<T> y, complex<T> * out, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) out[idx] = x[idx] + y;
}


template<typename T>
complex<T>* cvector_addf_const(complex<T> * x, complex<T> y, int n)
{
    complex<T> * output = find_memory<complex<T>>(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(complex<T>));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    cvector_add_const_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);    
    return output;
}


template<typename T>
void cvector_r_addf_const(complex<T> * x, complex<T> y, complex<T> * output, int n)
{    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    cvector_add_const_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);     
}


template<typename T>
__global__ void cvector_sub_const_device(complex<T> * x, complex<T> y, complex<T> * out, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) out[idx] = x[idx] - y;
}

template<typename T>
complex<T>* cvector_subf_const(complex<T> * x, complex<T> y, int n)
{  
    complex<T> * output = find_memory<complex<T>>(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(complex<T>));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    cvector_sub_const_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);        
    return output;
}


template<typename T>
void cvector_r_subf_const(complex<T> * x, complex<T> y, complex<T> *output, int n)
{  
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    cvector_sub_const_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);          
}


template<typename T>
__global__ void cvector_mul_const_device(complex<T> * x, complex<T> y, complex<T> * out, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) out[idx] = x[idx] * y;
}

template<typename T>
complex<T>* cvector_mulf_const(complex<T> * x, complex<T> y, int n)
{
    complex<T> * output = find_memory<complex<T>>(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(complex<T>));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    cvector_mul_const_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);     
    return output;
}

template<typename T>
void cvector_r_mulf_const(complex<T> * x, complex<T> y, complex<T> * output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    cvector_mul_const_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);        
}

template<typename T>
__global__ void cvector_div_const_device(complex<T> * x, complex<T> y, complex<T> * out, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) out[idx] = x[idx] / y;
}

template<typename T>
complex<T>* cvector_divf_const(complex<T> * x, complex<T> y, int n)
{
    complex<T> * output = find_memory<complex<T>>(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(complex<T>));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    cvector_div_const_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);        
    return output;
}

template<typename T>
void cvector_r_divf_const(complex<T> * x, complex<T> y, complex<T> *output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    cvector_div_const_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);        
}

template<typename T>
__global__ void cvector_powf_const_device(complex<T> * a, complex<T> b, complex<T> * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = pow(a[idx],b);
}

template<typename T>
complex<T>* cvector_powf_const(complex<T> * x, complex<T> y, int n)
{
    complex<T> * p = find_memory<complex<T>>(n);
    if(p == NULL)
        cudaMalloc((void**)&p,sizeof(complex<T>)*n);
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    cvector_powf_const_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,p,n);    
    return p;
}

template<typename T>
void cvector_r_powf_const(complex<T> * x, complex<T> y, complex<T> *p, int n)
{    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    cvector_powf_const_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,p,n);    
}

/////////////////////////////////
// const/scalar
/////////////////////////////////
template<typename T>
__global__ void cvector_add_scalar_device(complex<T> * x, complex<T> * y, complex<T> * out, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) out[idx] = x[idx] + y[0];
}


template<typename T>
complex<T>* cvector_addf_scalar(complex<T> * x, complex<T> * y, int n)
{    
    complex<T> * output = find_memory<complex<T>>(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(complex<T>));
    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    cvector_add_scalar_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);        
    return output;
}

template<typename T>
void cvector_r_addf_scalar(complex<T> * x, complex<T> * y, complex<T> *output, int n)
{      
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    cvector_add_scalar_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);            
}

template<typename T>
__global__ void cvector_sub_scalar_device(complex<T> * x, complex<T> * y, complex<T> * out, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) out[idx] = x[idx] - y[0];
}

template<typename T>
complex<T>* cvector_subf_scalar(complex<T> * x, complex<T> * y, int n)
{    
    complex<T> * output = find_memory<complex<T>>(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(complex<T>));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    cvector_sub_scalar_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);        
    return output;
}

template<typename T>
void cvector_r_subf_scalar(complex<T> * x, complex<T> * y, complex<T> *output, int n)
{    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    cvector_sub_scalar_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);        
}

template<typename T>
__global__ void cvector_mul_scalar_device(complex<T> * x, complex<T> * y, complex<T> * out, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) out[idx] = x[idx] * y[0];
}


template<typename T>
complex<T>* cvector_mulf_scalar(complex<T> * x, complex<T> * y, int n)
{    
    complex<T> * output = find_memory<complex<T>>(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(complex<T>));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    cvector_mul_scalar_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);        
    return output;
}

template<typename T>
void cvector_r_mulf_scalar(complex<T> * x, complex<T> * y, complex<T> * output, int n)
{        
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    cvector_mul_scalar_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);            
}

template<typename T>
__global__ void cvector_div_scalar_device(complex<T> * x, complex<T> * y, complex<T> * out, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) out[idx] = x[idx] / y[0];
}


template<typename T>
complex<T>* cvector_divf_scalar(complex<T> * x, complex<T> * y, int n)
{    
    complex<T> * output = find_memory<complex<T>>(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(complex<T>));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    cvector_div_scalar_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);        
    return output;
}

template<typename T>
void cvector_r_divf_scalar(complex<T> * x, complex<T> * y, complex<T> *output, int n)
{    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    cvector_div_scalar_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);        
}


template<typename T>
__global__ void cvector_powf_scalar_device(complex<T> * x, complex<T> * y, complex<T> * out, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) out[idx] = pow(x[idx],y[0]);
}


template<typename T>
complex<T>* cvector_powf_scalar(complex<T> * x, complex<T> * y, int n)
{    
    complex<T> * output = find_memory<complex<T>>(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(complex<T>));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    cvector_powf_scalar_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);        
    return output;
}

template<typename T>
void cvector_r_powf_scalar(complex<T> * x, complex<T> * y, complex<T> *output, int n)
{    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    cvector_powf_scalar_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);        
}

template<typename T>
__global__ void cvector_setrowf_device(complex<T> * dst, complex<T> * src, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) dst[idx] = src[idx];
}

template<typename T>
void cvector_setrowf(complex<T> * dst, int dst_row, complex<T> * src, int row_src, size_t n)
{    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);    
    cvector_setrowf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(dst+dst_row*n,src+row_src*n,n);        
}

template<typename T>
void cvector_copyrowf(complex<T> * dst, int dst_row, complex<T> * src, int row_src, int n) {
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);    
    cvector_setrowf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(dst+dst_row*n,src+row_src*n,n);        
}


template<typename T>
__global__ void cvector_add_rowf_device(complex<T> * x,complex<T> * y, complex<T> * out, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) out[idx] = x[idx] + y[idx];
}


template<typename T>
void cvector_addf_row(complex<T> * x, int row, complex<T> * y, int row_src, size_t n)
{    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);    
    cvector_add_rowf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x+row,y+row_src*n,x+row,n);        
}

template<typename T>
__global__ void cvector_sub_rowf_device(complex<T> * x, complex<T> * y, complex<T> * out, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) out[idx] = x[idx] - y[idx];
}


template<typename T>
void cvector_subf_row(complex<T> * x, int row, complex<T> * y, int row_src, size_t n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);    
    cvector_sub_rowf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x+row,y+row_src,x+row,n);        
}


template<typename T>
__global__ void cvector_mul_rowf_device(complex<T> * x, complex<T> * y, complex<T> * out, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) out[idx] = x[idx] * y[idx];    
}


template<typename T>
void cvector_mulf_row(complex<T> * x, int row, complex<T> * y, int row_src, size_t n)
{    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);    
    cvector_mul_rowf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x+row,y+row_src,x+row,n);            

}

template<typename T>
__global__ void cvector_div_rowf_device(complex<T> * x, complex<T> * y, complex<T> * out, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) out[idx] = x[idx] / y[idx];
}


template<typename T>
void cvector_divf_row(complex<T> * x,int row, complex<T> * y, int row_src, size_t n)
{   
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);    
    cvector_div_rowf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x+row,y+row_src,x+row,n);            
}



template<typename T> __global__ void cvector_realf_device(complex<T> * a, T * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = a[idx].real();
}


template<typename T>
void cvector_realf(complex<T> * x, T * y, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    cvector_realf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,n);
}


template<typename T> __global__ void cvector_imagf_device(complex<T> * a, T * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = a[idx].imag();
}

template<typename T>
void cvector_imagf(complex<T> * x, T * y, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    cvector_imagf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,n);
}

template<typename T> __global__ void cvector_argf_device(complex<T> * a, T * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = arg(a[idx]);
}

template<typename T>
void cvector_argf(complex<T> * devPtr,  T* x, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    cvector_argf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,x,n);
}

template<typename T> __global__ void cvector_conjf_device(complex<T> * a, complex<T> * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = conj(a[idx]);
}

template<typename T>
void cvector_conjf(complex<T> * devPtr, complex<T> * x, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    cvector_conjf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,x,n);
}

template<typename T> __global__ void cvector_normf_device(complex<T> * a, complex<T> * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = norm(a[idx]);
}

template<typename T>
void cvector_normf(complex<T> * devPtr, complex<T> * x, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    cvector_normf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,x,n);
}


template<typename T> __global__ void cvector_projf_device(complex<T> * a, complex<T> * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = proj(a[idx]);
}

template<typename T>
void cvector_projf(complex<T> * devPtr, complex<T> * x, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    cvector_projf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,x,n);
}

/* does not exit in std::cuda::complex yet
template<typename T> __global__ void cvector_polarf_device(T * r, T * t, complex<T> * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = polar(r[idx],t[idx]);
}

template<typename T>
void cvector_polar(T * rho, T *theta, complex<T> * x, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    cvector_polarf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(rho,theta,x,n);
}
*/