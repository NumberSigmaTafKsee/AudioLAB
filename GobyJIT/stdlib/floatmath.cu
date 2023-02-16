////////////////////////////////////////////////////////////////////////
// vector
////////////////////////////////////////////////////////////////////////


#include <cassert>
#include "cuda_runtime.h"
#include "math_constants.h"
#include <map>
#include <iostream>
#include <vector>
#include <map>
#include <string>
#include "vector_float.h"


cudaStream_t get_cuda_stream();
cudaStream_t random_stream();

// memory cache
std::multimap<int,float*> cuda_memory;

void add_memory(int length, float * f) {
    cuda_memory.insert(std::pair<int,float*>(length,f));
}
void return_memory(int length, float *fp) {        
    cuda_memory.insert(std::pair<int,float*>(length,fp));    
}
float* find_memory(int length) {            
    typename std::multimap<int,float*>::iterator i = cuda_memory.find(length);
    if(i == cuda_memory.end()) return NULL;                
    cuda_zero(i->second,length);
    cuda_memory.erase(i);
    return i->second;
}
void clear_cache() {
    typename std::multimap<int,float*>::iterator i = cuda_memory.begin();
    while(i != cuda_memory.end()) {
        cudaFree(i->second);       
        i++; 
    }
    cuda_memory.clear();    
}


/* what we want is this 
__global__ void vector_program(float * x, float * out, int N) {
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    out[idx] = eval(x)
}

////////////////////////////////////
Vector X;
Matrix M;
Vector output;

f(x) = M*X;
////////////////////////////////////

=>
    int row = threadIdx.x + blockIdx.x * blockDim.x;     
    int col = threadIdx.y + blockIdx.y * blockDim.y;     
    if(row < m && col < n)
        out[col] = M[row*n + col] * X[col]

*/

/*
std::map<std::string,vector_kernel1> vector_map1;
std::map<std::string,vector_kernel2> vector_map2;
std::map<std::string,vector_kernel3> vector_map3;
std::map<std::string,vector_kernel4> vector_map4;

void register_vector_kernel1(const char * name, vector_kernel1 kernel) {
    // should assert if already exist?
    vector_map1[name] = kernel;
}
void register_vector_kernel2(const char * name, vector_kernel2 kernel) {
    vector_map2[name] = kernel;
}
void register_vector_kernel3(const char * name, vector_kernel3 kernel) {
    vector_map3[name] = kernel;
}
void register_vector_kernel4(const char * name, vector_kernel4 kernel) {
    vector_map4[name] = kernel;
}
float* execute_vector_kernel1(const char * name, float * x, int n) {
    typename std::map<std::string,vector_kernel1>::iterator i = vector_map1.find(name);
    // assert or return NULL?
    if(i == vector_map1.end()) return NULL;
    return (i->second)(x,n);
}
float* execute_vector_kernel2(const char * name, float * x, float * y, int n) {
    typename std::map<std::string,vector_kernel2>::iterator i = vector_map2.find(name);
    // assert or return NULL?
    if(i == vector_map2.end()) return NULL;
    return (i->second)(x,y,n);
}
float* execute_vector_kernel3(const char * name, float * x, float * y, float * z, int n) {
    typename std::map<std::string,vector_kernel3>::iterator i = vector_map3.find(name);
    // assert or return NULL?
    if(i == vector_map3.end()) return NULL;
    return (i->second)(x,y,z,n);
}
float* execute_vector_kernel4(const char * name, float * x, float * y, float * z, float * w, int n) {
    typename std::map<std::string,vector_kernel4>::iterator i = vector_map4.find(name);
    // assert or return NULL?
    if(i == vector_map4.end()) return NULL;
    return (i->second)(x,y,z,w,n);
}
*/


__global__ void vector_dummy(float * x, float * out, int N) {
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) out[idx] =x[idx];
}


void calcSize(int N,int * gridSize, int * blockSize) {
    
    //int minGridSize = 0;    
    //cudaOccupancyMaxPotentialBlockSize(&minGridSize, blockSize, vector_dummy, 0, N); 
    //*gridSize = (N + *blockSize - 1) / *blockSize; 
    
    
    *blockSize=1024;    
    *gridSize=(N+*blockSize)/ *blockSize;    
}

__global__ void cuda_memcpy_device(float * dst, float * src, size_t N) {
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) dst[idx] = src[idx];
}

void cuda_memcpy(float * dst, float * src, int n)
{    
    //int gridSize,blockSize;
    //calcSize(n,&gridSize,&blockSize);        
    //cuda_memcpy_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(dst,src,n);        
    //cudaDeviceSynchronize();
    cudaMemcpyAsync(dst,src,n*sizeof(float),cudaMemcpyDeviceToDevice, get_cuda_stream());
}

__global__ void cuda_zero_device(float * dst, size_t N) {
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) dst[idx] = 0;
}

void cuda_zero(float * dst, int n)
{    
    /*
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    cuda_zero_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(dst,n);        
    */    
    cudaMemsetAsync(dst,0,n*sizeof(float), get_cuda_stream());    
}



/// making this parallel is alot harder than you would think it should be.
/// I need atomics so I can add all at once.
__global__ void vector_sum_device(float * x,float *out, int N)
{
    
    //int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    //if(idx < N) *out += x[idx];
    for(int i = 0; i < N; i++) 
    {
        *out += x[i];
        __syncthreads();
    }
}

float vector_sumf(float * x, int n)
{        
    /*
    float * output = find_memory(1);
    if(output == NULL) cudaMalloc((void**)&output,sizeof(float));    
    float zero = 0.0f;
    cudaMemcpyAsync(output,&zero,sizeof(float),cudaMemcpyHostToDevice);
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);    
    vector_sum_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);        
    float o=0;
    cudaMemcpyAsync(&o,output,sizeof(float),cudaMemcpyDeviceToHost);
    return_memory(1,output);
    */
    float *tmp = (float*)calloc(n,sizeof(float));
    cudaMemcpyAsync(tmp,x,n*sizeof(float),cudaMemcpyDeviceToHost, get_cuda_stream());
    float o = 0.0;
    for(size_t i = 0; i < n; i++) o += tmp[i];
    free(tmp);    
    return o;
}

__global__ void vector_prod_device(float * x,float *out, int N)
{
    for(int i = 0; i < N; i++) {
        *out *= x[i];
        __syncthreads();
    }
}

float vector_prodf(float * x, int n)
{
    /*
    float * output = find_memory(1);
    if(output == NULL) cudaMalloc((void**)&output,sizeof(float));    
    float zero = 1.0f;
    cudaMemcpyAsync(output,&zero,sizeof(float),cudaMemcpyHostToDevice);
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);    
    vector_prod_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);        
    float o=0;
    cudaMemcpyAsync(&o,output,sizeof(float),cudaMemcpyDeviceToHost);
    return_memory(1,output);
    */
    float *tmp = (float*)calloc(n,sizeof(float));
    cudaMemcpyAsync(tmp,x,n*sizeof(float),cudaMemcpyDeviceToHost, get_cuda_stream());
    float o = 1.0;
    for(size_t i = 0; i < n; i++) o *= tmp[i];
    free(tmp);    
    return o;
}



__global__ void vector_addf_device(float * x, float * y, float * out, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) out[idx] = x[idx] + y[idx];
}

float* vector_addf(float * x, float * y, int n)
{
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));
    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);    
    vector_addf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);    
    return output;
}


void vector_r_addf(float * x, float * y, float * output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);    
    vector_addf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);    
}



__global__ void vector_subf_device(float * x, float * y, float * out, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) out[idx] = x[idx] - y[idx];
}

float* vector_subf(float * x, float * y, int n)
{
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));
    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);    
    vector_subf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);    
    return output;
}


void vector_r_subf(float * x, float * y, float * output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);    
    vector_subf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);    
}



__global__ void vector_mulf_device(float * x, float * y, float * out, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) out[idx] = x[idx] * y[idx];
}

float* vector_mulf(float * x, float * y, int n)
{
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));
    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);    
    vector_mulf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);    
    return output;
}

void vector_r_mulf(float * x, float * y, float * output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);    
    vector_mulf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);    
}

__global__ void vector_divf_device(float * x, float * y, float * out, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) out[idx] = x[idx] / y[idx];
}

float* vector_divf(float * x, float * y, int n)
{
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);    
    vector_divf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);    
    return output;
}
void vector_r_divf(float * x, float * y, float * output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);    
    vector_divf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);    
}

__global__ void vector_modf_device(float * x, float * y, float * out, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) out[idx] = fmodf(x[idx],y[idx]);
}

float* vector_modf(float * x, float * y, int n)
{
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_modf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);    
    return output;
}
void vector_r_modf(float * x, float * y, float * output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_modf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);    
}



__global__ void vector_acosf_device(float * in, float * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n) out[idx] = acosf(in[idx]);
}

float* vector_acosf(float * devPtr, int n)
{
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);    
    vector_acosf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n);    
    return output;
}
void vector_r_acosf(float * devPtr, float * output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);    
    vector_acosf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n);        
}

__global__ void vector_acoshf_device(float * in, float * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n) out[idx] = acoshf(in[idx]);
}

float* vector_acoshf(float * devPtr, int n)
{
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));
    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);    
    vector_acoshf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n);
    return output;
}

void vector_r_acoshf(float * devPtr, float * output, int n)
{    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);    
    vector_acoshf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n); 
}

__global__ void vector_asinhf_device(float * in, float * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = asinhf(in[idx]);
}

float* vector_asinhf(float * devPtr, int n)
{
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));
    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);    
    vector_asinhf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n);
    return output;
}
void vector_r_asinhf(float * devPtr, float * output, int n)
{    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);    
    vector_asinhf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n); 
}


__global__ void vector_asinf_device(float * in, float * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = asinf(in[idx]);
}

float* vector_asinf(float * devPtr, int n)
{
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));
    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);    
    vector_asinf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n);
    return output;
}

void vector_r_asinf(float * devPtr, float * output, int n)
{    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);    
    vector_asinf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n);
}


__global__ void vector_atan2f_device(float * a, float * b, float * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = atan2f(a[idx],b[idx]);
}

float* vector_atan2f(float * a, float * b, int n)
{
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));
    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);    
    vector_atan2f_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(a,b,output,n);
    return output;
}

void vector_r_atan2f(float * a, float * b, float * output, int n)
{    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);    
    vector_atan2f_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(a,b,output,n); 
}



__global__ void vector_atanf_device(float * in, float * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = atanf(in[idx]);
}

float* vector_atanf(float * devPtr, int n)
{
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));
    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);    
    vector_atanf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n);
    return output;
}

void vector_r_atanf(float * devPtr, float * output, int n)
{    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);    
    vector_atanf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n);    
}



__global__ void vector_atanhf_device(float * in, float * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = atanhf(in[idx]);
}

float* vector_atanhf(float * devPtr, int n)
{
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);    
    vector_atanhf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n);
    return output;
}

void vector_r_atanhf(float * devPtr, float * output, int n)
{    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);    
    vector_atanhf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n); 
}



__global__ void vector_ceilf_device(float * in, float * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = ceilf(in[idx]);
}

float* vector_ceilf(float * devPtr, int n)
{
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));   
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);    
    vector_ceilf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n);
    return output;
}

void vector_r_ceilf(float * devPtr, float * output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);    
    vector_ceilf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n);
}



__global__ void vector_cosf_device(float * in, float * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = cosf(in[idx]);
}

float* vector_cosf(float * devPtr, int n)
{
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_cosf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n);
    return output;
}
void vector_r_cosf(float * devPtr, float * output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_cosf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n);    
}


__global__ void vector_coshf_device(float * in, float * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = coshf(in[idx]);
}

float* vector_coshf(float * devPtr, int n)
{
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));
    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_coshf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n);
    return output;
}
void vector_r_coshf(float * devPtr, float * output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_coshf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n);
}



__global__ void vector_exp10f_device(float * in, float * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = exp10f(in[idx]);
}

float* vector_exp10f(float * devPtr, int n)
{
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_exp10f_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n);
    return output;
}
void vector_r_exp10f(float * devPtr, float * output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_exp10f_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n);
}


__global__ void vector_exp2f_device(float * in, float * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = exp2f(in[idx]);
}

float* vector_exp2f(float * devPtr, int n)
{
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));
    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_exp2f_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n);
    return output;
}
void vector_r_exp2f(float * devPtr, float * output, int n)
{   
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_exp2f_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n); 
}



__global__ void vector_expf_device(float * in, float * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = expf(in[idx]);
}

float* vector_expf(float * devPtr, int n)
{
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));
    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_expf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n);
    return output;
}
void vector_r_expf(float * devPtr, float * output, int n)
{   
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_expf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n); 
}

__global__ void vector_expm1f_device(float * in, float * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = expm1f(in[idx]);
}

float* vector_expm1f(float * devPtr, int n)
{
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_expm1f_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n);
    return output;
}
void vector_r_expm1f(float * devPtr, float * output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_expm1f_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n);
}


__global__ void vector_fabsf_device(float * in, float * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = fabsf(in[idx]);
}

float* vector_fabsf(float * devPtr, int n)
{
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_fabsf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n);
    return output;
}
void vector_r_fabsf(float * devPtr, float * output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_fabsf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n);
}

__global__ void vector_floorf_device(float * a,float * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = floorf(a[idx]);
}

float* vector_floorf(float * devPtr, int n)
{
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));
    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_floorf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n);
    return output;
}
void vector_r_floorf(float * devPtr, float * output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_floorf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n);    
}



__global__ void vector_fmaxf_device(float * a,float * b, float * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = fmaxf(a[idx],b[idx]);
}

float* vector_fmaxf(float * x, float * y, int n)
{
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));
    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_fmaxf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);
    return output;
}

void vector_r_fmaxf(float * x, float * y, float * output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_fmaxf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);
}


__global__ void vector_fminf_device(float * a,float * b, float * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = fminf(a[idx],b[idx]);
}

float* vector_fminf(float * x, float * y, int n)
{
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_fminf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);
    return output;
}

void vector_r_fminf(float * x, float * y, float * output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_fminf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);
}

__global__ void vector_fmodf_device(float * a,float * b, float * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = fmodf(a[idx],b[idx]);
}

float* vector_fmodf(float * x, float * y, int n)
{
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_fmodf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);
    return output;
}
void vector_r_fmodf(float * x, float * y, float * output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_fmodf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);
}

__global__ void vector_log10f_device(float * a, float * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = log10f(a[idx]);
}

float* vector_log10f(float * x, int n)
{
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_log10f_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);
    return output;
}
void vector_r_log10f(float * x, float *output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_log10f_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);
}


__global__ void vector_log1pf_device(float * a, float * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = log1pf(a[idx]);
}


float* vector_log1pf(float * x, int n)
{
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));        
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_log1pf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);
    return output;
}

void vector_r_log1pf(float * x, float *output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_log1pf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);
}


__global__ void vector_log2f_device(float * a, float * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = log2f(a[idx]);
}

float* vector_log2f(float * x, int n)
{
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_log2f_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);
    return output;
}
void vector_r_log2f(float * x, float * output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_log2f_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);
}


__global__ void vector_logbf_device(float * a, float * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = logbf(a[idx]);
}

float* vector_logbf(float * x, int n)
{
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_logbf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);
    return output;
}
void vector_r_logbf(float * x, float *output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_logbf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);
}



__global__ void vector_powf_device(float * a, float * b, float * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = powf(a[idx],b[idx]);
}

float* vector_powf(float * x, float * y, int n)
{
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_powf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);
    return output;
}
void vector_r_powf(float * x, float * y, float *output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_powf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);
}




__global__ void vector_rsqrtf_device(float * a, float * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = rsqrtf(a[idx]);
}

float* vector_rsqrtf(float * x, int n)
{
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_rsqrtf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);
    return output;
}

void vector_r_rsqrtf(float * x, float * output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_rsqrtf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);
}



__global__ void vector_sinf_device(float * a, float * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = sinf(a[idx]);
}

float* vector_sinf(float * x, int n)
{
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_sinf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);
    return output;
}
void vector_r_sinf(float * x, float *output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_sinf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);
}


__global__ void vector_sinhf_device(float * a, float * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = sinhf(a[idx]);
}

float* vector_sinhf(float * x, int n)
{
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));
    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_sinhf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);
    return output;
}
void vector_r_sinhf(float * x, float *output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_sinhf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);    
}

__global__ void vector_sqrtf_device(float * a, float * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = sqrtf(a[idx]);
}

float* vector_sqrtf(float * x, int n)
{
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_sqrtf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);
    return output;
}

void vector_r_sqrtf(float * x, float *output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_sqrtf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);
}


__global__ void vector_tanf_device(float * a, float * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = tanf(a[idx]);
}

float* vector_tanf(float * x, int n)
{
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));
    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_tanf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);
    return output;
}
void vector_r_tanf(float * x, float *output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_tanf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);
}


__global__ void vector_tanhf_device(float * a, float * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = tanhf(a[idx]);
}

float* vector_tanhf(float * x, int n)
{
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_tanhf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);
    return output;
}

void vector_r_tanhf(float * x, float *output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_tanhf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);
}




__global__ void vector_softmax_device(float * x,float *out, float sum, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) out[idx] = expf(x[idx]) / sum;
}

float* vector_softmaxf(float * x, int n)
{
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));    
    float sum = vector_sumf(x,n);
    assert(sum != 0.0);
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);
    vector_softmax_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,sum,n);        
    return output;
}
void vector_r_softmaxf(float * x, float *output, int n)
{
    float sum = vector_sumf(x,n);
    assert(sum != 0.0);
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);
    vector_softmax_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,sum,n);            
}


__global__ void vector_sigmoid_device(float * x, float * out, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) out[idx] = 1.0 / (1.0 + expf(-x[idx]));
}

float* vector_sigmoidf(float * x, int n)
{
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);            
    vector_sigmoid_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);        
    return output;
}
void vector_r_sigmoidf(float * x, float *output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);            
    vector_sigmoid_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);        
}


__global__ void vector_sigmoid_grad_device(float * x, float * out, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) out[idx] = x[idx] * (1.0 - x[idx]);
}

float* vector_sigmoid_gradf(float * x, int n)
{
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_sigmoid_grad_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);        
    return output;
}


void vector_r_sigmoid_gradf(float * x, float *output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_sigmoid_grad_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);        
}




__global__ void vector_tanh_grad_device(float * x, float * out, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) out[idx] = 1.0 - (x[idx]*x[idx]);
}

float* vector_tanh_gradf(float * x, int n)
{
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_tanh_grad_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);        
    return output;
}
void vector_r_tanh_gradf(float * x, float *output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_tanh_grad_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);        
}


__global__ void vector_relu_device(float * x, float * out, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) 
    {
        if(x[idx] < 0) out[idx] = 0.0f;
        else out[idx] = x[idx]; 
    }
}

float* vector_reluf(float * x, int n)
{
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_relu_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);        
    return output;
}

void vector_r_reluf(float * x, float *output, int n)
{    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_relu_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);        
}


__global__ void vector_relu_grad_device(float * x, float * out, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) 
    {
        if(x[idx] > 0) out[idx] = 1.0;
        else out[idx] = 0.0f;
    }
}

float* vector_relu_gradf(float * x, int n)
{
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_relu_grad_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);        
    return output;
}
void vector_r_relu_gradf(float * x, float *output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_relu_grad_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);        
}



__global__ void vector_add_const_device(float * x, float y, float * out, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) out[idx] = x[idx] + y;
}

float* vector_addf_const(float * x, float y, int n)
{
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_add_const_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);    
    return output;
}

void vector_r_addf_const(float * x, float y, float * output, int n)
{    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_add_const_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);     
}


__global__ void vector_sub_const_device(float * x, float y, float * out, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) out[idx] = x[idx] - y;
}

float* vector_subf_const(float * x, float y, int n)
{  
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_sub_const_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);        
    return output;
}


void vector_r_subf_const(float * x, float y, float *output, int n)
{  
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_sub_const_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);          
}


__global__ void vector_mul_const_device(float * x, float y, float * out, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) out[idx] = x[idx] * y;
}

float* vector_mulf_const(float * x, float y, int n)
{
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_mul_const_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);     
    return output;
}
void vector_r_mulf_const(float * x, float y, float * output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_mul_const_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);        
}

__global__ void vector_div_const_device(float * x, float y, float * out, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) out[idx] = x[idx] / y;
}

float* vector_divf_const(float * x, float y, int n)
{
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_div_const_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);        
    return output;
}
void vector_r_divf_const(float * x, float y, float *output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_div_const_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);        
}


__global__ void vector_mod_const_device(float * x, float y, float * out, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) out[idx] = fmodf(x[idx],y);
}

float* vector_modf_const(float * x, float y, int n)
{    
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_mod_const_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);        
    return output;
}

void vector_r_modf_const(float * x, float y, float *output, int n)
{    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_mod_const_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);        
}



__global__ void vector_atan2f_const_device(float * a, float b, float * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = atan2f(a[idx],b);
}

float* vector_atan2f_const(float * a, float  b, int n)
{   
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_atan2f_const_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(a,b,output,n);    
    return output;
}

void vector_r_atan2f_const(float * a, float  b, float * output, int n)
{   
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_atan2f_const_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(a,b,output,n);    
}





__global__ void vector_fmaxf_const_device(float * a,float b, float * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = fmaxf(a[idx],b);
}

float* vector_fmaxf_const(float * x, float y, int n)
{
    float * p = find_memory(n);
    if(p == NULL)
        cudaMalloc((void**)&p,sizeof(float)*n);
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_fmaxf_const_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,p,n);    
    return p;
}

void vector_r_fmaxf_const(float * x, float y, float *p, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_fmaxf_const_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,p,n);        
}


__global__ void vector_fminf_const_device(float * a,float b, float * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = fminf(a[idx],b);
}

float* vector_fminf_const(float * x, float y, int n)
{
    float * p = find_memory(n);
    if(p == NULL)
        cudaMalloc((void**)&p,sizeof(float)*n);
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_fminf_const_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,p,n);    
    return p;
}
void vector_r_fminf_const(float * x, float y, float *p, int n)
{    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_fminf_const_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,p,n);        
}

__global__ void vector_fmodf_const_device(float * a,float b, float * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = fmodf(a[idx],b);
}

float* vector_fmodf_const(float * x, float y, int n)
{
    float * p = find_memory(n);
    if(p == NULL)
        cudaMalloc((void**)&p,sizeof(float)*n);

    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_fmodf_const_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,p,n);    
    return p;
}

void vector_r_fmodf_const(float * x, float y, float * p, int n)
{    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_fmodf_const_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,p,n);       
}



__global__ void vector_powf_const_device(float * a, float b, float * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = powf(a[idx],b);
}

float* vector_powf_const(float * x, float y, int n)
{
    float * p = find_memory(n);
    if(p == NULL)
        cudaMalloc((void**)&p,sizeof(float)*n);
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_powf_const_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,p,n);    
    return p;
}
void vector_r_powf_const(float * x, float y, float *p, int n)
{    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_powf_const_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,p,n);    
}






/////////////////////////////////
// const/scalar
/////////////////////////////////
__global__ void vector_add_scalar_device(float * x, float * y, float * out, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) out[idx] = x[idx] + y[0];
}



float* vector_addf_scalar(float * x, float * y, int n)
{    
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));
    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_add_scalar_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);        
    return output;
}
void vector_r_addf_scalar(float * x, float * y, float *output, int n)
{      
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_add_scalar_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);            
}




__global__ void vector_sub_scalar_device(float * x, float * y, float * out, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) out[idx] = x[idx] - y[0];
}

float* vector_subf_scalar(float * x, float * y, int n)
{    
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_sub_scalar_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);        
    return output;
}
void vector_r_subf_scalar(float * x, float * y, float *output, int n)
{    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_sub_scalar_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);        
}


__global__ void vector_mul_scalar_device(float * x, float * y, float * out, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) out[idx] = x[idx] * y[0];
}



float* vector_mulf_scalar(float * x, float * y, int n)
{    
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_mul_scalar_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);        
    return output;
}
void vector_r_mulf_scalar(float * x, float * y, float * output, int n)
{        
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_mul_scalar_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);            
}


__global__ void vector_div_scalar_device(float * x, float * y, float * out, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N && y[0] != 0.0f) out[idx] = x[idx] / y[0];
}

float* vector_divf_scalar(float * x, float * y, int n)
{    
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_div_scalar_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);        
    return output;
}
void vector_r_divf_scalar(float * x, float * y, float *output, int n)
{    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_div_scalar_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);        
}



__global__ void vector_mod_scalar_device(float * x, float * y, float * out, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) out[idx] = fmodf(x[idx],y[0]);
}

float* vector_modf_scalar(float * x, float * y, int n)
{    
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));
    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_mod_scalar_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);        
    return output;
}

void vector_r_modf_scalar(float * x, float * y, float *output, int n)
{    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_mod_scalar_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);        
}

float* vector_fmodf_scalar(float * x, float * y, int n)
{    
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));
    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_mod_scalar_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);        
    return output;
}

void vector_r_fmodf_scalar(float * x, float * y, float *output, int n)
{    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_mod_scalar_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);        
}



__global__ void vector_fmaxf_scalar_device(float * x, float * y, float * out, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) out[idx] = fmaxf(x[idx],y[0]);
}

float* vector_fmaxf_scalar(float * x, float * y, int n)
{    
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_fmaxf_scalar_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);        
    return output;
}
void vector_r_fmaxf_scalar(float * x, float * y, float *output, int n)
{    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_fmaxf_scalar_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);        
}


__global__ void vector_fminf_scalar_device(float * x, float * y, float * out, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) out[idx] = fminf(x[idx],y[0]);
}

float* vector_fminf_scalar(float * x, float * y, int n)
{    
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_fminf_scalar_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);        
    return output;
}
void vector_r_fminf_scalar(float * x, float * y, float *output, int n)
{    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_fminf_scalar_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);        
}

__global__ void vector_powf_scalar_device(float * x, float * y, float * out, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) out[idx] = powf(x[idx],y[0]);
}

float* vector_powf_scalar(float * x, float * y, int n)
{    
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_powf_scalar_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);        
    return output;
}
void vector_r_powf_scalar(float * x, float * y, float *output, int n)
{    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_powf_scalar_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);        
}

__global__ void vector_atan2f_scalar_device(float * x, float * y, float * out, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) out[idx] = atan2f(x[idx],y[0]);
}

float* vector_atan2f_scalar(float * x, float * y, int n)
{    
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_atan2f_scalar_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);        
    return output;
}
void vector_r_atan2f_scalar(float * x, float * y, float *output, int n)
{    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_atan2f_scalar_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);        
}

__global__ void vector_fdimf_scalar_device(float * x, float * y, float * out, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) out[idx] = fdimf(x[idx],y[0]);
}

float* vector_fdimf_scalar(float * x, float * y, int n)
{    
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_fdimf_scalar_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);        
    return output;
}
void vector_r_fdimf_scalar(float * x, float * y, float *output, int n)
{    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_fdimf_scalar_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);        
}

__global__ void vector_fdividef_scalar_device(float * x, float * y, float * out, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N && y[0] != 0.0) out[idx] = fdividef(x[idx],y[0]);
}

float* vector_fdividef_scalar(float * x, float * y, int n)
{    
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_fmaxf_scalar_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);        
    return output;
}
void vector_r_fdividef_scalar(float * x, float * y, float *output, int n)
{    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_fdividef_scalar_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);        
}

__global__ void vector_remainderf_scalar_device(float * x, float * y, float * out, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) out[idx] = remainderf(x[idx],y[0]);
}

float* vector_remainderf_scalar(float * x, float * y, int n)
{    
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_remainderf_scalar_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);        
    return output;
}
void vector_r_remainderf_scalar(float * x, float * y, float *output, int n)
{    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_remainderf_scalar_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);        
}

__global__ void vector_hypot_scalar_device(float * x, float * y, float * out, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) out[idx] = hypot(x[idx],y[0]);
}

float* vector_hypotf_scalar(float * x, float * y, int n)
{    
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_hypot_scalar_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);        
    return output;
}
void vector_r_hypotf_scalar(float * x, float * y, float *output, int n)
{    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_hypot_scalar_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);        
}


__global__ void vector_rhypot_scalar_device(float * x, float * y, float * out, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) out[idx] = rhypot(x[idx],y[0]);
}

float* vector_rhypotf_scalar(float * x, float * y, int n)
{    
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_rhypot_scalar_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);        
    return output;
}
void vector_r_rhypotf_scalar(float * x, float * y, float *output, int n)
{    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_rhypot_scalar_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);        
}





__global__ void vector_setrowf_device(float * dst, float * src, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) dst[idx] = src[idx];
}

void vector_setrowf(float * dst, int dst_row, float * src, int row_src, size_t n)
{    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);    
    vector_setrowf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(dst+dst_row*n,src+row_src*n,n);        
}

void vector_copyrowf(float * dst, int dst_row, float * src, int row_src, int n) {
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);    
    vector_setrowf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(dst+dst_row*n,src+row_src*n,n);        
}

__global__ void vector_add_rowf_device(float * x,float * y, float * out, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) out[idx] = x[idx] + y[idx];
}

void vector_addf_row(float * x, int row, float * y, int row_src, size_t n)
{    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);    
    vector_add_rowf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x+row,y+row_src*n,x+row,n);        
}

__global__ void vector_sub_rowf_device(float * x, float * y, float * out, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) out[idx] = x[idx] - y[idx];
}

void vector_subf_row(float * x, int row, float * y, int row_src, size_t n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);    
    vector_sub_rowf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x+row,y+row_src,x+row,n);        
}

__global__ void vector_mul_rowf_device(float * x, float * y, float * out, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) out[idx] = x[idx] * y[idx];    
}

void vector_mulf_row(float * x, int row, float * y, int row_src, size_t n)
{    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);    
    vector_mul_rowf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x+row,y+row_src,x+row,n);            

}


__global__ void vector_div_rowf_device(float * x, float * y, float * out, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) out[idx] = x[idx] / y[idx];
}

void vector_divf_row(float * x,int row, float * y, int row_src, size_t n)
{   
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);    
    vector_div_rowf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x+row,y+row_src,x+row,n);            
}

__global__ void vector_mod_rowf_device(float * x, float * y, float * out, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < N) out[idx] = fmod(x[idx],y[idx]);
}

void vector_modf_row(float * x, int row, float * y, int row_src, size_t n)
{    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);    
    vector_mod_rowf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x+row,y+row_src,x+row,n);            
}




// vector.cu
__global__ void vector_cbrtf_device(float * in, float * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = atanf(in[idx]);
}

float* vector_cbrtf(float * devPtr, int n)
{
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);    
    vector_cbrtf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n);
    return output;
}
void vector_r_cbrtf(float * devPtr, float * output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);    
    vector_cbrtf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n);
}


__global__ void vector_copysignf_device(float * x, float * y, float * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = copysignf(x[idx],y[idx]);
}

float* vector_copysignf(float * X, float *Y, int n)
{
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_copysignf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(X,Y,output,n);
    return output;
}
void vector_r_copysignf(float * X, float *Y, float *output, int n)
{    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_copysignf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(X,Y,output,n); 
}

__global__ void vector_cospif_device(float * in, float * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = cospif(in[idx]);
}

float* vector_cospif(float * devPtr, int n)
{
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_cospif_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n);
    return output;
}
void vector_r_cospif(float * devPtr, float * output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_cospif_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n);
}


__global__ void vector_cyl_bessel_i0f_device(float * in, float * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = cyl_bessel_i0f(in[idx]);
}

float* vector_cyl_bessel_i0f(float * devPtr, int n)
{
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));
    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_cyl_bessel_i0f_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n);
    return output;
}
void vector_r_cyl_bessel_i0f(float * devPtr, float * output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_cyl_bessel_i0f_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n);
}

__global__ void vector_cyl_bessel_i1f_device(float * in, float * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = cyl_bessel_i1f(in[idx]);
}

float* vector_cyl_bessel_i1f(float * devPtr, int n)
{
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_cyl_bessel_i1f_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n);
    return output;
}
void vector_r_cyl_bessel_i1f(float * devPtr, float * output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_cyl_bessel_i1f_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n);
}

__global__ void vector_erfcf_device(float * in, float * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = erfcf(in[idx]);
}

float* vector_erfcf(float * devPtr, int n)
{
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_erfcf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n);
    return output;
}
void vector_r_erfcf(float * devPtr, float * output, int n)
{ 
    int gridSize,blockSize;    
    calcSize(n,&gridSize,&blockSize);        
    vector_erfcf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n);
}



__global__ void vector_erfcinvf_device(float * in, float * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = erfcinvf(in[idx]);
}

float* vector_erfcinvf(float * devPtr, int n)
{
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));
    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_erfcinvf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n);
    return output;
}
void vector_r_erfcinvf(float * devPtr, float * output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_erfcinvf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n);
}


__global__ void vector_erfcxf_device(float * in, float * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = erfcxf(in[idx]);
}

float* vector_erfcxf(float * devPtr, int n)
{
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_erfcxf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n);
    return output;
}
void vector_r_erfcxf(float * devPtr, float * output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_erfcxf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n);
}

__global__ void vector_erff_device(float * in, float * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = erff(in[idx]);
}

float* vector_erff(float * devPtr, int n)
{
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));
    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_erff_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n);
    return output;
}
void vector_r_erff(float * devPtr, float * output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_erff_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n);
}


__global__ void vector_erfinvf_device(float * in, float * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = erfinvf(in[idx]);
}

float* vector_erfinvf(float * devPtr, int n)
{
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_erfinvf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n);
    return output;
}
void vector_r_erfinvf(float * devPtr, float * output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_erfinvf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(devPtr,output,n);
}


__global__ void vector_fdimf_device(float * a, float * b, float * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = fdimf(a[idx],b[idx]);
}

float* vector_fdimf(float * a, float * b, int n)
{
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_fdimf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(a,b,output,n);
    return output;
}
void vector_r_fdimf(float * a, float * b, float * output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_fdimf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(a,b,output,n);
}

__global__ void vector_fdividef_device(float * a, float * b, float * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = fdividef(a[idx],b[idx]);
}

float* vector_fdividef(float * a, float * b, int n)
{
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_fdividef_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(a,b,output,n);
    return output;
}
void vector_r_fdividef(float * a, float * b, float * output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_fdividef_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(a,b,output,n);
}


__global__ void vector_fmaf_device(float * a, float * b, float * c, float * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = fmaf(a[idx],b[idx],c[idx]);
}

float* vector_fmaf(float * x, float * y, float * z, int n)
{
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_fmaf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,z,output,n);
    return output;
}
void vector_r_fmaf(float * x, float * y, float * z, float *output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_fmaf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,z,output,n);
}


__global__ void vector_hypotf_device(float * a,float * b, float * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = hypotf(a[idx],b[idx]);
}

float* vector_hypotf(float * x, float * y, int n)
{
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));
    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_hypotf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);
    return output;
}
void vector_r_hypotf(float * x, float * y, float * output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_hypotf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);
}


__global__ void vector_ilogbf_device(float * a, float * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = ilogbf(a[idx]);
}

float* vector_ilogbf(float * x, int n)
{
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_ilogbf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);
    return output;
}
void vector_r_ilogbf(float * x, float *output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_ilogbf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);
}



__global__ void vector_j0f_device(float * a, float * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = j0f(a[idx]);
}

float* vector_j0f(float * x, int n)
{
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_j0f_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);
    return output;
}
void vector_r_j0f(float * x, float *output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_j0f_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);
}



__global__ void vector_j1f_device(float * a, float * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = j1f(a[idx]);
}

float* vector_j1f(float * x, int n)
{
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_j1f_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);
    return output;
}
void vector_r_j1f(float * x, float *output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_j1f_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);
}

__global__ void vector_jnf_device(float * a, int N, float * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = jnf(a[idx],N);
}

float* vector_jnf(float * x, int M, int n)
{
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_jnf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,M,output,n);
    return output;
}
void vector_r_jnf(float * x, float * output, int M, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_jnf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,M,output,n);
}

__global__ void vector_ldexpf_device(float * a, int exp, float * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = ldexpf(a[idx],exp);
}

float* vector_ldexpf(float * x, int exp, int n)
{
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_ldexpf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,exp,output,n);
    return output;
}
void vector_r_ldexpf(float * x, float * output, int exp, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_ldexpf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,exp,output,n);
}



__global__ void vector_lgammaf_device(float * a, float * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = lgammaf(a[idx]);
}

float* vector_lgammaf(float * x, int n)
{
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_lgammaf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);
    return output;
}
void vector_r_lgammaf(float * x, float *output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_lgammaf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);
}



__global__ void vector_llrintf_device(float * a, long long * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = llrintf(a[idx]);
}

long long* vector_llrintf(float * x, int n)
{
    long long * p;
    cudaMalloc((void**)&p,sizeof(long long)*n);
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_llrintf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,p,n);
    return p;
}

__global__ void vector_llroundf_device(float * a, long long * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = llroundf(a[idx]);
}

long long* vector_llroundf(float * x, int n)
{
    long long * p;
    cudaMalloc((void**)&p,sizeof(long long)*n);
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_llroundf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,p,n);
    return p;
}

__global__ void vector_lrintf_device(float * a, long * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = lrintf(a[idx]);
}

long * vector_lrintf(float * x, int n)
{
    long * p;
    cudaMalloc((void**)&p,sizeof(long )*n);
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_lrintf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,p,n);
    return p;
}

__global__ void vector_lroundf_device(float * a, long * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = lroundf(a[idx]);
}

long * vector_lroundf(float * x, int n)
{
    long * p;
    cudaMalloc((void**)&p,sizeof(long )*n);
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_lroundf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,p,n);
    return p;
}

__global__ void vector_nearbyintf_device(float * a, float * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = lroundf(a[idx]);
}

float* vector_nearbyintf(float * x, int n)
{
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_nearbyintf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);
    return output;
}
void vector_r_nearbyintf(float * x, float *output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_nearbyintf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);
}



__global__ void vector_norm3df_device(float * a, float * b, float * c, float* out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = norm3df(a[idx],b[idx],c[idx]);
}

float* vector_norm3df(float * x, float * y, float * z, int n)
{
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));   
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_norm3df_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,z,output,n);
    return output;
}
void vector_r_norm3df(float * x, float * y, float * z, float * output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_norm3df_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,z,output,n);
}


__global__ void vector_norm4df_device(float * a, float * b, float * c, float * d, float* out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = norm4df(a[idx],b[idx],c[idx],d[idx]);
}

float* vector_norm4df(float * x, float * y, float * z, float * q, int n)
{
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_norm4df_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,z,q,output,n);
    return output;
}
void vector_r_norm4df(float * x, float * y, float * z, float * q, float * output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_norm4df_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,z,q,output,n);
}


__global__ void vector_normcdff_device(float * a, float * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = normcdff(a[idx]);
}

float* vector_normcdff(float * x, int n)
{
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_normcdff_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);
    return output;
}
void vector_r_normcdff(float * x, float * output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_normcdff_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);    
}


__global__ void vector_normcdfinvf_device(float * a, float * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = normcdfinvf(a[idx]);
}

float* vector_normcdfinvf(float * x, int n)
{
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_normcdfinvf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);
    return output;
}
void vector_r_normcdfinvf(float * x, float *output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_normcdfinvf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);
}



__global__ void vector_normf_device(int dim, const float * a, float * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = normf(dim,a);
}

float* vector_normf(int dim, float * x, int n)
{
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));   
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_normf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(dim,x,output,n);
    return output;
}
void vector_r_normf(int dim, float * x, float * output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_normf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(dim,x,output,n);
}

__global__ void vector_rcbrtf_device(float * a, float * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = rcbrtf(a[idx]);
}

float* vector_rcbrtf(float * x, int n)
{
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_rcbrtf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);
    return output;
}
void vector_r_rcbrtf(float * x, float *output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_rcbrtf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);
}



__global__ void vector_remainderf_device(float * a, float * b, float * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = remainderf(a[idx],b[idx]);
}


float* vector_remainderf(float * x, float * y, int n)
{
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_remainderf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);
    return output;
}
void vector_r_remainderf(float * x, float * y, float *output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_remainderf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);
}


__global__ void vector_rhypotf_device(float * a, float * b, float * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = rhypotf(a[idx],b[idx]);
}

float* vector_rhypotf(float * x, float * y, int n)
{
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_rhypotf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n);
    return output;
}
void vector_r_rhypotf(float * x, float * y, float *output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_rhypotf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,output,n); 
}





__global__ void vector_rnorm3df_device(float * a, float * b, float * c, float* out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = rnorm3df(a[idx],b[idx],c[idx]);
}

float* vector_rnorm3df(float * x, float * y, float * z, int n)
{
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_rnorm3df_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,z,output,n);
    return output;
}
void vector_r_rnorm3df(float * x, float * y, float * z, float * output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_rnorm3df_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,z,output,n);
}

__global__ void vector_rnorm4df_device(float * a, float * b, float * c, float * d, float* out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = rnorm4df(a[idx],b[idx],c[idx],d[idx]);
}


float* vector_rnorm4df(float * x, float * y, float * z, float * q, int n)
{
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_rnorm4df_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,z,q,output,n);
    return output;
}
void vector_r_rnorm4df(float * x, float * y, float * z, float * q, float *output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_rnorm4df_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,z,q,output,n);
}


__global__ void vector_rnormf_device(int dim, const float * a, float * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = rnormf(dim,a);
}

float* vector_rnormf(int dim, float * x, int n)
{
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_rnormf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(dim,x,output,n);
    return output;
}
void vector_r_rnormf(int dim, float * x, float *output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_rnormf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(dim,x,output,n);
}


__global__ void vector_scalblnf_device(float * a, long int N,float * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = scalblnf(a[idx],N);
}


float* vector_scalblnf(float * x, long int M, int n)
{
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_scalblnf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,M,output,n);
    return output;
}
void vector_r_scalblnf(float * x, long int M, float * output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_scalblnf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,M,output,n);
}


__global__ void vector_sinpif_device(float * a, float * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = sinpif(a[idx]);
}

float* vector_sinpif(float * x, int n)
{
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_sinpif_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);
    return output;
}
void vector_r_sinpif(float * x, float *output, int n)
{    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_sinpif_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n); 
}



__global__ void vector_tgammaf_device(float * a, float * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = tgammaf(a[idx]);
}

float* vector_tgammaf(float * x, int n)
{
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_tgammaf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);
    return output;
}
void vector_r_tgammaf(float * x, float * output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_tgammaf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);
}




__global__ void vector_truncf_device(float * a, float * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = truncf(a[idx]);
}

float* vector_truncf(float * x, int n)
{
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_truncf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);
    return output;
}

void vector_r_truncf(float * x, float *output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_truncf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);
}


__global__ void vector_y0f_device(float * a, float * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = y0f(a[idx]);
}

float* vector_y0f(float * x, int n)
{
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_y0f_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);
    return output;
}
void vector_r_y0f(float * x, float *output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_y0f_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);
}




__global__ void vector_y1f_device(float * a, float * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = y1f(a[idx]);
}

float* vector_y1f(float * x, int n)
{
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_y1f_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);
    return output;
}
void vector_r_y1f(float * x, float * output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_y1f_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,output,n);
}

__global__ void vector_ynf_device(int N, float * a, float * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = ynf(n,a[idx]);
}

float* vector_ynf(int M, float * x, int n)
{
    float * output = find_memory(n);
    if(output == NULL) cudaMalloc((void**)&output,n*sizeof(float));
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_ynf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(M,x,output,n);
    return output;
}
void vector_r_ynf(int M, float * x, float *output, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_ynf_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(M,x,output,n);
}


__global__ void vector_fdimf_const_device(float * a, float b, float * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = fdimf(a[idx],b);
}

float* vector_fdimf_const(float * a, float  b, int n)
{
    float * p = find_memory(n);    
    if(p == NULL)
        cudaMalloc((void**)&p,sizeof(float)*n);
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_fdimf_const_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(a,b,p,n);    
    return p;
}
void vector_r_fdimf_const(float * a, float  b, float * output, int n)
{
    float * p = find_memory(n);        
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_fdimf_const_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(a,b,output,n);        
}


__global__ void vector_fdividef_const_device(float * a, float b, float * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = fdividef(a[idx],b);
}

float* vector_fdividef_const(float * a, float b, int n)
{
    float * p = find_memory(n);
    if(p == NULL)
        cudaMalloc((void**)&p,sizeof(float)*n);
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_fdividef_const_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(a,b,p,n);    
    return p;
}
void vector_r_fdividef_const(float * a, float b, float *p, int n)
{    
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_fdividef_const_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(a,b,p,n);     
}


__global__ void vector_hypotf_const_device(float * a,float b, float * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = hypotf(a[idx],b);
}

float* vector_hypotf_const(float * x, float y, int n)
{
    float * p = find_memory(n);
    if(p == NULL)
        cudaMalloc((void**)&p,sizeof(float)*n);
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_hypotf_const_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,p,n);    
    return p;
}
void vector_r_hypotf_const(float * x, float y, float *p, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_hypotf_const_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,p,n);    
}



__global__ void vector_remainderf_const_device(float * a, float b, float * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = remainderf(a[idx],b);
}


float* vector_remainderf_const(float * x, float y, int n)
{
    float * p = find_memory(n);
    if(p == NULL)
        cudaMalloc((void**)&p,sizeof(float)*n);
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_remainderf_const_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,p,n);    
    return p;
}
void vector_r_remainderf_const(float * x, float y, float *p, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_remainderf_const_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,p,n);        
}

__global__ void vector_rhypotf_const_device(float * a, float b, float * out, int n)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;     
    if(idx < n)
        out[idx] = rhypotf(a[idx],b);
}

float* vector_rhypotf_const(float * x, float y, int n)
{
    float * p = find_memory(n);
    if(p == NULL)
        cudaMalloc((void**)&p,sizeof(float)*n);
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_rhypotf_const_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,p,n);    
    return p;
}
void vector_r_rhypotf_const(float * x, float y, float *p, int n)
{
    int gridSize,blockSize;
    calcSize(n,&gridSize,&blockSize);        
    vector_rhypotf_const_device<<<gridSize,blockSize,0,get_cuda_stream()>>>(x,y,p,n);    
}
