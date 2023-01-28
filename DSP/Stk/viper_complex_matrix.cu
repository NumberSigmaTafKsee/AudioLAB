//////////////////////////////////////////////////////////////////////////
// matrix
//////////////////////////////////////////////////////////////////////////
#include <cassert>
#include "cuda_runtime.h"
#include "math_constants.h"
#include "viper_vector.hpp"
#include "matrix_complex.h"
#include "vector_complex.h"


cudaStream_t get_cuda_stream();


template void cmatrix_r_realf(complex<float> * x, float * y, int m,int n);
template void cmatrix_r_imagf(complex<float> * x, float * y, int m,int n);

template void cmatrix_r_argf(complex<float> * devPtr, float * x, int m,int n);
template void cmatrix_r_normf(complex<float> * devPtr, complex<float> * x, int m,int n);
template void cmatrix_r_conjf(complex<float> * devPtr, complex<float> * x, int m,int n);
template void cmatrix_r_projf(complex<float> * devPtr, complex<float> * x,int m,int n);
//template complex<float>* cmatrix_r_polar(float * r, float *theta,  int m,int n);

template void cmatrix_r_hadamardf( complex<float> * a, complex<float> * b, complex<float> *output, int M, int N, int K);
template void cmatrix_r_acosf(complex<float> * devPtr, complex<float> * outputs, int m,int n);
template void cmatrix_r_asinf(complex<float> * devPtr, complex<float> * outputs, int m,int n);
template void cmatrix_r_atanf(complex<float> * devPtr, complex<float> * outputs, int m,int n);
template void cmatrix_r_acoshf(complex<float> * devPtr, complex<float> * outputs, int m,int n);
template void cmatrix_r_asinhf(complex<float> * devPtr, complex<float> * outputs, int m,int n);
template void cmatrix_r_atanhf(complex<float> * devPtr, complex<float> * outputs, int m,int n);
template void cmatrix_r_cosf(complex<float> * devPtr, complex<float> * outputs, int m,int n);
template void cmatrix_r_sinf(complex<float> * devPtr, complex<float> * outputs, int m,int n);
template void cmatrix_r_tanf(complex<float> * devPtr, complex<float> * outputs, int m,int n);
template void cmatrix_r_coshf(complex<float> * devPtr, complex<float> * outputs, int m,int n);
template void cmatrix_r_sinhf(complex<float> * devPtr, complex<float> * outputs, int m,int n);
template void cmatrix_r_tanhf(complex<float> * devPtr, complex<float> * outputs, int m,int n);

template void cmatrix_2d_r_mulf(complex<float> * a, complex<float> * b, complex<float> *output, int M, int N);
template void cmatrix_2d_r_addf(complex<float> * a, complex<float> * b, complex<float> *output,int M, int N);
template void cmatrix_2d_r_subf(complex<float> * a, complex<float> * b, complex<float> *output,int M, int N);
template void cmatrix_2d_r_divf(complex<float> * a, complex<float> * b, complex<float> *output,int M, int N);


template void cmatrix_r_expf(complex<float> * devPtr, complex<float> * output, int m,int n);
template void cmatrix_r_logf(complex<float> * x, complex<float> *output, int m,int n);
template void cmatrix_r_log10f(complex<float> * x, complex<float> *output, int m,int n);
template void cmatrix_r_powf(complex<float> * x, complex<float> * y, complex<float> *output, int m,int n);
template void cmatrix_r_sqrtf(complex<float> * x, complex<float> *output, int m,int n);
template void cmatrix_r_sigmoidf(complex<float> * x, complex<float> *output, int m,int n);
template void cmatrix_r_sigmoid_gradf(complex<float> * x, complex<float> *output, int m,int n);
template void cmatrix_r_tanh_gradf(complex<float> * x, complex<float> *output, int m,int n);
template void cmatrix_r_addf_const(complex<float> * x, complex<float>  y, complex<float> *output, int m,int n);
template void cmatrix_r_subf_const(complex<float> * x, complex<float>  y, complex<float> *output,int m,int n);
template void cmatrix_r_mulf_const(complex<float> * x, complex<float>  y, complex<float> *output,int m,int n);
template void cmatrix_r_divf_const(complex<float> * x, complex<float>  y, complex<float> *output,int m,int n);
template void cmatrix_r_powf_const(complex<float> * x, complex<float> y, complex<float> *output,int m,int n);
template void cmatrix_r_addf_scalar(complex<float> * x, complex<float> * y, complex<float> *output,int m,int n);
template void cmatrix_r_subf_scalar(complex<float> * x, complex<float> * y, complex<float> *output,int m,int n);
template void cmatrix_r_mulf_scalar(complex<float> * x, complex<float> * y, complex<float> *output,int m,int n);
template void cmatrix_r_divf_scalar(complex<float> * x, complex<float> * y, complex<float> *output,int m,int n);
template void cmatrix_r_powf_scalar(complex<float> * x, complex<float> *y, complex<float> *output,int m,int n);

template void cmatrix_r_realf(complex<double> * x, double * y, int m,int n);
template void cmatrix_r_imagf(complex<double> * x, double * y, int m,int n);

template void cmatrix_r_argf(complex<double> * devPtr, double * x, int m,int n);
template void cmatrix_r_normf(complex<double> * devPtr, complex<double> * x, int m,int n);
template void cmatrix_r_conjf(complex<double> * devPtr, complex<double> * x, int m,int n);
template void cmatrix_r_projf(complex<double> * devPtr, complex<double> * x,int m,int n);
//template complex<double>* cmatrix_polar(double * r, double *theta,  int m,int n);

template void cmatrix_r_hadamardf( complex<double> * a, complex<double> * b, complex<double> *output, int M, int N, int K);
template<typename T>
void matrix_r_multiplyf(T * a, T * b, T * output, int M, int N, int K);
template<typename T>
void matrix_r_transposef(T * input, T *output, int M, int N);

template void cmatrix_2d_r_mulf(complex<double> * a, complex<double> * b, complex<double> *output, int M, int N);
template void cmatrix_2d_r_addf(complex<double> * a, complex<double> * b, complex<double> *output,int M, int N);
template void cmatrix_2d_r_subf(complex<double> * a, complex<double> * b, complex<double> *output,int M, int N);
template void cmatrix_2d_r_divf(complex<double> * a, complex<double> * b, complex<double> *output,int M, int N);


template void cmatrix_r_acosf(complex<double> * devPtr, complex<double> * outputs, int m,int n);
template void cmatrix_r_asinf(complex<double> * devPtr, complex<double> * outputs, int m,int n);
template void cmatrix_r_atanf(complex<double> * devPtr, complex<double> * outputs, int m,int n);
template void cmatrix_r_acoshf(complex<double> * devPtr, complex<double> * outputs, int m,int n);
template void cmatrix_r_asinhf(complex<double> * devPtr, complex<double> * outputs, int m,int n);
template void cmatrix_r_atanhf(complex<double> * devPtr, complex<double> * outputs, int m,int n);
template void cmatrix_r_cosf(complex<double> * devPtr, complex<double> * outputs, int m,int n);
template void cmatrix_r_sinf(complex<double> * devPtr, complex<double> * outputs, int m,int n);
template void cmatrix_r_tanf(complex<double> * devPtr, complex<double> * outputs, int m,int n);
template void cmatrix_r_coshf(complex<double> * devPtr, complex<double> * outputs, int m,int n);
template void cmatrix_r_sinhf(complex<double> * devPtr, complex<double> * outputs, int m,int n);
template void cmatrix_r_tanhf(complex<double> * devPtr, complex<double> * outputs, int m,int n);

template void cmatrix_r_expf(complex<double> * devPtr, complex<double> * output, int m,int n);
template void cmatrix_r_logf(complex<double> * x, complex<double> *output, int m,int n);
template void cmatrix_r_log10f(complex<double> * x, complex<double> *output, int m,int n);
template void cmatrix_r_powf(complex<double> * x, complex<double> * y, complex<double> *output, int m,int n);
template void cmatrix_r_sqrtf(complex<double> * x, complex<double> *output, int m,int n);
template void cmatrix_r_sigmoidf(complex<double> * x, complex<double> *output, int m,int n);
template void cmatrix_r_sigmoid_gradf(complex<double> * x, complex<double> *output, int m,int n);
template void cmatrix_r_tanh_gradf(complex<double> * x, complex<double> *output, int m,int n);
template void cmatrix_r_addf_const(complex<double> * x, complex<double>  y, complex<double> *output, int m,int n);
template void cmatrix_r_subf_const(complex<double> * x, complex<double>  y, complex<double> *output,int m,int n);
template void cmatrix_r_mulf_const(complex<double> * x, complex<double>  y, complex<double> *output,int m,int n);
template void cmatrix_r_divf_const(complex<double> * x, complex<double>  y, complex<double> *output,int m,int n);
template void cmatrix_r_powf_const(complex<double> * x, complex<double> y, complex<double> *output,int m,int n);
template void cmatrix_r_addf_scalar(complex<double> * x, complex<double> * y, complex<double> *output,int m,int n);
template void cmatrix_r_subf_scalar(complex<double> * x, complex<double> * y, complex<double> *output,int m,int n);
template void cmatrix_r_mulf_scalar(complex<double> * x, complex<double> * y, complex<double> *output,int m,int n);
template void cmatrix_r_divf_scalar(complex<double> * x, complex<double> * y, complex<double> *output,int m,int n);
template void cmatrix_r_powf_scalar(complex<double> * x, complex<double> *y, complex<double> *output,int m,int n);


static const int BLOCK_SIZE=16;

template<typename T> __global__ void cmatrix_gpu_2d_addf(complex<T> *a, complex<T> *b, complex<T> *c, int m, int n)
{ 
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;    
    if( col < n && row < m) 
    {        
        c[row * n + col] = a[row * n + col] + b[row*n + col];
    }
} 

template<typename T> void cmatrix_2d_r_addf(complex<T> * a, complex<T> * b, complex<T> *output, int M, int N)
{
       
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);
    cmatrix_gpu_2d_addf<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(a,b,output,M,N);
}

template<typename T> __global__ void cmatrix_gpu_2d_mulf(complex<T> *a, complex<T> *b, complex<T> *c, int m, int n)
{ 
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;    
    if( col < n && row < m) 
    {        
        c[row * n + col] = a[row * n + col] * b[row*n + col];
    }
} 


template<typename T> void cmatrix_2d_r_mulf( complex<T> * a, complex<T> * b, complex<T> *output, int M, int N)
{
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);
    cmatrix_gpu_2d_mulf<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(a,b,output,M,N);
}


template<typename T> __global__ void cmatrix_gpu_2d_subf(complex<T> *a, complex<T> *b, complex<T> *c, int m, int n)
{ 
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;    
    if( col < n && row < m) 
    {        
        c[row * n + col] = a[row * n + col] - b[row*n + col];
    }
} 


template<typename T>
void cmatrix_2d_r_subf( complex<T> * a, complex<T> * b, complex<T> *output, int M, int N)
{
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);
    cmatrix_gpu_2d_subf<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(a,b,output,M,N);
}

template<typename T> __global__ void cmatrix_gpu_2d_divf( complex<T> *a, complex<T> *b, complex<T> *c, int m, int n)
{ 
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;    
    if( col < n && row < m) 
    {        
        c[row * n + col] = a[row * n + col] / b[row*n + col];
    }
} 


template<typename T>
void cmatrix_2d_r_divf( complex<T> * a, complex<T> * b, complex<T> *output, int M, int N)
{
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);
    cmatrix_gpu_2d_divf<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(a,b,output,M,N);
}



template<typename T> __global__ void cmatrix_gpu_matrix_addf( complex<T> *a, complex<T> *b, complex<T> *c, int m, int n, int k)
{ 
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;    
    if( col < k && row < m) 
    {        
        c[row * n + col] = a[row * n + col] + b[row*n + col];
    }
} 

template<typename T> void cmatrix_r_addf( complex<T> * a, complex<T> * b, complex<T> *output, int M, int N, int K)
{
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (K + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);
    gpu_matrix_addf<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(a,b,output,M,N,K);
}

template<typename T> __global__ void cmatrix_gpu_matrix_subf( complex<T> *a, complex<T> *b, complex<T> *c, int m, int n, int k)
{ 
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;    
    if( col < k && row < m) 
    {                
        c[row * n + col] = a[row * n + col] - b[row*n + col];
    }
} 


template<typename T> void cmatrix_r_subf(complex<T> * a, complex<T> * b, complex<T> *output, int M, int N, int K)
{  
    
   // if(N < 1024) BLOCK_SIZE=N;
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (K + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);
    gpu_matrix_subf<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(a,b,output,M,N,K); 
}




template<typename T>
__global__ void gpu_cmatrix_realf(complex<T> *a ,T *r, int m, int n)
{ 
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;    
    if( col < n && row < m) 
    {        
        r[row * n + col] = a[row * n + col].real();
    }
} 


template<typename T>
void cmatrix_r_realf(complex<T> * a, T *output, int M, int N)
{
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);
    gpu_cmatrix_realf<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(a,output,M,N);
}


template<typename T>
__global__ void gpu_cmatrix_imagf(complex<T> *a ,T *r, int m, int n)
{ 
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;    
    if( col < n && row < m) 
    {        
        r[row * n + col] = a[row * n + col].imag();
    }
} 


template<typename T>
void cmatrix_r_imagf(complex<T> * a, T *output, int M, int N)
{
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);
    gpu_cmatrix_imagf<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(a,output,M,N);
}



template<typename T>
__global__ void gpu_cmatrix_hadamardf(complex<T> *a,complex<T> *b, complex<T> *c, int m, int n, int k)
{ 
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;    
    if( col < k && row < m) 
    {                
        c[row * n + col] = a[row * n + col] * b[row * k + col];
    }
} 

template<typename T>
void cmatrix_r_hadamardf(complex<T> * a, complex<T> * b, complex<T> *output, int M, int N, int K)
{   
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (K + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);
    gpu_cmatrix_hadamardf<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(a,b,output,M,N,K);
}

template<typename T>
__global__ void gpu_cmatrix_multiplyf(complex<T> *a,complex<T> *b, complex<T> *c, int m, int n, int k)
{ 
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    complex<T> sum = 0;    
    
    if( col < k && row < m) 
    {     
        for(int i = 0; i < n; i++) 
        {            
            sum += a[row*n+i]*b[i*k+col];
        }        
        c[row*n + col] = sum;
    }

} 




template<typename T>
void cmatrix_r_multiplyf(complex<T> * a, complex<T> * b, complex<T> * output, int M, int N, int K)
{
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (K + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_rows,grid_cols);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_cmatrix_multiplyf<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(a,b,output,M,N,K);
}


template<typename T>
__global__ void gpu_cmatrix_transposef(complex<T>* mat_in, complex<T>* mat_out, unsigned int rows, unsigned int cols) 
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


template<typename T>
void cmatrix_r_transposef(complex<T> * input, complex<T> *output, int M, int N)
{
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);
    gpu_cmatrix_transposef<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(input,output,M,N);
}


template<typename T>
__global__ void gpu_cmatrix_acosf(complex<T>* x, complex<T>* c, int m, int n) 
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (col < n && row < m) 
    {        
        c[row*n + col] = acos(x[row*n + col]);
    }
}



template<typename T>
void cmatrix_r_acosf(complex<T> * input, complex<T> *output, int M, int N)
{
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_cmatrix_acosf<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(input,output,M,N);
}


template<typename T>
__global__ void gpu_cmatrix_asinf(complex<T>* x, complex<T>* c, int m, int n) 
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (col < n && row < m) 
    {        
        c[row*n + col] = asin(x[row*n + col]);
    }
}



template<typename T>
void cmatrix_r_asinf(complex<T> * input, complex<T> * output,  int M, int N)
{    
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_cmatrix_asinf<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(input,output,M,N);
}

template<typename T>
__global__ void gpu_cmatrix_atanf(complex<T>* x, complex<T>* c, int m, int n) 
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (col < n && row < m) 
    {        
        c[row*n + col] = atan(x[row*n + col]);
    }
}

template<typename T>
void cmatrix_r_atanf(complex<T> * input, complex<T> * output, int M, int N)
{    
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_cmatrix_atanf<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(input,output,M,N);    
}




template<typename T>
__global__ void gpu_cmatrix_cosf(complex<T>* x, complex<T>* c, int m, int n) 
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (col < n && row < m) 
    {        
        c[row*n + col] = cos(x[row*n + col]);
    }
}



template<typename T>
void cmatrix_r_cosf(complex<T> * input, complex<T> * output, int M, int N)
{   
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_cmatrix_cosf<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(input,output,M,N); 
}


template<typename T>
__global__ void gpu_cmatrix_sinf(complex<T>* x, complex<T>* c, int m, int n) 
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (col < n && row < m) 
    {        
        c[row*n + col] = sin(x[row*n + col]);
    }
}



template<typename T>
void cmatrix_r_sinf(complex<T> * input, complex<T> * output, int M, int N)
{   
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_cmatrix_sinf<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(input,output,M,N); 
}


template<typename T>
__global__ void gpu_cmatrix_tanf(complex<T>* x, complex<T>* c, int m, int n) 
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (col < n && row < m) 
    {        
        c[row*n + col] = tan(x[row*n + col]);
    }
}


template<typename T>
void cmatrix_r_tanf(complex<T> * input, complex<T> * output, int M, int N)
{
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_cmatrix_tanf<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(input,output,M,N);
}


template<typename T>
__global__ void gpu_cmatrix_coshf(complex<T>* x, complex<T>* c, int m, int n) 
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (col < n && row < m) 
    {        
        c[row*n + col] = cosh(x[row*n + col]);
    }
}




template<typename T>
void cmatrix_r_coshf(complex<T> * input, complex<T> * output, int M, int N)
{    
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_cmatrix_coshf<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(input,output,M,N); 
}


template<typename T>
__global__ void gpu_cmatrix_sinhf(complex<T>* x, complex<T>* c, int m, int n) 
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (col < n && row < m) 
    {        
        c[row*n + col] = sinh(x[row*n + col]);
    }
}



template<typename T>
void cmatrix_r_sinhf(complex<T> * input, complex<T> *output, int M, int N)
{
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_cmatrix_sinhf<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(input,output,M,N);
}


template<typename T>
__global__ void gpu_cmatrix_tanhf(complex<T>* x, complex<T>* c, int m, int n) 
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (col < n && row < m) 
    {        
        c[row*n + col] = tanh(x[row*n + col]);
    }
}



template<typename T>
void cmatrix_r_tanhf(complex<T> * input, complex<T> * output, int M, int N)
{
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_cmatrix_tanhf<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(input,output,M,N);
}

template<typename T>
__global__ void gpu_cmatrix_acoshf(complex<T>* x, complex<T>* c, int m, int n) 
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (col < n && row < m) 
    {        
        c[row*n + col] = acosh(x[row*n + col]);
    }
}



template<typename T>
void cmatrix_r_acoshf(complex<T> * input, complex<T> * output, int M, int N)
{    
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_cmatrix_acoshf<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(input,output,M,N); 
}


template<typename T>
__global__ void gpu_cmatrix_asinhf(complex<T>* x, complex<T>* c, int m, int n) 
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (col < n && row < m) 
    {        
        c[row*n + col] = asinh(x[row*n + col]);
    }
}



template<typename T>
void cmatrix_r_asinhf(complex<T> * input, complex<T> * output, int M, int N)
{   
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_cmatrix_asinhf<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(input,output,M,N); 
}


template<typename T>
__global__ void gpu_cmatrix_atanhf(complex<T>* x, complex<T>* c, int m, int n) 
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (col < n && row < m) 
    {        
        c[row*n + col] = atanh(x[row*n + col]);
    }
}



template<typename T>
void cmatrix_r_atanhf(complex<T> * input, complex<T> *output, int M, int N)
{  
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_cmatrix_atanhf<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(input,output,M,N); 
}





template<typename T>
__global__ void gpu_cmatrix_expf(complex<T>* x, complex<T>* c, int m, int n) 
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    int idx = row*n + col;
    if (col < n && row < m)     {        
        c[idx] = exp(x[idx]);
    }    
}



template<typename T>
void cmatrix_r_expf(complex<T> * input, complex<T> *output, int M, int N)
{
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_cmatrix_expf<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(input,output,M,N);
}


template<typename T>
__global__ void gpu_cmatrix_logf(complex<T>* x, complex<T>* c, int m, int n) 
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    int idx = row*n + col;
    if (col < n && row < m)     {        
        c[idx] = log(x[idx]);
    }    
}



template<typename T>
void cmatrix_r_logf(complex<T> * input, complex<T> *output, int M, int N)
{
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_cmatrix_logf<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(input,output,M,N);
}


template<typename T>
__global__ void gpu_cmatrix_log10f(complex<T>* x, complex<T>* c, int m, int n) 
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    int idx = row*n + col;
    if (col < n && row < m)     {        
        c[idx] = log10(x[idx]);
    }    
}



template<typename T>
void cmatrix_r_log10f(complex<T> * input, complex<T> *output, int M, int N)
{
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_cmatrix_log10f<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(input,output,M,N);
}


template<typename T>
__global__ void gpu_cmatrix_powf(complex<T>* x, complex<T> * y, complex<T>* c, int m, int n) 
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (col < n && row < m) 
    {        
        c[row*n + col] = pow(x[row*n + col],y[row*n + col]);
    }
}


template<typename T>
void cmatrix_r_powf(complex<T> * x, complex<T> *y, complex<T> *output, int M, int N)
{
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_cmatrix_powf<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(x,y,output,M,N);
}


template<typename T>
__global__ void gpu_cmatrix_sqrtf(complex<T>* x, complex<T>* c, int m, int n) 
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (col < n && row < m) 
    {        
        c[row*n + col] = sqrt(x[row*n + col]);
    }
}



template<typename T>
void cmatrix_r_sqrtf(complex<T> * input, complex<T> *output, int M, int N)
{
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_cmatrix_sqrtf<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(input,output,M,N);
}


template<typename T>
__global__ void cmatrix_sigmoid_device(complex<T> * x, complex<T> * out, int m, int n)
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    int idx = row*n + col; 
    if (col < n && row < m) 
    {        
        complex<T> r = complex<T>(1.0,0.0) / ( complex<T>(1.0,0.0) + exp(-x[idx]));
        out[idx] = r;
    }
}



template<typename T>
void cmatrix_r_sigmoidf(complex<T> * x, complex<T> *output, int M, int N)
{
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);       
    cmatrix_sigmoid_device<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(x,output,M,N);        
}


template<typename T>
__global__ void cmatrix_sigmoid_grad_device(complex<T> * x, complex<T> * out, int m, int n)
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    int idx = row*n + col; 
    if (col < n && row < m) 
    {        
        complex<T> r = x[idx] * (complex<T>(1.0,0.0) - x[idx]);
        out[idx] = r;
    }    
}



template<typename T>
void cmatrix_r_sigmoid_gradf(complex<T> * x, complex<T> *output,int M, int N)
{
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    cmatrix_sigmoid_grad_device<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(x,output,M,N);        
}

template<typename T>
__global__ void cmatrix_tanh_grad_device(complex<T> * x, complex<T> * out, int m, int n)
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    int idx = row*n + col; 
    if (col < n && row < m) 
    {        
        out[idx] = complex<T>(1.0,0.0) - (x[idx] *x[idx]);
    }    
}


template<typename T>
void cmatrix_r_tanh_gradf(complex<T> * x, complex<T> *output, int M, int N)
{
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    cmatrix_tanh_grad_device<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(x,output,M,N);        
}





template<typename T>
__global__ void cmatrix_add_const_device(complex<T> * x, complex<T> y, complex<T> * out, int m, int n)
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
void cmatrix_r_addf_const(complex<T> * x, complex<T> y, complex<T> *output, int M, int N)
{    
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    cmatrix_add_const_device<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(x,y,output,M,N);     
}


template<typename T>
__global__ void cmatrix_sub_const_device(complex<T> * x, complex<T> y, complex<T> * out, int m, int n)
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
void cmatrix_r_subf_const(complex<T> * x, complex<T> y, complex<T> *output, int M, int N)
{
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    cmatrix_sub_const_device<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(x,y,output,M,N);    
}


template<typename T>
__global__ void cmatrix_mul_const_device(complex<T> * x, complex<T> y, complex<T> * out, int m, int n)
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
void cmatrix_r_mulf_const(complex<T> * x, complex<T> y, complex<T> *output,int M, int N)
{
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    cmatrix_mul_const_device<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(x,y,output,M,N);    
}

template<typename T>
__global__ void cmatrix_div_const_device(complex<T> * x, complex<T> y, complex<T> * out, int m, int n)
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    int idx = row*n + col; 
    if (col < n && row < m) 
    {        
        out[idx] = x[idx] / y;
    }
}


template<typename T>
void cmatrix_r_divf_const(complex<T> * x, complex<T> y, complex<T> *output, int M, int N)
{    
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    cmatrix_div_const_device<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(x,y,output,M,N);    
}


template<typename T>
__global__ void cmatrix_powf_const_device(complex<T> * x, complex<T> y, complex<T> * out, int m, int n)
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    int idx = row*n + col; 
    if (col < n && row < m) 
    {        
        out[idx] = pow(x[idx],y);
    }
}


template<typename T>
void cmatrix_r_powf_const(complex<T> * x, complex<T> y, complex<T> *output, int M, int N)
{
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    cmatrix_powf_const_device<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(x,y,output,M,N);    
}


template<typename T>
__global__ void cmatrix_add_scalar_device(complex<T> * x, complex<T> *y, complex<T> * out, int m, int n)
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
void cmatrix_r_addf_scalar(complex<T> * x, complex<T> *y, complex<T> *output, int M, int N)
{
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    cmatrix_add_scalar_device<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(x,y,output,M,N);    
}


template<typename T>
__global__ void cmatrix_subf_scalar_device(complex<T> * x, complex<T> *y, complex<T> * out, int m, int n)
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
void cmatrix_r_subf_scalar(complex<T> * x, complex<T> *y, complex<T> *output, int M, int N)
{
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    cmatrix_add_scalar_device<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(x,y,output,M,N);    
}


template<typename T>
__global__ void cmatrix_mul_scalar_device(complex<T> * x, complex<T> *y, complex<T> * out, int m, int n)
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
void cmatrix_r_mulf_scalar(complex<T> * x, complex<T> *y, complex<T> *output, int M, int N)
{
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    cmatrix_mul_scalar_device<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(x,y,output,M,N);    
}


template<typename T>
__global__ void cmatrix_div_scalar_device(complex<T> * x, complex<T> *y, complex<T> * out, int m, int n)
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    int idx = row*n + col; 
    if (col < n && row < m) 
    {        
        out[idx] = x[idx] / y[0];        
    }
}


template<typename T>
void cmatrix_r_divf_scalar(complex<T> * x, complex<T> *y, complex<T> *output, int M, int N)
{
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    cmatrix_div_scalar_device<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(x,y,output,M,N);    
}



template<typename T>
__global__ void cmatrix_pow_scalar_device(complex<T> * x, complex<T> *y, complex<T> * out, int m, int n)
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    int idx = row*n + col; 
    if (col < n && row < m) 
    {            
        out[idx] = pow(x[idx],y[0]);        
    }
}


template<typename T>
void cmatrix_r_powf_scalar(complex<T> * x, complex<T> *y, complex<T> *output, int M, int N)
{
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    cmatrix_pow_scalar_device<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(x,y,output,M,N);    
}


template<typename T>
__global__ void gpu_cmatrix_normf(const complex<T>* x, complex<T>* c, int m, int n) 
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (col < n && row < m) 
    {        
        c[row*n + col] = norm(x[row*n+col]);
    }
}


template<typename T>
void cmatrix_r_normf(complex<T> *input, complex<T> *output, int M, int N)
{
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_cmatrix_normf<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(input,output,M,N);
}




template<typename T> __global__ void gpu_cmatrix_argf_device(complex<T> * a, T * out, int m, int n)
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (col < n && row < m) 
    {        
        out[row*n + col] = arg(a[row*n + col]);
    }       
}

template<typename T> 
void cmatrix_r_argf(complex<T> * devPtr,  T* x, int M, int N)
{
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_cmatrix_argf_device<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(devPtr,x,M,N);

}

template<typename T> __global__ void gpu_cmatrix_conjf_device(complex<T> * a, complex<T> * out, int m, int n)
{
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (col < n && row < m) 
    {        
        out[row*n + col] = conj(a[row*n + col]);
    }       
}

template<typename T>
void cmatrix_r_conjf(complex<T> * devPtr, complex<T> * x,int M, int N)
{
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_cmatrix_conjf_device<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(devPtr,x,M,N);

}



template<typename T> __global__ void gpu_cmatrix_projf_device(complex<T> * a, complex<T> * out, int m, int n)
{        
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (col < n && row < m) 
    {        
        out[row*n + col] = proj(a[row*n + col]);
    }       
}

template<typename T>
void cmatrix_r_projf(complex<T> * devPtr, complex<T> * x, int M, int N)
{
    unsigned int grid_rows = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);    
    gpu_cmatrix_projf_device<<<dimGrid,dimBlock,0,get_cuda_stream()>>>(devPtr,x,M,N);

}