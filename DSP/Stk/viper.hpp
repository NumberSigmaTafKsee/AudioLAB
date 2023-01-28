#pragma once

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <cmath>
#include <vector>
#include <random>
#include <chrono>
#include <memory>
#include <map>
#include <complex>
#include <stdexcept>
#include <algorithm>

#include <pthread.h>

#include "cuda_runtime.h"
#include "cublas_v2.h"
#include "curand.h"

// libcu++
// cuda::std::sample_vector

#include "viper_vector.hpp"
#include "viper_matrix.hpp"



static const char *cublasGetErrorString(cublasStatus_t error)
{
    switch (error)
    {
        case CUBLAS_STATUS_SUCCESS:
            return "CUBLAS_STATUS_SUCCESS";

        case CUBLAS_STATUS_NOT_INITIALIZED:
            return "CUBLAS_STATUS_NOT_INITIALIZED";

        case CUBLAS_STATUS_ALLOC_FAILED:
            return "CUBLAS_STATUS_ALLOC_FAILED";

        case CUBLAS_STATUS_INVALID_VALUE:
            return "CUBLAS_STATUS_INVALID_VALUE";

        case CUBLAS_STATUS_ARCH_MISMATCH:
            return "CUBLAS_STATUS_ARCH_MISMATCH";

        case CUBLAS_STATUS_MAPPING_ERROR:
            return "CUBLAS_STATUS_MAPPING_ERROR";

        case CUBLAS_STATUS_EXECUTION_FAILED:
            return "CUBLAS_STATUS_EXECUTION_FAILED";

        case CUBLAS_STATUS_INTERNAL_ERROR:
            return "CUBLAS_STATUS_INTERNAL_ERROR";
    }

    return "<unknown>";
}
//#define ASSERTOK(status) if(status != cudaSuccess) throw std::runtime_error(cublasGetErrorString(status));        
//(assert(status == CUBLAS_STATUS_SUCCESS))
template<typename T> void ASSERTOK(T x) {

}
template<> void ASSERTOK<cublasStatus_t>(cublasStatus_t status)
{
    if(status != CUBLAS_STATUS_SUCCESS)
    throw std::runtime_error(cublasGetErrorString(status));        
}
template<> void ASSERTOK<cudaError_t>(cudaError_t status)
{
    if(status != cudaSuccess)
        throw std::runtime_error(cudaGetErrorString(status));        
}

typedef int array_index;

pthread_mutex_t mutex;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// cuda/cublas
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
struct CudaStream
{
    cudaStream_t stream;

    CudaStream() {
        cudaStreamCreate(&stream);
    }
    ~CudaStream() {
        cudaStreamDestroy(stream);
    }
    void reset() {
        cudaStreamCreate(&stream);
    }
};



struct CublasPointerMode 
{
    cublasPointerMode_t pointer_mode;
};

struct CublasAtomicsMode 
{
    cublasAtomicsMode_t atomics_mode;
};

struct CublasMathMode 
{
    cublasMath_t math_mode;
};



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Cublas 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
struct Cublas
{
    cublasStatus_t   status;
    cublasHandle_t   handle;

    Cublas()
    {        
        status = cublasCreate(&handle);
        ASSERTOK(status);
        cudaSetDeviceFlags(cudaDeviceMapHost);
        pthread_mutexattr_t attr;
        pthread_mutexattr_init(&attr);
        pthread_mutex_init(&mutex,&attr);        
    }
    ~Cublas() 
    {
        if(handle) cublasDestroy(handle);
    }    

    void reset() {
        status = cublasCreate(&handle);
        ASSERTOK(status);
        cudaSetDeviceFlags(cudaDeviceMapHost);
    }
    int GetVersion() 
    {
        int v = -1;
        status = cublasGetVersion(handle,&v);
        ASSERTOK(status);
        return v;
    }

    const char* GetStatusName()
    {
        const char * r  = cublasGetStatusName(status);
        return r;
    }
    
    void SetWorkspace(void * workspace, size_t workspace_size)
    {
        status = cublasSetWorkspace(handle,workspace,workspace_size);
        ASSERTOK(status);
    }

    void SetStream(const CudaStream& stream)
    {
        status = cublasSetStream(handle,stream.stream);
        ASSERTOK(status);
    }

    void GetStream(CudaStream & stream)
    {
        status = cublasGetStream(handle,&stream.stream);
        ASSERTOK(status);
    }

    void SetPointerMode(CublasPointerMode &p)
    {
        status = cublasSetPointerMode(handle,p.pointer_mode);
        ASSERTOK(status);
    }
    void GetPointerMode(CublasPointerMode & p)
    {
        status = cublasGetPointerMode(handle,&p.pointer_mode);
        ASSERTOK(status);
    }
    void SetAtomicsMode(CublasAtomicsMode & a)
    {
        status = cublasSetAtomicsMode(handle,a.atomics_mode);
        ASSERTOK(status);
    }
    void GetAtomicsMode(CublasAtomicsMode & a)
    {
        status = cublasGetAtomicsMode(handle,&a.atomics_mode);
        ASSERTOK(status);
    }
    void SetMathMode(CublasMathMode & m)
    {
        status = cublasSetMathMode(handle,m.math_mode);
        ASSERTOK(status);
    }
    void GetMathMode(CublasMathMode & m)
    {
        status = cublasGetMathMode(handle,&m.math_mode);
        ASSERTOK(status);
    }
    void SetSmCountTarget(int countTarget)
    {
        status = cublasSetSmCountTarget(handle,countTarget);
        ASSERTOK(status);
    }
    int GetSmCountTarget()
    {
        int sm = -1;
        status = cublasGetSmCountTarget(handle,&sm);
        ASSERTOK(status);
        return sm;
    }

    void LoggerConfigure(int logIsOn,int logToStdOut, int logToStdErr, const char * filename)
    {
        status = cublasLoggerConfigure(logIsOn, logToStdOut, logToStdErr, filename);
    }
};


struct CuRand 
{
    // curand 
    curandGenerator_t gen;

    CuRand(unsigned long long seed=0) {
        if(seed == 0)
            curandSetPseudoRandomGeneratorSeed(gen,time(NULL));
        else 
            curandSetPseudoRandomGeneratorSeed(gen,seed);
        curandCreateGenerator(&gen,CURAND_RNG_PSEUDO_DEFAULT);
    }
    ~CuRand() {
        curandDestroyGenerator(gen);
    }

    template<typename T>
    void curand_uniform(T * dev, size_t size) {                
        curandGenerateUniform(gen,dev,size*sizeof(T));                
    }
    template<typename T>
    void curand_normal(T * dev, size_t size, T mean, T stddev) {                
        curandGenerateNormal(gen,dev,size*sizeof(T), mean, stddev);                
    }
    template<typename T>
    void curand_lognormal(T *dev, size_t size, T mean, T stddev) {                
        curandGenerateLogNormal(gen,dev,size*sizeof(T), mean, stddev);                
    }    
};

extern Cublas *cublas;
cudaStream_t get_cuda_stream();
void set_stream(int streamid);
int get_stream();
cudaStream_t random_stream();


int current_stream = 0;
CudaStream cuda_streams[16];

cudaStream_t get_cuda_stream() {
    return cuda_streams[current_stream].stream;
}


void set_stream(int streamid) {
    assert(streamid >= 0 && streamid < 16);
    current_stream = streamid;    
    cublas->SetStream(cuda_streams[current_stream]);    
}
int get_stream() { 
    return current_stream; 
}


void Memcpy(void * dst, const void *src, size_t count, enum cudaMemcpyKind kind)
{
    //if(kind == cudaMemcpyHostToHost) memcpy(dst,src,count);
    //assert(kind != cudaMemcpyHostToHost);
    cudaMemcpyAsync(dst,src,count,kind,get_cuda_stream());
    //cudaMemcpy(dst,src,count,kind);
} 


cudaStream_t random_stream() {
    return cuda_streams[rand() % 16].stream;
}

struct SchemaMemory 
{
    void * ptr;
    size_t length;

    SchemaMemory() {
        ptr = NULL;
        length = 0;
    }
    SchemaMemory(size_t l) {
        ptr = calloc(length=l,sizeof(uint8_t));
        assert(ptr != NULL);
    }
    ~SchemaMemory() {
        if(ptr) free(ptr);
    }

    void allocate(size_t s) {
        if(ptr) free(ptr);
        ptr = calloc(s,sizeof(uint8_t));
        length = s;
    }

    SchemaMemory& operator = (const SchemaMemory & s) {
        if(ptr) free(ptr);
        length = s.length;
        ptr = calloc(length,sizeof(uint8_t));
        memcpy(ptr,s.ptr,length);
        return *this;
    }

    void zero() {
        assert(ptr != NULL);
        memset(ptr,0x00,length);
    }
    void fillui8(uint8_t v) {
        assert(ptr != NULL);
        for(size_t i = 0; i < length; i++)
            ((uint8_t*)ptr)[i] = v;
    }
    void fillui16(uint16_t v) {
        assert(ptr != NULL);
        for(size_t i = 0; i < length/sizeof(uint16_t); i++)
            ((uint16_t*)ptr)[i] = v;
    }
    void fillui32(uint32_t v) {
        assert(ptr != NULL);
        for(size_t i = 0; i < length/sizeof(uint32_t); i++)
            ((uint32_t*)ptr)[i] = v;
    }
    void fillui64(uint64_t v) {
        assert(ptr != NULL);
        for(size_t i = 0; i < length/sizeof(uint64_t); i++)
            ((uint64_t*)ptr)[i] = v;
    }
    void filli8(int8_t v) {
        assert(ptr != NULL);
        for(size_t i = 0; i < length; i++)
            ((int8_t*)ptr)[i] = v;
    }
    void filli16(int16_t v) {
        assert(ptr != NULL);
        for(size_t i = 0; i < length/sizeof(int16_t); i++)
            ((int16_t*)ptr)[i] = v;
    }
    void filli32(int32_t v) {
        assert(ptr != NULL);
        for(size_t i = 0; i < length/sizeof(int32_t); i++)
            ((int32_t*)ptr)[i] = v;
    }
    void filli64(int64_t v) {
        assert(ptr != NULL);
        for(size_t i = 0; i < length/sizeof(int64_t); i++)
            ((int64_t*)ptr)[i] = v;
    }
    void fillf32(float v) {
        assert(ptr != NULL);
        for(size_t i = 0; i < length/sizeof(float); i++)
            ((float*)ptr)[i] = v;
    }
    void fillf64(double v) {
        assert(ptr != NULL);
        for(size_t i = 0; i < length/sizeof(double); i++)
            ((double*)ptr)[i] = v;
    }
    /*
    void fill(T _Complex v) {
        assert(ptr != NULL);
        for(size_t i = 0; i < length/sizeof(T _Complex); i++)
            ((T _Complex*)ptr)[i] = v;
    }
    void fill(double _Complex v) {
        assert(ptr != NULL);
        for(size_t i = 0; i < length/sizeof(double _Complex); i++)
            ((double _Complex*)ptr)[i] = v;
    }
    */

    float&      f32(size_t pos) { return ((float*)ptr)[pos]; }
    double&     f64(size_t pos) { return ((double*)ptr)[pos]; }
    int8_t&     i8(size_t pos) { return ((int8_t*)ptr)[pos]; }
    uint8_t&    ui8(size_t pos) { return ((uint8_t*)ptr)[pos]; }
    int16_t&    i16(size_t pos) { return ((int16_t*)ptr)[pos]; }
    uint16_t&   ui16(size_t pos) { return ((uint16_t*)ptr)[pos]; }
    int32_t&    i32(size_t pos) { return ((int32_t*)ptr)[pos]; }
    uint32_t&   ui32(size_t pos) { return ((uint32_t*)ptr)[pos]; }
    int64_t&    i64(size_t pos) { return ((int64_t*)ptr)[pos]; }
    uint64_t&   ui64(size_t pos) { return ((uint64_t*)ptr)[pos]; }
    //T _Complex& cf32(size_t pos) { return ((T _Complex*)ptr)[pos]; }
    //double _Complex& cf64(size_t pos) { return ((double _Complex*)ptr)[pos]; }


    void flip(size_t bit) {
        size_t pos = bit/8;
        size_t m   = bit%8;
        assert(ptr != NULL);
        uint8_t *p = (uint8_t*)ptr;
        uint8_t bits = p[pos] & (1 << m);
        if(bits) p[pos] &= ~m;
        else p[pos] |= m;
    }
    size_t num_bits() const { return length*8; }
    void set(size_t bit) {
        size_t pos = bit/8;
        size_t m   = bit%8;
        assert(ptr != NULL);
        uint8_t *p = (uint8_t*)ptr;
        p[pos] |= 1 << m;        
    }
    int8_t get(size_t bit) {
        size_t pos = bit/8;
        size_t m   = bit%8;
        assert(ptr != NULL);
        uint8_t *p = (uint8_t*)ptr;
        return p[pos] & (1 << m);
    }
    void set_range(size_t start, size_t end) {
        for(size_t i = start; i < end; i++)
            set(i);
    }
};

// memory cache
std::multimap<int,std::pair<void*,void*>> host_memory;


void clear_memory()
{
    pthread_mutex_lock(&mutex);
    
    for(auto i = host_memory.begin(); i != host_memory.end(); i++)
    {
        std::pair<void*,void*> t = i->second;
        //cudaFree(t.first);
        cudaFree(t.second);
    }
    host_memory.clear();
    
    //cudaDeviceReset();    
    //host_memory.clear();
    //cublas->reset();
    //for(size_t i = 0; i < 16; i++) cuda_streams[i].reset();
    pthread_mutex_unlock(&mutex);    
}

void return_host(int length, const std::pair<void*,void*> &fp) {        
    pthread_mutex_lock(&mutex);
    host_memory.insert(std::pair<int,std::pair<void*,void*>>(length,fp));    
    pthread_mutex_unlock(&mutex);
}


std::pair<void*,void*> find_host(int length) {   
    pthread_mutex_lock(&mutex);                
    typename std::multimap<int,std::pair<void*,void*>>::iterator i = host_memory.find(length);
    
    if(i == host_memory.end()) 
    {
        i = std::find_if(host_memory.begin(),host_memory.end(),[&](std::pair<int,std::pair<void*,void*>> p){ return p.first  >= length;});
        if(i == host_memory.end()) {
            pthread_mutex_unlock(&mutex);
            return std::pair<void*,void*>(nullptr,nullptr);                
        }
    }
    //cuda_zero(i->second.first,length);
    //cudaMemsetAsync(i->second.second,0,length,get_cuda_stream());
    host_memory.erase(i);
    pthread_mutex_unlock(&mutex);
    return (i->second);
}


template<typename T>
struct DevPtr 
{
    void * dev;        
    void * host;
    int    len;
    

    DevPtr(size_t s, size_t e = sizeof(T)) {                 
        host = NULL;
        len   = s*e;
        std::pair<void*,void*> p = find_host(len);        
        dev = (p.first);
        host= (p.second);
        if(dev == nullptr) {            
            cudaError_t err;        
            err = cudaHostAlloc(&host,len,cudaHostAllocMapped);         
            if(err != cudaSuccess) throw std::runtime_error(cudaGetErrorString(err));        
            assert(host != NULL);                        
            err = cudaHostGetDevicePointer((void**)&dev,(void*)host,0);
            if(err != cudaSuccess) throw std::runtime_error(cudaGetErrorString(err));                        
            //err = cudaMalloc(&dev,n*sizeof(T));
            //if(err != cudaSuccess) throw std::runtime_error(cudaGetErrorString(err));        
            assert(dev != NULL);                                    
        }                            
        zero();
    }
    ~DevPtr() {                         
        return_host(len,std::pair<void*,void*>(dev,host));        
    }
    size_t size() const { return len; }

    void zero() { 
        memset(host,0,size());  
    }
    void copy(DevPtr * p) {
        assert(p != nullptr);
        assert(p->size() == size());
        assert(host != nullptr);
        assert(p->host != nullptr);
        Memcpy(host,p->host,size(),cudaMemcpyHostToHost);        
    }

    int8_t& int8(size_t i) { return ((int8_t*)host)[i]; }
    uint8_t& uint8(size_t i) { return ((uint8_t*)host)[i]; }
    int16_t& int16(size_t i) { return ((int16_t*)host)[i]; }
    uint16_t& uint16(size_t i) { return ((uint16_t*)host)[i]; }
    int32_t& int32(size_t i) { return ((int32_t*)host)[i]; }
    uint32_t& uint32(size_t i) { return ((uint32_t*)host)[i]; }
    int64_t& int64(size_t i) { return ((int64_t*)host)[i]; }
    uint64_t& uint64(size_t i) { return ((uint64_t*)host)[i]; }
    T& f32(size_t i) { return ((T*)host)[i]; }
    double& f64(size_t i) { return ((double*)host)[i]; }

    // bitset
    // bitget
    // upload_device
    // download_host
    // bitflip
    // fill
    // copy
    // sort
    
    
};


template<typename T>
inline void rand(T * devPtr, T *host, int size, T min, T max) {
    typedef std::chrono::high_resolution_clock myclock;
    myclock::time_point beginning = myclock::now();
    myclock::duration d = myclock::now() - beginning;
    unsigned seed = d.count();
    std::default_random_engine generator(seed);
    std::uniform_real_distribution<double> distribution(min,max);                
    for(size_t i = 0; i < size; i++) host[i] = distribution(generator);
    Memcpy(devPtr,host,size*sizeof(T),cudaMemcpyHostToDevice);                
}


namespace Viper
{
    template<typename T>
    struct Vector
    {
    public:

        void init(size_t n) {        
            N    = n;
            dev_ptr  = { new DevPtr<T>(N), [](DevPtr<T>* p){ delete p; } };
            dev  = (T*)dev_ptr.get()->dev;
            host = (T*)dev_ptr.get()->host;            
            //zero();
        }

        T * dev;
        T * host;            
        std::shared_ptr<DevPtr<T>> dev_ptr;
        size_t   N;
        
        
    public:

        Vector() {
            dev = host = nullptr;
            N = 0;    
        }
        Vector(size_t n) {
            init(n);        
        }
        Vector(int n, const std::initializer_list<T> & input) {
            init(n);
            std::vector<T> tmp(input.begin(),input.end());
            Memcpy(host,tmp.data(),n*sizeof(T),cudaMemcpyHostToHost);
        }
        Vector(int n, const std::vector<T> & tmp) {
            init(n);        
            Memcpy(host,tmp.data(),n*sizeof(T),cudaMemcpyHostToHost);
        }    
        Vector(const Vector<T> & v) {
            *this = v;
            //init(v.N);
            //copy(v);
        }        
        Vector(T * p, size_t n) {
            N = n;
            dev_ptr  = { new DevPtr<T>(N), [](DevPtr<T>* p){ delete p; } };
            dev  = (T*)dev_ptr.get()->dev;
            host = (T*)dev_ptr.get()->host;
            zero();
        }
        ~Vector() = default;

        Vector& operator = (const Vector<T> & v) {
            if(dev == v.dev) return *this;                
            dev_ptr.reset();
            dev_ptr  = v.dev_ptr;
            N = v.N;
            host = v.host;            
            dev  = v.dev;
            return *this;
        }

        size_t size() const { return N; }    
        
        void download_host()     {                  
            Memcpy(host,dev,size()*sizeof(T),cudaMemcpyDeviceToHost);            
        }
        void upload_device()     {
            Memcpy(dev,host,size()*sizeof(T),cudaMemcpyHostToDevice);            
        }    
        void zero()    {            
            dev_ptr->zero();
        }
        void ones()    {
            fill(1.0f);
        }    
        void randu() {
            rand<T>(dev,host,size(),0.0f,1.0f);
        }        
        void random(T min, T max) {
            rand(dev,host,size(),min,max);
        }    
        void fill(const T val)     {
            for(size_t i = 0; i < size(); i++) host[i] = val;            
        }

        T& operator[](array_index pos)    {        
            assert(pos < size());
            while(pos < 0) pos += size();
            return host[pos];        
        }
        T __getitem__(array_index pos) { return (*this)[pos]; }
        void __setitem__(array_index pos, T val) { while(pos < 0) pos += N; host[pos] = val; }

        Vector<T> operator - ();
        Vector<T> operator + (const Vector<T> & a);
        Vector<T> operator - (const Vector<T> & a);
        Vector<T> operator * (const Vector<T> & a);

        Vector<T> operator + (const T a);
        Vector<T> operator - (const T a);
        Vector<T> operator * (const T a);
        Vector<T> operator / (const T a);
        
        virtual void copy(const Vector<T> & v) {
            resize(v.size());
            dev_ptr.get()->copy(v.dev_ptr.get());
        }
        void resize(size_t n) {
            N = n;
            dev_ptr.reset();        
            init(n);
        }
                
        void print()  {
            // if you use the index operator or modify it on the CPU make sure to synchronize before using it...
            download_host();  
            std::cout << "vector[" << N << "]" << std::endl;
            for(size_t w = 0; w < size(); w++) {
                std::cout << host[w] << ",";
            }
            std::cout << std::endl;        
        }

        Vector<T> clone(const Vector<T> & v) { return Vector<T>(v); }
        Vector<T> eval() { return Vector<T>(*this); }

        void abs() {
            vector_r_fabsf(dev,dev,N);
        }
        void exp() {
            vector_r_expf(dev,dev,N);
        }    
        void log2() {
            vector_r_log2f(dev,dev,N);
        }        
        void log10() {
            vector_r_log10f(dev,dev,N);
        }
        void pow(T p) {
            vector_r_powf_const(dev,p,dev,N);
        }
        void pow(const Vector<T> & p) {
            vector_r_powf(dev,p.dev,dev,N);
        }
        void sqrt() {
            vector_r_sqrtf(dev,dev,N);
        }
        void rsqrt() {
            vector_r_rsqrtf(dev,dev,N);
        }        
        void sin() {
            vector_r_sinf(dev,dev,N);        
        }
        void cos() {
            vector_r_cosf(dev,dev,N);        
        }
        void tan() {
            vector_r_tanf(dev,dev,N);        
        }
        void asin() {
            vector_r_asinf(dev,dev,N);        
        }
        void acos() {
            vector_r_acosf(dev,dev,N);        
        }
        void atan() {
            vector_r_atanf(dev,dev,N);        
        }

        void sinh() {
            vector_r_sinhf(dev,dev,N);        
        }
        void cosh() {
            vector_r_coshf(dev,dev,N);        
        }
        void tanh() {
            vector_r_tanhf(dev,dev,N);        
        }
        void asinh() {
            vector_r_asinhf(dev,dev,N);        
        }
        void acosh() {
            vector_r_acoshf(dev,dev,N);        
        }
        void atanh() {
            vector_r_atanhf(dev,dev,N);        
        }    
        void atan2(T v) {
            vector_r_atan2f_const(dev,v,dev,N);        
        }

        
        void atan2(const Vector<T> & v) {
            vector_r_atan2f(dev,v.dev,dev,N);        
        }
        void sigmoid() {
            vector_r_sigmoidf(dev,dev,N);
        }
        void sigmoid_deriv() {
            vector_r_sigmoid_gradf(dev,dev,N);
        }
        void relu() {
            vector_r_reluf(dev,dev,N);
        }
        void relu_deriv() {
            vector_r_relu_gradf(dev,dev,N);
        }
        void softmax() {
            vector_r_softmaxf(dev,dev,N);
        }
        void tanh_deriv() {
            vector_r_tanh_gradf(dev,dev,N);
        }

        
        void cbrt() {
            vector_r_cbrtf(dev,dev,N);        
        }
        void ceil() {
            vector_r_ceilf(dev,dev,N);        
        }
        /*
        void copysign(const T val) {
            vector_r_copysignf_const(dev,val,dev,N);        
        }
        void copysign(const Vector<T>& val) {
            vector_r_copysignf_const(dev,val.dev,dev,N);        
        }*/
        void cospi() {
            vector_r_cospif(dev,dev,N);        
        }
        void cyl_bessel_i0() {
            vector_r_cyl_bessel_i0f(dev,dev,N);
        }
        void cyl_bessel_i1() {
            vector_r_cyl_bessel_i1f(dev,dev,N);
        }
        
        void erfc() {
            vector_r_erfcf(dev,dev,N);
        }
        void erfcx() {
            vector_r_erfcxf(dev,dev,N);
        }
        void erfcinv() {
            vector_r_erfcinvf(dev,dev,N);
        }
        void erf() {
            vector_r_erff(dev,dev,N);
        }
        void erfinv() {
            vector_r_erfinvf(dev,dev,N);
        }
        void exp10() {
            vector_r_exp10f(dev,dev,N);
        }
        
        void exp2() {
            vector_r_exp2f(dev,dev,N);
        }
        void expm1() {
            vector_r_expm1f(dev,dev,N);
        }
        void fabs() {
            vector_r_fabsf(dev,dev,N);
        }    
        void fdim(const Vector<T> & b) {
            vector_r_fdimf(dev,b.dev,dev,N);
        }
        void fmod(const Vector<T> & b) {
            vector_r_fmodf(dev,b.dev,dev,N);
        }
        void hypot(const Vector<T> & b) {
            vector_r_hypotf(dev,b.dev,dev,N);
        }
        void ilogb() {
            vector_r_ilogbf(dev,dev,N);
        }
        void j0() {
            vector_r_j0f(dev,dev,N);
        }

        
        void j1() {
            vector_r_j1f(dev,dev,N);
        }
        void jn(int n) {
            vector_r_jnf(dev,dev,N,n);
        }
        void lgamma() {
            vector_r_lgammaf(dev,dev,N);
        }
        void log1p() {
            vector_r_log1pf(dev,dev,N);
        }
        void logb() {
            vector_r_logbf(dev,dev,N);
        }
        void norm3d(const Vector<T> & a, const Vector<T> & b, const Vector<T> & c) {
            vector_r_norm3df(a.dev,b.dev,c.dev,dev,N);
        }
        void norm4d(const Vector<T> & a, const Vector<T> & b, const Vector<T> & c, const Vector<T> & d) {
            vector_r_norm4df(a.dev,b.dev,c.dev,d.dev,dev,N);
        }
        
        void normcdf() {
            vector_r_normcdff(dev,dev,N);
        }
        void normcdfinv() {
            vector_r_normcdfinvf(dev,dev,N);
        }
        void norm() {
            vector_r_normf(1,dev,dev,N);
        }
        void rcbrt() {
            vector_r_rcbrtf(dev,dev,N);        
        }
        
        void remainder(const Vector<T> & b) {
            vector_r_remainderf(dev,b.dev,dev,N);
        }
        void rhypot(const Vector<T> & b) {
            vector_r_rhypotf(dev,b.dev,dev,N);
        }
        /*
        void rint() {
            vector_r_rintf(dev,dev,N);        
        }
        */
        void rnorm3d(const Vector<T> & a, const Vector<T> & b, const Vector<T> & c) {
            vector_r_rnorm3df(a.dev,b.dev,c.dev,dev,N);
        }
        void rnorm4d(const Vector<T> & a, const Vector<T> & b, const Vector<T> & c, const Vector<T> & d) {
            vector_r_rnorm4df(a.dev,b.dev,c.dev,d.dev,dev,N);
        }
        /*
        void round() {
            vector_r_roundf(dev,dev,N);        
        }
        */
        
        void rnorm() {
            vector_r_rnormf(1,dev,dev,N);
        }
        void tgamma() {
            vector_r_tgammaf(dev,dev,N);
        }
        void trunc() {
            vector_r_truncf(dev,dev,N);
        }
        void y0() {
            vector_r_y0f(dev,dev,N);
        }
        void y1() {
            vector_r_y1f(dev,dev,N);
        }
        void yn(int n) {
            vector_r_ynf(n,dev,dev,N);
        }
    };

    int max_index(const Vector<float> &v) {
        int r=0;
        cublasStatus_t     status; 
        status = cublasIsamax(cublas->handle,v.N,v.dev,1,&r);
        ASSERTOK(status);
        return r-1;
    }
    int min_index(const Vector<float> & v) {
        int r=0;
        cublasStatus_t     status; 
            status = cublasIsamin(cublas->handle,v.N,v.dev,1,&r);
            ASSERTOK(status);
        return r-1;
    }    
    float nrm2(const Vector<float> & v) {
        float r = 0;
        cublasStatus_t     status; 
        status = cublasSnrm2(cublas->handle,v.N,v.dev,1,&r);
        ASSERTOK(status);
        return r;
    }
    Vector<float> scale(const Vector<float> & v, float x) {
        Vector<float> r(v.N);
        r.copy(v);
        cublasStatus_t     status; 
        status = cublasSscal(cublas->handle,v.N,&x,r.dev,1);
        ASSERTOK(status);
        return r;
    }
    void swap(Vector<float> & a, Vector<float> & b) {
        cublasStatus_t     status; 
        status = cublasSswap(cublas->handle,a.N,a.dev,1,b.dev,1);
        ASSERTOK(status);
    }
    int max_index(const Vector<double> &v) {
        int r=0;
        cublasStatus_t     status; 
        status = cublasIdamax(cublas->handle,v.N,v.dev,1,&r);
        ASSERTOK(status);
        return r-1;
    }
    int min_index(const Vector<double> & v) {
        int r=0;
        cublasStatus_t     status; 
            status = cublasIdamin(cublas->handle,v.N,v.dev,1,&r);
            ASSERTOK(status);
        return r-1;
    }    
    double nrm2(const Vector<double> & v) {
        double r = 0;
        cublasStatus_t     status; 
        status = cublasDnrm2(cublas->handle,v.N,v.dev,1,&r);
        ASSERTOK(status);
        return r;
    }
    Vector<double> scale(const Vector<double> & v, double x) {
        Vector<double> r(v.N);
        r.copy(v);
        cublasStatus_t     status; 
        status = cublasDscal(cublas->handle,v.N,&x,r.dev,1);
        ASSERTOK(status);
        return r;
    }
    void swap(Vector<double> & a, Vector<double> & b) {
        cublasStatus_t     status; 
        status = cublasDswap(cublas->handle,a.N,a.dev,1,b.dev,1);
        ASSERTOK(status);
    }

    void copy(const Vector<float> & a, Vector<float> & r) {       
        if(a.N != r.N) r.resize(a.N);
        cublasStatus_t     status; 
        status = cublasScopy(cublas->handle,a.N,a.dev,1,r.dev,1);
        ASSERTOK(status);            
    }
    void copy(const Vector<double> & a, Vector<double> & r) {       
        if(a.N != r.N) r.resize(a.N);
        cublasStatus_t     status; 
        status = cublasDcopy(cublas->handle,a.N,a.dev,1,r.dev,1);
        ASSERTOK(status);            
    }

    template<typename T>
    Vector<T> shift(const Vector<T> &a, int count)
    {
        Vector<T> r(a);
        count = count % a.N;
        if(count < 0)
        {
            for(size_t i = count; i < 0; i++) {
                T x = r[0];
                
                for(size_t j = 0; j < r.size()-1; j++)
                    r[j] = r[j+1];
                r[j+1] = x;
            }
        }
        else
        {
            for(size_t i = 0; i < count; i++) {
                T x = r[r.size()-1];
                
                for(size_t j = r.size()-2; j >0; j++)
                    r[j+1] = r[j];
                r[0] = x;
            }
        }
    }

    // I want this to be Vector<complex<T>>
    template<typename Type>
    struct ComplexVector
    {
        using T = complex<Type>;

        void init(size_t n) {        
            N    = n;
            dev_ptr  = { new DevPtr<T>(N), [](DevPtr<T>* p){ delete p; } };
            dev  = (T*)dev_ptr.get()->dev;
            host = (T*)dev_ptr.get()->host;            
            //zero();
        }

        T * dev;
        T * host;            
        std::shared_ptr<DevPtr<T>> dev_ptr;
        size_t   N;
        
        
    public:

        ComplexVector() {
            dev = host = nullptr;
            N = 0;    
        }
        ComplexVector(size_t n) {
            init(n);        
        }
        ComplexVector(int n, const std::initializer_list<T> & input) {
            init(n);
            std::vector<T> tmp(input.begin(),input.end());
            Memcpy(host,tmp.data(),n*sizeof(T),cudaMemcpyHostToHost);
        }
        ComplexVector(int n, const std::vector<T> & tmp) {
            init(n);        
            Memcpy(host,tmp.data(),n*sizeof(T),cudaMemcpyHostToHost);
        }    
        ComplexVector(const ComplexVector<Type> & v) {            
            init(v.N);
            copy(v);
        }        
        ComplexVector(T * p, size_t n) {
            N = n;
            dev_ptr  = { new DevPtr<T>(N), [](DevPtr<T>* p){ delete p; } };
            dev  = (T*)dev_ptr.get()->dev;
            host = (T*)dev_ptr.get()->host;
            //zero();
        }
        ~ComplexVector() = default;

        ComplexVector<Type>& operator = (const ComplexVector<Type> & v) {
            if(dev == v.dev) return *this;                
            dev_ptr.reset();
            dev_ptr  = v.dev_ptr;
            N = v.N;
            host = v.host;
            dev  = v.dev;
            return *this;
        }

        size_t size() const { return N; }    
        
        void download_host()     {
            Memcpy(host,dev,size()*sizeof(T),cudaMemcpyDeviceToHost);
        }
        void upload_device()     {
            Memcpy(dev,host,size()*sizeof(T),cudaMemcpyHostToDevice);
        }    
        void zero()    {            
            dev_ptr->zero();
        }       
        void ones()    {
            fill(T(1.0f,0.0f));
        }    
        
        void fill(const T val)     {
            for(size_t i = 0; i < size(); i++) host[i] = val;            
        }

        T& operator[](array_index pos)    {        
            assert(pos < size());
            while(pos < 0) pos += size();
            return host[pos];        
        }
        T    __getitem__(array_index pos) { return (*this)[pos]; }
        void __setitem__(array_index pos, T val) { while(pos < 0) pos += N; host[pos] = val; }

        
        ComplexVector<Type> operator - ();
        ComplexVector<Type> operator + (const ComplexVector<Type> & a);
        ComplexVector<Type> operator - (const ComplexVector<Type> & a);
        ComplexVector<Type> operator * (const ComplexVector<Type> & a);

        ComplexVector<Type> operator + (const complex<Type> a);
        ComplexVector<Type> operator - (const complex<Type> a);
        ComplexVector<Type> operator * (const complex<Type> a);
        ComplexVector<Type> operator / (const complex<Type> a);
        
        virtual void copy(const ComplexVector<Type> & v) {
            resize(v.size());
            dev_ptr.get()->copy(v.dev_ptr.get());
        }
        void resize(size_t n) {
            N = n;
            dev_ptr.reset();        
            init(n);
        }
                
        void print()  {
            download_host();               
            std::cout << "vector[" << N << "]=";
            for(size_t w = 0; w < size(); w++) {
                std::cout << host[w] << ",";
            }
            std::cout << std::endl;        
        }

        ComplexVector<Type> clone(const ComplexVector<Type> & v) { return ComplexVector<Type>(v); }
        ComplexVector<Type> eval() { return ComplexVector<Type>(*this); }

        Vector<Type> real() {
            Vector<Type> r(size());
            cvector_realf(dev,r.dev,size());
            return r;
        } 
        Vector<Type> imag() {
            Vector<Type> r(size());
            cvector_imagf(dev,r.dev,size());
            return r;
        } 

        Vector<Type> abs() {
            Vector<Type> r(size());
            cvector_absf(dev,r.dev,size());
            return r;
        } 
        Vector<Type> arg() {
            Vector<Type> r(size());
            cvector_argf(dev,r.dev,size());
            return r;
        } 
        ComplexVector<Type> norm() {
            Vector<Type> r(size());
            cvector_normf(dev,r.dev,size());
            return r;
        } 
        ComplexVector<Type> conj() {
            Vector<Type> r(size());
            cvector_conjf(dev,r.dev,size());
            return r;
        } 
        ComplexVector<Type> proj() {
            Vector<Type> r(size());
            cvector_projf(dev,r.dev,size());
            return r;
        } 

        ComplexVector<Type> exp() {
            Vector<Type> r(size());
            cvector_expf(dev,r.dev,size());
            return r;
        } 
        ComplexVector<Type> log() {
            Vector<Type> r(size());
            cvector_logf(dev,r.dev,size());
            return r;
        } 
        ComplexVector<Type> log10() {
            Vector<Type> r(size());
            cvector_log10f(dev,r.dev,size());
            return r;        
        } 
        ComplexVector<Type> cos() {
            Vector<Type> r(size());
            cvector_cosf(dev,r.dev,size());
            return r;
        } 
        ComplexVector<Type> sin() {
            Vector<Type> r(size());
            cvector_sinf(dev,r.dev,size());
            return r;
        } 
        ComplexVector<Type> tan() {
            Vector<Type> r(size());
            cvector_tanf(dev,r.dev,size());
            return r;
        } 
        ComplexVector<Type> acos() {
            Vector<Type> r(size());
            cvector_acosf(dev,r.dev,size());
            return r;
        } 
        ComplexVector<Type> asin() {
            Vector<Type> r(size());
            cvector_asinf(dev,r.dev,size());
            return r;
        } 
        ComplexVector<Type> atan() {
            Vector<Type> r(size());
            cvector_atanf(dev,r.dev,size());
            return r;
        } 
        ComplexVector<Type> cosh() {
            Vector<Type> r(size());
            cvector_coshf(dev,r.dev,size());
            return r;
        } 
        ComplexVector<Type> sinh() {
            Vector<Type> r(size());
            cvector_sinhf(dev,r.dev,size());
            return r;
        } 
        ComplexVector<Type> tanh() {
            Vector<Type> r(size());
            cvector_tanhf(dev,r.dev,size());
            return r;
        } 
        ComplexVector<Type> acosh() {
            Vector<Type> r(size());
            cvector_acoshf(dev,r.dev,size());
            return r;
        } 
        ComplexVector<Type> asinh() {
            Vector<Type> r(size());
            cvector_asinhf(dev,r.dev,size());
            return r;
        } 
        ComplexVector<Type> atanh() {
            Vector<Type> r(size());
            cvector_atanhf(dev,r.dev,size());
            return r;
        } 
        ComplexVector<Type> sqrt() {
            Vector<Type> r(size());
            cvector_sqrtf(dev,r.dev,size());
            return r;
        } 
        ComplexVector<Type> pow(const ComplexVector<T> & c) {
            Vector<Type> r(size());
            cvector_powf(dev,c.dev,r.dev,size());
            return r;
        } 
        ComplexVector<Type> pow(const complex<Type> & c) {
            Vector<Type> r(size());
            cvector_powf_scalar(dev,c,r.dev,size());
            return r;
        } 
        ComplexVector<Type> sigmoid() {
            Vector<Type> r(size());
            cvector_sigmoidf(dev,r.dev,size());
            return r;
        } 
        ComplexVector<Type> sigmoid_grad() {
            Vector<Type> r(size());
            cvector_sigmoid_gradf(dev,r.dev,size());
            return r;
        } 
        ComplexVector<Type> tanh_grad() {
            Vector<Type> r(size());
            cvector_tanh_gradf(dev,r.dev,size());
            return r;
        } 
        ComplexVector<Type> hadmard(const ComplexVector<Type> & c) {
            Vector<Type> r(size());
            cvector_hadamardf(dev,c.dev,r.dev,size());
            return r;
        } 
        
    };

    template<typename T>
    ComplexVector<T> shift(const ComplexVector<T> &a, int count)
    {
        ComplexVector<T> r(a);
        count = count % a.N;
        if(count < 0)
        {
            for(size_t i = count; i < 0; i++) {
                complex<T> x = r[0];
                
                for(size_t j = r.size()-2; j >0; j++)
                    r[j+1] = r[j];
                r[j+1] = x;
            }
        }
        else
        {
            for(size_t i = 0; i < count; i++) {
                complex<T> x = r[r.size()-1];
                
                for(size_t j = r.size()-2; j > 0; j++)
                    r[j] = r[j+1];
                r[0] = x;
            }
        }
    }

    // this is for lua
    template<typename T>
    struct VectorView {
        T * host;
        int     r,c;

        VectorView(T *h, int row, int cols) {        
            host = h;
            r = row;
            c = cols;
        }
        VectorView(const VectorView & v) {
            host = v.host;
            r    = v.r;
            c    = v.c;
        }
        VectorView& operator = (const VectorView & v) {
            host = v.host;
            r    = v.r;
            c    = v.c;
            return *this;
        }

        T& operator[](array_index pos) {
            while(pos < 0) pos += c;
            return host[r*c+pos];
        }
        T __getitem(array_index pos) {
            while(pos < 0) pos += c;
            return host[r*c+pos];
        }
        void  __setitem(array_index pos, T val) {
            while(pos < 0) pos += c;
            host[r*c+pos] = val;                
        }
    };

    
    template<typename T>
    struct Matrix 
    {
    public:

        void init(size_t m, size_t n) {        
            M = m;
            N = n;
            dev_ptr  = { new DevPtr<T>(M*N), [](DevPtr<T>* p){ delete p; } };
            dev  = (T*)dev_ptr.get()->dev;
            host = (T*)dev_ptr.get()->host;    
            //zero();
        }

        T * dev;
        T * host;        
        std::shared_ptr<DevPtr<T>> dev_ptr;
        size_t   M,N;

    public:

        Matrix() {
            dev = host = nullptr;
            M = N = 0;
        }
        Matrix(size_t i, size_t j) {
            init(i,j);           
        }        
        Matrix(size_t i, size_t j, const std::initializer_list<T> & in) {
            init(i,j);
            std::vector<T> v(in.begin(),in.end());            
            Memcpy(host,v.data(), size()*sizeof(T),cudaMemcpyHostToHost);            
        }
        Matrix(size_t i, size_t j, const std::vector<T> & v) {
            init(i,j);            
            Memcpy(host,v.data(), size()*sizeof(T),cudaMemcpyHostToHost);            
        }    
        
        Matrix(const Matrix<T> & m) {        
            *this = m;
        } 
        
        ~Matrix() = default;
            

        Matrix<T>& operator = (const Matrix<T> & m) {              
            M = m.M;
            N = m.N;
            host = m.host;
            dev  = m.dev;        
            dev_ptr.reset();
            dev_ptr = m.dev_ptr;                
            return *this;
        }

        
        size_t size() const { return M*N; }
        size_t rows() const { return M; }
        size_t cols() const { return N; }


        void resize(int i, int j) {
            if(i == M && j == N) return;
            dev_ptr.reset();
            init(i,j);
        }
        
        Matrix<T> t();

        
        Matrix<T> transpose(Matrix<T> & a){
            return a.t();
        }

        void copy(const Matrix<T> & a) {   
            /*     
            if(a.M != M || a.N != N) resize(a.M,a.N);
            cublasStatus_t     status; 
            status = cublasScopy(cublas->handle,M*N,a.dev,1,dev,1);
            ASSERTOK(status);            
            */
            if(M != a.M || N != a.N) resize(a.M,a.N);
            dev_ptr.get()->copy(a.dev_ptr.get());            
        }
        
        Matrix<T> operator + (const Matrix<T> & b);
        Matrix<T> operator - (const Matrix<T> & b);    
        Matrix<T> operator * (const Matrix<T> & m);
        
        Matrix<T> operator -();

        
        Matrix<T> operator + (T v) {
            Matrix<T> r(M,N);
            matrix_r_addf_const(dev,v,r.dev,M,N);
            //vector_r_addf_const(dev,v,r.dev,M*N);
            return r;
        }
        Matrix<T> operator - (T v) {
            Matrix<T> r(M,N);
            matrix_r_subf_const(dev,v,r.dev,M,N);
            //vector_r_subf_const(dev,v,r.dev,M*N);
            return r;
        }    
        Matrix<T> operator * (T v) {
            Matrix<T> r(M,N); 
            // dont know why it needs this
            //r.zero();       
            matrix_r_mulf_const(dev,v,r.dev,M,N);          
            //vector_r_mulf_const(dev,v,r.dev,M*N);        
            return r;
        }
        Matrix<T> operator / (T v) {
            Matrix<T> r(M,N);
            matrix_r_divf_const(dev,v,r.dev,M,N);
            //vector_r_divf_const(dev,v,r.dev,M*N);
            return r;
        }
        Matrix<T> operator % (T v) {
            Matrix<T> r(M,N);
            matrix_r_modf_const(dev,v,r.dev,M,N);
            //vector_r_modf_const(dev,v,r.dev,M*N);
            return r;
        } 
        
        
        void hadamard(const Matrix<T> & b) {
            assert(M == b.M && N == b.N);
            vector_r_mulf(dev,b.dev,dev,M*N);                
        }
        
        
        void download_host()     {
            Memcpy(host,dev,size()*sizeof(T),cudaMemcpyDeviceToHost);            
        }
        void upload_device()     {
            Memcpy(dev,host,size()*sizeof(T),cudaMemcpyHostToDevice);            
        }    
        void zero()    {            
            dev_ptr->zero();
        }
        void ones()    {
            fill(1.0f);
        }    
        void randu() {
            rand<T>(dev,host,size(),0.0f,1.0f);
        }
        void random(T min, T max) {
            rand(dev,host,size(),min,max);
        }    
        void fill(const T val)     {        
            for(size_t i = 0; i < size(); i++) host[i] = val;                    
        }

        T& operator()(int i, int j) {
            while(i < 0) i+= M;
            while(j < 0) j+= N;
            return host[i*N+j];
        }
        T& operator[](array_index pos)    {        
            assert(pos < size());
            while(pos < 0) pos += size();
            return host[pos];        
        }
        
        size_t index(int r, int c) { 
            while(r < 0) r += M;
            while(c < 0) c += N;
            return r*N + c;
        }
        
        T sum();
        
        Matrix<T> eval();
        Matrix<T> eval() const;
        
        // there is a problem I do not know what it is
        // the first is for tanh and sigmoid
        // the second if for relu
        // there is some problem with the Tmath and the Tmatrix 
        void addToEachRow(const Matrix<T> &b, int row=0) {
            assert(N == b.N);
            T alpha=1.0;        
            for(size_t i = 0; i < M; i++)
                //vector_addf_row(dev,(i*N),b.dev,row,N);        
                cublasSaxpy(cublas->handle,N,&alpha,dev+i*N,1,b.dev+row*N,1);
        }
        
        
        T get(int row, int col) {
            return host[row*N + col];
        }
        void set(int row, int col, T val) {
            while(row < 0) row += M;
            while(col < 0) col += N;
            host[row*N + col] = val;

        }

        Matrix<T> get_row(int row) {
            while(row < 0) row += M;     
            assert(row < M);   
            Matrix<T> r(1,cols());                    
            memcpy(r.host,host+row*N,cols()*sizeof(T));
            return r;
        }
        void set_row(const Matrix<T> & m, int dst_row=0, int src_row=0) {
            while(dst_row < 0) dst_row += M;
            while(dst_row < 0) dst_row += M;
            while(src_row < 0) src_row += M;                    
            memcpy(host + dst_row*N, m.host + src_row*m.N, m.N*sizeof(T));
        }
        void set_row(const Vector<T> & v, int dst_row=0) {
            assert(N == v.N);       
            while(dst_row < 0) dst_row += M;             
            memcpy(host + dst_row*N, v.host, v.N*sizeof(T));
        }

        void print() {        
            download_host();        
            for(size_t i = 0; i < M; i++) {
                for(size_t j = 0; j < N; j++) 
                {
                    std::cout << (*this)(i,j);
                    if(j < (N-1)) std::cout << ",";
                }            
                std::cout << std::endl;
            }        
        }
        void print_dims() const {
            std::cout << "Matrix<T>(" << M << "," << N << ")" << std::endl;
        }

        void identity()  {     
            // identity only makes sense on square matrix.   
            assert(M == N);
            size_t c = 0;
            //download_host();
            zero();
            for(size_t i = 0; i < M; i++) {
                host[i*N + c++] = 1;
            }            
            //upload_device();
        }    


        void abs() {
            matrix_r_fabsf(dev,dev,M,N);
            //vector_r_fabsf(dev,dev,M*N);
        }
        void exp() {
            matrix_r_expf(dev,dev,M,N);
            //vector_r_expf(dev,dev,M*N);
        }    
        void log2() {
            matrix_r_log2f(dev,dev,M,N);
            //vector_r_log2f(dev,dev,M*N);
        }
        void log10() {
            matrix_r_log10f(dev,dev,M,N);
            //vector_r_log10f(dev,dev,M*N);
        }
        void pow(T p) {
            matrix_r_powf_const(dev,p,dev,M,N);
            //vector_r_powf_const(dev,p,dev,M*N);
        }
        void pow(const Matrix<T> & p) {
            matrix_r_powf(dev,p.dev,dev,M,N);
            //vector_r_powf(dev,p.dev,dev,M*N);
        }
        void sqrt() {
            matrix_r_sqrtf(dev,dev,M,N);
            //vector_r_sqrtf(dev,dev,M*N);
        }
        void rsqrt() {
            matrix_r_rsqrtf(dev,dev,M,N);
            //vector_r_rsqrtf(dev,dev,M*N);
        }
        void sin() {
            matrix_r_sinf(dev,dev,M,N);
            //vector_r_sinf(dev,dev,M*N);
        }
        void cos() {
            matrix_r_cosf(dev,dev,M,N);
            //vector_r_cosf(dev,dev,M*N);
        }
        void tan() {
            matrix_r_tanf(dev,dev,M,N);
            //vector_r_tanf(dev,dev,M*N);
        }
        void asin() {
            matrix_r_asinf(dev,dev,M,N);
            //vector_r_asinf(dev,dev,M*N);
        }
        void acos() {
            matrix_r_acosf(dev,dev,M,N);
            //vector_r_acosf(dev,dev,M*N);
        }
        void atan() {
            matrix_r_atanf(dev,dev,M,N);
            //vector_r_atanf(dev,dev,M*N);
        }
        void sinh() {
            matrix_r_sinhf(dev,dev,M,N);
            //vector_r_sinhf(dev,dev,M*N);
        }
        void cosh() {
            matrix_r_coshf(dev,dev,M,N);
            //vector_r_coshf(dev,dev,M*N);
        }
        void tanh() {
            matrix_r_tanhf(dev,dev,M,N);
            //vector_r_tanhf(dev,dev,M*N);
        }
        void asinh() {
            matrix_r_asinhf(dev,dev,M,N);
            //vector_r_asinhf(dev,dev,M*N);
        }
        void acosh() {
            matrix_r_acoshf(dev,dev,M,N);
            //vector_r_acoshf(dev,dev,M*N);
        }
        void atanh() {
            matrix_r_atanhf(dev,dev,M,N);
            //vector_r_atanhf(dev,dev,M*N);
        }    
        void atan2(T v) {
            matrix_r_atan2f_const(dev,v,dev,M,N);
            //vector_r_atan2f_const(dev,v,dev,M*N);
        }
        void atan2(const Matrix<T> & v) {
            matrix_r_atan2f(dev,v.dev,dev,M,N);
            //vector_r_atan2f(dev,v.dev,dev,M*N);
        }
        void sigmoid() {
            matrix_r_sigmoidf(dev,dev,M,N);
            ///vector_r_sigmoidf(dev,dev,M*N);
        }
        void sigmoid_deriv() {
            matrix_r_sigmoid_gradf(dev,dev,M,N);
            //vector_r_sigmoid_gradf(dev,dev,M*N);
        }
        void relu() {
            matrix_r_reluf(dev,dev,M,N);
            //vector_r_reluf(dev,dev,M*N);               
        }
        void relu_deriv() {
            matrix_r_relu_gradf(dev,dev,M,N);
            //vector_r_relu_gradf(dev,dev,M*N);                    
        }
        void softmax() {
            matrix_r_softmaxf(dev,dev,M,N);
            //vector_r_softmaxf(dev,dev,M*N);
        }
        void tanh_deriv() {
            matrix_r_tanh_gradf(dev,dev,M,N);
            //vector_r_tanh_gradf(dev,dev,M*N);
        }
        void cbrt() {
            matrix_r_cbrtf(dev,dev,M,N);        
            //vector_r_cbrtf(dev,dev,M*N);        
        }
        void ceil() {
            matrix_r_ceilf(dev,dev,M,N);        
            //vector_r_ceilf(dev,dev,M*N);        
        }
        /*
        void copysign(const T val) {
            vector_r_cbrtf_const(dev,val,dev,M*N);        
        }
        void copysign(const Vector<T>& val) {
            vector_r_cbrtf_const(dev,val.dev,dev,M*N);        
        }
        */
        void cospi() {
            matrix_r_cospif(dev,dev,M,N);        
            //vector_r_cospif(dev,dev,M*N);        
        }
        void cyl_bessel_i0() {
            matrix_r_cyl_bessel_i0f(dev,dev,M,N);
            //vector_r_cyl_bessel_i0f(dev,dev,M*N);
        }
        void cyl_bessel_i1() {
            matrix_r_cyl_bessel_i1f(dev,dev,M,N);
            //vector_r_cyl_bessel_i1f(dev,dev,M*N);
        }
        void erfc() {
            matrix_r_erfcf(dev,dev,M,N);
            //vector_r_erfcf(dev,dev,M*N);
        }
        void erfcx() {
            matrix_r_erfcxf(dev,dev,M,N);
            //vector_r_erfcxf(dev,dev,M*N);
        }
        void erfcinv() {
            matrix_r_erfcinvf(dev,dev,M,N);
            //vector_r_erfcinvf(dev,dev,M*N);
        }
        void erf() {
            matrix_r_erff(dev,dev,M,N);
            //vector_r_erff(dev,dev,M*N);
        }
        void erfinv() {
            matrix_r_erfinvf(dev,dev,M,N);
            //vector_r_erfinvf(dev,dev,M*N);
        }
        void exp10() {
            matrix_r_exp10f(dev,dev,M,N);
            //vector_r_exp10f(dev,dev,M*N);
        }
        void exp2() {
            matrix_r_exp2f(dev,dev,M,N);
            //vector_r_exp2f(dev,dev,M*N);
        }
        void expm1() {
            matrix_r_expm1f(dev,dev,M,N);
            //vector_r_expm1f(dev,dev,M*N);
        }
        void fabs() {
            matrix_r_fabsf(dev,dev,M,N);
            //vector_r_fabsf(dev,dev,M*N);
        }
        void fdim(const Matrix<T> & b) {
            matrix_r_fdimf(dev,b.dev,dev,M,N);
            //vector_r_fdimf(dev,b.dev,dev,M*N);
        }
        void fmod(const Matrix<T> & b) {
            matrix_r_fmodf(dev,b.dev,dev,M,N);
            //vector_r_fmodf(dev,b.dev,dev,M*N);
        }
        void hypot(const Matrix<T> & b) {
            matrix_r_hypotf(dev,b.dev,dev,M,N);
            //vector_r_hypotf(dev,b.dev,dev,M*N);
        }
        void ilogb() {
            matrix_r_ilogbf(dev,dev,M,N);
            //vector_r_ilogbf(dev,dev,M*N);
        }
        void j0() {
            matrix_r_j0f(dev,dev,M,N);
            //vector_r_j0f(dev,dev,M*N);
        }
        void j1() {
            matrix_r_j1f(dev,dev,M*N);
            //vector_r_j1f(dev,dev,M*N);
        }
        void jn(int n) {
            matrix_r_jnf(dev,dev,M,N,n);
            //vector_r_jnf(dev,dev,M*N,n);
        }
        void lgamma() {
            matrix_r_lgammaf(dev,dev,M,N);
            //vector_r_lgammaf(dev,dev,M*N);
        }
        void log1p() {
            matrix_r_log1pf(dev,dev,M,N);
            //vector_r_log1pf(dev,dev,M*N);
        }
        void logb() {
            matrix_r_logbf(dev,dev,M,N);
            //vector_r_logbf(dev,dev,M*N);
        }
        void norm3d(const Matrix<T> & a, const Matrix<T> & b, const Matrix<T> & c) {
            matrix_r_norm3df(a.dev,b.dev,c.dev,dev,M,N);
            //vector_r_norm3df(a.dev,b.dev,c.dev,dev,M*N);
        }
        void norm4d(const Matrix<T> & a, const Matrix<T> & b, const Matrix<T> & c, const Matrix<T> & d) {
            matrix_r_norm4df(a.dev,b.dev,c.dev,d.dev,dev,M,N);
            //vector_r_norm4df(a.dev,b.dev,c.dev,d.dev,dev,M*N);
        }
        void normcdf() {
            matrix_r_normcdff(dev,dev,M,N);
            //vector_r_normcdff(dev,dev,M*N);
        }
        void normcdfinv() {
            matrix_r_normcdfinvf(dev,dev,M,N);
            //vector_r_normcdfinvf(dev,dev,M*N);
        }
        void norm() {
            matrix_r_normf(1,dev,dev,M,N);
            //vector_r_normf(1,dev,dev,M*N);
        }
        void rcbrt() {
            matrix_r_rcbrtf(dev,dev,M,N);        
            //vector_r_rcbrtf(dev,dev,M*N);        
        }
        void remainder(const Matrix<T>& b) {
            matrix_r_remainderf(dev,b.dev,dev,M,N);
            //vector_r_remainderf(dev,b.dev,dev,M*N);
        }
        void rhypot(const Matrix<T> & b) {
            matrix_r_rhypotf(dev,b.dev,dev,M,N);
            //vector_r_rhypotf(dev,b.dev,dev,M*N);
        }
        /*
        void rint() {
            vector_r_rintf(dev,dev,M*N);        
        }
        */
        void rnorm3d(const Matrix<T> & a, const Matrix<T> & b, const Matrix<T> & c) {
            matrix_r_rnorm3df(a.dev,b.dev,c.dev,dev,M,N);
            //vector_r_rnorm3df(a.dev,b.dev,c.dev,dev,M*N);
        }
        void rnorm4d(const Matrix<T> & a, const Matrix<T> & b, const Matrix<T> & c, const Matrix<T> & d) {
            matrix_r_rnorm4df(a.dev,b.dev,c.dev,d.dev,dev,M,N);
            //vector_r_rnorm4df(a.dev,b.dev,c.dev,d.dev,dev,M*N);
        }
        /*
        void round() {
            vector_r_round(dev,dev,M*N);        
        }
        */
        void rnorm() {
            matrix_r_rnormf(1,dev,dev,M,N);
            //vector_r_rnormf(1,dev,dev,M*N);
        }
        void tgamma() {
            matrix_r_tgammaf(dev,dev,M,N);
            //vector_r_tgammaf(dev,dev,M*N);
        }
        void trunc() {
            matrix_r_truncf(dev,dev,M,N);
            //vector_r_truncf(dev,dev,M*N);
        }
        void y0() {
            matrix_r_y0f(dev,dev,M,N);
            //vector_r_y0f(dev,dev,M*N);
        }
        void y1() {
            matrix_r_y1f(dev,dev,M,N);
            //vector_r_y1f(dev,dev,M*N);
        }
        void yn(int n) {
            matrix_r_ynf(n,dev,dev,M,N);
            //vector_r_ynf(n,dev,dev,M*N);
        }

    };

    template<typename Type>
    struct ComplexMatrix 
    {
    public:

        using T = complex<Type>;

        void init(size_t m, size_t n) {        
            M = m;
            N = n;
            dev_ptr  = { new DevPtr<T>(M*N), [](DevPtr<T>* p){ delete p; } };
            dev  = (T*)dev_ptr.get()->dev;
            host = (T*)dev_ptr.get()->host;    
            //zero();
        }

        T * dev;
        T * host;        
        std::shared_ptr<DevPtr<T>> dev_ptr;
        size_t   M,N;

    public:

        ComplexMatrix() {
            dev = host = nullptr;
            M = N = 0;
        }
        ComplexMatrix(size_t i, size_t j) {
            init(i,j);           
        }
        ComplexMatrix(T * p, size_t m, size_t n) {
            M = m;
            N = n;        
            dev_ptr  = {new DevPtr<T>(dev,M*N), [](DevPtr<T>* p){delete p;} };
            dev  = (T*)dev_ptr.get()->dev;
            host = (T*)dev_ptr.get()->host;
            ///zero();
        }
        ComplexMatrix(size_t i, size_t j, const std::initializer_list<T> & in) {
            init(i,j);
            std::vector<T> v(in.begin(),in.end());
            Memcpy(host,v.data(), size()*sizeof(T),cudaMemcpyHostToHost);
        }
        ComplexMatrix(size_t i, size_t j, const std::vector<T> & v) {
            init(i,j);
            Memcpy(host,v.data(), size()*sizeof(T),cudaMemcpyHostToHost);
        }    
        ComplexMatrix(const ComplexMatrix<T> & m) {        
            //init(m.M,m.N);
            //copy(m);        
            *this = m;
        }    
        ~ComplexMatrix() = default;
            

        ComplexMatrix<T>& operator = (const ComplexMatrix<T> & m) {              
            M = m.M;
            N = m.N;
            host = m.host;
            dev  = m.dev;        
            dev_ptr.reset();
            dev_ptr = m.dev_ptr;                
            return *this;
        }

        
        size_t size() const { return M*N; }
        size_t rows() const { return M; }
        size_t cols() const { return N; }


        void resize(int i, int j) {
            if(i == M && j == N) return;
            dev_ptr.reset();
            init(i,j);
        }
        
        ComplexMatrix<T> t();

        
        ComplexMatrix<T> transpose(ComplexMatrix<T> & a){
            return a.t();
        }

        void copy(const ComplexMatrix<T> & a) {   
            /*     
            if(a.M != M || a.N != N) resize(a.M,a.N);
            cublasStatus_t     status; 
            status = cublasScopy(cublas->handle,M*N,a.dev,1,dev,1);
            ASSERTOK(status);            
            */
            resize(a.M,a.N);
            dev_ptr.get()->copy(a.dev_ptr.get());
        }
        
        ComplexMatrix<T> operator + (const ComplexMatrix<T> & b);
        ComplexMatrix<T> operator - (const ComplexMatrix<T> & b);    
        ComplexMatrix<T> operator * (const ComplexMatrix<T> & m);
        
        ComplexMatrix<T> operator -();

        
        ComplexMatrix<T> operator + (T v) {
            ComplexMatrix<T> r(M,N);
            matrix_r_addf_const(dev,v,r.dev,M,N);
            //vector_r_addf_const(dev,v,r.dev,M*N);
            return r;
        }
        ComplexMatrix<T> operator - (T v) {
            ComplexMatrix<T> r(M,N);
            cmatrix_r_subf_const(dev,v,r.dev,M,N);
            //vector_r_subf_const(dev,v,r.dev,M*N);
            return r;
        }    
        ComplexMatrix<T> operator * (T v) {
            ComplexMatrix<T> r(M,N); 
            // dont know why it needs this
            //r.zero();       
            cmatrix_r_mulf_const(dev,v,r.dev,M,N);          
            //vector_r_mulf_const(dev,v,r.dev,M*N);        
            return r;
        }
        ComplexMatrix<T> operator / (T v) {
            Matrix<T> r(M,N);
            cmatrix_r_divf_const(dev,v,r.dev,M,N);
            //vector_r_divf_const(dev,v,r.dev,M*N);
            return r;
        }
        
        
        void hadamard(const ComplexMatrix<T> & b) {
            assert(M == b.M && N == b.N);
            vector_r_mulf(dev,b.dev,dev,M*N);                
        }
        
        void download_host()     {
            Memcpy(host,dev,size()*sizeof(T),cudaMemcpyDeviceToHost);
        }
        void upload_device()     {
            Memcpy(dev,host,size()*sizeof(T),cudaMemcpyHostToDevice);
        }    
        void zero()    {
            //cudaMemsetAsync(dev,0x00,size()*sizeof(T));
            dev_ptr->zero();
        }
        void ones()    {
            fill(T(1.0,0.0));
        }            
        void fill(const T val)     {        
            for(size_t i = 0; i < size(); i++) host[i] = val;                    
        }

        T& operator()(int i, int j) {
            while(i < 0) i+= M;
            while(j < 0) j+= N;
            return host[i*N+j];
        }
        T& operator[](array_index pos)    {        
            assert(pos < size());
            while(pos < 0) pos += size();
            return host[pos];        
        }
        
        size_t index(int r, int c) { 
            while(r < 0) r += M;
            while(c < 0) c += N;
            return r*N + c;
        }
      
        T sum();
        
        ComplexMatrix<T> eval();

        ComplexMatrix<Type> real() {
            ComplexMatrix<Type> r(M,N);
            cmatrix_realf(dev,r.dev,M,N);
            return r;
        } 
        ComplexMatrix<Type> imag() {
            ComplexMatrix<Type> r(M,N);
            cmatrix_r_imagf(dev,r.dev,M,N);
            return r;
        } 

        ComplexMatrix<Type> abs() {
            ComplexMatrix<Type>  r(M,N);
            cmatrix_r_absf(dev,r.dev,M,N);
            return r;
        } 
        ComplexMatrix<Type>  arg() {
            ComplexMatrix<Type>  r(M,N);
            matrix_r_argf(dev,r.dev,M,N);
            return r;
        } 
        ComplexMatrix<Type> norm() {
            ComplexMatrix<Type>  r(M,N);
            cmatrix_r_normf(dev,r.dev,M,N);
            return r;
        } 
        ComplexMatrix<Type> conj() {
            ComplexMatrix<Type>  r(M,N);
            cmatrix_r_conjf(dev,r.dev,M,N);
            return r;
        } 
        ComplexMatrix<Type> proj() {
            ComplexMatrix<Type>  r(M,N);
            cmatrix_r_projf(dev,r.dev,M,N);
            return r;
        } 
        
        // there is a problem I do not know what it is
        // the first is for tanh and sigmoid
        // the second if for relu
        // there is some problem with the Tmath and the Tmatrix 
        void addToEachRow(const ComplexMatrix<Type> &b, int row=0);
        
        
        T get(int row, int col) {
            return host[row*N + col];
        }
        void set(int row, int col, T val) {
            while(row < 0) row += M;
            while(col < 0) col += N;
            host[row*N + col] = val;

        }

        ComplexMatrix<T> get_row(int row) {
            while(row < 0) row += M;     
            assert(row < M);   
            ComplexMatrix<T> r(1,cols());        
            Memcpy(r.dev,dev+row*N,cols()*sizeof(T),cudaMemcpyDeviceToDevice);
            return r;
        }
        void set_row(const ComplexMatrix<T> & m, int dst_row=0, int src_row=0) {
            while(dst_row < 0) dst_row += M;
            while(dst_row < 0) dst_row += M;
            while(src_row < 0) src_row += M;        
            Memcpy(dev + dst_row*N, m.dev + src_row*m.N, m.N*sizeof(T), cudaMemcpyDeviceToDevice);
        }
        void set_row(const ComplexVector<T> & v, int dst_row=0) {
            assert(N == v.N);       
            while(dst_row < 0) dst_row += M; 
            Memcpy(dev + dst_row*N, v.dev, v.N*sizeof(T), cudaMemcpyDeviceToDevice);
        }

        void print() {        
            download_host();        
            for(size_t i = 0; i < M; i++) {
                for(size_t j = 0; j < N; j++) 
                {
                    std::cout << (*this)(i,j);
                    if(j < (N-1)) std::cout << ",";
                }            
                std::cout << std::endl;
            }        
        }
        void print_dims() const {
            std::cout << "Matrix<T>(" << M << "," << N << ")" << std::endl;
        }

        void identity()  {     
            // identity only makes sense on square matrix.   
            assert(M == N);
            size_t c = 0;
            //download_host();
            zero();
            for(size_t i = 0; i < M; i++) {
                host[i*N + c++] = 1;
            }            
            //upload_device();
        }    
    };

    
    template<> void ComplexMatrix<float>::addToEachRow(const ComplexMatrix<float> &b, int row) {
            assert(N == b.N);
            cuComplex alpha={1.0,0.0};        
            for(size_t i = 0; i < M; i++)                
                cublasCaxpy(cublas->handle,N,&alpha,(cuComplex*)(dev+i*N),1,(cuComplex*)(b.dev+row*N),1);
        }
    template<> void ComplexMatrix<double>::addToEachRow(const ComplexMatrix<double> &b, int row) {
            assert(N == b.N);
            cuDoubleComplex alpha={1.0,0.0};        
            for(size_t i = 0; i < M; i++)                
                cublasZaxpy(cublas->handle,N,&alpha,(cuDoubleComplex*)(dev+i*N),1,(cuDoubleComplex*)(b.dev+row*N),1);
        }
    
    template<> Matrix<float> Matrix<float>::t() {                    
            Matrix<float> r(N,M);                
            int m = M;
            int n = N;    
            float alpha = 1.0;
            float beta = 1.0;
            cublasStatus_t status = cublasSgeam(cublas->handle,CUBLAS_OP_T,CUBLAS_OP_N,m,n,&alpha,dev,n,&beta,dev,m,r.dev,m);    
            ASSERTOK(status);
            return r;            
        }
    template<> Matrix<double> Matrix<double>::t() {                    
            Matrix<double> r(N,M);                
            int m = M;
            int n = N;    
            double alpha = 1.0;
            double beta = 1.0;
            cublasStatus_t status = cublasDgeam(cublas->handle,CUBLAS_OP_T,CUBLAS_OP_N,m,n,&alpha,dev,n,&beta,dev,m,r.dev,m);    
            ASSERTOK(status);
            return r;            
        }
    
    template<> float Matrix<float>::sum() {
            float r = 0;
            cublasStatus_t     status; 
            status = cublasSasum(cublas->handle,M*N,dev,1,&r);
            ASSERTOK(status);
            return r;
        }        
    template<> double Matrix<double>::sum() {
            double r = 0;
            cublasStatus_t     status; 
            status = cublasDasum(cublas->handle,M*N,dev,1,&r);
            ASSERTOK(status);
            return r;
        }                
    /*
    void dgmm(Matrix<float> & m, const Vector<float> & diagonal) {
            cublasStatus_t status = cublasSdgmm(cublas->handle,CUBLAS_SIDE_RIGHT,m.N,m.M,m.dev,diagonal.N,diagonal.dev,1,m.dev,N);
            ASSERTOK(status);
        }
    void dgmm(Matrix<double> & m, const Vector<double> & diagonal) {
            cublasStatus_t status = cublasSdgmm(cublas->handle,CUBLAS_SIDE_RIGHT,m.N,m.M,m.dev,diagonal.N,diagonal.dev,1,m.dev,N);
            ASSERTOK(status);
        }        
    */
    void copy(const Matrix<float> & a, Matrix<float> & r) {       
        if(a.M != r.M || a.N != r.N) r.resize(a.M,a.N);
        //cublasStatus_t     status; 
        //status = cublasScopy(cublas->handle,a.M*a.N,a.dev,1,r.dev,1);
        //ASSERTOK(status);            
        //cudaMemcpy2DAsync(r.dev,r.N*sizeof(float),a.dev,a.N*sizeof(float),a.N*sizeof(float),r.M,cudaMemcpyDeviceToDevice,get_cuda_stream());
        cudaMemcpy2D(r.dev,r.N*sizeof(float),a.dev,a.N*sizeof(float),a.N*sizeof(float),r.M,cudaMemcpyDeviceToDevice);
    }
    void copy(const Matrix<double> & a, Matrix<double> & r) {       
        if(a.M != r.M || a.N != r.N) r.resize(a.M,a.N);
        //cublasStatus_t     status; 
        //status = cublasDcopy(cublas->handle,a.M*a.N,a.dev,1,r.dev,1);
        //ASSERTOK(status);            
        //cudaMemcpy2DAsync(r.dev,r.N*sizeof(double),a.dev,a.N*sizeof(double),a.N*sizeof(double),r.M,cudaMemcpyDeviceToDevice,get_cuda_stream());
        cudaMemcpy2D(r.dev,r.N*sizeof(double),a.dev,a.N*sizeof(double),a.N*sizeof(double),r.M,cudaMemcpyDeviceToDevice);
    }
    void copy(const ComplexMatrix<float> & a, ComplexMatrix<float> & r) {       
        if(a.M != r.M || a.N != r.N) r.resize(a.M,a.N);
        //cublasStatus_t     status; 
        //status = cublasScopy(cublas->handle,a.M*a.N,a.dev,1,r.dev,1);
        //ASSERTOK(status);            
        //cudaMemcpy2DAsync(r.dev,r.N*sizeof(float),a.dev,a.N*sizeof(float),a.N*sizeof(float),r.M,cudaMemcpyDeviceToDevice,get_cuda_stream());
        cudaMemcpy2D(r.dev,r.N*sizeof(complex<float>),a.dev,a.N*sizeof(complex<float>),a.N*sizeof(complex<float>),r.M,cudaMemcpyDeviceToDevice);
    }
    void copy(const ComplexMatrix<double> & a, ComplexMatrix<double> & r) {       
        if(a.M != r.M || a.N != r.N) r.resize(a.M,a.N);
        //cublasStatus_t     status; 
        //status = cublasDcopy(cublas->handle,a.M*a.N,a.dev,1,r.dev,1);
        //ASSERTOK(status);            
        //cudaMemcpy2DAsync(r.dev,r.N*sizeof(double),a.dev,a.N*sizeof(double),a.N*sizeof(double),r.M,cudaMemcpyDeviceToDevice,get_cuda_stream());
        cudaMemcpy2D(r.dev,r.N*sizeof(complex<double>),a.dev,a.N*sizeof(complex<double>),a.N*sizeof(complex<double>),r.M,cudaMemcpyDeviceToDevice);
    }
    
    
    template<typename T>
    inline Matrix<T> CopyMatrix(const Matrix<T> & m) {        
        Matrix<T> r(m.M,m.N);
        r.copy(m);        
        //cudaMemcpy2DAsync(r.dev,r.N*sizeof(T),m.dev,m.N*sizeof(T),m.N*sizeof(T),r.M,cudaMemcpyDeviceToDevice,get_cuda_stream());
        return r;
    }

    template<typename T>
    inline Vector<T> CopyVector(const Vector<T> & m) {
        Vector<T> r(m.N);
        r.copy(m);
        //Memcpy(r.dev,m.dev,r.size()*sizeof(T),cudaMemcpyDeviceToDevice);
        return r;
    }
    template<typename T>
    inline ComplexMatrix<T> CopyComplexMatrix(const ComplexMatrix<T> & m) {        
        ComplexMatrix<T> r(m.M,m.N);
        r.copy(m);        
        //cudaMemcpy2DAsync(r.dev,r.N*sizeof(T),m.dev,m.N*sizeof(T),m.N*sizeof(T),r.M,cudaMemcpyDeviceToDevice,get_cuda_stream());
        return r;
    }

    template<typename T>
    inline ComplexVector<T> CopyComplexVector(const ComplexVector<T> & m) {
        ComplexVector<T> r(m.N);
        r.copy(m);
        //Memcpy(r.dev,m.dev,r.size()*sizeof(T),cudaMemcpyDeviceToDevice);
        return r;
    }

    template<typename T>
    inline Vector<T> Vector<T>::operator - () {
        Vector<T> r(N);
        r.copy(*this);
        vector_r_mulf_const(r.dev,T(-1.0),r.dev,N);
        return r;
    }

    
    template<typename T>
    inline Vector<T> Vector<T>::operator + (const Vector<T> & a) {
        assert(size() == a.size());
        Vector<T> r = CopyVector(a);
        cublasStatus_t     status; 
        T alpha = 1.0;
        status = cublasSaxpy(cublas->handle, N,&alpha,dev,1,r.dev,1);
        ASSERTOK(status);
        return r;
    }
    
    template<>
    inline Vector<float> Vector<float>::operator + (const Vector<float> & a) {
        assert(size() == a.size());
        Vector<float> r = CopyVector(a);
        cublasStatus_t     status; 
        float alpha = 1.0;
        status = cublasSaxpy(cublas->handle, N,&alpha,dev,1,r.dev,1);
        ASSERTOK(status);
        return r;
    }
    template<>
    inline Vector<double> Vector<double>::operator + (const Vector<double> & a) {
        assert(size() == a.size());
        Vector<double> r = CopyVector(a);
        cublasStatus_t     status; 
        double alpha = 1.0;
        status = cublasDaxpy(cublas->handle, N,&alpha,dev,1,r.dev,1);
        ASSERTOK(status);
        return r;
    }


    template<typename T>
    inline Vector<T> Vector<T>::operator - (const Vector<T> & a) {
        assert(size() == a.size());
        Vector<T> r = CopyVector(*this);
        cublasStatus_t     status; 
        T alpha = T(-1.0);
        status = cublasSaxpy(cublas->handle,N,&alpha,a.dev,1,r.dev,1);
        ASSERTOK(status);
        return r;
    }


    template<>
    inline Vector<float> Vector<float>::operator - (const Vector<float> & a) {
        assert(size() == a.size());
        Vector<float> r = CopyVector(*this);
        cublasStatus_t     status; 
        float alpha = -1.0f;
        status = cublasSaxpy(cublas->handle,N,&alpha,a.dev,1,r.dev,1);
        ASSERTOK(status);
        return r;
    }

    template<>
    inline Vector<double> Vector<double>::operator - (const Vector<double> & a) {
        assert(size() == a.size());
        Vector<double> r = CopyVector(*this);
        cublasStatus_t     status; 
        double alpha = -1.0;
        status = cublasDaxpy(cublas->handle,N,&alpha,a.dev,1,r.dev,1);
        ASSERTOK(status);
        return r;
    }

    
    template<typename T>
    inline Vector<T> Vector<T>::operator * (const Vector<T> & a) {
        assert(size() == a.size());
        Vector<T> r = CopyVector(a);
        cublasStatus_t     status;         
        status = cublasSdgmm(cublas->handle,CUBLAS_SIDE_RIGHT,1,N,dev,1,a.dev,1,r.dev,1);
        ASSERTOK(status);
        return r;
    }

    template<>
    inline Vector<float> Vector<float>::operator * (const Vector<float> & a) {
        assert(size() == a.size());
        Vector<float> r = CopyVector(a);
        cublasStatus_t     status;         
        status = cublasSdgmm(cublas->handle,CUBLAS_SIDE_RIGHT,1,N,dev,1,a.dev,1,r.dev,1);
        ASSERTOK(status);
        return r;
    }

    template<>
    inline Vector<double> Vector<double>::operator * (const Vector<double> & a) {
        assert(size() == a.size());
        Vector<double> r = CopyVector(a);
        cublasStatus_t     status;         
        status = cublasDdgmm(cublas->handle,CUBLAS_SIDE_RIGHT,1,N,dev,1,a.dev,1,r.dev,1);
        ASSERTOK(status);
        return r;
    }

    template<typename T>
    inline Vector<T> Vector<T>::operator + (const T a) {    
        Vector<T> r(size());
        vector_r_addf_const(dev,a,r.dev,size());
        return r;
    }

    template<typename T>
    inline Vector<T> Vector<T>::operator - (const T a) {    
        Vector<T> r(size());
        vector_r_subf_const(dev,a,r.dev,size());
        return r;
    }

    template<typename T>
    inline Vector<T> Vector<T>::operator * (const T a) {    
        Vector<T> r(size());    
        vector_r_mulf_const(dev,a,r.dev,r.size());
        return r;
    }

    template<typename T>
    inline Vector<T> Vector<T>::operator / (const T a) {    
        Vector<T> r(size());
        vector_r_divf_const(dev,a,r.dev,r.size());
        return r;
    }

    
    template<typename T>
    inline void add(const Vector<T> & a, const Vector<T> & b, Vector<T> & r)
    {    
        assert(b.size() == a.size());        
        cublasStatus_t     status; 
        T alpha = 1.0;
        r.resize(a.N);
        r.copy(b);
        status = cublasSaxpy(cublas->handle, a.N,&alpha,a.dev,1,r.dev,1);
        ASSERTOK(status);    
    }
    template<>
    inline void add(const Vector<float> & a, const Vector<float> & b, Vector<float> & r)
    {    
        assert(b.size() == a.size());        
        cublasStatus_t     status; 
        float alpha = 1.0;
        r.resize(a.N);
        r.copy(b);
        status = cublasSaxpy(cublas->handle, a.N,&alpha,a.dev,1,r.dev,1);
        ASSERTOK(status);    
    }


    template<>
    inline void add(const Vector<double> & a, const Vector<double> & b, Vector<double> & r)
    {    
        assert(b.size() == a.size());        
        cublasStatus_t     status; 
        double alpha = 1.0;
        r.resize(a.N);
        r.copy(b);
        status = cublasDaxpy(cublas->handle, a.N,&alpha,a.dev,1,r.dev,1);
        ASSERTOK(status);    
    }

    
    template<typename T>
    inline void sub(const Vector<T> & a, const Vector<T> & b, Vector<T> & r)
    {
        assert(b.size() == a.size());        
        cublasStatus_t     status; 
        T alpha = -1.0;
        r.resize(a.N);
        r.copy(b);
        status = cublasSaxpy(cublas->handle,a.N,&alpha,a.dev,1,r.dev,1);
        ASSERTOK(status);
    }
        
    template<>
    inline void sub(const Vector<float> & a, const Vector<float> & b, Vector<float> & r)
    {
        assert(b.size() == a.size());        
        cublasStatus_t     status; 
        float alpha = -1.0;
        r.resize(a.N);
        r.copy(b);
        status = cublasSaxpy(cublas->handle,a.N,&alpha,a.dev,1,r.dev,1);
        ASSERTOK(status);
    }

    
    template<>
    inline void sub(const Vector<double> & a, const Vector<double> & b, Vector<double> & r)
    {
        assert(b.size() == a.size());        
        cublasStatus_t     status; 
        double alpha = -1.0;
        r.resize(a.N);
        r.copy(b);
        status = cublasDaxpy(cublas->handle,a.N,&alpha,a.dev,1,r.dev,1);
        ASSERTOK(status);
    }

    
    template<typename T>
    inline void mul(const Vector<T> & a, const Vector<T> & b, Vector<T> & r)
    {
        assert(b.size() == a.size());
        cublasStatus_t     status;         
        r.resize(a.N);    
        status = cublasSdgmm(cublas->handle,CUBLAS_SIDE_RIGHT,1,a.N,b.dev,1,a.dev,1,r.dev,1);
        ASSERTOK(status);
    }
    
    template<>
    inline void mul(const Vector<float> & a, const Vector<float> & b, Vector<float> & r)
    {
        assert(b.size() == a.size());
        cublasStatus_t     status;         
        r.resize(a.N);    
        status = cublasSdgmm(cublas->handle,CUBLAS_SIDE_RIGHT,1,a.N,b.dev,1,a.dev,1,r.dev,1);
        ASSERTOK(status);
    }
    template<>
    inline void mul(const Vector<double> & a, const Vector<double> & b, Vector<double> & r)
    {
        assert(b.size() == a.size());
        cublasStatus_t     status;         
        r.resize(a.N);    
        status = cublasDdgmm(cublas->handle,CUBLAS_SIDE_RIGHT,1,a.N,b.dev,1,a.dev,1,r.dev,1);
        ASSERTOK(status);
    }

    
    template<typename T>
    inline void div(const Vector<T> & a, const Vector<T> & b, Vector<T> & r)
    {
        assert(b.size() == a.size());
        cublasStatus_t     status;         
        r.resize(a.N);    
        status = cublasSdgmm(cublas->handle,CUBLAS_SIDE_RIGHT,1,a.N,b.dev,1,a.dev,1,r.dev,1);
        ASSERTOK(status);
    }
    
    template<>
    inline void div(const Vector<float> & a, const Vector<float> & b, Vector<float> & r)
    {
        assert(b.size() == a.size());
        cublasStatus_t     status;         
        r.resize(a.N);    
        status = cublasSdgmm(cublas->handle,CUBLAS_SIDE_RIGHT,1,a.N,b.dev,1,a.dev,1,r.dev,1);
        ASSERTOK(status);
    }

    template<typename T>
    inline Vector<T> operator + (const Vector<T> & a, const Vector<T> & b) {
        Vector<T> r(a.N);    
        add(a,b,r);
        return r;
    }

    template<typename T>
    inline Vector<T> operator - (const Vector<T> & a, const Vector<T> & b) {
        Vector<T> r(a.N);    
        sub(a,b,r);
        return r;
    }

    template<typename T>
    inline Vector<T> operator * (const Vector<T> & a, const Vector<T> & b) {
        Vector<T> r(a.N);    
        mul(a,b,r);
        return r;
    }


    template<typename T>
    inline Vector<T> operator + (const Vector<T> & a, const T v)
    {
        Vector<T> r(a.N);    
        vector_r_addf_const(a.dev,v,r.dev,a.N);
        return r;
    }

    template<typename T>
    inline Vector<T> operator + (const T v, const Vector<T> & a)
    {    
        Vector<T> r(a.N);
        vector_r_addf_const(a.dev,v,r.dev,a.N);
        return r;
    }

    template<typename T>
    inline Vector<T> operator - (const Vector<T> & a, const T v)
    {    
        Vector<T> r(a.N);    
        vector_r_subf_const(a.dev,v,r.dev,a.N);
        return r;
    }

    template<typename T>
    inline Vector<T> operator - (const T v, const Vector<T> & a)
    {    
        Vector<T> r(a.N);        
        vector_r_subf_const(a.dev,v,r.dev,a.N);
        return r;
    }

    template<typename T>
    inline Vector<T> operator * (const Vector<T> & a, const T v)
    {    
        Vector<T> r(a.N);    
        vector_r_mulf_const(a.dev,v,r.dev,a.N);
        //cublasSscal(cublas->handle,a.N,&v,r.dev,1);
        return r;
    }

    template<typename T>
    inline Vector<T> operator * (const T v, const Vector<T> & a)
    {    
        Vector<T> r(a.N);    
        vector_r_mulf_const(a.dev,v,r.dev,a.N);
        //cublasSscal(cublas->handle,a.N,&v,r.dev,1);
        return r;
    }

    template<typename T>
    inline Vector<T> operator / (const Vector<T> & a, const T v)
    {    
        Vector<T> r(a.N);    
        vector_r_divf_const(a.dev,v,r.dev,a.N);    
        return r;
    }

    template<typename T>
    inline Vector<T> operator % (const Vector<T> & a, const T v)
    {    
        Vector<T> r(a.N);        
        vector_r_modf_const(r.dev,v,r.dev,a.N);
        return r;
    }

    


    template<typename T>
    inline Vector<T> abs(const Vector<T> & m) {    
        Vector<T> r = CopyVector(m);
        r.abs();
        return r;    
    }

    template<typename T>
    inline Vector<T> exp(const Vector<T> & m) {
        Vector<T> r = CopyVector(m);
        r.exp(); 
        return r;   
    }    

    template<typename T>
    inline Vector<T> log2(const Vector<T> & m) {
        Vector<T> r = CopyVector(m);
        r.log2();
        return r;
    }

    template<typename T>
    inline Vector<T> log10(const Vector<T> & m) {
        Vector<T> r = CopyVector(m);
        r.log10();
        return r;
    }

    template<typename T>
    inline Vector<T> pow(const Vector<T> & m,T p) {
        Vector<T> r = CopyVector(m);
        r.pow(p);
        return r;
    }

    template<typename T>
    inline Vector<T> pow(const Vector<T> & m,const Vector<T> & p) {
        Vector<T> r = CopyVector(m);
        r.pow(p);
        return r;
    }

    template<typename T>
    inline Vector<T> sqrt(const Vector<T> & m) {
        Vector<T> r = CopyVector(m);
        r.sqrt();
        return r;
    }

    template<typename T>
    inline Vector<T> rsqrt(const Vector<T> & m) {
        Vector<T> r = CopyVector(m);
        r.rsqrt();
        return r;
    }

    template<typename T>
    inline Vector<T> sin(const Vector<T> & m) {
        Vector<T> r = CopyVector(m);
        r.sin();
        return r;
    }

    template<typename T>
    inline Vector<T> cos(const Vector<T> & m) {
        Vector<T> r = CopyVector(m);
        r.cos();
        return r;
    }

    template<typename T>
    inline Vector<T> tan(const Vector<T> & m) {
        Vector<T> r = CopyVector(m);
        r.tan();
        return r;
    }

    template<typename T>
    inline Vector<T> asin(const Vector<T> & m) {
        Vector<T> r = CopyVector(m);
        r.asin();
        return r;
    }

    template<typename T>
    inline Vector<T> acos(const Vector<T> & m) {
        Vector<T> r = CopyVector(m);
        r.acos();
        return r;
    }

    template<typename T>
    inline Vector<T> atan(const Vector<T> & m) {
        Vector<T> r = CopyVector(m);
        r.atan();
        return r;
    }

    template<typename T>
    inline Vector<T> sinh(const Vector<T> & m) {
        Vector<T> r = CopyVector(m);
        r.sinh();
        return r;
    }

    template<typename T>
    inline Vector<T> cosh(const Vector<T> & m) {
        Vector<T> r = CopyVector(m);
        r.cosh();
        return r;
    }

    template<typename T>
    inline Vector<T> tanh(const Vector<T> & m) {
        Vector<T> r = CopyVector(m);
        r.tanh();
        return r;
    }

    template<typename T>
    inline Vector<T> asinh(const Vector<T> & m) {
        Vector<T> r = CopyVector(m);
        r.asinh();
        return r;
    }

    template<typename T>
    inline Vector<T> acosh(const Vector<T> & m) {
        Vector<T> r = CopyVector(m);
        r.acosh();
        return r;
    }

    template<typename T>
    inline Vector<T> atanh(const Vector<T> & m) {
        Vector<T> r = CopyVector(m);
        r.atanh();
        return r;
    }    

    template<typename T>
    inline Vector<T> atan2(const Vector<T> & m,T v) {
        Vector<T> r = CopyVector(m);
        r.atan2(v);
        return r;
    }

    template<typename T>
    inline Vector<T> atan2(const Vector<T> & m,const Vector<T> & v) {
        Vector<T> r = CopyVector(m);
        r.atan2(v);
        return r;
    }

    template<typename T>
    inline Vector<T> sigmoid(const Vector<T> & m) {
        Vector<T> r = CopyVector(m);
        r.sigmoid();
        return r;
    }

    template<typename T>
    inline Vector<T> sigmoid_deriv(const Vector<T> & m) {
        Vector<T> r = CopyVector(m);
        r.sigmoid_deriv();
        return r;
    }

    template<typename T>
    inline Vector<T> relu(const Vector<T> & m) {
        Vector<T> r = CopyVector(m);
        r.relu();
        return r;
    }

    template<typename T>
    inline Vector<T> relu_deriv(const Vector<T> & m) {
        Vector<T> r = CopyVector(m);
        r.relu_deriv();
        return r;
    }

    template<typename T>
    inline Vector<T> softmax(const Vector<T> & m) {
        Vector<T> r = CopyVector(m);
        r.softmax();
        return r;
    }

    template<typename T>
    inline Vector<T> tanh_deriv(const Vector<T> & m) {
        Vector<T> r = CopyVector(m);
        r.tanh_deriv();
        return r;
    }


    template<typename T>
    inline Vector<T> cbrt(const Vector<T> & m) {
        Vector<T> r = CopyVector(m);
        r.cbrt();
        return r;
    }

    template<typename T>
    inline Vector<T> ceil(const Vector<T> & m) {
        Vector<T> r = CopyVector(m);
        r.ceil();
        return r;
    }
    /*
    inline Vector<T> copysign(const Vector<T> & m, T x) {
        Vector<T> r = CopyVector(m);
        r.copysign(x);
        return r;
    }
    inline Vector<T> copysign(const Vector<T> & m, const Vector<T> & x) {
        Vector<T> r = CopyVector(m);
        r.copysign(x);
        return r;
    }
    */

    template<typename T>
    inline Vector<T> cospi(const Vector<T> & m) {
        Vector<T> r = CopyVector(m);
        r.cospi();
        return r;
    }

    template<typename T>
    inline Vector<T> cyl_bessel_i0f(const Vector<T> & m) {
        Vector<T> r = CopyVector(m);
        r.cyl_bessel_i0();
        return r;
    }

    template<typename T>
    inline Vector<T> cyl_bessel_i1f(const Vector<T> & m) {
        Vector<T> r = CopyVector(m);
        r.cyl_bessel_i1();
        return r;
    }

    template<typename T>
    inline Vector<T> erfc(const Vector<T> & m) {
        Vector<T> r = CopyVector(m);
        r.erfc();
        return r;
    }

    template<typename T>
    inline Vector<T> erfcx(const Vector<T> & m) {
        Vector<T> r = CopyVector(m);
        r.erfcx();
        return r;
    }

    template<typename T>
    inline Vector<T> erfcinv(const Vector<T> & m) {
        Vector<T> r = CopyVector(m);
        r.erfcinv();
        return r;
    }

    template<typename T>
    inline Vector<T> erf(const Vector<T> & m) {
        Vector<T> r = CopyVector(m);
        r.erf();
        return r;
    }

    template<typename T>
    inline Vector<T> erfinv(const Vector<T> & m) {
        Vector<T> r = CopyVector(m);
        r.erfinv();
        return r;
    }

    template<typename T>
    inline Vector<T> exp10(const Vector<T> & m) {
        Vector<T> r = CopyVector(m);
        r.exp10();
        return r;
    }

    template<typename T>
    inline Vector<T> exp2(const Vector<T> & m) {
        Vector<T> r = CopyVector(m);
        r.exp2();
        return r;
    }

    template<typename T>
    inline Vector<T> expm1(const Vector<T> & m) {
        Vector<T> r = CopyVector(m);
        r.expm1();
        return r;
    }

    template<typename T>
    inline Vector<T> fabs(const Vector<T> & m) {
        Vector<T> r = CopyVector(m);
        r.fabs();
        return r;
    }

    template<typename T>
    inline Vector<T> fdim(const Vector<T> & m, const Vector<T> & b) {
        Vector<T> r = CopyVector(m);
        r.fdim(b);
        return r;
    }

    template<typename T>
    inline Vector<T> fmod(const Vector<T> & m, const Vector<T> & b) {
        Vector<T> r = CopyVector(m);
        r.fmod(b);
        return r;
    }

    template<typename T>
    inline Vector<T> hypot(const Vector<T> & m, const Vector<T> & b) {
        Vector<T> r = CopyVector(m);
        r.hypot(b);
        return r;
    }

    template<typename T>
    inline Vector<T> ilogb(const Vector<T> & m) {
        Vector<T> r = CopyVector(m);
        r.ilogb();
        return r;
    }

    template<typename T>
    inline Vector<T> j0(const Vector<T> & m) {
        Vector<T> r = CopyVector(m);
        r.j0();
        return r;
    }

    template<typename T>
    inline Vector<T> j1(const Vector<T> & m) {
        Vector<T> r = CopyVector(m);
        r.j1();
        return r;
    }

    template<typename T>
    inline Vector<T> jn(const Vector<T> & m, int n) {
        Vector<T> r = CopyVector(m);
        r.jn(n);
        return r;
    }

    template<typename T>
    inline Vector<T> lgamma(const Vector<T> & m) {
        Vector<T> r = CopyVector(m);
        r.lgamma();
        return r;
    }

    template<typename T>
    inline Vector<T> log1p(const Vector<T> & m) {
        Vector<T> r = CopyVector(m);
        r.log1p();
        return r;
    }

    template<typename T>
    inline Vector<T> logb(const Vector<T> & m) {
        Vector<T> r = CopyVector(m);
        r.logb();
        return r;
    }

    template<typename T>
    inline Vector<T> norm3d(const Vector<T> & m,const Vector<T> & a, const Vector<T> & b, const Vector<T> & c) {
        Vector<T> r = CopyVector(m);
        r.norm3d(a,b,c);
        return r;
    }

    template<typename T>
    inline Vector<T> norm4d(const Vector<T> & m,const Vector<T> & a, const Vector<T> & b, const Vector<T> & c, const Vector<T> & d) {
        Vector<T> r = CopyVector(m);
        r.norm4d(a,b,c,d);
        return r;
    }

    template<typename T>
    inline Vector<T> normcdf(const Vector<T> & m) {
        Vector<T> r = CopyVector(m);
        r.normcdf();
        return r;
    }

    template<typename T>
    inline Vector<T> normcdfinv(const Vector<T> & m) {
        Vector<T> r = CopyVector(m);
        r.normcdfinv();
        return r;
    }

    template<typename T>
    inline Vector<T> rcbrt(const Vector<T> & m) {
        Vector<T> r = CopyVector(m);
        r.rcbrt();
        return r;
    }

    template<typename T>
    inline Vector<T> norm(const Vector<T> & m) {
        Vector<T> r = CopyVector(m);
        r.norm();
        return r;
    }

    template<typename T>
    inline Vector<T> rhypot(const Vector<T> & m, const Vector<T> & b) {
        Vector<T> r = CopyVector(m);
        r.rhypot(b);
        return r;
    }
    /*
    inline Vector<T> rint(const Vector<T> & m) {
        Vector<T> r = CopyVector(m);
        r.rint();
        return r;
    }*/

    template<typename T>
    inline Vector<T> rnorm3d(const Vector<T> & m,const Vector<T> & a, const Vector<T> & b, const Vector<T> & c) {
        Vector<T> r = CopyVector(m);
        r.rnorm3d(a,b,c);
        return r;
    }

    template<typename T>
    inline Vector<T> rnorm4d(const Vector<T> & m,const Vector<T> & a, const Vector<T> & b, const Vector<T> & c, const Vector<T> & d) {
        Vector<T> r = CopyVector(m);
        r.rnorm4d(a,b,c,d);
        return r;
    }
    /*
    inline Vector<T> round(const Vector<T> & m) {
        Vector<T> r = CopyVector(m);
        r.round();
        return r;
    }*/

    template<typename T>
    inline Vector<T> rnorm(const Vector<T> & m, const Vector<T> & b) {
        Vector<T> r = CopyVector(m);
        r.rnorm();
        return r;
    }

    template<typename T>
    inline Vector<T> tgamma(const Vector<T> & m) {
        Vector<T> r = CopyVector(m);
        r.tgamma();
        return r;
    }

    template<typename T>
    inline Vector<T> trunc(const Vector<T> & m) {
        Vector<T> r = CopyVector(m);
        r.tgamma();
        return r;
    }

    template<typename T>
    inline Vector<T> y0(const Vector<T> & m) {
        Vector<T> r = CopyVector(m);
        r.y0();
        return r;
    }

    template<typename T>
    inline Vector<T> y1(const Vector<T> & m) {
        Vector<T> r = CopyVector(m);
        r.y1();
        return r;
    }

    template<typename T>
    inline Vector<T> yn(const Vector<T> & m, int n) {
        Vector<T> r = CopyVector(m);
        r.yn(n);
        return r;
    }


    template<typename Type>
    inline Vector<Type> real(const ComplexVector<Type> & m) {    
        Vector<Type> r(m.size());
        r = m.real();
        return r;    
    }
    template<typename Type>
    inline Vector<Type> imag(const ComplexVector<Type> & m) {    
        Vector<Type> r(m.size());
        r = m.imag();
        return r;    
    }
    template<typename Type>
    inline Vector<Type> abs(const ComplexVector<Type> & m) {    
        Vector<Type> r(m.size());
        r = m.abs();
        return r;    
    }
    template<typename Type>
    inline Vector<Type> arg(const ComplexVector<Type> & m) {    
        Vector<Type> r(m.size());
        r = m.arg();
        return r;    
    }
    template<typename Type>
    inline ComplexVector<Type> norm(const ComplexVector<Type> & m) {    
        ComplexVector<Type> r(m.size());
        r = m.norm();
        return r;    
    }
    template<typename Type>
    inline ComplexVector<Type> conj(const ComplexVector<Type> & m) {    
        ComplexVector<Type> r(m.size());
        r = m.conj();
        return r;    
    }
    template<typename Type>
    inline ComplexVector<Type> proj(const ComplexVector<Type> & m) {
        ComplexVector<Type> r(m.size());
        r = m.proj();
        return r;    
    }
    template<typename Type>
    inline ComplexVector<Type> polar(const Vector<Type> & rho, const Vector<Type> & theta) {    
        ComplexVector<Type> r(rho.size());
        r =cvector_polarf(rho.dev,theta.dev,r.dev,r.size());
        return r;   
    }
    template<typename Type>
    inline ComplexVector<Type> exp(const ComplexVector<Type> & m) {    
        ComplexVector<Type> r(m.size());
        r = m.exp();
        return r;    
    }
    template<typename Type>
    inline ComplexVector<Type> log(const ComplexVector<Type> & m) {    
        ComplexVector<Type> r(m.size());
        r = m.log();
        return r;    
    }
    template<typename Type>
    inline ComplexVector<Type> log10(const ComplexVector<Type> & m) {    
        ComplexVector<Type> r(m.size());
        r = m.log10();
        return r;    
    }
    template<typename Type>
    inline ComplexVector<Type> pow(const ComplexVector<Type> & a, const ComplexVector<Type> & b) {    
        ComplexVector<Type> r(a.size());
        r = a.pow(b);
        return r;    
    }
    template<typename Type>
    inline ComplexVector<Type> pow(const ComplexVector<Type> & a, const complex<Type> & b) {    
        ComplexVector<Type> r(a.size());
        r = a.pow(b);
        return r;    
    }
    template<typename Type>
    inline ComplexVector<Type> sqrt(const ComplexVector<Type> & m) {    
        ComplexVector<Type> r(m.size());
        r = m.sqrt();
        return r;    
    }
    template<typename Type>
    inline ComplexVector<Type> sin(const ComplexVector<Type> & m) {    
        ComplexVector<Type> r(m.size());
        r = m.sin();
        return r;    
    }
    template<typename Type>
    inline ComplexVector<Type> cos(const ComplexVector<Type> & m) {    
        ComplexVector<Type> r(m.size());
        r = m.cos();
        return r;    
    }
    template<typename Type>
    inline ComplexVector<Type> tan(const ComplexVector<Type> & m) {    
        ComplexVector<Type> r(m.size());
        r = m.tan();
        return r;    
    }
    template<typename Type>
    inline ComplexVector<Type> asin(const ComplexVector<Type> & m) {    
        ComplexVector<Type> r(m.size());
        r = m.asin();
        return r;    
    }
    template<typename Type>
    inline ComplexVector<Type> acos(const ComplexVector<Type> & m) {    
        ComplexVector<Type> r(m.size());
        r = m.acos();
        return r;    
    }
    template<typename Type>
    inline ComplexVector<Type> atan(const ComplexVector<Type> & m) {    
        ComplexVector<Type> r(m.size());
        r = m.atan();
        return r;    
    }
    template<typename Type>
    inline ComplexVector<Type> sinh(const ComplexVector<Type> & m) {    
        ComplexVector<Type> r(m.size());
        r = m.sinh();
        return r;    
    }
    template<typename Type>
    inline ComplexVector<Type> cosh(const ComplexVector<Type> & m) {    
        ComplexVector<Type> r(m.size());
        r = m.cosh();
        return r;    
    }
    template<typename Type>
    inline ComplexVector<Type> tanh(const ComplexVector<Type> & m) {    
        ComplexVector<Type> r(m.size());
        r = m.tanh();
        return r;    
    }
    template<typename Type>
    inline ComplexVector<Type> asinh(const ComplexVector<Type> & m) {    
        ComplexVector<Type> r(m.size());
        r = m.asinh();
        return r;    
    }
    template<typename Type>
    inline ComplexVector<Type> acosh(const ComplexVector<Type> & m) {    
        ComplexVector<Type> r(m.size());
        r = m.acosh();
        return r;    
    }
    template<typename Type>
    inline ComplexVector<Type> atanh(const ComplexVector<Type> & m) {    
        ComplexVector<Type> r(m.size());
        r = m.atanh();
        return r;    
    }

    template<typename Type>
    inline ComplexVector<Type> sigmoid(const ComplexVector<Type> & m) {
        ComplexVector<Type> r(m.size());
        r.sigmoid();
        return r;
    }

    template<typename Type>
    inline ComplexVector<Type> sigmoid_deriv(const ComplexVector<Type> & m) {
        ComplexVector<Type> r(m.size());
        r.sigmoid_grad();
        return r;
    }

    template<typename Type>
    inline ComplexVector<Type> tanh_deriv(const ComplexVector<Type> & m) {
        ComplexVector<Type> r (m.size());
        r.tanh_grad();
        return r;
    }
    template<typename Type>
    inline ComplexVector<Type> hadamard(const ComplexVector<Type> & a, const ComplexVector<Type> & b) {
        ComplexVector<Type> r;
        r = a;
        r = r.hadamard(b);
        return r;
    }

    template<typename Type>
    inline ComplexVector<Type> ComplexVector<Type>::operator - () {
        ComplexVector<Type> r(N);
        r.copy(*this);
        cvector_r_mulf_const((cuComplex*)r.dev,complex<Type>(-1.0,0.0),r.dev,N);
        return r;
    }

    template<>
    inline ComplexVector<float> ComplexVector<float>::operator + (const ComplexVector<float> & a) {
        assert(size() == a.size());
        ComplexVector<float> r = CopyComplexVector(a);
        cublasStatus_t     status; 
        complex<float> alpha(1.0,0.0);
        status = cublasCaxpy(cublas->handle, N,(cuComplex*)&alpha,(cuComplex*)dev,1,(cuComplex*)r.dev,1);
        ASSERTOK(status);
        return r;
    }
    template<>
    inline ComplexVector<double> ComplexVector<double>::operator + (const ComplexVector<double> & a) {
        assert(size() == a.size());
        ComplexVector<double> r = CopyComplexVector(a);
        cublasStatus_t     status; 
        complex<double> alpha(1.0,0.0);
        status = cublasZaxpy(cublas->handle, N,(cuDoubleComplex*)&alpha,(cuDoubleComplex*)dev,1,(cuDoubleComplex*)r.dev,1);
        ASSERTOK(status);
        return r;
    }    

    template<>
    inline ComplexVector<float> ComplexVector<float>::operator - (const ComplexVector<float> & a) {
        assert(size() == a.size());
        ComplexVector<float> r = CopyComplexVector(*this);
        cublasStatus_t     status; 
        complex<float> alpha(-1.0f,0.0f);
        status = cublasCaxpy(cublas->handle,N,(cuComplex*)&alpha,(cuComplex*)a.dev,1,(cuComplex*)r.dev,1);
        ASSERTOK(status);
        return r;
    }

    template<>
    inline ComplexVector<double> ComplexVector<double>::operator - (const ComplexVector<double> & a) {
        assert(size() == a.size());
        ComplexVector<double> r = CopyComplexVector(*this);
        cublasStatus_t     status; 
        complex<double> alpha(-1.0f,0.0f);
        status = cublasZaxpy(cublas->handle,N,(cuDoubleComplex*)&alpha,(cuDoubleComplex*)a.dev,1,(cuDoubleComplex*)r.dev,1);
        ASSERTOK(status);
        return r;
    }
    
    template<>
    inline ComplexVector<float> ComplexVector<float>::operator * (const ComplexVector<float> & a) {
        assert(size() == a.size());
        ComplexVector<float> r = CopyComplexVector(a);
        cublasStatus_t     status;         
        status = cublasCdgmm(cublas->handle,CUBLAS_SIDE_RIGHT,1,N,(cuComplex*)dev,1,(cuComplex*)a.dev,1,(cuComplex*)r.dev,1);
        ASSERTOK(status);
        return r;
    }

    template<>
    inline ComplexVector<double> ComplexVector<double>::operator * (const ComplexVector<double> & a) {
        assert(size() == a.size());
        ComplexVector<double> r = CopyComplexVector(a);
        cublasStatus_t     status;         
        status = cublasZdgmm(cublas->handle,CUBLAS_SIDE_RIGHT,1,N,(cuDoubleComplex*)dev,1,(cuDoubleComplex*)a.dev,1,(cuDoubleComplex*)r.dev,1);
        ASSERTOK(status);
        return r;
    }
    
    template<typename Type>
    inline ComplexVector<Type> ComplexVector<Type>::operator - (const complex<Type> a) {    
        ComplexVector<Type> r(size());
        cvector_r_subf_const(dev,a,r.dev,size());
        return r;
    }

    template<typename Type>
    inline ComplexVector<Type> ComplexVector<Type>::operator * (const complex<Type> a) {    
        ComplexVector<Type> r(size());    
        cvector_r_mulf_const(dev,a,r.dev,r.size());
        return r;
    }

    template<typename Type>
    inline ComplexVector<Type> ComplexVector<Type>::operator / (const complex<Type> a) {    
        ComplexVector<Type> r(size());
        cvector_r_divf_const(dev,a,r.dev,r.size());
        return r;
    }

    template<typename Type>
    inline ComplexVector<Type> operator + (const ComplexVector<Type> & a, const ComplexVector<Type> & b) {
        ComplexVector<Type> r(a.N);    
        add(a,b,r);
        return r;
    }

    template<typename Type>
    inline ComplexVector<Type> operator - (const ComplexVector<Type> & a, const ComplexVector<Type> & b) {
        ComplexVector<Type> r(a.N);    
        sub(a,b,r);
        return r;
    }

    template<typename Type>
    inline ComplexVector<Type> operator * (const ComplexVector<Type> & a, const ComplexVector<Type>& b) {
        ComplexVector<Type> r(a.N);    
        mul(a,b,r);
        return r;
    }

    template<typename Type>
    inline void add(const ComplexVector<Type> & a, const ComplexVector<Type> & b, ComplexVector<Type> & r)
    {    
        assert(b.size() == a.size());        
        cublasStatus_t     status; 
        complex<Type> alpha(1.0,0.0);
        r.resize(a.N);
        r.copy(b);
        status = cublasCaxpy(cublas->handle, a.N,&alpha,a.dev,1,r.dev,1);
        ASSERTOK(status);    
    }
    template<>
    inline void add(const ComplexVector<float> & a, const ComplexVector<float> & b, ComplexVector<float> & r)
    {    
        assert(b.size() == a.size());        
        cublasStatus_t     status; 
        complex<float> alpha(1.0,0.0);
        r.resize(a.N);
        r.copy(b);
        status = cublasCaxpy(cublas->handle, a.N, (cuComplex*)&alpha, (cuComplex*)a.dev,1,(cuComplex*)r.dev,1);
        ASSERTOK(status);    
    }

    template<>
    inline void add(const ComplexVector<double> & a, const ComplexVector<double> & b, ComplexVector<double> & r)
    {    
        assert(b.size() == a.size());        
        cublasStatus_t     status; 
        complex<double> alpha(1.0,0.0);
        r.resize(a.N);
        r.copy(b);
        status = cublasZaxpy(cublas->handle, a.N,(cuDoubleComplex*)&alpha,(cuDoubleComplex*)a.dev,1,(cuDoubleComplex*)r.dev,1);
        ASSERTOK(status);    
    }

    template<typename Type>
    inline void sub(const ComplexVector<Type> & a, const ComplexVector<Type> & b, ComplexVector<Type> & r)
    {
        assert(b.size() == a.size());        
        cublasStatus_t     status; 
        complex<Type> alpha(-1.0,0.0);
        r.resize(a.N);
        r.copy(b);
        status = cublasSaxpy(cublas->handle,a.N,&alpha,a.dev,1,r.dev,1);
        ASSERTOK(status);
    }
    template<>
    inline void sub(const ComplexVector<float> & a, const ComplexVector<float> & b, ComplexVector<float> & r)
    {
        assert(b.size() == a.size());        
        cublasStatus_t     status; 
        complex<float> alpha(-1.0,0.0);
        r.resize(a.N);
        r.copy(b);
        status = cublasCaxpy(cublas->handle,a.N,(cuComplex*)&alpha,(cuComplex*)a.dev,1,(cuComplex*)r.dev,1);
        ASSERTOK(status);
    }
    
    template<>
    inline void sub(const ComplexVector<double> & a, const ComplexVector<double> & b, ComplexVector<double> & r)
    {
        assert(b.size() == a.size());        
        cublasStatus_t     status; 
        complex<double> alpha(-1.0,0.0);
        r.resize(a.N);
        r.copy(b);
        status = cublasZaxpy(cublas->handle,a.N,(cuDoubleComplex*)&alpha,(cuDoubleComplex*)a.dev,1,(cuDoubleComplex*)r.dev,1);
        ASSERTOK(status);
    }

    template<typename Type>
    inline void mul(const ComplexVector<Type> & a, const ComplexVector<Type> & b, ComplexVector<Type> & r)
    {
        assert(b.size() == a.size());
        cublasStatus_t     status;         
        r.resize(a.N);    
        status = cublasSdgmm(cublas->handle,CUBLAS_SIDE_RIGHT,1,a.N,b.dev,1,a.dev,1,r.dev,1);
        ASSERTOK(status);
    }
    
    template<>
    inline void mul(const ComplexVector<float> & a, const ComplexVector<float> & b, ComplexVector<float> & r)
    {
        assert(b.size() == a.size());
        cublasStatus_t     status;         
        r.resize(a.N);    
        status = cublasCdgmm(cublas->handle,CUBLAS_SIDE_RIGHT,1,a.N,(cuComplex*)b.dev,1,(cuComplex*)a.dev,1,(cuComplex*)r.dev,1);
        ASSERTOK(status);
    }
    template<>
    inline void mul(const ComplexVector<double> & a, const ComplexVector<double> & b, ComplexVector<double> & r)
    {
        assert(b.size() == a.size());
        cublasStatus_t     status;         
        r.resize(a.N);    
        status = cublasZdgmm(cublas->handle,CUBLAS_SIDE_RIGHT,1,a.N,(cuDoubleComplex*)b.dev,1,(cuDoubleComplex*)a.dev,1,(cuDoubleComplex*)r.dev,1);
        ASSERTOK(status);
    }


    template<typename Type>
    inline void div(const ComplexVector<Type> & a, const ComplexVector<Type> & b, ComplexVector<Type> & r)
    {
        assert(b.size() == a.size());
        cublasStatus_t     status;         
        r.resize(a.N);    
        status = cublasSdgmm(cublas->handle,CUBLAS_SIDE_RIGHT,1,a.N,b.dev,1,a.dev,1,r.dev,1);
        ASSERTOK(status);
    }
    
    template<>
    inline void div(const Vector<double> & a, const Vector<double> & b, Vector<double> & r)
    {
        assert(b.size() == a.size());
        cublasStatus_t     status;         
        r.resize(a.N);    
        status = cublasDdgmm(cublas->handle,CUBLAS_SIDE_RIGHT,1,a.N,b.dev,1,a.dev,1,r.dev,1);
        ASSERTOK(status);
    }

    template<>
    inline void div(const ComplexVector<float> & a, const ComplexVector<float> & b, ComplexVector<float> & r)
    {
        assert(b.size() == a.size());
        cublasStatus_t     status;         
        r.resize(a.N);    
        status = cublasCdgmm(cublas->handle,CUBLAS_SIDE_RIGHT,1,a.N,(cuComplex*)b.dev,1,(cuComplex*)a.dev,1,(cuComplex*)r.dev,1);
        ASSERTOK(status);
    }
    template<>
    inline void div(const ComplexVector<double> & a, const ComplexVector<double> & b, ComplexVector<double> & r)
    {
        assert(b.size() == a.size());
        cublasStatus_t     status;         
        r.resize(a.N);    
        status = cublasZdgmm(cublas->handle,CUBLAS_SIDE_RIGHT,1,a.N,(cuDoubleComplex*)b.dev,1,(cuDoubleComplex*)a.dev,1,(cuDoubleComplex*)r.dev,1);
        ASSERTOK(status);
    }

    template<typename Type>
    inline ComplexVector<Type> operator + (const complex<Type> v, const ComplexVector<Type> & a)
    {    
        ComplexVector<Type> r(a.N);
        cvector_r_addf_const(a.dev,v,r.dev,a.N);
        return r;
    }


    template<typename Type>
    inline ComplexVector<Type> operator - (const ComplexVector<Type> & a, const complex<Type> v)
    {    
        ComplexVector<Type> r(a.N);    
        cvector_r_subf_const(a.dev,v,r.dev,a.N);
        return r;
    }

    template<typename Type>
    inline ComplexVector<Type> operator - (const complex<Type> v, const ComplexVector<Type> & a)
    {    
        ComplexVector<Type> r(a.N);        
        cvector_r_subf_const(a.dev,v,r.dev,a.N);
        return r;
    }

    template<typename Type>
    inline ComplexVector<Type> operator * (const ComplexVector<Type> & a, const complex<Type> v)
    {    
        ComplexVector<Type> r(a.N);    
        cvector_r_mulf_const(a.dev,v,r.dev,a.N);
        //cublasSscal(cublas->handle,a.N,&v,r.dev,1);
        return r;
    }

    template<typename Type>
    inline ComplexVector<Type> operator * (const complex<Type> v, const ComplexVector<Type> & a)
    {    
        ComplexVector<Type> r(a.N);    
        cvector_r_mulf_const(a.dev,v,r.dev,a.N);
        //cublasSscal(cublas->handle,a.N,&v,r.dev,1);
        return r;
    }

    template<typename Type>
    inline ComplexVector<Type> operator / (const ComplexVector<Type> & a, const Type v)
    {    
        ComplexVector<Type> r(a.N);    
        cvector_r_divf_const(a.dev,v,r.dev,a.N);    
        return r;
    }


    template<typename Type>
    inline Matrix<Type> hadamard(const Matrix<Type> & a, const Matrix<Type> &b)
    {
        Matrix<Type> r = CopyMatrix(a);    
        r.hadamard(b);
        return r;
    }

    template<typename Type>
    inline void hadamard_fast(Matrix<Type> & r, const Matrix<Type> & a, const Matrix<Type> &b)
    {           
        assert(a.M == b.M && a.N == b.N && a.M == r.M && a.N == r.N);          
        vector_r_mulf(a.dev,b.dev,r.dev,a.M*a.N);                
    }

    
    template<typename T>
    inline void add(const Matrix<T> & a, const Matrix<T> & b, Matrix<T> & r)
    {
        assert(a.M == b.M && a.N == b.N);                
        int m = a.M;
        int k = a.N;    
        int n = a.N;
        T alpha = 1.0;
        T beta  = 1.0;
        r.resize(a.M,b.N);            
        cublasStatus_t status = cublasSgeam(cublas->handle,CUBLAS_OP_N,CUBLAS_OP_N,n,m,&alpha,b.dev,n,&beta,a.dev,n,r.dev,n);
        ASSERTOK(status);        
    }
    
    
    template<>
    inline void add(const Matrix<float> & a, const Matrix<float> & b, Matrix<float> & r)
    {
        assert(a.M == b.M && a.N == b.N);                
        int m = a.M;
        int k = a.N;    
        int n = a.N;
        float alpha = 1.0f;
        float beta  = 1.0f;
        r.resize(a.M,b.N);            
        cublasStatus_t status = cublasSgeam(cublas->handle,CUBLAS_OP_N,CUBLAS_OP_N,n,m,&alpha,b.dev,n,&beta,a.dev,n,r.dev,n);
        ASSERTOK(status);        
    }

    template<>
    inline void add(const Matrix<double> & a, const Matrix<double> & b, Matrix<double> & r)
    {
        assert(a.M == b.M && a.N == b.N);                
        int m = a.M;
        int k = a.N;    
        int n = a.N;
        double alpha = 1.0;
        double beta  = 1.0;
        r.resize(a.M,b.N);            
        cublasStatus_t status = cublasDgeam(cublas->handle,CUBLAS_OP_N,CUBLAS_OP_N,n,m,&alpha,b.dev,n,&beta,a.dev,n,r.dev,n);
        ASSERTOK(status);        
    }

    
    template<typename T>
    inline void sub(const Matrix<T> & a, const Matrix<T> & b, Matrix<T> & r)
    {
        assert(a.M == b.M && a.N == b.N);                
        int m = a.M;
        int k = a.N;    
        int n = a.N;
        T alpha = -1.0;
        T beta  = 1.0;        
        r.resize(a.M,b.N);    
        cublasStatus_t status = cublasSgeam(cublas->handle,CUBLAS_OP_N,CUBLAS_OP_N,n,m,&alpha,b.dev,n,&beta,a.dev,n,r.dev,n);
        ASSERTOK(status);        
    }
    
    template<>
    inline void sub(const Matrix<float> & a, const Matrix<float> & b, Matrix<float> & r)
    {
        assert(a.M == b.M && a.N == b.N);                
        int m = a.M;
        int k = a.N;    
        int n = a.N;
        float alpha = -1.0f;
        float beta  = 1.0f;        
        r.resize(a.M,b.N);    
        cublasStatus_t status = cublasSgeam(cublas->handle,CUBLAS_OP_N,CUBLAS_OP_N,n,m,&alpha,b.dev,n,&beta,a.dev,n,r.dev,n);
        ASSERTOK(status);        
    }

    template<>
    inline void sub(const Matrix<double> & a, const Matrix<double> & b, Matrix<double> & r)
    {
        assert(a.M == b.M && a.N == b.N);                
        int m = a.M;
        int k = a.N;    
        int n = a.N;
        double alpha = -1.0;
        double beta  = 1.0;        
        r.resize(a.M,b.N);    
        cublasStatus_t status = cublasDgeam(cublas->handle,CUBLAS_OP_N,CUBLAS_OP_N,n,m,&alpha,b.dev,n,&beta,a.dev,n,r.dev,n);
        ASSERTOK(status);        
    }

    
    template<typename T>
    inline void mul(const Matrix<T> & a, const Matrix<T> & b, Matrix<T> & r)
    {
        assert(a.N == b.M);             
        cublasOperation_t ta = CUBLAS_OP_N;
        cublasOperation_t tb = CUBLAS_OP_N;          
        int m = a.M;
        int k = a.N;    
        int n = b.N;
        T alpha=1.0;
        T beta=0.0;
        r.resize(a.M,b.N);    
        cublasStatus_t status = cublasSgemm(cublas->handle,ta,tb,n,m,k,&alpha,b.dev,n,a.dev,k,&beta,r.dev,n);                        
        ASSERTOK(status);    
    }
    

    template<>
    inline void mul(const Matrix<float> & a, const Matrix<float> & b, Matrix<float> & r)
    {
        assert(a.N == b.M);             
        cublasOperation_t ta = CUBLAS_OP_N;
        cublasOperation_t tb = CUBLAS_OP_N;          
        int m = a.M;
        int k = a.N;    
        int n = b.N;
        float alpha=1.0f;
        float beta=0.0f;
        r.resize(a.M,b.N);    
        cublasStatus_t status = cublasSgemm(cublas->handle,ta,tb,n,m,k,&alpha,b.dev,n,a.dev,k,&beta,r.dev,n);                        
        ASSERTOK(status);    
    }

    template<>
    inline void mul(const Matrix<double> & a, const Matrix<double> & b, Matrix<double> & r)
    {
        assert(a.N == b.M);             
        cublasOperation_t ta = CUBLAS_OP_N;
        cublasOperation_t tb = CUBLAS_OP_N;          
        int m = a.M;
        int k = a.N;    
        int n = b.N;
        double alpha=1.0;
        double beta=0.0;
        r.resize(a.M,b.N);    
        cublasStatus_t status = cublasDgemm(cublas->handle,ta,tb,n,m,k,&alpha,b.dev,n,a.dev,k,&beta,r.dev,n);                        
        ASSERTOK(status);    
    }

    template<typename T>
    inline Matrix<T> operator + (const Matrix<T> & a, const Matrix<T> & b)
    {
        assert(a.M == b.M && a.N == b.N);
        Matrix<T> r(a.M,a.N);    
        add(a,b,r);
        return r;
    }

    template<typename T>
    inline void add_matrix(Matrix<T> & r, const Matrix<T> & a, const Matrix<T> & b)
    {
        assert(a.M == b.M && a.N == b.N);        
        add(a,b,r);        
    }

    template<typename T>
    inline Matrix<T> operator - (const Matrix<T> & a, const Matrix<T> & b)
    {
        assert(a.M == b.M && a.N == b.N);
        Matrix<T> r(a.M,a.N);    
        sub(a,b,r);
        return r;
    }

    template<typename T>
    inline void sub_matrix(Matrix<T> & r, const Matrix<T> & a, const Matrix<T> & b)
    {
        assert(a.M == b.M && a.N == b.N);        
        sub(a,b,r);        
    }

    template<typename T>
    inline Matrix<T> operator * (const Matrix<T> & a, const Matrix<T> & b)
    {        
        Matrix<T> r(a.M,a.N);    
        mul(a,b,r);
        return r;
    }

    template<typename T>
    inline void mul_matrix(Matrix<T> & r, const Matrix<T> & a, const Matrix<T> & b)
    {                
        mul(a,b,r);        
    }
    template<typename T>
    inline Matrix<T> operator + (const Matrix<T> & a, const T v)
    {    
        Matrix<T> r(a.M,a.N);        
        matrix_r_addf_const(a.dev,v,r.dev,a.M,a.N);    
        return r;
    }

    template<typename T>
    inline Matrix<T> operator + (const T v, const Matrix<T> & a)
    {
        Matrix<T> r(a.M,a.N);    
        matrix_r_addf_const(a.dev,v,r.dev,a.M,a.N);    
        return r;
    }

    template<typename T>
    inline Matrix<T> operator - (const Matrix<T> & a, const T v)
    {
        Matrix<T> r(a.M,a.N);    
        matrix_r_subf_const(a.dev,v,r.dev,a.M,a.N);
        return r;
    }

    template<typename T>
    inline Matrix<T> operator - (const T v, const Matrix<T> & a)
    {    
        Matrix<T> r = CopyMatrix(a);
        r = -r;
        matrix_r_addf_const<T>(r.dev,v,r.dev,a.M,a.N);
        return r;
    }

    template<typename T>
    inline Matrix<T> operator * (const Matrix<T> & a, const T v)
    {
        Matrix<T> r = CopyMatrix(a);        
        matrix_r_mulf_const(a.dev,v,r.dev,a.M,a.N);
        return r;
    }

    template<typename T>
    inline void mul_const (const Matrix<T> & a, const T v, Matrix<T> & b)
    {        
        matrix_r_mulf_const(a.dev,v,b.dev,a.M,a.N);        
    }

    template<typename T>
    inline Matrix<T> operator * (const T v, const Matrix<T> & a)
    {    
        Matrix<T> r = CopyMatrix(a);    
        matrix_r_mulf_const(a.dev,v,r.dev,a.M,a.N);
        return r;
    }

    template<typename T>
    inline Matrix<T> operator / (const Matrix<T> & a, T v)
    {    
        Matrix<T> r = CopyMatrix(a);    
        matrix_r_divf_const(a.dev,v,r.dev,a.M,a.N);
        return r;
    }

    template<typename T>
    inline Matrix<T> operator / (const T a, const Matrix<T> & b)
    {    
        Matrix<T> r = CopyMatrix(b);    
        matrix_r_divf_matrix_const(b.dev,a,r.dev,b.M,b.N);
        return r;
    }

    template<typename T>
    inline Matrix<T> operator % (const Matrix<T> & a, const T v)
    {    
        Matrix<T> r(a.M,a.N);        
        matrix_r_modf_const(a.dev,v,r.dev,a.M,a.N);    
        return r;
    }


    template<typename T>
    Matrix<T> Matrix<T>::operator + (const Matrix<T> & b) {        
        Matrix<T> r(M,N);                
        add(*this,b,r);
        return r;
    }

    template<typename T>
    Matrix<T> Matrix<T>::operator - (const Matrix<T> & b) {                 
        Matrix<T> r(M,N);                
        sub(*this,b,r);
        return r;
    }

    template<typename T>
    Matrix<T> Matrix<T>::operator * (const Matrix<T> & m) 
    {   
        Matrix<T> r(M,N);      
        mul(*this,m,r);
        return r;
    }


    template<typename T>
    inline Matrix<T> abs(const Matrix<T> & m) {    
        Matrix<T> r = CopyMatrix(m);
        r.abs();
        return r;
    }

    template<typename T>
    inline Matrix<T> exp(const Matrix<T> & m) {
        Matrix<T> r = CopyMatrix(m);
        r.exp(); 
        return r;   
    }    

    template<typename T>
    inline Matrix<T> log2(const Matrix<T> & m) {
        Matrix<T> r = CopyMatrix(m);
        r.log2();
        return r;
    }

    template<typename T>
    inline Matrix<T> log10(const Matrix<T> & m) {
        Matrix<T> r = CopyMatrix(m);
        r.log10();
        return r;
    }

    template<typename T>
    inline Matrix<T> pow(const Matrix<T> & m,T p) {
        Matrix<T> r = CopyMatrix(m);
        r.pow(p);
        return r;
    }
    
    template<typename T>
    inline Matrix<T> pow(const Matrix<T> & m,const Matrix<T> & p) {
        Matrix<T> r = CopyMatrix(m);
        r.pow(p);
        return r;
    }

    template<typename T>
    inline Matrix<T> sqrt(const Matrix<T> & m) {
        Matrix<T> r = CopyMatrix(m);
        r.sqrt();
        return r;
    }

    template<typename T>
    inline Matrix<T> rsqrt(const Matrix<T> & m) {
        Matrix<T> r = CopyMatrix(m);
        r.rsqrt();
        return r;
    }

    template<typename T>
    inline Matrix<T> sin(const Matrix<T> & m) {
        Matrix<T> r = CopyMatrix(m);
        r.sin();
        return r;
    }

    template<typename T>
    inline Matrix<T> cos(const Matrix<T> & m) {
        Matrix<T> r = CopyMatrix(m);
        r.cos();
        return r;
    }

    template<typename T>
    inline Matrix<T> tan(const Matrix<T> & m) {
        Matrix<T> r = CopyMatrix(m);
        r.tan();
        return r;
    }

    template<typename T>
    inline Matrix<T> asin(const Matrix<T> & m) {
        Matrix<T> r = CopyMatrix(m);
        r.asin();
        return r;
    }


    template<typename T>
    inline Matrix<T> acos(const Matrix<T> & m) {
        Matrix<T> r = CopyMatrix(m);
        r.acos();
        return r;
    }

    template<typename T>
    inline Matrix<T> atan(const Matrix<T> & m) {
        Matrix<T> r = CopyMatrix(m);
        r.atan();
        return r;
    }

    template<typename T>
    inline Matrix<T> sinh(const Matrix<T> & m) {
        Matrix<T> r = CopyMatrix(m);
        r.sinh();
        return r;
    }

    template<typename T>
    inline Matrix<T> cosh(const Matrix<T> & m) {
        Matrix<T> r = CopyMatrix(m);
        r.cosh();
        return r;
    }

    template<typename T>
    inline Matrix<T> tanh(const Matrix<T> & m) {
        Matrix<T> r = CopyMatrix(m);
        r.tanh();
        return r;
    }

    template<typename T>
    inline Matrix<T> asinh(const Matrix<T> & m) {
        Matrix<T> r = CopyMatrix(m);
        r.asinh();
        return r;
    }

    template<typename T>
    inline Matrix<T> acosh(const Matrix<T> & m) {
        Matrix<T> r = CopyMatrix(m);
        r.acosh();
        return r;
    }

    template<typename T>
    inline Matrix<T> atanh(const Matrix<T> & m) {
        Matrix<T> r = CopyMatrix(m);
        r.atanh();
        return r;
    }    

    template<typename T>
    inline Matrix<T> atan2(const Matrix<T> & m,T v) {
        Matrix<T> r = CopyMatrix(m);
        r.atan2(v);
        return r;
    }

    template<typename T>
    inline Matrix<T> atan2(const Matrix<T> & m,const Matrix<T> & v) {
        Matrix<T> r = CopyMatrix(m);
        r.atan2(v);
        return r;
    }

    template<typename T>
    inline Matrix<T> sigmoid(const Matrix<T> & m) {
        Matrix<T> r = CopyMatrix(m);
        r.sigmoid();
        return r;
    }

    template<typename T>
    inline Matrix<T> sigmoid_deriv(const Matrix<T> & m) {
        Matrix<T> r = CopyMatrix(m);
        r.sigmoid_deriv();
        return r;
    }

    template<typename T>
    inline Matrix<T> relu(const Matrix<T> & m) {
        Matrix<T> r = CopyMatrix(m);
        r.relu();
        return r;
    }

    template<typename T>
    inline Matrix<T> relu_deriv(const Matrix<T> & m) {
        Matrix<T> r = CopyMatrix(m);
        r.relu_deriv();
        return r;
    }

    template<typename T>
    inline Matrix<T> softmax(const Matrix<T> & m) {
        Matrix<T> r = CopyMatrix(m);
        r.softmax();
        return r;
    }

    template<typename T>
    inline Matrix<T> tanh_deriv(const Matrix<T> & m) {
        Matrix<T> r = CopyMatrix(m);
        r.tanh_deriv();
        return r;
    }

    template<typename T>
    inline Matrix<T> cbrt(const Matrix<T> & m) {
        Matrix<T> r = CopyMatrix(m);
        r.cbrt();
        return r;
    }

    template<typename T>
    inline Matrix<T> ceil(const Matrix<T> & m) {
        Matrix<T> r = CopyMatrix(m);
        r.ceil();
        return r;
    }
    /*
    inline Matrix<T> copysign(const Matrix<T> & m, T x) {
        Matrix<T> r = CopyMatrix(m);
        r.copysign(x);
        return r;
    }
    inline Matrix<T> copysign(const Matrix<T> & m, const Matrix<T> & x) {
        Matrix<T> r = CopyMatrix(m);
        r.copysign(x);
        return r;
    }
    */

    template<typename T>
    inline Matrix<T> cospi(const Matrix<T> & m) {
        Matrix<T> r = CopyMatrix(m);
        r.cospi();
        return r;
    }

    template<typename T>
    inline Matrix<T> cyl_bessel_i0f(const Matrix<T> & m) {
        Matrix<T> r = CopyMatrix(m);
        r.cyl_bessel_i0();
        return r;
    }

    template<typename T>
    inline Matrix<T> cyl_bessel_i1f(const Matrix<T> & m) {
        Matrix<T> r = CopyMatrix(m);
        r.cyl_bessel_i1();
        return r;
    }

    template<typename T>
    inline Matrix<T> erfcf(const Matrix<T> & m) {
        Matrix<T> r = CopyMatrix(m);
        r.erfc();
        return r;
    }

    template<typename T>
    inline Matrix<T> erfcx(const Matrix<T> & m) {
        Matrix<T> r = CopyMatrix(m);
        r.erfcx();
        return r;
    }

    template<typename T>
    inline Matrix<T> erfcinv(const Matrix<T> & m) {
        Matrix<T> r = CopyMatrix(m);
        r.erfcinv();
        return r;
    }

    template<typename T>
    inline Matrix<T> erf(const Matrix<T> & m) {
        Matrix<T> r = CopyMatrix(m);
        r.erf();
        return r;
    }

    template<typename T>
    inline Matrix<T> erfinv(const Matrix<T> & m) {
        Matrix<T> r = CopyMatrix(m);
        r.erfinv();
        return r;
    }

    template<typename T>
    inline Matrix<T> exp10(const Matrix<T> & m) {
        Matrix<T> r = CopyMatrix(m);
        r.exp10();
        return r;
    }

    template<typename T>
    inline Matrix<T> exp2(const Matrix<T> & m) {
        Matrix<T> r = CopyMatrix(m);
        r.exp2();
        return r;
    }

    template<typename T>
    inline Matrix<T> expm1(const Matrix<T> & m) {
        Matrix<T> r = CopyMatrix(m);
        r.expm1();
        return r;
    }

    template<typename T>
    inline Matrix<T> fabs(const Matrix<T> & m) {
        Matrix<T> r = CopyMatrix(m);
        r.fabs();
        return r;
    }

    template<typename T>
    inline Matrix<T> fdim(const Matrix<T> & m, const Matrix<T> & b) {
        Matrix<T> r = CopyMatrix(m);
        r.fdim(b);
        return r;
    }

    template<typename T>
    inline Matrix<T> fmod(const Matrix<T> & m, const Matrix<T> & b) {
        Matrix<T> r = CopyMatrix(m);
        r.fmod(b);
        return r;
    }

    template<typename T>
    inline Matrix<T> hypot(const Matrix<T> & m, const Matrix<T> & b) {
        Matrix<T> r = CopyMatrix(m);
        r.hypot(b);
        return r;
    }

    template<typename T>
    inline Matrix<T> ilogb(const Matrix<T> & m) {
        Matrix<T> r = CopyMatrix(m);
        r.ilogb();
        return r;
    }

    template<typename T>
    inline Matrix<T> j0(const Matrix<T> & m) {
        Matrix<T> r = CopyMatrix(m);
        r.j0();
        return r;
    }

    template<typename T>
    inline Matrix<T> j1(const Matrix<T> & m) {
        Matrix<T> r = CopyMatrix(m);
        r.j1();
        return r;
    }

    template<typename T>
    inline Matrix<T> jn(const Matrix<T> & m, int n) {
        Matrix<T> r = CopyMatrix(m);
        r.jn(n);
        return r;
    }

    template<typename T>
    inline Matrix<T> lgamma(const Matrix<T> & m) {
        Matrix<T> r = CopyMatrix(m);
        r.lgamma();
        return r;
    }

    template<typename T>
    inline Matrix<T> log1p(const Matrix<T> & m) {
        Matrix<T> r = CopyMatrix(m);
        r.log1p();
        return r;
    }

    template<typename T>
    inline Matrix<T> logb(const Matrix<T> & m) {
        Matrix<T> r = CopyMatrix(m);
        r.logb();
        return r;
    }

    template<typename T>
    inline Matrix<T> norm3d(const Matrix<T> & m,const Matrix<T> & a, const Matrix<T> & b, const Matrix<T> & c) {
        Matrix<T> r = CopyMatrix(m);
        r.norm3d(a,b,c);
        return r;
    }

    template<typename T>
    inline Matrix<T> norm4d(const Matrix<T> & m,const Matrix<T> & a, const Matrix<T> & b, const Matrix<T> & c, const Matrix<T> & d) {
        Matrix<T> r = CopyMatrix(m);
        r.norm4d(a,b,c,d);
        return r;
    }

    template<typename T>
    inline Matrix<T> normcdf(const Matrix<T> & m) {
        Matrix<T> r = CopyMatrix(m);
        r.normcdf();
        return r;
    }

    template<typename T>
    inline Matrix<T> normcdfinv(const Matrix<T> & m) {
        Matrix<T> r = CopyMatrix(m);
        r.normcdfinv();
        return r;
    }

    template<typename T>
    inline Matrix<T> norm(const Matrix<T> & m) {
        Matrix<T> r = CopyMatrix(m);
        r.norm();
        return r;
    }

    template<typename T>
    inline Matrix<T> rcbrt(const Matrix<T> & m) {
        Matrix<T> r = CopyMatrix(m);
        r.rcbrt();
        return r;
    }

    template<typename T>
    inline Matrix<T> rhypot(const Matrix<T> & m, const Matrix<T> & b) {
        Matrix<T> r = CopyMatrix(m);
        r.rhypot(b);
        return r;
    }

    template<typename T>
    inline Matrix<T> rnorm3d(const Matrix<T> & m,const Matrix<T> & a, const Matrix<T> & b, const Matrix<T> & c) {
        Matrix<T> r = CopyMatrix(m);
        r.rnorm3d(a,b,c);
        return r;
    }

    template<typename T>
    inline Matrix<T> rnorm4d(const Matrix<T> & m,const Matrix<T> & a, const Matrix<T> & b, const Matrix<T> & c, const Matrix<T> & d) {
        Matrix<T> r = CopyMatrix(m);
        r.rnorm4d(a,b,c,d);
        return r;
    }

    template<typename T>
    inline Matrix<T> rnorm(const Matrix<T> & m, const Matrix<T> & b) {
        Matrix<T> r = CopyMatrix(m);
        r.rnorm();
        return r;
    }

    template<typename T>
    inline Matrix<T> tgamma(const Matrix<T> & m) {
        Matrix<T> r = CopyMatrix(m);
        r.tgamma();
        return r;
    }

    template<typename T>
    inline Matrix<T> trunc(const Matrix<T> & m) {
        Matrix<T> r = CopyMatrix(m);
        r.tgamma();
        return r;
    }

    template<typename T>
    inline Matrix<T> y0(const Matrix<T> & m) {
        Matrix<T> r = CopyMatrix(m);
        r.y0();
        return r;
    }

    template<typename T>
    inline Matrix<T> y1(const Matrix<T> & m) {
        Matrix<T> r = CopyMatrix(m);
        r.y1();
        return r;
    }

    template<typename T>
    inline Matrix<T> yn(const Matrix<T> & m, int n) {
        Matrix<T> r = CopyMatrix(m);
        r.yn(n);
        return r;
    }

    
    template<typename T>
    inline T sum(const Matrix<T> & m) {
        T r = 0;
        cublasStatus_t     status; 
        status = cublasSasum(cublas->handle,m.M*m.N,m.dev,1,&r);
        ASSERTOK(status);
        return r;
    }
    
    template<>
    inline double sum(const Matrix<double> & m) {
        double r = 0;
        cublasStatus_t     status; 
        status = cublasDasum(cublas->handle,m.M*m.N,m.dev,1,&r);
        ASSERTOK(status);
        return r;
    }    
    template<>
    inline float sum(const Matrix<float> & v) {
        float r = 0;
        cublasStatus_t     status; 
        status = cublasSasum(cublas->handle,v.N,v.dev,1,&r);
        ASSERTOK(status);
        return r;
    }

    template<typename T>
    inline Vector<T> matvec(const Matrix<T> & a, const Vector<T> & v,bool transa, T alpha=1.0, T beta=0.0) {    
        Vector<T> R(v.size());
        int m = a.M;
        int n = a.N;
        cublasOperation_t ta = CUBLAS_OP_N;
        cublasStatus_t status = cublasSgemv(cublas->handle,ta,n,m,&alpha,a.dev,n,a.dev,1,&beta,R.dev,1);
        ASSERTOK(status);
        return R;
    }
    
    template<>
    inline Vector<float> matvec(const Matrix<float> & a, const Vector<float> & v,bool transa, float alpha, float beta) {    
        Vector<float> R(v.size());
        int m = a.M;
        int n = a.N;
        cublasOperation_t ta = CUBLAS_OP_N;
        cublasStatus_t status = cublasSgemv(cublas->handle,ta,n,m,&alpha,a.dev,n,a.dev,1,&beta,R.dev,1);
        ASSERTOK(status);
        return R;
    }

    template<>
    inline Vector<double> matvec(const Matrix<double> & a, const Vector<double> & v,bool transa, double alpha, double beta) {    
        Vector<double> R(v.size());
        int m = a.M;
        int n = a.N;
        cublasOperation_t ta = CUBLAS_OP_N;
        cublasStatus_t status = cublasDgemv(cublas->handle,ta,n,m,&alpha,a.dev,n,a.dev,1,&beta,R.dev,1);
        ASSERTOK(status);
        return R;
    }

    template<typename T>
    inline Matrix<T> Matrix<T>::operator -() {
        Matrix<T> r = CopyMatrix(*this);
        matrix_r_mulf_const(r.dev,(T)-1.0,r.dev,r.M,r.N);    
        return r;
    }

    template<typename T>
    inline Matrix<T> Matrix<T>::eval() {
        Matrix<T> r = CopyMatrix(*this);
        return r;
    }
    template<typename T>
    inline Matrix<T> Matrix<T>::eval() const {
        Matrix<T> r = CopyMatrix(*this);
        return r;
    }
    

    // C^T = B^T * A^T
    // assumes a and b are in row major order        
    template<typename T>
    inline Matrix<T> matmul(const Matrix<T> & a, const Matrix<T> & b, T alpha=1.0, T beta=0.0) {    
        assert(a.N == b.M);     
        Matrix<T> C(a.M,b.N);
        cublasOperation_t ta = CUBLAS_OP_N;
        cublasOperation_t tb = CUBLAS_OP_N;      
        int m = a.M;
        int k = a.N;
        int n = b.N;    
        cublasStatus_t status = cublasSgemm(cublas->handle,ta,tb,n,m,k,&alpha,b.dev,n,a.dev,k,&beta,C.dev,n);                        
        ASSERTOK(status);
        return C;        
    }    
    template<>
    inline Matrix<float> matmul(const Matrix<float> & a, const Matrix<float> & b, float alpha, float beta) {    
        assert(a.N == b.M);     
        Matrix<float> C(a.M,b.N);
        cublasOperation_t ta = CUBLAS_OP_N;
        cublasOperation_t tb = CUBLAS_OP_N;      
        int m = a.M;
        int k = a.N;
        int n = b.N;    
        cublasStatus_t status = cublasSgemm(cublas->handle,ta,tb,n,m,k,&alpha,b.dev,n,a.dev,k,&beta,C.dev,n);                        
        ASSERTOK(status);
        return C;        
    }

    template<>
    inline Matrix<double> matmul(const Matrix<double> & a, const Matrix<double> & b, double alpha, double beta) {    
        assert(a.N == b.M);     
        Matrix<double> C(a.M,b.N);
        cublasOperation_t ta = CUBLAS_OP_N;
        cublasOperation_t tb = CUBLAS_OP_N;      
        int m = a.M;
        int k = a.N;
        int n = b.N;    
        cublasStatus_t status = cublasDgemm(cublas->handle,ta,tb,n,m,k,&alpha,b.dev,n,a.dev,k,&beta,C.dev,n);                        
        ASSERTOK(status);
        return C;        
    }
    /*
    // C^T = A*B
    // A and B are in column major order
    template<typename T>
    inline Matrix<T> matmulTT(const Matrix<T> & a,const Matrix<T> & b, T alpha=1.0, T beta=0.0) {         
        assert(a.N == b.M);    
        cublasOperation_t ta = CUBLAS_OP_T;
        cublasOperation_t tb = CUBLAS_OP_T;      
        Matrix<T> C(a.N,b.M);            
        int m = a.M;
        int k = a.N;
        int n = b.N; 
        cublasStatus_t status = cublasSgemm(cublas->handle,ta,tb,n,m,k,&alpha,b.dev,n,a.dev,k,&beta,C.dev,n);                                            
        ASSERTOK(status);
        return C;    
    }
    // C^T = B^T*A
    // A is in column major order
    template<typename T>
    inline Matrix<T> matmulTN(const Matrix<T> & a,const Matrix<T> & b, T alpha=1.0, T beta=0.0) {         
        assert(a.N == b.N);    
        cublasOperation_t ta = CUBLAS_OP_N;
        cublasOperation_t tb = CUBLAS_OP_T;               
        Matrix<T> C(a.N,b.N); 
        int m = a.M;
        int k = a.N;
        int n = b.N; 
        cublasStatus_t status = cublasSgemm(cublas->handle,ta,tb,n,k,m,&alpha,b.dev,n,a.dev,k,&beta,C.dev,n); 
        ASSERTOK(status);
        return C;        
    }
    // C^T = B*A^T
    // A is in column major order
    template<typename T>
    inline Matrix<T> matmulNT(const Matrix<T> & a,const Matrix<T> & b, T alpha=1.0, T beta=0.0) {         
        assert(a.M == b.M);    
        cublasOperation_t ta = CUBLAS_OP_T;
        cublasOperation_t tb = CUBLAS_OP_N;            
        Matrix<T> C(a.M,b.M);        
        int m = a.M;
        int k = a.N;
        int n = b.N; 
        cublasStatus_t status = cublasSgemm(cublas->handle,ta,tb,m,n,k,&alpha,b.dev,n,a.dev,k,&beta,C.dev,n);                        
        ASSERTOK(status);
        return C;    
    }
    */

    template<typename T>
    std::ostream& operator << (std::ostream& o, Matrix<T> & m)
    {
        m.print();
        return o;
    }
}