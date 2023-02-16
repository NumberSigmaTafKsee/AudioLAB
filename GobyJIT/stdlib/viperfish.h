// NVRTC kernels
// Expressions exprtk
// 
#ifndef __VIPERFISH_H
#define __VIPERFISH_H

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
#include <ccomplex>

#include "cuda_runtime.h"
#include "cublas_v2.h"
#include "curand.h"

#include "vector_float.h"
//#include "matrix_float.h"
//#include "cube_float.h"
//#include "field_float.h"
//#include "array_float.h"


#define ASSERTOK(status) (assert(status == CUBLAS_STATUS_SUCCESS))

typedef int array_index;



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
    }
    ~Cublas() 
    {
        if(handle) cublasDestroy(handle);
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
    void curand_uniform(float * dev, size_t size) {                
        curandGenerateUniform(gen,dev,size*sizeof(float));                
    }
    void curand_normal(float * dev, size_t size, float mean, float stddev) {                
        curandGenerateNormal(gen,dev,size*sizeof(float), mean, stddev);                
    }
    void curand_lognormal(float *dev, size_t size, float mean, float stddev) {                
        curandGenerateLogNormal(gen,dev,size*sizeof(float), mean, stddev);                
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
    void fill(float _Complex v) {
        assert(ptr != NULL);
        for(size_t i = 0; i < length/sizeof(float _Complex); i++)
            ((float _Complex*)ptr)[i] = v;
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
    //float _Complex& cf32(size_t pos) { return ((float _Complex*)ptr)[pos]; }
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

struct DevPtr 
{
    void * dev;
    void * host;
    int    len;
    int    elems;
    DevPtr(void * d, size_t s, size_t e = sizeof(float)) { 
        dev = d; 
        len = s;
        elems = e;         
        host = calloc(s,e);
    }
    ~DevPtr() {                 
        return_memory(len,(float*)dev);                 
        free(host);
    }
    size_t size() const { return len*elems; }
    void memcpy_host() { cudaMemcpy(host,dev,size(),cudaMemcpyDeviceToHost); }
    void memcpy_device() { cudaMemcpy(dev,host,size(),cudaMemcpyHostToDevice); }
    void memcpy_device(void *ptr) { cudaMemcpy(dev,ptr,size(),cudaMemcpyHostToDevice); }
    

    int8_t& int8(size_t i) { return ((int8_t*)host)[i]; }
    uint8_t& uint8(size_t i) { return ((uint8_t*)host)[i]; }
    int16_t& int16(size_t i) { return ((int16_t*)host)[i]; }
    uint16_t& uint16(size_t i) { return ((uint16_t*)host)[i]; }
    int32_t& int32(size_t i) { return ((int32_t*)host)[i]; }
    uint32_t& uint32(size_t i) { return ((uint32_t*)host)[i]; }
    int64_t& int64(size_t i) { return ((int64_t*)host)[i]; }
    uint64_t& uint64(size_t i) { return ((uint64_t*)host)[i]; }
    float& f32(size_t i) { return ((float*)host)[i]; }
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

inline void rand(float * devPtr, float *host, int size, float min, float max) {
    typedef std::chrono::high_resolution_clock myclock;
    myclock::time_point beginning = myclock::now();
    myclock::duration d = myclock::now() - beginning;
    unsigned seed = d.count();
    std::default_random_engine generator(seed);
    std::uniform_real_distribution<double> distribution(min,max);                
    for(size_t i = 0; i < size; i++) host[i] = distribution(generator);
    cudaMemcpy(devPtr,host,size*sizeof(float),cudaMemcpyHostToDevice);            
}



struct Vector 
{
public:

    void init(size_t n) {
        dev = find_memory(n);
        if(dev == NULL) {
            cudaMalloc((void**)&dev,n*sizeof(float));
            assert(dev != NULL);
        }
        N    = n;
        host = (float*)calloc(n,sizeof(float));
        host_ptr = { host, free };
        dev_ptr  = { new DevPtr(dev,N), [](DevPtr* p){ delete p; } };
        //zero();
    }

    float * dev;
    float * host;
    std::shared_ptr<float> host_ptr;
    std::shared_ptr<DevPtr> dev_ptr;
    size_t   N;
    
    
public:

    Vector() {
        dev = host = nullptr;
        N = 0;    
    }
    Vector(size_t n) {
        init(n);        
    }
    Vector(int n, const std::initializer_list<float> & input) {
        init(n);
        std::vector<float> tmp(input.begin(),input.end());
        cudaMemcpy(dev,tmp.data(),n*sizeof(float),cudaMemcpyHostToDevice);
    }
    Vector(int n, const std::vector<float> & tmp) {
        init(n);        
        cudaMemcpy(dev,tmp.data(),n*sizeof(float),cudaMemcpyHostToDevice);
    }    
    Vector(const Vector & v) {
        //*this = v;
        init(v.N);
        copy(v);
    }        
    Vector(float * p, size_t n) {
        N = n;
        host = (float*)calloc(n,sizeof(float));
        host_ptr = { host, free };
        dev = p;
        dev_ptr  = { new DevPtr(dev,N), [](DevPtr* p){ delete p; } };
    }
    ~Vector() = default;

    Vector& operator = (const Vector & v) {
        if(dev == v.dev) return *this;        
        host_ptr.reset();
        host_ptr = v.host_ptr;
        dev_ptr.reset();
        dev_ptr  = v.dev_ptr;
        N = v.N;
        host = v.host;
        dev  = v.dev;
        return *this;
    }

    size_t size() const { return N; }    
    void download_host()     {
        cudaMemcpy(host,dev,size()*sizeof(float),cudaMemcpyDeviceToHost);
    }
    void upload_device()     {
        cudaMemcpy(dev,host,size()*sizeof(float),cudaMemcpyHostToDevice);
    }    
    void zero()    {
        cudaMemsetAsync(dev,0x00,size()*sizeof(float));
    }
    void ones()    {
        fill(1.0f);
    }    
    void randu() {
        rand(dev,host,size(),0.0f,1.0f);
    }
    void random(float min, float max) {
        rand(dev,host,size(),min,max);
    }    
    void fill(const float val)     {
        for(size_t i = 0; i < size(); i++) host[i] = val;
        cudaMemcpy(dev,host,size()*sizeof(float),cudaMemcpyHostToDevice);
    }

    float& operator[](array_index pos)    {        
        assert(pos < size());
        while(pos < 0) pos += size();
        return host[pos];        
    }
    float __getitem(array_index pos) { return (*this)[pos]; }
    void __setitem(array_index pos, float val) { while(pos < 0) pos += N; host[pos] = val; }

    Vector operator - ();
    Vector operator + (const Vector & a);
    Vector operator - (const Vector & a);
    Vector operator * (const Vector & a);

    
    float dot(const Vector & a) {
        float r;
        cublasStatus_t     status; 
        status = cublasSdot(cublas->handle,N,dev,1,a.dev,1,&r);
        ASSERTOK(status);
        return r;
    }
    void resize(size_t n) {
        N = n;
        dev_ptr.reset();
        host_ptr.reset();
        init(n);
    }
    void copy(const Vector & a) {
        if(a.size() != size()) resize(a.size());
        cublasStatus_t     status; 
        status = cublasScopy(cublas->handle,N,a.dev,1,dev,1);
        ASSERTOK(status);        
    }
    int max_index() {
       int r=0;
       cublasStatus_t     status; 
       status = cublasIsamax(cublas->handle,N,dev,1,&r);
       ASSERTOK(status);
       return r-1;
    }
    int min_index() {
       int r=0;
       cublasStatus_t     status; 
        status = cublasIsamin(cublas->handle,N,dev,1,&r);
        ASSERTOK(status);
       return r-1;
    }
    float sum() {
        float r = 0;
        cublasStatus_t     status; 
        status = cublasSasum(cublas->handle,N,dev,1,&r);
        ASSERTOK(status);
        return r;
    }
    float nrm2() {
        float r = 0;
        cublasStatus_t     status; 
        status = cublasSnrm2(cublas->handle,N,dev,1,&r);
        ASSERTOK(status);
        return r;
    }
    Vector scale(float x) {
        Vector r(N);
        r.copy(*this);
        cublasStatus_t     status; 
        status = cublasSscal(cublas->handle,N,&x,r.dev,1);
        ASSERTOK(status);
        return r;
    }
    void swap(Vector & v) {
        cublasStatus_t     status; 
        status = cublasSswap(cublas->handle,N,dev,1,v.dev,1);
        ASSERTOK(status);
    }

    void print()  {
        download_host();               
        std::cout << "vector[" << N << "]" << std::endl;
        for(size_t w = 0; w < size(); w++) {
            std::cout << host[w] << ",";
        }
        std::cout << std::endl;        
    }

    Vector clone(const Vector & v) { return Vector(v); }
    Vector eval() { return Vector(*this); }

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
    void pow(float p) {
        vector_r_powf_const(dev,p,dev,N);
    }
    void pow(const Vector & p) {
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
    void atan2(float v) {
        vector_r_atan2f_const(dev,v,dev,N);        
    }
    void atan2(const Vector & v) {
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
    void copysign(const float val) {
        vector_r_copysignf_const(dev,val,dev,N);        
    }
    void copysign(const Vector& val) {
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
    void fdim(const Vector & b) {
        vector_r_fdimf(dev,b.dev,dev,N);
    }
    void fmod(const Vector & b) {
        vector_r_fmodf(dev,b.dev,dev,N);
    }
    void hypot(const Vector & b) {
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
    void norm3d(const Vector & a, const Vector & b, const Vector & c) {
        vector_r_norm3df(a.dev,b.dev,c.dev,dev,N);
    }
    void norm4d(const Vector & a, const Vector & b, const Vector & c, const Vector & d) {
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
    void remainder(const Vector & b) {
        vector_r_remainderf(dev,b.dev,dev,N);
    }
    void rhypot(const Vector & b) {
        vector_r_rhypotf(dev,b.dev,dev,N);
    }
    /*
    void rint() {
        vector_r_rintf(dev,dev,N);        
    }
    */
    void rnorm3d(const Vector & a, const Vector & b, const Vector & c) {
        vector_r_rnorm3df(a.dev,b.dev,c.dev,dev,N);
    }
    void rnorm4d(const Vector & a, const Vector & b, const Vector & c, const Vector & d) {
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


struct VectorView {
    float * host;
    int     r,c;

    VectorView(float *h, int row, int cols) {        
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

    float& operator[](array_index pos) {
        while(pos < 0) pos += c;
        return host[r*c+pos];
    }
    float __getitem(array_index pos) {
        while(pos < 0) pos += c;
        return host[r*c+pos];
    }
    void  __setitem(array_index pos, float val) {
        while(pos < 0) pos += c;
        host[r*c+pos] = val;                
    }
};


struct Matrix 
{
public:

    void init(size_t m, size_t n) {        
        dev = find_memory(m*n);
        if(dev == NULL) {
            cudaMalloc((void**)&dev,m*n*sizeof(float));
            assert(dev != NULL);
        }        
        M    = m;
        N    = n;        
        host = (float*)calloc(m*n,sizeof(float));
        assert(host != NULL);
        host_ptr = { host, free };
        dev_ptr  = {new DevPtr(dev,M*N), [](DevPtr* p){delete p;} };
        //zero();
    }

    float * dev;
    float * host;
    std::shared_ptr<float> host_ptr;
    std::shared_ptr<DevPtr> dev_ptr;
    size_t   M,N;

public:

    Matrix() {
        dev = host = nullptr;
        M = N = 0;
    }
    Matrix(size_t i, size_t j) {
        init(i,j);           
    }
    Matrix(float * p, size_t m, size_t n) {
        M = m;
        N = n;
        host = (float*)calloc(m*n,sizeof(float));
        host_ptr = { host, free };
        dev = p;
        dev_ptr  = {new DevPtr(dev,M*N), [](DevPtr* p){delete p;} };
    }
    Matrix(size_t i, size_t j, const std::initializer_list<float> & in) {
        init(i,j);
        std::vector<float> v(in.begin(),in.end());
        cudaMemcpy(dev,v.data(), size()*sizeof(float),cudaMemcpyHostToDevice);
    }
    Matrix(size_t i, size_t j, const std::vector<float> & v) {
        init(i,j);
        cudaMemcpy(dev,v.data(), size()*sizeof(float),cudaMemcpyHostToDevice);
    }    
    Matrix(const Matrix & m) {        
        //init(m.M,m.N);
        //copy(m);        
        *this = m;
    }    
    ~Matrix() = default;
        

    Matrix& operator = (const Matrix & m) {              
        M = m.M;
        N = m.N;
        host = m.host;
        dev  = m.dev;
        host_ptr.reset();
        host_ptr = m.host_ptr;
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
        host_ptr.reset();
        init(i,j);
    }
    
    Matrix t() {                    
        Matrix r(N,M);                
        int m = M;
        int n = N;    
        float alpha = 1.0;
        float beta = 1.0;
        cublasStatus_t status = cublasSgeam(cublas->handle,CUBLAS_OP_T,CUBLAS_OP_N,m,n,&alpha,dev,n,&beta,dev,m,r.dev,m);    
        ASSERTOK(status);
        return r;            
    }

    Matrix transpose(Matrix & a){
        return a.t();
    }

    void copy(const Matrix & a) {        
        if(a.M != M || a.N != N) resize(a.M,a.N);
        cublasStatus_t     status; 
        status = cublasScopy(cublas->handle,M*N,a.dev,1,dev,1);
        ASSERTOK(status);            
    }
    
    Matrix operator + (const Matrix & b);
    Matrix operator - (const Matrix & b);    
    Matrix operator * (const Matrix & m);
    
    Matrix operator -();

    
    Matrix operator + (float v) {
        Matrix r(M,N);
        //matrix_r_addf_const(dev,v,r.dev,M,N);
        vector_r_addf_const(dev,v,r.dev,M*N);
        return r;
    }
    Matrix operator - (float v) {
        Matrix r(M,N);
        //matrix_r_subf_const(dev,v,r.dev,M,N);
        vector_r_subf_const(dev,v,r.dev,M*N);
        return r;
    }    
    Matrix operator * (float v) {
        Matrix r(M,N); 
        // dont know why it needs this
        //r.zero();       
        //matrix_r_mulf_const(dev,v,r.dev,M,N);          
        vector_r_mulf_const(dev,v,r.dev,M*N);        
        return r;
    }
    Matrix operator / (float v) {
        Matrix r(M,N);
        //matrix_r_divf_const(dev,v,r.dev,M,N);
        vector_r_divf_const(dev,v,r.dev,M*N);
        return r;
    }
    Matrix operator % (float v) {
        Matrix r(M,N);
        //matrix_r_modf_const(dev,v,r.dev,M,N);
        vector_r_modf_const(dev,v,r.dev,M*N);
        return r;
    } 
    
    
    void hadamard(const Matrix & b) {
        assert(M == b.M && N == b.N);
        vector_r_mulf(dev,b.dev,dev,M*N);                
    }
    void dgmm(const Vector & diagonal) {
        cublasStatus_t status = cublasSdgmm(cublas->handle,CUBLAS_SIDE_RIGHT,N,M,dev,N,diagonal.dev,1,dev,N);
        ASSERTOK(status);
    }
    void download_host()     {
        cudaMemcpy(host,dev,size()*sizeof(float),cudaMemcpyDeviceToHost);
    }
    void upload_device()     {
        cudaMemcpy(dev,host,size()*sizeof(float),cudaMemcpyHostToDevice);
    }    
    void zero()    {
        cudaMemsetAsync(dev,0x00,size()*sizeof(float));
    }
    void ones()    {
        fill(1.0f);
    }    
    void randu() {
        rand(dev,host,size(),0.0f,1.0f);
    }
    void random(float min, float max) {
        rand(dev,host,size(),min,max);
    }    
    void fill(const float val)     {
        for(size_t i = 0; i < size(); i++) host[i] = val;
        cudaMemcpy(dev,host,size()*sizeof(float),cudaMemcpyHostToDevice);
    }

    float& operator()(int i, int j) {
        while(i < 0) i+= M;
        while(j < 0) j+= N;
        return host[i*N+j];
    }
    float& operator[](array_index pos)    {        
        assert(pos < size());
        while(pos < 0) pos += size();
        return host[pos];        
    }
    
    VectorView __getitem(array_index row) {        
        return VectorView(host,row,N);
    }
    void __setitem(array_index pos, float val) { while(pos < 0) pos += N; host[pos] = val; }

    size_t index(int r, int c) { 
        while(r < 0) r += M;
        while(c < 0) c += N;
        return r*N + c;
    }
    void abs() {
        //matrix_r_fabsf(dev,dev,M,N);
        vector_r_fabsf(dev,dev,M*N);
    }
    void exp() {
        //matrix_r_expf(dev,dev,M,N);
        vector_r_expf(dev,dev,M*N);
    }    
    void log2() {
        //matrix_r_log2f(dev,dev,M,N);
        vector_r_log2f(dev,dev,M*N);
    }
    void log10() {
        //matrix_r_log10f(dev,dev,M,N);
        vector_r_log10f(dev,dev,M*N);
    }
    void pow(float p) {
        //matrix_r_powf_const(dev,p,dev,M,N);
        vector_r_powf_const(dev,p,dev,M*N);
    }
    void pow(const Matrix & p) {
        //matrix_r_powf(dev,p.dev,dev,M,N);
        vector_r_powf(dev,p.dev,dev,M*N);
    }
    void sqrt() {
        //matrix_r_sqrtf(dev,dev,M,N);
        vector_r_sqrtf(dev,dev,M*N);
    }
    void rsqrt() {
        //matrix_r_rsqrtf(dev,dev,M,N);
        vector_r_rsqrtf(dev,dev,M*N);
    }
    void sin() {
        //matrix_r_sinf(dev,dev,M,N);
        vector_r_sinf(dev,dev,M*N);
    }
    void cos() {
        //matrix_r_cosf(dev,dev,M,N);
        vector_r_cosf(dev,dev,M*N);
    }
    void tan() {
        //matrix_r_tanf(dev,dev,M,N);
        vector_r_tanf(dev,dev,M*N);
    }
    void asin() {
        //matrix_r_asinf(dev,dev,M,N);
        vector_r_asinf(dev,dev,M*N);
    }
    void acos() {
        //matrix_r_acosf(dev,dev,M,N);
        vector_r_acosf(dev,dev,M*N);
    }
    void atan() {
        //matrix_r_atanf(dev,dev,M,N);
        vector_r_atanf(dev,dev,M*N);
    }
    void sinh() {
        //matrix_r_sinhf(dev,dev,M,N);
        vector_r_sinhf(dev,dev,M*N);
    }
    void cosh() {
        //matrix_r_coshf(dev,dev,M,N);
        vector_r_coshf(dev,dev,M*N);
    }
    void tanh() {
        //matrix_r_tanhf(dev,dev,M,N);
        vector_r_tanhf(dev,dev,M*N);
    }
    void asinh() {
        //matrix_r_asinhf(dev,dev,M,N);
        vector_r_asinhf(dev,dev,M*N);
    }
    void acosh() {
        //matrix_r_acoshf(dev,dev,M,N);
        vector_r_acoshf(dev,dev,M*N);
    }
    void atanh() {
        //matrix_r_atanhf(dev,dev,M,N);
        vector_r_atanhf(dev,dev,M*N);
    }    
    void atan2(float v) {
        //matrix_r_atan2f_const(dev,v,dev,M,N);
        vector_r_atan2f_const(dev,v,dev,M*N);
    }
    void atan2(const Matrix & v) {
        //matrix_r_atan2f(dev,v.dev,dev,M,N);
        vector_r_atan2f(dev,v.dev,dev,M*N);
    }
    void sigmoid() {
        //matrix_r_sigmoidf(dev,dev,M,N);
        vector_r_sigmoidf(dev,dev,M*N);
    }
    void sigmoid_deriv() {
        //matrix_r_sigmoid_gradf(dev,dev,M,N);
        vector_r_sigmoid_gradf(dev,dev,M*N);
    }
    void relu() {
        //matrix_r_reluf(dev,dev,M,N);
        vector_r_reluf(dev,dev,M*N);        
    }
    void relu_deriv() {
        //matrix_r_relu_gradf(dev,dev,M,N);
        vector_r_relu_gradf(dev,dev,M*N);        
    }
    void softmax() {
        //matrix_r_softmaxf(dev,dev,M,N);
        vector_r_softmaxf(dev,dev,M*N);
    }
    void tanh_deriv() {
        //matrix_r_tanh_gradf(dev,dev,M,N);
        vector_r_tanh_gradf(dev,dev,M*N);
    }
    float sum() {
        float r = 0;
        cublasStatus_t     status; 
        status = cublasSasum(cublas->handle,M*N,dev,1,&r);
        ASSERTOK(status);
        return r;
    }
    Matrix eval();
    
    void addToEachRow(const Matrix &b, int row=0) {
        assert(N == b.N);
        float alpha=1.0;        
        for(size_t i = 0; i < M; i++)
            //vector_addf_row(dev,(i*N),b.dev,row,N);        
            cublasSaxpy(cublas->handle,N,&alpha,dev+i*N,1,b.dev+row*N,1);
    }

    
    float get(int row, int col) {
        return host[row*N + col];
    }
    void set(int row, int col, float val) {
        while(row < 0) row += M;
        while(col < 0) col += N;
        host[row*N + col] = val;

    }

    Matrix get_row(int row) {
        while(row < 0) row += M;     
        assert(row < M);   
        Matrix r(1,cols());        
        cudaMemcpy(r.dev,dev+row*N,cols()*sizeof(float),cudaMemcpyDeviceToDevice);
        return r;
    }
    void set_row(const Matrix & m, int dst_row=0, int src_row=0) {
        while(dst_row < 0) dst_row += M;
        while(dst_row < 0) dst_row += M;
        while(src_row < 0) src_row += M;        
        cudaMemcpy(dev + dst_row*N, m.dev + src_row*m.N, m.N*sizeof(float), cudaMemcpyDeviceToDevice);
    }
    void set_row(const Vector & v, int dst_row=0) {
        assert(N == v.N);       
        while(dst_row < 0) dst_row += M; 
        cudaMemcpy(dev + dst_row*N, v.dev, v.N*sizeof(float), cudaMemcpyDeviceToDevice);
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
        std::cout << "Matrix(" << M << "," << N << ")" << std::endl;
    }

    void identity()  {     
        // identity only makes sense on square matrix.   
        assert(M == N);
        size_t c = 0;
        download_host();
        zero();
        for(size_t i = 0; i < M; i++) {
            host[i*N + c++] = 1;
        }            
        upload_device();
    }    


    void cbrt() {
        vector_r_cbrtf(dev,dev,M*N);        
    }
    void ceil() {
        vector_r_ceilf(dev,dev,M*N);        
    }
    /*
    void copysign(const float val) {
        vector_r_cbrtf_const(dev,val,dev,M*N);        
    }
    void copysign(const Vector& val) {
        vector_r_cbrtf_const(dev,val.dev,dev,M*N);        
    }
    */
    void cospi() {
        vector_r_cospif(dev,dev,M*N);        
    }
    void cyl_bessel_i0() {
        vector_r_cyl_bessel_i0f(dev,dev,M*N);
    }
    void cyl_bessel_i1() {
        vector_r_cyl_bessel_i1f(dev,dev,M*N);
    }
    void erfc() {
        vector_r_erfcf(dev,dev,M*N);
    }
    void erfcx() {
        vector_r_erfcxf(dev,dev,M*N);
    }
    void erfcinv() {
        vector_r_erfcinvf(dev,dev,M*N);
    }
    void erf() {
        vector_r_erff(dev,dev,M*N);
    }
    void erfinv() {
        vector_r_erfinvf(dev,dev,M*N);
    }
    void exp10() {
        vector_r_exp10f(dev,dev,M*N);
    }
    void exp2() {
        vector_r_exp2f(dev,dev,M*N);
    }
    void expm1() {
        vector_r_expm1f(dev,dev,M*N);
    }
    void fabs() {
        vector_r_fabsf(dev,dev,M*N);
    }
    void fdim(const Matrix & b) {
        vector_r_fdimf(dev,b.dev,dev,M*N);
    }
    void fmod(const Matrix & b) {
        vector_r_fmodf(dev,b.dev,dev,M*N);
    }
    void hypot(const Matrix & b) {
        vector_r_hypotf(dev,b.dev,dev,M*N);
    }
    void ilogb() {
        vector_r_ilogbf(dev,dev,M*N);
    }
    void j0() {
        vector_r_j0f(dev,dev,M*N);
    }
    void j1() {
        vector_r_j1f(dev,dev,M*N);
    }
    void jn(int n) {
        vector_r_jnf(dev,dev,M*N,n);
    }
    void lgamma() {
        vector_r_lgammaf(dev,dev,M*N);
    }
    void log1p() {
        vector_r_log1pf(dev,dev,M*N);
    }
    void logb() {
        vector_r_logbf(dev,dev,M*N);
    }
    void norm3d(const Matrix & a, const Matrix & b, const Matrix & c) {
        vector_r_norm3df(a.dev,b.dev,c.dev,dev,M*N);
    }
    void norm4d(const Matrix & a, const Matrix & b, const Matrix & c, const Matrix & d) {
        vector_r_norm4df(a.dev,b.dev,c.dev,d.dev,dev,M*N);
    }
    void normcdf() {
        vector_r_normcdff(dev,dev,M*N);
    }
    void normcdfinv() {
        vector_r_normcdfinvf(dev,dev,M*N);
    }
    void norm() {
        vector_r_normf(1,dev,dev,M*N);
    }
    void rcbrt() {
        vector_r_rcbrtf(dev,dev,M*N);        
    }
    void remainder(const Matrix& b) {
            vector_r_remainderf(dev,b.dev,dev,M*N);
    }
    void rhypot(const Matrix & b) {
        vector_r_rhypotf(dev,b.dev,dev,M*N);
    }
    /*
    void rint() {
        vector_r_rintf(dev,dev,M*N);        
    }
    */
    void rnorm3d(const Matrix & a, const Matrix & b, const Matrix & c) {
        vector_r_rnorm3df(a.dev,b.dev,c.dev,dev,M*N);
    }
    void rnorm4d(const Matrix & a, const Matrix & b, const Matrix & c, const Matrix & d) {
        vector_r_rnorm4df(a.dev,b.dev,c.dev,d.dev,dev,M*N);
    }
    /*
    void round() {
        vector_r_round(dev,dev,M*N);        
    }
    */
    void rnorm() {
        vector_r_rnormf(1,dev,dev,M*N);
    }
    void tgamma() {
        vector_r_tgammaf(dev,dev,M*N);
    }
    void trunc() {
        vector_r_truncf(dev,dev,M*N);
    }
    void y0() {
        vector_r_y0f(dev,dev,M*N);
    }
    void y1() {
        vector_r_y1f(dev,dev,M*N);
    }
    void yn(int n) {
        vector_r_ynf(n,dev,dev,M*N);
    }

};

inline Matrix CopyMatrix(const Matrix & m) {
    Matrix r(m.M,m.N);
    r.copy(m);
    return r;
}
inline Vector CopyVector(const Vector & m) {
    Vector r(m.N);
    r.copy(m);
    return r;
}

inline Vector Vector::operator - () {
    Vector r(N);
    r.copy(*this);
    vector_r_mulf_const(r.dev,-1.0,r.dev,N);
    return r;
}
inline Vector Vector::operator + (const Vector & a) {
    assert(size() == a.size());
    Vector r = CopyVector(a);
    cublasStatus_t     status; 
    float alpha = 1.0;
    status = cublasSaxpy(cublas->handle, N,&alpha,dev,1,r.dev,1);
    ASSERTOK(status);
    return r;
}
inline Vector Vector::operator - (const Vector & a) {
    assert(size() == a.size());
    Vector r = CopyVector(a);
    cublasStatus_t     status; 
    float alpha = -1.0;
    status = cublasSaxpy(cublas->handle,N,&alpha,a.dev,1,r.dev,1);
    ASSERTOK(status);
    return r;
}
inline Vector Vector::operator * (const Vector & a) {
    assert(size() == a.size());
    Vector r = CopyVector(a);
    cublasStatus_t     status;         
    status = cublasSdgmm(cublas->handle,CUBLAS_SIDE_RIGHT,1,N,dev,1,a.dev,1,r.dev,1);
    ASSERTOK(status);
    return r;
}
inline Matrix Matrix::operator -() {
    Matrix r = CopyMatrix(*this);
    //matrix_r_mulf_const(r.dev,-1.0,r.dev,r.M,r.N);
    vector_r_mulf_const(r.dev,-1.0,r.dev,r.M*r.N);
    return r;
}

inline Matrix Matrix::eval() {
    Matrix r = CopyMatrix(*this);
    return r;
}

inline Vector matvec(const Matrix & a, const Vector & v,bool transa, float alpha=1.0, float beta=0.0) {    
    Vector R(v.size());
    int m = a.M;
    int n = a.N;
    cublasOperation_t ta = CUBLAS_OP_N;
    cublasStatus_t status = cublasSgemv(cublas->handle,ta,n,m,&alpha,a.dev,n,a.dev,1,&beta,R.dev,1);
    ASSERTOK(status);
    return R;
}

// C^T = B^T * A^T
// assumes a and b are in row major order
inline Matrix matmul(const Matrix & a, const Matrix & b, float alpha=1.0, float beta=0.0) {    
    assert(a.N == b.M);     
    Matrix C(a.M,b.N);
    cublasOperation_t ta = CUBLAS_OP_N;
    cublasOperation_t tb = CUBLAS_OP_N;      
    int m = a.M;
    int k = a.N;
    int n = b.N;    
    cublasStatus_t status = cublasSgemm(cublas->handle,ta,tb,n,m,k,&alpha,b.dev,n,a.dev,k,&beta,C.dev,n);                        
    ASSERTOK(status);
    return C;        
}

// C^T = A*B
// A and B are in column major order
inline Matrix matmulTT(const Matrix & a,const Matrix & b, float alpha=1.0, float beta=0.0) {         
    assert(a.N == b.M);    
    cublasOperation_t ta = CUBLAS_OP_T;
    cublasOperation_t tb = CUBLAS_OP_T;      
    Matrix C(a.N,b.M);            
    int m = a.M;
    int k = a.N;
    int n = b.N; 
    cublasStatus_t status = cublasSgemm(cublas->handle,ta,tb,n,m,k,&alpha,b.dev,n,a.dev,k,&beta,C.dev,n);                                            
    ASSERTOK(status);
    return C;    
}

// C^T = B^T*A
// A is in column major order
inline Matrix matmulTN(const Matrix & a,const Matrix & b, float alpha=1.0, float beta=0.0) {         
    assert(a.N == b.N);    
    cublasOperation_t ta = CUBLAS_OP_N;
    cublasOperation_t tb = CUBLAS_OP_T;               
    Matrix C(a.N,b.N);        
    int m = a.M;
    int k = a.N;
    int n = b.N; 
    cublasStatus_t status = cublasSgemm(cublas->handle,ta,tb,n,k,m,&alpha,b.dev,n,a.dev,k,&beta,C.dev,n); 
    ASSERTOK(status);
    return C;        
}
// C^T = B*A^T
// A is in column major order
inline Matrix matmulNT(const Matrix & a,const Matrix & b, float alpha=1.0, float beta=0.0) {         
    assert(a.M == b.M);    
    cublasOperation_t ta = CUBLAS_OP_T;
    cublasOperation_t tb = CUBLAS_OP_N;            
    Matrix C(a.M,b.M);        
    int m = a.M;
    int k = a.N;
    int n = b.N; 
    cublasStatus_t status = cublasSgemm(cublas->handle,ta,tb,m,n,k,&alpha,b.dev,n,a.dev,k,&beta,C.dev,n);                        
    ASSERTOK(status);
    return C;    
}


inline void add(const Vector & a, const Vector & b, Vector & r)
{
    assert(b.size() == a.size());        
    cublasStatus_t     status; 
    float alpha = 1.0;
    r.resize(a.N);
    r.copy(b);
    status = cublasSaxpy(cublas->handle, a.N,&alpha,a.dev,1,r.dev,1);
    ASSERTOK(status);
}
inline void sub(const Vector & a, const Vector & b, Vector & r)
{
    assert(b.size() == a.size());        
    cublasStatus_t     status; 
    float alpha = -1.0;
    r.resize(a.N);
    r.copy(b);
    status = cublasSaxpy(cublas->handle,a.N,&alpha,a.dev,1,r.dev,1);
    ASSERTOK(status);
}
inline void mul(const Vector & a, const Vector & b, Vector & r)
{
    assert(b.size() == a.size());
    cublasStatus_t     status;         
    r.resize(a.N);    
    status = cublasSdgmm(cublas->handle,CUBLAS_SIDE_RIGHT,1,a.N,b.dev,1,a.dev,1,r.dev,1);
    ASSERTOK(status);
}
inline Vector operator + (const Vector & a, const Vector & b) {
    Vector r(a.N);    
    add(a,b,r);
    return r;
}
inline Vector operator - (const Vector & a, const Vector & b) {
    Vector r(a.N);    
    sub(a,b,r);
    return r;
}
inline Vector operator * (const Vector & a, const Vector & b) {
    Vector r(a.N);    
    mul(a,b,r);
    return r;
}

inline Vector operator + (const Vector & a, const float v)
{
    Vector r(a.N);    
    vector_r_addf_const(a.dev,v,r.dev,a.N);
    return r;
}
inline Vector operator + (const float v, const Vector & a)
{    
    Vector r(a.N);
    vector_r_addf_const(a.dev,v,r.dev,a.N);
    return r;
}
inline Vector operator - (const Vector & a, const float v)
{    
    Vector r(a.N);    
    vector_r_subf_const(a.dev,v,r.dev,a.N);
    return r;
}
inline Vector operator - (const float v, const Vector & a)
{    
    Vector r(a.N);    
    r = -r;
    vector_r_addf_const(r.dev,v,r.dev,a.N);
    return r;
}
inline Vector operator * (const Vector & a, const float v)
{    
    Vector r(a.N);    
    //vector_r_mulf_const(a.dev,v,r.dev,a.N);
    cublasSscal(cublas->handle,a.N,&v,r.dev,1);
    return r;
}
inline Vector operator * (const float v, const Vector & a)
{    
    Vector r(a.N);    
    //vector_r_mulf_const(a.dev,v,r.dev,a.N);
    cublasSscal(cublas->handle,a.N,&v,r.dev,1);
    return r;
}
/*
inline Vector operator / (const Vector & a, const float v)
{    
    Vector r(a.N);    
    //vector_r_mulf_const(a.dev,v,r.dev,a.N);
    assert(v != 0);
    float q = 1.0/v;
    cublasSscal(cublas->handle,a.N,&q,r.dev,1);
    return r;
}
inline Vector operator % (const Vector & a, const float v)
{    
    Vector r(a.N);        
    vector_r_modf_const(r.dev,v,r.dev,a.N);
    return r;
}
*/
inline Matrix hadamard(const Matrix & a, const Matrix &b)
{
    Matrix r = CopyMatrix(a);    
    r.hadamard(b);
    return r;
}

inline void add(const Matrix & a, const Matrix & b, Matrix & r)
{
    assert(a.M == b.M && a.N == b.N);                
    int m = a.M;
    int k = a.N;    
    int n = a.N;
    float alpha = 1.0;
    float beta  = 1.0;
    r.resize(a.M,b.N);            
    cublasStatus_t status = cublasSgeam(cublas->handle,CUBLAS_OP_N,CUBLAS_OP_N,n,m,&alpha,b.dev,n,&beta,a.dev,n,r.dev,n);
    ASSERTOK(status);        
}
inline void sub(const Matrix & a, const Matrix & b, Matrix & r)
{
    assert(a.M == b.M && a.N == b.N);                
    int m = a.M;
    int k = a.N;    
    int n = a.N;
    float alpha = -1.0;
    float beta  = 1.0;        
    r.resize(a.M,b.N);    
    cublasStatus_t status = cublasSgeam(cublas->handle,CUBLAS_OP_N,CUBLAS_OP_N,n,m,&alpha,b.dev,n,&beta,a.dev,n,r.dev,n);
    ASSERTOK(status);        
}
inline void mul(const Matrix & a, const Matrix & b, Matrix & r)
{
    assert(a.N == b.M);             
    cublasOperation_t ta = CUBLAS_OP_N;
    cublasOperation_t tb = CUBLAS_OP_N;          
    int m = a.M;
    int k = a.N;    
    int n = b.N;
    float alpha=1.0;
    float beta=0.0;
    r.resize(a.M,b.N);    
    cublasStatus_t status = cublasSgemm(cublas->handle,ta,tb,n,m,k,&alpha,b.dev,n,a.dev,k,&beta,r.dev,n);                        
    ASSERTOK(status);    
}

inline Matrix operator + (const Matrix & a, const Matrix & b)
{
    assert(a.M == b.M && a.N == b.N);
    Matrix r(a.M,a.N);    
    add(a,b,r);
    return r;
}
inline Matrix operator - (const Matrix & a, const Matrix & b)
{
    assert(a.M == b.M && a.N == b.N);
    Matrix r(a.M,a.N);    
    sub(a,b,r);
    return r;
}
inline Matrix operator * (const Matrix & a, const Matrix & b)
{        
    Matrix r(a.M,a.N);    
    mul(a,b,r);
    return r;
}
inline Matrix operator + (const Matrix & a, const float v)
{    
    Matrix r(a.M,a.N);        
    //matrix_r_addf_const(a.dev,v,r.dev,a.M,a.N);
    vector_r_addf_const(a.dev,v,r.dev,a.M*a.N);
    return r;
}
inline Matrix operator + (const float v, const Matrix & a)
{
    Matrix r(a.M,a.N);    
    //matrix_r_addf_const(a.dev,v,r.dev,a.M,a.N);
    vector_r_addf_const(a.dev,v,r.dev,a.M*a.N);
    return r;
}
inline Matrix operator - (const Matrix & a, const float v)
{
    Matrix r(a.M,a.N);    
    vector_r_subf_const(a.dev,v,r.dev,a.M*a.N);
    return r;
}
inline Matrix operator - (const float v, const Matrix & a)
{    
    Matrix r = CopyMatrix(a);
    r = -r;
    vector_r_addf_const(r.dev,v,r.dev,a.M*a.N);
    return r;
}
inline Matrix operator * (const Matrix & a, const float v)
{
    Matrix r = CopyMatrix(a);    
    cublasSscal(cublas->handle,a.M*a.N,&v,r.dev,1);
    //vector_r_mulf_const(a.dev,v,r.dev,a.M*a.N);
    return r;
}
inline Matrix operator * (const float v, const Matrix & a)
{    
    Matrix r = CopyMatrix(a);
    cublasSscal(cublas->handle,a.M*a.N,&v,r.dev,1);
    //vector_r_mulf_const(a.dev,v,r.dev,a.M*a.N);
    return r;
}
inline Matrix operator / (const Matrix & a, float v)
{    
    Matrix r = CopyMatrix(a);
    assert(v != 0);
    v= 1.0/v;
    cublasSscal(cublas->handle,a.M*a.N,&v,r.dev,1);
    //vector_r_divf_const(a.dev,v,r.dev,a.M*a.N);
    return r;
}
inline Matrix operator % (const Matrix & a, const float v)
{    
    Matrix r(a.M,a.N);        
    //matrix_r_modf_const(a.dev,v,r.dev,a.M,a.N);
    vector_r_modf_const(a.dev,v,r.dev,a.M*a.N);
    return r;
}


Matrix Matrix::operator + (const Matrix & b) {        
    Matrix r(M,N);                
    add(*this,b,r);
    return r;
}

Matrix Matrix::operator - (const Matrix & b) {                 
    Matrix r(M,N);                
    sub(*this,b,r);
    return r;
}
Matrix Matrix::operator * (const Matrix & m) 
{   
    Matrix r(M,N);      
    mul(*this,m,r);
    return r;
}


inline Matrix abs(const Matrix & m) {    
    Matrix r = CopyMatrix(m);
    r.abs();
    return r;
}
inline Matrix exp(const Matrix & m) {
    Matrix r = CopyMatrix(m);
    r.exp(); 
    return r;   
}    
inline Matrix log2(const Matrix & m) {
    Matrix r = CopyMatrix(m);
    r.log2();
    return r;
}
inline Matrix log10(const Matrix & m) {
    Matrix r = CopyMatrix(m);
    r.log10();
    return r;
}
inline Matrix pow(const Matrix & m,float p) {
    Matrix r = CopyMatrix(m);
    r.pow(p);
    return r;
}
inline Matrix pow(const Matrix & m,const Matrix & p) {
    Matrix r = CopyMatrix(m);
    r.pow(p);
    return r;
}
inline Matrix sqrt(const Matrix & m) {
    Matrix r = CopyMatrix(m);
    r.sqrt();
    return r;
}
inline Matrix rsqrt(const Matrix & m) {
    Matrix r = CopyMatrix(m);
    r.rsqrt();
    return r;
}
inline Matrix sin(const Matrix & m) {
    Matrix r = CopyMatrix(m);
    r.sin();
    return r;
}
inline Matrix cos(const Matrix & m) {
    Matrix r = CopyMatrix(m);
    r.cos();
    return r;
}
inline Matrix tan(const Matrix & m) {
    Matrix r = CopyMatrix(m);
    r.tan();
    return r;
}
inline Matrix asin(const Matrix & m) {
    Matrix r = CopyMatrix(m);
    r.asin();
    return r;
}
inline Matrix acos(const Matrix & m) {
    Matrix r = CopyMatrix(m);
    r.acos();
    return r;
}
inline Matrix atan(const Matrix & m) {
    Matrix r = CopyMatrix(m);
    r.atan();
    return r;
}
inline Matrix sinh(const Matrix & m) {
    Matrix r = CopyMatrix(m);
    r.sinh();
    return r;
}
inline Matrix cosh(const Matrix & m) {
    Matrix r = CopyMatrix(m);
    r.cosh();
    return r;
}
inline Matrix tanh(const Matrix & m) {
    Matrix r = CopyMatrix(m);
    r.tanh();
    return r;
}
inline Matrix asinh(const Matrix & m) {
    Matrix r = CopyMatrix(m);
    r.asinh();
    return r;
}
inline Matrix acosh(const Matrix & m) {
    Matrix r = CopyMatrix(m);
    r.acosh();
    return r;
}
inline Matrix atanh(const Matrix & m) {
    Matrix r = CopyMatrix(m);
    r.atanh();
    return r;
}    
inline Matrix atan2(const Matrix & m,float v) {
    Matrix r = CopyMatrix(m);
    r.atan2(v);
    return r;
}
inline Matrix atan2(const Matrix & m,const Matrix & v) {
    Matrix r = CopyMatrix(m);
    r.atan2(v);
    return r;
}
inline Matrix sigmoid(const Matrix & m) {
    Matrix r = CopyMatrix(m);
    r.sigmoid();
    return r;
}
inline Matrix sigmoid_deriv(const Matrix & m) {
    Matrix r = CopyMatrix(m);
    r.sigmoid_deriv();
    return r;
}
inline Matrix relu(const Matrix & m) {
    Matrix r = CopyMatrix(m);
    r.relu();
    return r;
}
inline Matrix relu_deriv(const Matrix & m) {
    Matrix r = CopyMatrix(m);
    r.relu_deriv();
    return r;
}
inline Matrix softmax(const Matrix & m) {
    Matrix r = CopyMatrix(m);
    r.softmax();
    return r;
}
inline Matrix tanh_deriv(const Matrix & m) {
    Matrix r = CopyMatrix(m);
    r.tanh_deriv();
    return r;
}
inline Matrix cbrt(const Matrix & m) {
    Matrix r = CopyMatrix(m);
    r.cbrt();
    return r;
}
inline Matrix ceil(const Matrix & m) {
    Matrix r = CopyMatrix(m);
    r.ceil();
    return r;
}
/*
inline Matrix copysign(const Matrix & m, float x) {
    Matrix r = CopyMatrix(m);
    r.copysign(x);
    return r;
}

inline Matrix copysign(const Matrix & m, const Matrix & x) {
    Matrix r = CopyMatrix(m);
    r.copysign(x);
    return r;
}
*/
inline Matrix cospi(const Matrix & m) {
    Matrix r = CopyMatrix(m);
    r.cospi();
    return r;
}
inline Matrix cyl_bessel_i0f(const Matrix & m) {
    Matrix r = CopyMatrix(m);
    r.cyl_bessel_i0();
    return r;
}
inline Matrix cyl_bessel_i1f(const Matrix & m) {
    Matrix r = CopyMatrix(m);
    r.cyl_bessel_i1();
    return r;
}
inline Matrix erfcf(const Matrix & m) {
    Matrix r = CopyMatrix(m);
    r.erfc();
    return r;
}
inline Matrix erfcx(const Matrix & m) {
    Matrix r = CopyMatrix(m);
    r.erfcx();
    return r;
}
inline Matrix erfcinv(const Matrix & m) {
    Matrix r = CopyMatrix(m);
    r.erfcinv();
    return r;
}
inline Matrix erf(const Matrix & m) {
    Matrix r = CopyMatrix(m);
    r.erf();
    return r;
}
inline Matrix erfinv(const Matrix & m) {
    Matrix r = CopyMatrix(m);
    r.erfinv();
    return r;
}
inline Matrix exp10(const Matrix & m) {
    Matrix r = CopyMatrix(m);
    r.exp10();
    return r;
}
inline Matrix exp2(const Matrix & m) {
    Matrix r = CopyMatrix(m);
    r.exp2();
    return r;
}
inline Matrix expm1(const Matrix & m) {
    Matrix r = CopyMatrix(m);
    r.expm1();
    return r;
}
inline Matrix fabs(const Matrix & m) {
    Matrix r = CopyMatrix(m);
    r.fabs();
    return r;
}
inline Matrix fdim(const Matrix & m, const Matrix & b) {
    Matrix r = CopyMatrix(m);
    r.fdim(b);
    return r;
}
inline Matrix fmod(const Matrix & m, const Matrix & b) {
    Matrix r = CopyMatrix(m);
    r.fmod(b);
    return r;
}
inline Matrix hypot(const Matrix & m, const Matrix & b) {
    Matrix r = CopyMatrix(m);
    r.hypot(b);
    return r;
}
inline Matrix ilogb(const Matrix & m) {
    Matrix r = CopyMatrix(m);
    r.ilogb();
    return r;
}
inline Matrix j0(const Matrix & m) {
    Matrix r = CopyMatrix(m);
    r.j0();
    return r;
}
inline Matrix j1(const Matrix & m) {
    Matrix r = CopyMatrix(m);
    r.j1();
    return r;
}
inline Matrix jn(const Matrix & m, int n) {
    Matrix r = CopyMatrix(m);
    r.jn(n);
    return r;
}
inline Matrix lgamma(const Matrix & m) {
    Matrix r = CopyMatrix(m);
    r.lgamma();
    return r;
}
inline Matrix log1p(const Matrix & m) {
    Matrix r = CopyMatrix(m);
    r.log1p();
    return r;
}
inline Matrix logb(const Matrix & m) {
    Matrix r = CopyMatrix(m);
    r.logb();
    return r;
}
inline Matrix norm3d(const Matrix & m,const Matrix & a, const Matrix & b, const Matrix & c) {
    Matrix r = CopyMatrix(m);
    r.norm3d(a,b,c);
    return r;
}
inline Matrix norm4d(const Matrix & m,const Matrix & a, const Matrix & b, const Matrix & c, const Matrix & d) {
    Matrix r = CopyMatrix(m);
    r.norm4d(a,b,c,d);
    return r;
}
inline Matrix normcdf(const Matrix & m) {
    Matrix r = CopyMatrix(m);
    r.normcdf();
    return r;
}
inline Matrix normcdfinv(const Matrix & m) {
    Matrix r = CopyMatrix(m);
    r.normcdfinv();
    return r;
}
inline Matrix norm(const Matrix & m) {
    Matrix r = CopyMatrix(m);
    r.norm();
    return r;
}
inline Matrix rcbrt(const Matrix & m) {
    Matrix r = CopyMatrix(m);
    r.rcbrt();
    return r;
}
inline Matrix rhypot(const Matrix & m, const Matrix & b) {
    Matrix r = CopyMatrix(m);
    r.rhypot(b);
    return r;
}
inline Matrix rnorm3d(const Matrix & m,const Matrix & a, const Matrix & b, const Matrix & c) {
    Matrix r = CopyMatrix(m);
    r.rnorm3d(a,b,c);
    return r;
}
inline Matrix rnorm4d(const Matrix & m,const Matrix & a, const Matrix & b, const Matrix & c, const Matrix & d) {
    Matrix r = CopyMatrix(m);
    r.rnorm4d(a,b,c,d);
    return r;
}
inline Matrix rnorm(const Matrix & m, const Matrix & b) {
    Matrix r = CopyMatrix(m);
    r.rnorm();
    return r;
}
inline Matrix tgamma(const Matrix & m) {
    Matrix r = CopyMatrix(m);
    r.tgamma();
    return r;
}
inline Matrix trunc(const Matrix & m) {
    Matrix r = CopyMatrix(m);
    r.tgamma();
    return r;
}
inline Matrix y0(const Matrix & m) {
    Matrix r = CopyMatrix(m);
    r.y0();
    return r;
}
inline Matrix y1(const Matrix & m) {
    Matrix r = CopyMatrix(m);
    r.y1();
    return r;
}
inline Matrix yn(const Matrix & m, int n) {
    Matrix r = CopyMatrix(m);
    r.yn(n);
    return r;
}
inline float sum(const Matrix & m) {
    float r = 0;
    cublasStatus_t     status; 
    status = cublasSasum(cublas->handle,m.M*m.N,m.dev,1,&r);
    ASSERTOK(status);
    return r;
}
inline float sum(const Vector & v) {
    float r = 0;
    cublasStatus_t     status; 
    status = cublasSasum(cublas->handle,v.N,v.dev,1,&r);
    ASSERTOK(status);
    return r;
}

inline Vector abs(const Vector & m) {    
    Vector r = CopyVector(m);
    r.abs();
    return r;
}
inline Vector exp(const Vector & m) {
    Vector r = CopyVector(m);
    r.exp(); 
    return r;   
}    
inline Vector log2(const Vector & m) {
    Vector r = CopyVector(m);
    r.log2();
    return r;
}
inline Vector log10(const Vector & m) {
    Vector r = CopyVector(m);
    r.log10();
    return r;
}
inline Vector pow(const Vector & m,float p) {
    Vector r = CopyVector(m);
    r.pow(p);
    return r;
}
inline Vector pow(const Vector & m,const Vector & p) {
    Vector r = CopyVector(m);
    r.pow(p);
    return r;
}
inline Vector sqrt(const Vector & m) {
    Vector r = CopyVector(m);
    r.sqrt();
    return r;
}
inline Vector rsqrt(const Vector & m) {
    Vector r = CopyVector(m);
    r.rsqrt();
    return r;
}
inline Vector sin(const Vector & m) {
    Vector r = CopyVector(m);
    r.sin();
    return r;
}
inline Vector cos(const Vector & m) {
    Vector r = CopyVector(m);
    r.cos();
    return r;
}
inline Vector tan(const Vector & m) {
    Vector r = CopyVector(m);
    r.tan();
    return r;
}
inline Vector asin(const Vector & m) {
    Vector r = CopyVector(m);
    r.asin();
    return r;
}
inline Vector acos(const Vector & m) {
    Vector r = CopyVector(m);
    r.acos();
    return r;
}
inline Vector atan(const Vector & m) {
    Vector r = CopyVector(m);
    r.atan();
    return r;
}
inline Vector sinh(const Vector & m) {
    Vector r = CopyVector(m);
    r.sinh();
    return r;
}
inline Vector cosh(const Vector & m) {
    Vector r = CopyVector(m);
    r.cosh();
    return r;
}
inline Vector tanh(const Vector & m) {
    Vector r = CopyVector(m);
    r.tanh();
    return r;
}
inline Vector asinh(const Vector & m) {
    Vector r = CopyVector(m);
    r.asinh();
    return r;
}
inline Vector acosh(const Vector & m) {
    Vector r = CopyVector(m);
    r.acosh();
    return r;
}
inline Vector atanh(const Vector & m) {
    Vector r = CopyVector(m);
    r.atanh();
    return r;
}    
inline Vector atan2(const Vector & m,float v) {
    Vector r = CopyVector(m);
    r.atan2(v);
    return r;
}
inline Vector atan2(const Vector & m,const Vector & v) {
    Vector r = CopyVector(m);
    r.atan2(v);
    return r;
}
inline Vector sigmoid(const Vector & m) {
    Vector r = CopyVector(m);
    r.sigmoid();
    return r;
}
inline Vector sigmoid_deriv(const Vector & m) {
    Vector r = CopyVector(m);
    r.sigmoid_deriv();
    return r;
}
inline Vector relu(const Vector & m) {
    Vector r = CopyVector(m);
    r.relu();
    return r;
}
inline Vector relu_deriv(const Vector & m) {
    Vector r = CopyVector(m);
    r.relu_deriv();
    return r;
}
inline Vector softmax(const Vector & m) {
    Vector r = CopyVector(m);
    r.softmax();
    return r;
}
inline Vector tanh_deriv(const Vector & m) {
    Vector r = CopyVector(m);
    r.tanh_deriv();
    return r;
}

inline Vector cbrt(const Vector & m) {
    Vector r = CopyVector(m);
    r.cbrt();
    return r;
}
inline Vector ceil(const Vector & m) {
    Vector r = CopyVector(m);
    r.ceil();
    return r;
}
/*
inline Vector copysign(const Vector & m, float x) {
    Vector r = CopyVector(m);
    r.copysign(x);
    return r;
}
inline Vector copysign(const Vector & m, const Vector & x) {
    Vector r = CopyVector(m);
    r.copysign(x);
    return r;
}
*/
inline Vector cospi(const Vector & m) {
    Vector r = CopyVector(m);
    r.cospi();
    return r;
}
inline Vector cyl_bessel_i0f(const Vector & m) {
    Vector r = CopyVector(m);
    r.cyl_bessel_i0();
    return r;
}
inline Vector cyl_bessel_i1f(const Vector & m) {
    Vector r = CopyVector(m);
    r.cyl_bessel_i1();
    return r;
}
inline Vector erfc(const Vector & m) {
    Vector r = CopyVector(m);
    r.erfc();
    return r;
}
inline Vector erfcx(const Vector & m) {
    Vector r = CopyVector(m);
    r.erfcx();
    return r;
}
inline Vector erfcinv(const Vector & m) {
    Vector r = CopyVector(m);
    r.erfcinv();
    return r;
}
inline Vector erf(const Vector & m) {
    Vector r = CopyVector(m);
    r.erf();
    return r;
}
inline Vector erfinv(const Vector & m) {
    Vector r = CopyVector(m);
    r.erfinv();
    return r;
}
inline Vector exp10(const Vector & m) {
    Vector r = CopyVector(m);
    r.exp10();
    return r;
}
inline Vector exp2(const Vector & m) {
    Vector r = CopyVector(m);
    r.exp2();
    return r;
}
inline Vector expm1(const Vector & m) {
    Vector r = CopyVector(m);
    r.expm1();
    return r;
}
inline Vector fabs(const Vector & m) {
    Vector r = CopyVector(m);
    r.fabs();
    return r;
}
inline Vector fdim(const Vector & m, const Vector & b) {
    Vector r = CopyVector(m);
    r.fdim(b);
    return r;
}
inline Vector fmod(const Vector & m, const Vector & b) {
    Vector r = CopyVector(m);
    r.fmod(b);
    return r;
}
inline Vector hypot(const Vector & m, const Vector & b) {
    Vector r = CopyVector(m);
    r.hypot(b);
    return r;
}
inline Vector ilogb(const Vector & m) {
    Vector r = CopyVector(m);
    r.ilogb();
    return r;
}
inline Vector j0(const Vector & m) {
    Vector r = CopyVector(m);
    r.j0();
    return r;
}
inline Vector j1(const Vector & m) {
    Vector r = CopyVector(m);
    r.j1();
    return r;
}
inline Vector jn(const Vector & m, int n) {
    Vector r = CopyVector(m);
    r.jn(n);
    return r;
}
inline Vector lgamma(const Vector & m) {
    Vector r = CopyVector(m);
    r.lgamma();
    return r;
}
inline Vector log1p(const Vector & m) {
    Vector r = CopyVector(m);
    r.log1p();
    return r;
}
inline Vector logb(const Vector & m) {
    Vector r = CopyVector(m);
    r.logb();
    return r;
}
inline Vector norm3d(const Vector & m,const Vector & a, const Vector & b, const Vector & c) {
    Vector r = CopyVector(m);
    r.norm3d(a,b,c);
    return r;
}
inline Vector norm4d(const Vector & m,const Vector & a, const Vector & b, const Vector & c, const Vector & d) {
    Vector r = CopyVector(m);
    r.norm4d(a,b,c,d);
    return r;
}
inline Vector normcdf(const Vector & m) {
    Vector r = CopyVector(m);
    r.normcdf();
    return r;
}
inline Vector normcdfinv(const Vector & m) {
    Vector r = CopyVector(m);
    r.normcdfinv();
    return r;
}
inline Vector rcbrt(const Vector & m) {
    Vector r = CopyVector(m);
    r.rcbrt();
    return r;
}
inline Vector norm(const Vector & m) {
    Vector r = CopyVector(m);
    r.norm();
    return r;
}
inline Vector rhypot(const Vector & m, const Vector & b) {
    Vector r = CopyVector(m);
    r.rhypot(b);
    return r;
}
/*
inline Vector rint(const Vector & m) {
    Vector r = CopyVector(m);
    r.rint();
    return r;
}*/
inline Vector rnorm3d(const Vector & m,const Vector & a, const Vector & b, const Vector & c) {
    Vector r = CopyVector(m);
    r.rnorm3d(a,b,c);
    return r;
}
inline Vector rnorm4d(const Vector & m,const Vector & a, const Vector & b, const Vector & c, const Vector & d) {
    Vector r = CopyVector(m);
    r.rnorm4d(a,b,c,d);
    return r;
}
/*
inline Vector round(const Vector & m) {
    Vector r = CopyVector(m);
    r.round();
    return r;
}*/
inline Vector rnorm(const Vector & m, const Vector & b) {
    Vector r = CopyVector(m);
    r.rnorm();
    return r;
}
inline Vector tgamma(const Vector & m) {
    Vector r = CopyVector(m);
    r.tgamma();
    return r;
}
inline Vector trunc(const Vector & m) {
    Vector r = CopyVector(m);
    r.tgamma();
    return r;
}
inline Vector y0(const Vector & m) {
    Vector r = CopyVector(m);
    r.y0();
    return r;
}
inline Vector y1(const Vector & m) {
    Vector r = CopyVector(m);
    r.y1();
    return r;
}
inline Vector yn(const Vector & m, int n) {
    Vector r = CopyVector(m);
    r.yn(n);
    return r;
}


    

#endif
