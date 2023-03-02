#pragma once

#include "Std.h"
#include <ccomplex>
#include <complex>
#include "cppmkl/cppmkl_allocator.h"
#include "cppmkl/cppmkl_vml.h"
#include "cppmkl/cppmkl_cblas.h"


namespace MKL {

    template<typename T>
    struct Vector
    {
        std::vector<T, cppmkl::cppmkl_allocator<T> > vector;
        
        Vector() = default;
        Vector(size_t n) { assert(n > 0); vector.resize(n); }
        Vector(const Vector<T> & v) { vector = v.vector; }
        Vector(const std::vector<T> & v) {
            vector.resize(v.size());
            memcpy(vector.data(),v.data(),v.size()*sizeof(T));
        }

        Vector<T>& operator = (const Vector<T> & v)
        {
            vector = v.vector;
            return *this;
        }

        T min() { return *std::min_element(vector.begin(),vector.end()); }
        T max() { return *std::max_element(vector.begin(),vector.end()); }
        size_t min_index() { return std::distance(vector.begin(), std::min_element(vector.begin(),vector.end())); }
        size_t max_index() { return std::distance(vector.begin(), std::max_element(vector.begin(),vector.end())); }
        
        T& operator[](size_t i) { return vector[i]; }
        const T& operator[](size_t i) const { return vector[i]; }
        T __getitem(size_t i) { return vector[i]; }
        void __setitem(size_t i, T x) { vector[i] = x; }

        size_t size() const { return vector.size(); }
        T* data() { return vector.data(); }
        void resize(size_t i) { vector.resize(i); }

        T front() const { return vector.front(); }
        T back() const  { return vector.back();  }

        void push_back(T x) { vector.push_back(x); }
        T    pop_back() { T x = vector.back(); vector.pop_back(); return x; }

        Vector<T>& operator +=  (const Vector<T> & v) { 
            cppmkl::vadd(vector,v.vector,vector);
            return *this;
        }
        Vector<T>& operator -=  (const Vector<T> & v) { 
            cppmkl::vsub(vector,v.vector,vector);
            return *this;
        }
        Vector<T>& operator *=  (const Vector<T> & v) { 
            cppmkl::vmul(vector,v.vector,vector);
            return *this;
        }
        Vector<T>& operator /=  (const Vector<T> & v) { 
            cppmkl::vdiv(vector,v.vector,vector);
            return *this;
        }
        Vector<T>& operator %=  (const Vector<T> & v) { 
            cppmkl::vmodf(vector,v.vector,vector);
            return *this;
        }

        Vector<T> operator + (const Vector<T> & v) { 
            Vector<T> r(size());            
            cppmkl::vadd(vector,v.vector,r.vector);
            return r;
        }
        Vector<T> operator - (const Vector<T> & v) { 
            Vector<T> r(size());            
            cppmkl::vsub(vector,v.vector,r.vector);
            return r;            
        }
        Vector<T> operator * (const Vector<T> & v) { 
            Vector<T> r(size());            
            cppmkl::vmul(vector,v.vector,r.vector);
            return r;
        }
        Vector<T> operator / (const Vector<T> & v) { 
            Vector<T> r(size());            
            cppmkl::vdiv(vector,v.vector,r.vector);
            return r;
        }
        Vector<T> operator % (const Vector<T> & v) { 
            Vector<T> r(size());            
            cppmkl::vmodf(vector,v.vector,r.vector);
            return r;
        }        

        void zero() {
            memset(vector.data(),0x00,size()*sizeof(T));
        }
        void fill(T x) {
            for(size_t i = 0; i < size(); i++) vector[i] = x;
        }
        void ones() {
            fill((T)1);
        }
        void random(T min = T(0), T max = T(1)) {
            Std::StdRandomUniform r;
            for(size_t i = 0; i < size(); i++) vector[i] = r.random(min,max);
        }

        Vector<T> slice(size_t start, size_t len) {
            Vector<T> x(len);
            memcpy(x.vector.data(),vector.data()+start,len*sizeof(T));
            return x;
        }
    };

    template<typename T>
    struct RealFFT1D
    {
        DFTI_DESCRIPTOR_HANDLE handle1;        
        size_t size;
        
        RealFFT1D(size_t size) {
            DFTI_CONFIG_VALUE prec;
            if(typeid(T) == typeid(float)) prec = DFTI_SINGLE;
            else prec = DFTI_DOUBLE;
            DftiCreateDescriptor(&handle1, prec, DFTI_REAL,  1, size );
            DftiSetValue(handle1, DFTI_PLACEMENT, DFTI_NOT_INPLACE); //Out of place FFT
            DftiSetValue(handle1, DFTI_BACKWARD_SCALE, 1.0f / size);
            DftiCommitDescriptor(handle1);            
            this->size = size;
        }
        ~RealFFT1D() {
            DftiFreeDescriptor(&handle1);            
        }

        void Forward( Vector<T> & input, Vector<std::complex<T>> & output) {
            output.resize(size);
            Vector<float> x(size*2);            
            DftiComputeForward(handle1, input.vector.data(),x.vector.data());
            memcpy(output.data(),x.data(), x.size()*sizeof(float));            
        }
        void Backward( Vector<std::complex<T>> & input, Vector<T> & output) {
            output.resize(size);
            Vector<float> x(size*2);            
            memcpy(x.data(),input.data(),x.size()*sizeof(float));
            DftiComputeBackward(handle1, x.vector.data(), output.vector.data());
        }                
    };

    template<typename T = float>
    struct ComplexFFT1D
    {
        DFTI_DESCRIPTOR_HANDLE handle1;        
        size_t size;
        
        ComplexFFT1D(size_t size) {
            DFTI_CONFIG_VALUE prec;
            if(typeid(T) == typeid(float)) prec = DFTI_SINGLE;
            else prec = DFTI_DOUBLE;
            DftiCreateDescriptor(&handle1, prec, DFTI_COMPLEX, 1, size );
            DftiSetValue(handle1, DFTI_PLACEMENT, DFTI_NOT_INPLACE); //Out of place FFT
            DftiCommitDescriptor(handle1);            
            this->size = size;
        }
        ~ComplexFFT1D() {
            DftiFreeDescriptor(&handle1);            
        }

        void Forward( Vector<std::complex<T>> & input, Vector<std::complex<T>> & output) {
            output.resize(size);
            DftiComputeForward(handle1, input.vector.data(),output.vector.data());
        }
        void Backward( Vector<std::complex<T>> & input, Vector<std::complex<T>> & output) {
            output.resize(size);
            DftiComputeBackward(handle1, input.vector.data(), output.vector.data());
        }        
    };

    template<typename T>
    Vector<T> sqr(Vector<T> & a) {
        Vector<T> r(a.size());                
        cppmkl::vsqr(a.vector,r.vector);
        return r;
    }
    template<typename T>
    Vector<T> abs(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vabs(a.vector,r.vector);
        return r;
    }
    template<typename T>
    Vector<T> inv(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vinv(a.vector,r.vector);
        return r;            
    }
    template<typename T>
    Vector<T> sqrt(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vsqrt(a.vector,r.vector);
        return r;
    }
    template<typename T>
    Vector<T> rsqrt(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vinvsqrt(a.vector,r.vector);
        return r;
    }
    template<typename T>
    Vector<T> cbrt(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vcbrt(a.vector,r.vector);
        return r;
    }
    template<typename T>
    Vector<T> rcbrt(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vinvcbrt(a.vector,r.vector);
        return r;
    }
    template<typename T>
    Vector<T> pow(const Vector<T> & a,const Vector<T> & b) {
        Vector<T> r(a.size());
        cppmkl::vpow(a.vector,b.vector,r.vector);
        return r;
    }
    template<typename T>
    Vector<T> pow2o3(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vpow2o3(a.vector,r.vector);
        return r;
    }
    template<typename T>
    Vector<T> pow3o2(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vpow3o2(a.vector,r.vector);
        return r;
    }
    template<typename T>
    Vector<T> pow(const Vector<T> & a,const T b) {
        Vector<T> r(a.size());
        cppmkl::vpowx(a.vector,b,r.vector);
        return r;
    }
    template<typename T>
    Vector<T> hypot(const Vector<T> & a,const Vector<T> & b) {
        Vector<T> r(a.size());
        cppmkl::vhypot(a.vector,b.vector,r.vector);
        return r;
    }
    template<typename T>
    Vector<T> exp(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vexp(a.vector,r.vector);
        return r;
    }
    template<typename T>
    Vector<T> exp2(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vexp2(a.vector,r.vector);
        return r;
    }
    template<typename T>
    Vector<T> exp10(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vexp10(a.vector,r.vector);
        return r;
    }
    template<typename T>
    Vector<T> expm1(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vexpm1(a.vector,r.vector);
        return r;
    }
    template<typename T>
    Vector<T> ln(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vln(a.vector,r.vector);
        return r;
    }
    template<typename T>
    Vector<T> log10(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vlog10(a.vector,r.vector);
        return r;
    }
    template<typename T>
    Vector<T> log2(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vlog2(a.vector,r.vector);
        return r;
    }
    template<typename T>
    Vector<T> logb(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vlogb(a.vector,r.vector);
        return r;
    }
    template<typename T>
    Vector<T> log1p(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vlog1p(a.vector,r.vector);
        return r;
    }
    template<typename T>
    Vector<T> cos(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vcos(a.vector,r.vector);
        return r;
    }
    template<typename T>
    Vector<T> sin(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vsin(a.vector,r.vector);
        return r;            
    }
    template<typename T>
    Vector<T> tan(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vtan(a.vector,r.vector);
        return r;
    }
    template<typename T>
    Vector<T> cosh(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vcosh(a.vector,r.vector);
        return r;
    }
    template<typename T>
    Vector<T> sinh(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vsinh(a.vector,r.vector);
        return r;            
    }
    template<typename T>
    Vector<T> tanh(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vtanh(a.vector,r.vector);
        return r;            
    }
    template<typename T>
    Vector<T> acos(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vacos(a.vector,r.vector);
        return r;
    }
    template<typename T>
    Vector<T> asin(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vasin(a.vector,r.vector);
        return r;            
    }
    template<typename T>
    Vector<T> atan(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vatan(a.vector,r.vector);
        return r;
    }
    template<typename T>
    Vector<T> atan2(const Vector<T> & a,const Vector<T> &n) {
        Vector<T> r(a.size());
        cppmkl::vatan2(a.vector,n.vector,r.vector);
        return r;
    }
    template<typename T>
    Vector<T> acosh(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vacosh(a.vector,r.vector);
        return r;
    }
    template<typename T>
    Vector<T> asinh(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vasinh(a.vector,r.vector);
        return r;            
    }
    template<typename T>
    Vector<T> atanh(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vatanh(a.vector,r.vector);
        return r;
    }        
    template<typename T>
    void sincos(const Vector<T> & a, Vector<T> & b, Vector<T> & r) {        
        cppmkl::vsincos(a.vector,b.vector,r.vector);        
    }
    template<typename T>
    Vector<T> erf(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::verf(a.vector,r.vector);
        return r;
    }
    template<typename T>
    Vector<T> erfinv(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::verfinv(a.vector,r.vector);
        return r;
    }
    template<typename T>
    Vector<T> erfc(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::verfc(a.vector,r.vector);
        return r;
    }
    template<typename T>
    Vector<T> cdfnorm(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vcdfnorm(a.vector,r.vector);
        return r;
    }    
    template<typename T>
    Vector<T> cdfnorminv(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vcdfnorm(a.vector,r.vector);
        return r;
    }    
    template<typename T>
    Vector<T> floor(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vfloor(a.vector,r.vector);
        return r;
    }    
    template<typename T>
    Vector<T> ceil(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vceil(a.vector,r.vector);
        return r;
    }    
    template<typename T>
    Vector<T> trunc(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vtrunc(a.vector,r.vector);
        return r;
    }    
    template<typename T>
    Vector<T> round(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vround(a.vector,r.vector);
        return r;        
    }    
    template<typename T>
    Vector<T> nearbyint(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vnearbyint(a.vector,r.vector);
        return r;
    }    
    template<typename T>
    Vector<T> rint(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vrint(a.vector,r.vector);
        return r;
    }    
    template<typename T>
    Vector<T> fmod(const Vector<T> & a, Vector<T> & b) {
        Vector<T> r(a.size());
        cppmkl::vmodf(a.vector,b.vector,r.vector);
        return r;
    }    
    /*
    template<typename T>
    Vector<T> mulbyconj(const Vector<T> & a, const Vector<T> & b) {
        Vector<T> r(a.size());
        cppmkl::vmulbyconj(a.vector,b.vector.data(),r.vector);
        return r;
    }    
    template<typename T>
    Vector<T> conj(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vconj(a.vector,r.vector);
        return r;
    }    
    template<typename T>
    Vector<T> arg(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::varg(a.vector,r.vector);
        return r;
    } 
    */   
    template<typename T>
    Vector<T> CIS(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vCIS(a.vector,r.vector);
        return r;
    }    
    template<typename T>
    Vector<T> cospi(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vcospi(a.vector,r.vector);
        return r;
    }    
    template<typename T>
    Vector<T> sinpi(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vsinpi(a.vector,r.vector);
        return r;
    }    
    template<typename T>
    Vector<T> tanpi(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vtanpi(a.vector,r.vector);
        return r;
    }    
    template<typename T>
    Vector<T> acospi(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vacospi(a.vector,r.vector);
        return r;
    }    
    template<typename T>
    Vector<T> asinpi(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vasinpi(a.vector,r.vector);
        return r;
    }    
    template<typename T>
    Vector<T> atanpi(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vatanpi(a.vector,r.vector);
        return r;
    }    
    template<typename T>
    Vector<T> atan2pi(const Vector<T> & a, Vector<T> & b) {
        Vector<T> r(a.size());
        cppmkl::vatan2pi(a.vector,b.vector,r.vector);
        return r;
    }    
    template<typename T>
    Vector<T> cosd(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vcosd(a.size(),a.vector.data(),r.vector.data());
        return r;
    }    
    template<typename T>
    Vector<T> sind(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vsind(a.size(),a.vector.data(),r.vector.data());
        return r;
    }    
    template<typename T>
    Vector<T> tand(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vtand(a.size(),a.vector.data(),r.vector.data());
        return r;
    }       
    template<typename T>
    Vector<T> lgamma(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vlgamma(a.vector,r.vector);
        return r;
    }       
    template<typename T>
    Vector<T> tgamma(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vtgamma(a.vector,r.vector);
        return r;
    }       
    template<typename T>
    Vector<T> expint1(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vexpint1(a.vector,r.vector);
        return r;
    }       
    template<typename T>
    Vector<T> copy(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::cblas_copy(a.size(),a.vector.data(),1,r.vector.data(),1);
        return r;
    }       

    template<typename T> T sum(const Vector<T> & a) {        
        return cppmkl::cblas_asum(a.size(), a.vector.data(),1);        
    }       

    template<typename T>
    Vector<T> add(const Vector<T> & a, const Vector<T> & b) {
        Vector<T> r(b);
        cppmkl::cblas_axpy(a.size(),1.0,a.vector.data(),1,r.vector.data(),1);
        return r;
    }       
    template<typename T>
    Vector<T> sub(const Vector<T> & a, const Vector<T> & b) {
        Vector<T> r(a);
        cppmkl::cblas_axpy(a.size(),-1.0,b.vector.data(),1,r.vector.data(),1);
        return r;
    }       
    template<typename T>
    T dot(const Vector<T> & a, Vector<T> & b) {
        return cppmkl::cblas_dot(a.vector,b.vector);
    }       
    template<typename T>
    T nrm2(const Vector<T> & a) {
        Vector<T> r(a);
        return cppmkl::cblas_nrm2(a.vector);        
    }       
    
    template<typename T>
    void scale(Vector<T> & x, T alpha) {
        cppmkl::cblas_scal(x.size(),alpha,x.vector.data(),1);
    }


    template<typename T>
    size_t min_index(const Matrix<T> & m) { return cppmkl::cblas_iamin(m.size(),m.matrix.vector.data(),1); }
    template<typename T>
    size_t max_index(const Matrix<T> & m) { return cppmkl::cblas_iamax(m.size(),m.matrix.vector.data(),1); }
    
    template<typename T>
    size_t min_index(const Vector<T> & v) { return cppmkl::cblas_iamin(v.size(),v.vector.data(),1); }
    template<typename T>
    size_t max_index(const Vector<T> & v) { return cppmkl::cblas_iamax(v.size(),v.vector.data(),1); }

    template<typename T>
    Vector<T> linspace(size_t n, T start, T inc=1) {
        Vector<T> r(n);
        for(size_t i = 0; i < n; i++) {
            r[i] = start + i*inc;
        }
        return r;
    }
    template<typename T>
    Vector<T> linspace(T start, T end, T inc=1) {
        size_t n = (end - start)/inc;
        Vector<T> r(n);
        for(size_t i = 0; i < n; i++) {
            r[i] = start + i*inc;
        }
        return r;
    }
}