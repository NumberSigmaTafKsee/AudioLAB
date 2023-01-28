#pragma once

#include <iostream>
#include <ccomplex>
#include <complex>
#include <vector>
#include <algorithm>
#include <cstring>
#include <chrono>
#include <random>
#include <functional>

#include "StdNoise.hpp"

#include "cppmkl/cppmkl_allocator.h"
#include "cppmkl/cppmkl_vml.h"
#include "cppmkl/cppmkl_cblas.h"
#include "cppmkl/matrix.h"

//#include "samples/Allocator.hpp"

namespace Casino::MKL
{
    
    template<typename T> using vector_base = std::vector<T, cppmkl::cppmkl_allocator<T> >;

    template<typename T>
    struct Vector : public vector_base<T>
    {                
        using vector_base<T>::size;
        using vector_base<T>::resize;
        using vector_base<T>::data;
        using vector_base<T>::push_back;
        using vector_base<T>::pop_back;
        using vector_base<T>::front;
        using vector_base<T>::back;
        using vector_base<T>::at;
        using vector_base<T>::operator [];
        using vector_base<T>::operator =;

        Vector() = default;
        Vector(size_t n) { assert(n > 0); resize(n); zero(); }
        Vector(const Vector<T> & v) { (*this) = v; }
        Vector(const std::vector<T> & v) {
            resize(v.size());
            memcpy(this->data(),v.data(),v.size()*sizeof(T));
        }

        T min() { return *std::min_element(this->begin(),this->end()); }
        T max() { return *std::max_element(this->begin(),this->end()); }

        size_t min_index() { return std::distance(this->begin(), std::min_element(this->begin(),this->end())); }
        size_t max_index() { return std::distance(this->begin(), std::max_element(this->begin(),this->end())); }
                
        Vector<T>& operator +=  (const Vector<T> & v) { 
            cppmkl::vadd(*this,v,*this);
            return *this;
        }
        Vector<T>& operator -=  (const Vector<T> & v) { 
            cppmkl::vsub(*this,v,*this);
            return *this;
        }
        Vector<T>& operator *=  (const Vector<T> & v) { 
            cppmkl::vmul(*this,v,*this);
            return *this;
        }
        Vector<T>& operator /=  (const Vector<T> & v) { 
            cppmkl::vdiv(*this,v,*this);
            return *this;        
	    }

        Vector<T> operator - () {
            Vector<T> r(*this);
            return (T)-1.0 * r;
        }

        Vector<T>& operator = (const T v)
        {
            fill(v);
            return *this;
        }
                
        void zero() {
            memset(data(),0x00,size()*sizeof(T));
        }
        void zeros() { zero(); }
        void fill(T x) {
            for(size_t i = 0; i < size(); i++) (*this)[i] = x;
        }
        void ones() {
            fill((T)1);
        }
        void random(T min = T(0), T max = T(1)) {
            Default noise;
            for(size_t i = 0; i < size(); i++) (*this)[i] = noise.random(min,max);
        }
        void randu(T min = T(0), T max = T(1)) { random(min,max); }
        void randn(T min = T(0), T max = T(1)) { random(min,max); }

        
        void clamp(T min = T(-1), T max = T(1)) {
            for(size_t i = 0; i < size(); i++)
            {
                if((*this)[i] < min) (*this)[i] = min;
                if((*this)[i] < max) (*this)[i] = max;
            }
        }

        void set_size(size_t n) { resize(n); }


        Vector<T> eval() {
            return Vector<T>(*this); 
        }
        Vector<T> slice(size_t start, size_t len) {
            Vector<T> x(len);
            memcpy(x.data(),data()+start,len*sizeof(T));
            return x;
        }

        void print() {
            std::cout << "Vector[" << size() << "]=";
            for(size_t i = 0; i < size()-1; i++) std::cout << (*this)[i] << ",";
            std::cout << (*this)[size()-1] << std::endl;
        }
        
        // eigen compatibility
        void setZero() { zero(); }
        void setOnes() { ones(); }
        void setRandom() { random(); }
        
        
    };

    
    template<typename T>
    std::ostream& operator << (std::ostream & o, const Vector<T> & v )
    {
        for(size_t i = 0; i < v.size(); i++)
        {               
            o << v[i] << ",";            
        }
        o << std::endl;
        return o;
    }
        
    template<typename T>
    struct Matrix : public Vector<T>
    {        
        size_t M;
        size_t N;

        Matrix() { M = N = 0; }
        Matrix(size_t m, size_t n) {
            resize(m,n);
            assert(M > 0);
            assert(N > 0);
        }
        Matrix(const Matrix<T> & m) {
            
            resize(m.M,m.N);
            memcpy(data(),m.data(),size()*sizeof(T));            
        }
        Matrix(T * ptr, size_t m, size_t n)
        {
            resize(m,n);            
            memcpy(data(),ptr,m*n*sizeof(T));
        }
        Matrix<T>& operator = (const T v) {
            fill(v);
            return *this;
        }

        Matrix<T>& operator = (const Matrix<T> & m) {            
            resize(m.M,m.N);
            memcpy(data(),m.data(),size()*sizeof(T));            
            return *this;
        }
        Matrix<T>& operator = (const Vector<T> & m) {            
            resize(1,m.size());
            memcpy(data(),m.data(),size()*sizeof(T));            
            return *this;
        }

        using Vector<T>::size;
        using Vector<T>::data;
        using Vector<T>::resize;
        using Vector<T>::at;
        using Vector<T>::operator [];

        size_t rows() const { return M; }
        size_t cols() const { return N; }
        
        Matrix<T> cwiseMax(T v) {
            Matrix<T> r(*this);
            for(size_t i = 0; i < rows(); i++)
                for(size_t j = 0;  j < cols(); j++)
                    if(r(i,j) < v) r(i,j) = v;
            return r;
        }
        Matrix<T> row(size_t m) { 
            Matrix<T> r(1,cols());
            for(size_t i = 0; i < cols(); i++) r(0,i) = (*this)(m,i);
            return r;
        }
        Matrix<T> col(size_t c) { 
            Matrix<T> r(cols(),1);
            for(size_t i = 0; i < rows(); i++) r(i,0) = (*this)(i,c);
            return r;
        }
        void row(size_t m, const Vector<T> & v)
        {
            for(size_t i = 0; i < cols(); i++) (*this)(m,i) = v[i];
        }
        void col(size_t n, const Vector<T> & v)
        {
            for(size_t i = 0; i < rows(); i++) (*this)(i,n) = v[i];
        }
        
        void resize(size_t r, size_t c) {
            M = r;
            N = c;
            resize(r*c);            
        }

        T& operator()(size_t i, size_t j) {             
            return (*this)[i*N + j]; }

        T  operator()(size_t i, size_t j) const {             
            return (*this)[i*N + j]; }
        
        
        std::ostream& operator << (std::ostream & o )
        {
            for(size_t i = 0; i < rows(); i++)
            {
                for(size_t j = 0; j < cols(); j++)
                    std::cout << (*this)(i,j) << ",";
                std::cout << std::endl;
            }
            return o;
        }

        Matrix<T> operator - () {
            Matrix<T> r(*this);
            return T(-1.0)*r;
        }
        Matrix<T> addToEachRow(Vector<T> & v) {
            Matrix<T> r(*this);
            for(size_t i = 0; i < M; i++)
            {
                for(size_t j = 0; j < N; j++)
                {
                    r(i,j) += v[j];
                }
            }
            return r;
        }
        Matrix<T> eval() {
            Matrix<T> r(*this);
            return r;
        }
        Matrix<T>& operator += (const T b)
        {
            Matrix<T> x(rows(),cols());
            x.fill(b);
            cppmkl::vadd(*this,x,*this);
            return *this;
        }
        Matrix<T>& operator -= (const T b)
        {
            Matrix<T> x(rows(),cols());
            x.fill(b);
            cppmkl::vsub(*this,x,*this);
            return *this;
        }
        Matrix<T>& operator *= (const T b)
        {
            Matrix<T> x(rows(),cols());
            x.fill(b);
            cppmkl::vmul(*this,x,*this);
            return *this;
        }
        Matrix<T>& operator /= (const T b)
        {
            Matrix<T> x(rows(),cols());
            x.fill(b);
            cppmkl::vdiv(*this,x,*this);
            return *this;
        }
        
        Matrix<T>& operator += (const Matrix<T> b)
        {
            assert(rows() == b.rows() && cols() == b.cols());            
            cppmkl::vadd(*this,b,*this);
            return *this;
        }
        Matrix<T>& operator -= (const Matrix<T> b)
        {
            assert(rows() == b.rows() && cols() == b.cols());            ;
            cppmkl::vsub(*this,b,*this);
            return *this;
        }
        Matrix<T>& operator *= (const Matrix<T> b)
        {
            assert(rows() == b.rows() && cols() == b.cols());            
            cppmkl::vmul(*this,b,*this);
            return *this;
        }
        Matrix<T>& operator /= (const Matrix<T> b)
        {
            assert(rows() == b.rows() && cols() == b.cols());            
            cppmkl::vdiv(*this,b,*this);
            return *this;
        }

        Matrix<T> hadamard(const Matrix<T> & m) {
            Matrix<T> r(rows(),cols());
            assert(rows() == m.rows() && cols() == m.cols());
            cppmkl::vmul(*this,m,r);
            return r;
        }
        
        Matrix<T> transpose() {
            Matrix<T> r(cols(),rows());
            for(size_t i = 0; i < rows(); i++)
                for(size_t j = 0; j < cols(); j++)
                    r(j,i) =(*this)[i*N + j];
            return r;
        }
        Matrix<T> t() {
            return transpose();
        }
        void zero() {
            memset(data(),0x00,size()*sizeof(T));
        }
        void zeros() { zero(); }

        void fill(T x) {
            for(size_t i = 0; i < size(); i++) (*this)[i] = x;
        }
        void ones() {
            fill((T)1);
        }
        void minus_ones() {
            fill((T)-1.0);
        }

        T& cell(size_t i, size_t j) { return (*this)(i,j); }
        void swap_cells(size_t i1, size_t j1, size_t i2, size_t j2)
        {
            T x = (*this)(i1,j1);
            (*this)(i1,j1) = (*this)(i2,j2);
            (*this)(i2,j2) = x;
        }
            
        
        void random(T min = T(0), T max = T(1)) {
            Default r;
            for(size_t i = 0; i < size(); i++) (*this)[i] = r.random(min,max);
        }
        void identity() {
            size_t x = 0;
            zero();
            for(size_t i = 0; i < rows(); i++)
            {
                (*this)[i*N + x++] = 1;
            }
        }
        Vector<T>& get_vector() { return *this; }

        
        void print() const {
            std::cout << "Matrix[" << M << "," << N << "]=";
            for(size_t i = 0; i < rows(); i++) 
            {
                for(size_t j = 0; j < cols(); j++) std::cout << (*this)(i,j) << ",";
                std::cout << std::endl;
            }
        }

        // eigen compatibility
        void setZero() { zero(); }
        void setOnes() { ones(); }
        void setRandom() { random(); }
        Matrix<T> matrix() { return eval(); }        
        void setIdentity() { identity(); }
    };

    template<typename T>
    std::ostream& operator << (std::ostream & o, const Matrix<T> & m )
    {
        for(size_t i = 0; i < m.rows(); i++)
        {
            for(size_t j = 0; j < m.cols(); j++)
                o << m(i,j) << ",";
            o << std::endl;
        }
        return o;
    }

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
            DftiComputeForward(handle1, input.data(),x.data());
            memcpy(output.data(),x.data(), x.size()*sizeof(float));            
        }
        void Backward( Vector<std::complex<T>> & input, Vector<T> & output) {
            output.resize(size);
            Vector<float> x(size*2);            
            memcpy(x.data(),input.data(),x.size()*sizeof(float));
            DftiComputeBackward(handle1, x.data(), output.data());
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
            DftiComputeForward(handle1, input.data(),output.data());
        }
        void Backward( Vector<std::complex<T>> & input, Vector<std::complex<T>> & output) {
            output.resize(size);
            DftiComputeBackward(handle1, input.data(), output.data());
        }        
    };


    template<typename T>
    Matrix<T> matmul(Matrix<T> & a, Matrix<T> & b)
    {
        Matrix<T> r = a * b;
        return r;
    }

    
    template<typename T>
    Vector<T> sqr(Vector<T> & a) {
        Vector<T> r(a.size());                
        cppmkl::vsqr(a,r);
        return r;
    }
    template<typename T>
    Vector<T> abs(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vabs(a,r);
        return r;
    }
    template<typename T>
    Vector<T> inv(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vinv(a,r);
        return r;            
    }
    template<typename T>
    Vector<T> sqrt(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vsqrt(a,r);
        return r;
    }
    template<typename T>
    Vector<T> rsqrt(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vinvsqrt(a,r);
        return r;
    }
    template<typename T>
    Vector<T> cbrt(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vcbrt(a,r);
        return r;
    }
    template<typename T>
    Vector<T> rcbrt(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vinvcbrt(a,r);
        return r;
    }
    template<typename T>
    Vector<T> pow(const Vector<T> & a,const Vector<T> & b) {
        Vector<T> r(a.size());
        cppmkl::vpow(a,b,r);
        return r;
    }
    template<typename T>
    Vector<T> pow2o3(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vpow2o3(a,r);
        return r;
    }
    template<typename T>
    Vector<T> pow3o2(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vpow3o2(a,r);
        return r;
    }
    template<typename T>
    Vector<T> pow(const Vector<T> & a,const T b) {
        Vector<T> r(a.size());
        cppmkl::vpowx(a,b,r);
        return r;
    }
    template<typename T>
    Vector<T> hypot(const Vector<T> & a,const Vector<T> & b) {
        Vector<T> r(a.size());
        cppmkl::vhypot(a,b,r);
        return r;
    }
    template<typename T>
    Vector<T> exp(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vexp(a,r);
        return r;
    }
    template<typename T>
    Vector<T> exp2(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vexp2(a,r);
        return r;
    }
    template<typename T>
    Vector<T> exp10(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vexp10(a,r);
        return r;
    }
    template<typename T>
    Vector<T> expm1(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vexpm1(a,r);
        return r;
    }
    template<typename T>
    Vector<T> ln(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vln(a,r);
        return r;
    }
    template<typename T>
    Vector<T> log10(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vlog10(a,r);
        return r;
    }
    template<typename T>
    Vector<T> log2(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vlog2(a,r);
        return r;
    }
    template<typename T>
    Vector<T> logb(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vlogb(a,r);
        return r;
    }
    template<typename T>
    Vector<T> log1p(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vlog1p(a,r);
        return r;
    }
    template<typename T>
    Vector<T> cos(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vcos(a,r);
        return r;
    }
    template<typename T>
    Vector<T> sin(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vsin(a,r);
        return r;            
    }
    template<typename T>
    Vector<T> tan(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vtan(a,r);
        return r;
    }
    template<typename T>
    Vector<T> cosh(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vcosh(a,r);
        return r;
    }
    template<typename T>
    Vector<T> sinh(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vsinh(a,r);
        return r;            
    }
    template<typename T>
    Vector<T> tanh(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vtanh(a,r);
        return r;            
    }
    template<typename T>
    Vector<T> acos(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vacos(a,r);
        return r;
    }
    template<typename T>
    Vector<T> asin(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vasin(a,r);
        return r;            
    }
    template<typename T>
    Vector<T> atan(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vatan(a,r);
        return r;
    }
    template<typename T>
    Vector<T> atan2(const Vector<T> & a,const Vector<T> &n) {
        Vector<T> r(a.size());
        cppmkl::vatan2(a,n,r);
        return r;
    }
    template<typename T>
    Vector<T> acosh(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vacosh(a,r);
        return r;
    }
    template<typename T>
    Vector<T> asinh(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vasinh(a,r);
        return r;            
    }
    template<typename T>
    Vector<T> atanh(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vatanh(a,r);
        return r;
    }        
    template<typename T>
    void sincos(const Vector<T> & a, Vector<T> & b, Vector<T> & r) {        
        cppmkl::vsincos(a,b,r);        
    }
    template<typename T>
    Vector<T> erf(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::verf(a,r);
        return r;
    }
    template<typename T>
    Vector<T> erfinv(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::verfinv(a,r);
        return r;
    }
    template<typename T>
    Vector<T> erfc(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::verfc(a,r);
        return r;
    }
    template<typename T>
    Vector<T> cdfnorm(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vcdfnorm(a,r);
        return r;
    }    
    template<typename T>
    Vector<T> cdfnorminv(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vcdfnorm(a,r);
        return r;
    }    
    template<typename T>
    Vector<T> floor(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vfloor(a,r);
        return r;
    }    
    template<typename T>
    Vector<T> ceil(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vceil(a,r);
        return r;
    }    
    template<typename T>
    Vector<T> trunc(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vtrunc(a,r);
        return r;
    }    
    template<typename T>
    Vector<T> round(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vround(a,r);
        return r;        
    }    
    template<typename T>
    Vector<T> nearbyint(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vnearbyint(a,r);
        return r;
    }    
    template<typename T>
    Vector<T> rint(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vrint(a,r);
        return r;
    }    
    template<typename T>
    Vector<T> fmod(const Vector<T> & a, Vector<T> & b) {
        Vector<T> r(a.size());
        cppmkl::vmodf(a,b,r);
        return r;
    }    
    
    template<typename T>
    Vector<T> CIS(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vCIS(a,r);
        return r;
    }    
    template<typename T>
    Vector<T> cospi(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vcospi(a,r);
        return r;
    }    
    template<typename T>
    Vector<T> sinpi(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vsinpi(a,r);
        return r;
    }    
    template<typename T>
    Vector<T> tanpi(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vtanpi(a,r);
        return r;
    }    
    template<typename T>
    Vector<T> acospi(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vacospi(a,r);
        return r;
    }    
    template<typename T>
    Vector<T> asinpi(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vasinpi(a,r);
        return r;
    }    
    template<typename T>
    Vector<T> atanpi(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vatanpi(a,r);
        return r;
    }    
    template<typename T>
    Vector<T> atan2pi(const Vector<T> & a, Vector<T> & b) {
        Vector<T> r(a.size());
        cppmkl::vatan2pi(a,b,r);
        return r;
    }    
    template<typename T>
    Vector<T> cosd(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vcosd(a.size(),a.data(),r.data());
        return r;
    }    
    template<typename T>
    Vector<T> sind(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vsind(a.size(),a.data(),r.data());
        return r;
    }    
    template<typename T>
    Vector<T> tand(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vtand(a.size(),a.data(),r.data());
        return r;
    }       
    template<typename T>
    Vector<T> lgamma(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vlgamma(a,r);
        return r;
    }       
    template<typename T>
    Vector<T> tgamma(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vtgamma(a,r);
        return r;
    }       
    template<typename T>
    Vector<T> expint1(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vexpint1(a,r);
        return r;
    }       
    template<typename T>
    Matrix<T> sqr(Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vsqr(a,r);
        return r;
    }
    template<typename T>
    Matrix<T> abs(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vabs(a,r);
        return r;
    }
    template<typename T>
    Matrix<T> inv(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vinv(a,r);
        return r;            
    }
    template<typename T>
    Matrix<T> sqrt(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vsqrt(a,r);
        return r;
    }
    template<typename T>
    Matrix<T> rsqrt(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vinvsqrt(a,r);
        return r;
    }
    template<typename T>
    Matrix<T> cbrt(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vcbrt(a,r);
        return r;
    }
    template<typename T>
    Matrix<T> rcbrt(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vinvcbrt(a,r);
        return r;
    }
    template<typename T>
    Matrix<T> pow(const Matrix<T> & a,const Matrix<T> & b) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vpow(a,b,r);
        return r;
    }
    template<typename T>
    Matrix<T> pow2o3(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vpow2o3(a,r);
        return r;
    }
    template<typename T>
    Matrix<T> pow3o2(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vpow3o2(a,r);
        return r;
    }
    template<typename T>
    Matrix<T> pow(const Matrix<T> & a,const T b) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vpowx(a,b,r);
        return r;
    }
    template<typename T>
    Matrix<T> hypot(const Matrix<T> & a,const Matrix<T> & b) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vhypot(a,b,r);
        return r;
    }
    template<typename T>
    Matrix<T> exp(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vexp(a,r);
        return r;
    }
    template<typename T>
    Matrix<T> exp2(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vexp2(a,r);
        return r;
    }
    template<typename T>
    Matrix<T> exp10(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vexp10(a,r);
        return r;
    }
    template<typename T>
    Matrix<T> expm1(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vexpm1(a,r);
        return r;
    }
    template<typename T>
    Matrix<T> ln(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vln(a,r);
        return r;
    }
    template<typename T>
    Matrix<T> log10(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vlog10(a,r);
        return r;
    }
    template<typename T>
    Matrix<T> log2(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vlog2(a,r);        
        return r;
    }
    template<typename T>
    Matrix<T> logb(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vlogb(a,r);
        return r;
    }
    template<typename T>
    Matrix<T> log1p(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vlog1p(a,r);
        return r;
    }
    template<typename T>
    Matrix<T> cos(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vcos(a,r);
        return r;
    }
    template<typename T>
    Matrix<T> sin(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vsin(a,r);
        return r;            
    }
    template<typename T>
    Matrix<T> tan(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vtan(a,r);
        return r;
    }
    template<typename T>
    Matrix<T> cosh(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vcosh(a,r);
        return r;
    }
    template<typename T>
    Matrix<T> sinh(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vsinh(a,r);
        return r;            
    }
    template<typename T>
    Matrix<T> tanh(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vtanh(a,r);
        return r;            
    }
    template<typename T>
    Matrix<T> acos(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vacos(a,r);
        return r;
    }
    template<typename T>
    Matrix<T> asin(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vasin(a,r);
        return r;            
    }
    template<typename T>
    Matrix<T> atan(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vatan(a,r);
        return r;
    }
    template<typename T>
    Matrix<T> atan2(const Matrix<T> & a,const Matrix<T> &n) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vatan2(a,n,r);
        return r;
    }
    template<typename T>
    Matrix<T> acosh(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vacosh(a,r);
        return r;
    }
    template<typename T>
    Matrix<T> asinh(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vasinh(a,r);
        return r;            
    }
    template<typename T>
    Matrix<T> atanh(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vatanh(a,r);
        return r;
    }        
    template<typename T>
    void sincos(const Matrix<T> & a, Matrix<T> & b, Matrix<T> & r) {        
        cppmkl::vsincos(a,b,r);
    }
    template<typename T>
    Matrix<T> erf(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::verf(a,r);
        return r;
    }
    template<typename T>
    Matrix<T> erfinv(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::verfinv(a,r);
        return r;
    }
    template<typename T>
    Matrix<T> erfc(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::verfc(a,r);
        return r;
    }
    template<typename T>
    Matrix<T> cdfnorm(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vcdfnorm(a,r);
        return r;
    }    
    template<typename T>
    Matrix<T> cdfnorminv(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vcdfnorm(a,r);
        return r;
    }    
    template<typename T>
    Matrix<T> floor(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vfloor(a,r);
        return r;
    }    
    template<typename T>
    Matrix<T> ceil(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vceil(a,r);
        return r;
    }    
    template<typename T>
    Matrix<T> trunc(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vtrunc(a,r);
        return r;
    }    
    template<typename T>
    Matrix<T> round(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vround(a,r);
        return r;        
    }    
    template<typename T>
    Matrix<T> nearbyint(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vnearbyint(a,r);
        return r;
    }    
    template<typename T>
    Matrix<T> rint(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vrint(a,r);
        return r;
    }    
    template<typename T>
    Matrix<T> fmod(const Matrix<T> & a, Matrix<T> & b) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vmodf(a,b,r);
        return r;
    }    
    template<typename T>
    Matrix<T> mulbyconj(const Matrix<T> & a, const Matrix<T> & b) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vmulbyconj(a,b,r);
        return r;
    }    
    template<typename T>
    Matrix<T> conj(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vconj(a,r);
        return r;
    }    
    template<typename T>
    Matrix<T> arg(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::varg(a,r);
        return r;
    }    
    template<typename T>
    Matrix<T> CIS(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vCIS(a,r);
        return r;
    }    
    template<typename T>
    Matrix<T> cospi(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vcospi(a,r);
        return r;
    }    
    template<typename T>
    Matrix<T> sinpi(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vsinpi(a,r);
        return r;
    }    
    template<typename T>
    Matrix<T> tanpi(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vtanpi(a,r);
        return r;
    }    
    template<typename T>
    Matrix<T> acospi(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vacospi(a,r);
        return r;
    }    
    template<typename T>
    Matrix<T> asinpi(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vasinpi(a,r);
        return r;
    }    
    template<typename T>
    Matrix<T> atanpi(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vatanpi(a,r);
        return r;
    }    
    template<typename T>
    Matrix<T> atan2pi(const Matrix<T> & a, Matrix<T> & b) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vatan2pi(a,b,r);
        return r;
    }    
    template<typename T>
    Matrix<T> cosd(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vcosd(a.size(),a.data(),r.data());
        return r;
    }    
    template<typename T>
    Matrix<T> sind(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vsind(a.size(),a.data(),r.data());
        return r;
    }    
    template<typename T>
    Matrix<T> tand(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vtand(a.size(),a.data(),r.data());
        return r;
    }       
    template<typename T>
    Matrix<T> lgamma(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vlgamma(a,r);
        return r;
    }       
    template<typename T>
    Matrix<T> tgamma(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vtgamma(a,r);
        return r;
    }       
    template<typename T>
    Matrix<T> expint1(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vexpint1(a,r);
        return r;
    }       

    template<typename T>
    Vector<T> copy(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::cblas_copy(a.size(),a.data(),1,r.data(),1);
        return r;
    }       

    template<typename T> T sum(const Vector<T> & a) {        
        return cppmkl::cblas_asum(a.size(), a.data(),1);        
    }       

    template<typename T>
    Vector<T> add(const Vector<T> & a, const Vector<T> & b) {
        Vector<T> r(b);
        cppmkl::cblas_axpy(a.size(),1.0,a.data(),1,r.data(),1);
        return r;
    }       
    template<typename T>
    Vector<T> sub(const Vector<T> & a, const Vector<T> & b) {
        Vector<T> r(a);
        cppmkl::cblas_axpy(a.size(),-1.0,b.data(),1,r.data(),1);
        return r;
    }       
    template<typename T>
    T dot(const Vector<T> & a, Vector<T> & b) {
        return cppmkl::cblas_dot(a,b);
    }       
    template<typename T>
    T nrm2(const Vector<T> & a) {
        Vector<T> r(a);
        return cppmkl::cblas_nrm2(a);        
    }       
    
    template<typename T>
    void scale(Vector<T> & x, T alpha) {
        cppmkl::cblas_scal(x.size(),alpha,x.data(),1);
    }


    template<typename T>
    size_t min_index(const Matrix<T> & m) { return cppmkl::cblas_iamin(m.size(),m.data(),1); }
    template<typename T>
    size_t max_index(const Matrix<T> & m) { return cppmkl::cblas_iamax(m.size(),m.data(),1); }
    
    template<typename T>
    size_t min_index(const Vector<T> & v) { return cppmkl::cblas_iamin(v.size(),v.data(),1); }
    template<typename T>
    size_t max_index(const Vector<T> & v) { return cppmkl::cblas_iamax(v.size(),v.data(),1); }

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

    template<typename T>
    Vector<T> operator * (T a, const Vector<T> & b) {
        Vector<T> x(b.size());
        x.fill(a);
        Vector<T> r(b.size());
        cppmkl::vmul(x,b,r);            
        return r;
    }
    template<typename T>
    Vector<T> operator / (T a, const Vector<T> & b) {
        Vector<T> x(b.size());
        x.fill(a);
        Vector<T> r(b.size());
        cppmkl::vdiv(x,b,r);            
        return r;
    }
    template<typename T>
    Vector<T> operator + (T a, const Vector<T> & b) {
        Vector<T> x(b.size());
        x.fill(a);
        Vector<T> r(b.size());
        cppmkl::vadd(x,b,r);            
        return r;
    }
    template<typename T>
    Vector<T> operator - (T a, const Vector<T> & b) {
        Vector<T> x(b.size());
        x.fill(a);
        Vector<T> r(b.size());
        cppmkl::vsub(x,b,r);            
        return r;
    }

    template<typename T>
    Vector<T> operator * (const Vector<T> & a, T b) {
        Vector<T> x(a.size());
        x.fill(b);
        Vector<T> r(a.size());
        cppmkl::vmul(a,x,r);            
        return r;
    }
    template<typename T>
    Vector<T> operator / (const Vector<T> & a , T b) {
        Vector<T> x(a.size());
        x.fill(b);
        Vector<T> r(a.size());
        cppmkl::vdiv(a,x,r);            
        return r;
    }
    template<typename T>
    Vector<T> operator + (const Vector<T> & a, T b) {
        Vector<T> x(a.size());
        x.fill(b);
        Vector<T> r(a.size());
        cppmkl::vadd(a,x,r);            
        return r;
    }
    template<typename T>
    Vector<T> operator - (const Vector<T> & a, T b) {
        Vector<T> x(a.size());
        x.fill(b);
        Vector<T> r(a.size());
        cppmkl::vsub(a,x,r);            
        return r;
    }
    
    template<typename T>
    Vector<T> operator - (const Vector<T> & a, const Vector<T> & b) {        
        Vector<T> r(a.size());
        cppmkl::vsub(a,b,r);            
        return r;
    }
    template<typename T>
    Vector<T> operator + (const Vector<T> & a, const Vector<T> & b) {        
        Vector<T> r(a.size());
        cppmkl::vadd(a,b,r);                    
        return r;
    }
    template<typename T>
    Vector<T> operator * (const Vector<T> & a, const Vector<T> & b) {        
        Vector<T> r(a.size());
        cppmkl::vmul(a,b,r);            
        return r;
    }
    template<typename T>
    Vector<T> operator / (const Vector<T> & a, const Vector<T> & b) {        
        Vector<T> r(a.size());
        cppmkl::vdiv(a,b,r);            
        return r;
    }

    template<typename T>
    Matrix<T> operator * (const Matrix<T> & a, const Matrix<T> & b) {
        assert(a.N == b.M);
        Matrix<T> r(a.rows(),b.cols());
        r.zero();
                        
        int m = a.rows();
        int n = b.cols();
        int k = a.cols();       
        
        cppmkl::cblas_gemm(CblasRowMajor,
                    CblasNoTrans,
                    CblasNoTrans,
                    m,
                    n,
                    k,
                    T(1.0),
                    a.data(),
                    k,
                    b.data(),
                    n,
                    T(0.0),
                    r.data(),
                    n);
        /*
        for(size_t i = 0; i < a.rows(); i++)
            for(size_t j = 0; j < b.cols(); j++)
            {
                float sum = 0;
                for(size_t k = 0; k < b.M; k++) sum += a(i,k) * b(k,j);
                r(i,j) = sum;
            }
        */
        return r;
    } 
    
    template<typename T>
    Matrix<T> operator + (const Matrix<T> & a, const Matrix<T> & b) {
        assert(a.rows() == b.rows() && a.cols() == b.cols());
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vadd(a,b,r);            
        return r;
    }
    template<typename T>
    Matrix<T> operator / (const Matrix<T> & a, const Matrix<T> & b) {
        assert(a.rows() == b.rows() && a.cols() == b.cols());
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vdiv(a,b,r);            
        return r;
    }        
    template<typename T>
    Matrix<T> operator - (const Matrix<T> & a, const Matrix<T> & b) {
        assert(a.rows() == b.rows() && a.cols() == b.cols());
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vsub(a,b,r);            
        return r;
    }        

    template<typename T>
    Matrix<T> operator * (T a, const Matrix<T> & b) {
        Matrix<T> x(b.M,b.N);
        x.fill(a);
        Matrix<T> r(b.M,b.N);
        cppmkl::vmul(x,b,r);
        return r;
    }
    
    template<typename T>
    Matrix<T> operator + (T a, const Matrix<T> & b) {
        Matrix<T> x(b.M,b.N);
        x.fill(a);
        Matrix<T> r(b.M,b.N);
        cppmkl::vadd(x,b,r);
        return r;
    }
    template<typename T>
    Matrix<T> operator - (T a, const Matrix<T> & b) {
        Matrix<T> x(b.M,b.N);
        x.fill(a);
        Matrix<T> r(b.M,b.N);
        cppmkl::vsub(x,b,r);
        return r;
    }
    template<typename T>
    Matrix<T> operator / (T a, const Matrix<T> & b) {
        Matrix<T> x(b.M,b.N);
        x.fill(a);
        Matrix<T> r(b.M,b.N);
        cppmkl::vdiv(x,b,r);
        return r;
    }
    template<typename T>
    Matrix<T> operator * (const Matrix<T> & a, T b) {
        Matrix<T> x(a.M,a.N);
        x.fill(b);
        Matrix<T> r(a.M,a.N);
        cppmkl::vmul(a,x,r);
        return r;
    }
    
    template<typename T>
    Matrix<T> operator + (const Matrix<T> & a, T b) {
        Matrix<T> x(a.M,a.N);
        x.fill(b);
        Matrix<T> r(a.M,a.N);
        cppmkl::vadd(a,x,r);
        return r;
    }
    template<typename T>
    Matrix<T> operator - (const Matrix<T> & a, T b) {
        Matrix<T> x(a.M,a.N);
        x.fill(b);
        Matrix<T> r(a.M,a.N);
        cppmkl::vsub(a,x,r);
        return r;
    }
    template<typename T>
    Matrix<T> operator / (const Matrix<T> & a, T b) {
        Matrix<T> x(a.M,a.N);
        x.fill(b);
        Matrix<T> r(a.M,a.N);
        cppmkl::vdiv(a,x,r);
        return r;
    }    
    template<typename T>
    Matrix<T> hadamard(const Matrix<T> & a, const Matrix<T> & b) {
        Matrix<T> r(a.rows(),a.cols());
        assert(a.rows() == b.rows() && a.cols() == b.cols());
        cppmkl::vmul(a,b,r);
        return r;
    }      
    template<typename T>
    Vector<T> operator *(const Vector<T> &a, const Matrix<T> &b) {
        Vector<T> r(a.size());        
        sgemv(CblasRowMajor,CblasNoTrans,b.rows(),b.cols(),1.0,b.data(),b.cols(),a.data(),1,1.0,r.data(),1);
        return r;
    }      
    template<typename T>
    Vector<T> operator *(const Matrix<T> &a, const Vector<T> &b) {
        Vector<T> r(b.size());        
        sgemv(CblasRowMajor,CblasNoTrans,a.rows(),a.cols(),1.0,a.data(),a.cols(),b.data(),1,1.0,r.data(),1);
        return r;
    }      
}