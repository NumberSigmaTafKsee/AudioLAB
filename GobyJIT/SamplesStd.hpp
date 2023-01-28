#pragma once

#include <vector>
#include <algorithm>
#include <cmath>
#include <algorithm>
#include <functional>
#include <iostream>
#include <complex>
#include <vector>
#include <new>
#include <chrono>
#include <random>
#include <cassert>

#include "cppmkl/cppmkl_allocator.h"
#include "cppmkl/cppmkl_vml.h"
#include "cppmkl/cppmkl_cblas.h"
#include "StdNoise.hpp"

namespace AudioDSP
{
    template<typename T>
    struct sample_vector : public std::vector<T, cppmkl::cppmkl_allocator<T>>
    {
        using base = std::vector<T,cppmkl::cppmkl_allocator<T>>;    
        using base::operator [];
        using base::fill;
        using base::pop_back;
        using base::push_front;
        using base::front;
        using base::back;
        using base::begin;
        using base::end;
        using base::size;
        using base::resize;
        using base::at;
        using base::data;
        using base::assign;
        using base::insert;        
        using base::erase;
        using base::swap;
        using base::clear;
        using base::emplace;
        using base::emplace_back;
        using base::get_allocator;

        sample_vector() = default;
        sample_vector(size_t n) : base(n) {}
        sample_vector(const base& b) : base(b) { }
        sample_vector(const sample_vector<T> & v) { *this = v; }

        sample_vector<T> operator + (const T x) {
            sample_vector<T> r(size());
            #pragma omp simd
            for(size_t i = 0; i < size(); i++) r[i] = (*this)[i] + x;
            return *this;
        }
        sample_vector<T> operator - (const T x) {
            sample_vector<T> r(size());
            #pragma omp simd
            for(size_t i = 0; i < size(); i++) r[i] = (*this)[i] - x;
            return *this;
        }
        sample_vector<T> operator * (const T x) {
            sample_vector<T> r(size());
            #pragma omp simd
            for(size_t i = 0; i < size(); i++) r[i] = (*this)[i] * x;
            return *this;
        }
        sample_vector<T> operator / (const T x) {
            sample_vector<T> r(size());
            #pragma omp simd
            for(size_t i = 0; i < size(); i++) r[i] = (*this)[i] / x;
            return *this;
        }
        sample_vector<T> operator % (const T x) {
            sample_vector<T> r(size());
            #pragma omp simd
            for(size_t i = 0; i < size(); i++) r[i] = std::fmod(*this[i],x);
            return *this;
        }
        sample_vector<T>& operator += (const T x) {
            #pragma omp simd
            for(size_t i = 0; i < size(); i++) (*this)[i] += x;
            return *this;
        }
        sample_vector<T>& operator -= (const T x) {
            #pragma omp simd
            for(size_t i = 0; i < size(); i++) (*this)[i] -= x;
            return *this;
        }
        sample_vector<T>& operator *= (const T x) {
            #pragma omp simd
            for(size_t i = 0; i < size(); i++) (*this)[i] *= x;
            return *this;
        }
        sample_vector<T>& operator /= (const T x) {
            #pragma omp simd
            for(size_t i = 0; i < size(); i++) (*this)[i] /= x;
            return *this;
        }
        sample_vector<T>& operator %= (const T x) {
            #pragma omp simd
            for(size_t i = 0; i < size(); i++) (*this)[i] = std::fmod(*this[i],x);
            return *this;
        }

        T min() { return *std::min_element(begin(),end()); }
        T max() { return *std::max_element(begin(),end()); }

        size_t min_index() { return std::distance(begin(), std::min_element(begin(),end())); }
        size_t max_index() { return std::distance(begin(), std::max_element(begin(),end())); }
        
        sample_vector<T>& operator +=  (const sample_vector<T> & v) { 
            cppmkl::vadd(*this,v,*this);
            return *this;
        }
        sample_vector<T>& operator -=  (const sample_vector<T> & v) { 
            cppmkl::vsub(*this,v,*this);
            return *this;
        }
        sample_vector<T>& operator *=  (const sample_vector<T> & v) { 
            cppmkl::vmul(*this,v,*this);
            return *this;
        }
        sample_vector<T>& operator /=  (const sample_vector<T> & v) { 
            cppmkl::vdiv(*this,v,*this);
            return *this;        
	    }
        sample_vector<T> operator + (const sample_vector<T> & v) { 
            sample_vector<T> r(size());            
            cppmkl::vadd(*this,v,r);
            return r;
        }
        sample_vector<T> operator - (const sample_vector<T> & v) { 
            sample_vector<T> r(size());            
            cppmkl::vsub(*this,v,r);
            return r;            
        }
        sample_vector<T> operator * (const sample_vector<T> & v) { 
            sample_vector<T> r(size());            
            cppmkl::vmul(*this,v,r);
            return r;
        }
        sample_vector<T> operator / (const sample_vector<T> & v) { 
            sample_vector<T> r(size());            
            cppmkl::vdiv(*this,v,r);
            return r;
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

        sample_vector<T> eval() {
            return sample_vector<T>(*this); 
        }
        sample_vector<T> slice(size_t start, size_t len) {
            sample_vector<T> x(len);
            memcpy(x.data(),data()+start,len*sizeof(T));
            return x;
        }

        void print() {
            std::cout << "sample_vector[" << size() << "]=";
            for(size_t i = 0; i < size()-1; i++) std::cout << (*this)[i] << ",";
            std::cout << (*this)[size()-1] << std::endl;
        }
        
        // eigen compatibility
        void setZero() { zero(); }
        void setOnes() { ones(); }
        void setRandom() { random(); }
        
    };

    template<typename T> struct sample_matrix;
    template<typename T>
    struct MatrixView 
    {
        sample_matrix<T> * matrix;
        size_t row;
        MatrixView(sample_matrix<T> * m, size_t r) {
            matrix = m;
            row = r;
        }

        T& operator[](size_t i);
        T  __getitem(size_t i);
        void __setitem(size_t i, T x);
    };
    template<typename T>
    struct sample_matrix : public sample_vector<T>
    {
        size_t M;
        size_t N;

        sample_matrix() { M = N = 0; }
        sample_matrix(size_t m, size_t n) {
            resize(m,n);
            assert(M > 0);
            assert(N > 0);
        }
        sample_matrix(const sample_matrix<T> & m) {
            
            resize(m.M,m.N);
            memcpy(data(),m.data(),size()*sizeof(T));            
        }
        sample_matrix(T * ptr, size_t m, size_t n)
        {
            resize(m,n);            
            memcpy(data(),ptr,m*n*sizeof(T));
        }
        sample_matrix<T>& operator = (const T v) {
            fill(v);
            return *this;
        }

        sample_matrix<T>& operator = (const sample_matrix<T> & m) {            
            resize(m.M,m.N);
            memcpy(data(),m.data(),size()*sizeof(T));            
            return *this;
        }
        sample_matrix<T>& operator = (const sample_vector<T> & m) {            
            resize(1,m.size());
            memcpy(data(),m.data(),size()*sizeof(T));            
            return *this;
        }

        using sample_vector<T>::size;
        using sample_vector<T>::data;
        using sample_vector<T>::resize;
        using sample_vector<T>::at;
        using sample_vector<T>::operator [];

        size_t rows() const { return M; }
        size_t cols() const { return N; }
        
        sample_matrix<T> cwiseMax(T v) {
            sample_matrix<T> r(*this);
            for(size_t i = 0; i < rows(); i++)
                for(size_t j = 0;  j < cols(); j++)
                    if(r(i,j) < v) r(i,j) = v;
            return r;
        }
        sample_matrix<T> row(size_t m) { 
            sample_matrix<T> r(1,cols());
            for(size_t i = 0; i < cols(); i++) r(0,i) = (*this)(m,i);
            return r;
        }
        sample_matrix<T> col(size_t c) { 
            sample_matrix<T> r(cols(),1);
            for(size_t i = 0; i < rows(); i++) r(i,0) = (*this)(i,c);
            return r;
        }
        void row(size_t m, const sample_vector<T> & v)
        {
            for(size_t i = 0; i < cols(); i++) (*this)(m,i) = v[i];
        }
        void col(size_t n, const sample_vector<T> & v)
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

        sample_matrix<T> operator - () {
            sample_matrix<T> r(*this);
            return T(-1.0)*r;
        }
        sample_matrix<T> addToEachRow(sample_vector<T> & v) {
            sample_matrix<T> r(*this);
            for(size_t i = 0; i < M; i++)
            {
                for(size_t j = 0; j < N; j++)
                {
                    r(i,j) += v[j];
                }
            }
            return r;
        }
        sample_matrix<T> eval() {
            sample_matrix<T> r(*this);
            return r;
        }
        sample_matrix<T>& operator += (const T b)
        {
            sample_matrix<T> x(rows(),cols());
            x.fill(b);
            cppmkl::vadd(*this,x,*this);
            return *this;
        }
        sample_matrix<T>& operator -= (const T b)
        {
            sample_matrix<T> x(rows(),cols());
            x.fill(b);
            cppmkl::vsub(*this,x,*this);
            return *this;
        }
        sample_matrix<T>& operator *= (const T b)
        {
            sample_matrix<T> x(rows(),cols());
            x.fill(b);
            cppmkl::vmul(*this,x,*this);
            return *this;
        }
        sample_matrix<T>& operator /= (const T b)
        {
            sample_matrix<T> x(rows(),cols());
            x.fill(b);
            cppmkl::vdiv(*this,x,*this);
            return *this;
        }
        
        sample_matrix<T>& operator += (const sample_matrix<T> b)
        {
            assert(rows() == b.rows() && cols() == b.cols());            
            cppmkl::vadd(*this,b,*this);
            return *this;
        }
        sample_matrix<T>& operator -= (const sample_matrix<T> b)
        {
            assert(rows() == b.rows() && cols() == b.cols());            ;
            cppmkl::vsub(*this,b,*this);
            return *this;
        }
        sample_matrix<T>& operator *= (const sample_matrix<T> b)
        {
            assert(rows() == b.rows() && cols() == b.cols());            
            cppmkl::vmul(*this,b,*this);
            return *this;
        }
        sample_matrix<T>& operator /= (const sample_matrix<T> b)
        {
            assert(rows() == b.rows() && cols() == b.cols());            
            cppmkl::vdiv(*this,b,*this);
            return *this;
        }

        sample_matrix<T> hadamard(const sample_matrix<T> & m) {
            sample_matrix<T> r(rows(),cols());
            assert(rows() == m.rows() && cols() == m.cols());
            cppmkl::vmul(*this,m,r);
            return r;
        }
        
        sample_matrix<T> transpose() {
            sample_matrix<T> r(cols(),rows());
            for(size_t i = 0; i < rows(); i++)
                for(size_t j = 0; j < cols(); j++)
                    r(j,i) =(*this)[i*N + j];
            return r;
        }
        sample_matrix<T> t() {
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
        sample_vector<T>& get_vector() { return *this; }

        
        void print() const {
            std::cout << "sample_matrix[" << M << "," << N << "]=";
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
        sample_matrix<T> matrix() { return eval(); }        
        void setIdentity() { identity(); }
    };


    template<typename T>
    using complex_vector = sample_vector<std::complex<T>>;

    template<typename T>
    using complex_matrix = sample_matrix<std::complex<T>>;    


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

        void Forward( sample_vector<T> & input, sample_vector<std::complex<T>> & output) {
            output.resize(size);
            sample_vector<float> x(size*2);            
            DftiComputeForward(handle1, input.data(),x.data());
            memcpy(output.data(),x.data(), x.size()*sizeof(float));            
        }
        void Backward( sample_vector<std::complex<T>> & input, sample_vector<T> & output) {
            output.resize(size);
            sample_vector<float> x(size*2);            
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

        void Forward( sample_vector<std::complex<T>> & input, sample_vector<std::complex<T>> & output) {
            output.resize(size);
            DftiComputeForward(handle1, input.data(),output.data());
        }
        void Backward( sample_vector<std::complex<T>> & input, sample_vector<std::complex<T>> & output) {
            output.resize(size);
            DftiComputeBackward(handle1, input.data(), output.data());
        }        
    };
 
    template<typename T>
    std::ostream& operator << (std::ostream & o, const sample_matrix<T> & m )
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
    T& MatrixView<T>::operator[](size_t i) {
        return matrix->matrix[row*matrix->N + i];
    }

    template<typename T>    
    T  MatrixView<T>::__getitem(size_t i) {
        return matrix->matrix[row*matrix->N + i];
    }

    template<typename T>    
    void MatrixView<T>::__setitem(size_t i, T v)
    {
        matrix->matrix[row*matrix->N + i] = v;
    }

    template<typename T>
    sample_matrix<T> matmul(sample_matrix<T> & a, sample_matrix<T> & b)
    {
        sample_matrix<T> r = a * b;
        return r;
    }

#ifdef NO_VML
    template<class T>
    sample_vector<T> cos(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::cos(v[i]);
        return r;
    }
    template<class T>
    sample_vector<T> sin(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::sin(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> tan(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::tan(v[i]);
        return r;
    }

    template<class T>
    sample_vector<T> acos(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::acos(v[i]);
        return r;
    }
    template<class T>
    sample_vector<T> asin(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::asin(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> atan(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::atan(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> atan2(const sample_vector<T> & v, const T value) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::atan2(v[i], value);
        return r;
    }    
    template<class T>
    sample_vector<T> atan2(const sample_vector<T> & v, const sample_vector<T> value) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::atan2(v[i], value[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> cosh(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::cosh(v[i]);
        return r;
    }
    template<class T>
    sample_vector<T> sinh(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::sinh(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> tanh(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::tanh(v[i]);
        return r;
    }

    template<class T>
    sample_vector<T> acosh(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::acosh(v[i]);
        return r;
    }
    template<class T>
    sample_vector<T> asinh(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::asinh(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> atanh(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::atanh(v[i]);
        return r;
    }    

    template<class T>
    sample_vector<T> exp(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::exp(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> log(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::log(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> log10(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::log10(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> exp2(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::exp2(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> expm1(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::expm1(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> ilogb(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::ilogb(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> log2(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::log2(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> log1p(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::log1p(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> logb(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::logb(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> scalbn(const sample_vector<T> & v, const sample_vector<int> & x) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::scalbn(v[i],x[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> scalbn(const sample_vector<T> & v, const int x) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::scalbn(v[i],x);
        return r;
    }    
    template<class T>
    sample_vector<T> scalbln(const sample_vector<T> & v, const sample_vector<long int> & x) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::scalbln(v[i],x[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> scalbln(const sample_vector<T> & v, const long int x) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::scalbln(v[i],x);
        return r;
    }    
    template<class T>
    sample_vector<T> pow(const sample_vector<T> & v, const sample_vector<T> & x) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::pow(v[i],x[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> pow(const sample_vector<T> & v, const T x) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::pow(v[i],x);
        return r;
    }    
    template<class T>
    sample_vector<T> pow(const T x, const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::pow(x,v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> sqrt(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::sqrt(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> cbrt(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::cbrt(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> hypot(const sample_vector<T> & v, const sample_vector<T> & x) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::hypot(v[i],x[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> hypot(const sample_vector<T> & v, const T x) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::hypot(v[i],x);
        return r;
    }    
    template<class T>
    sample_vector<T> hypot(const T x, const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::hypot(x,v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> erf(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::erf(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> erfc(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::erfc(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> tgamma(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::tgamma(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> lgamma(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::lgamma(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> ceil(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::ceil(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> floor(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::floor(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> fmod(const sample_vector<T> & v, const sample_vector<T> & x) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::fmod(v[i],x[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> fmod(const sample_vector<T> & v, const T x) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::fmod(v[i],x);
        return r;
    }    
    template<class T>
    sample_vector<T> fmod(const T x, const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::fmod(x,v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> trunc(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::trunc(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> round(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::round(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<long int> lround(const sample_vector<T> & v) {
        sample_vector<long int> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::lround(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<long long int> llround(const sample_vector<T> & v) {
        sample_vector<long long int> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::llround(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> nearbyint(const sample_vector<T> & v) {
        sample_vector<long long int> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::nearbyint(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> remainder(const sample_vector<T> & a, const sample_vector<T> & b) {
        sample_vector<long long int> r(a.size());
        #pragma omp simd
        for(size_t i = 0; i < a.size(); i++) r[i] = std::remainder(a[i],b[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> copysign(const sample_vector<T> & a, const sample_vector<T> & b) {
        sample_vector<long long int> r(a.size());
        #pragma omp simd
        for(size_t i = 0; i < a.size(); i++) r[i] = std::copysign(a[i],b[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> fdim(const sample_vector<T> & a, const sample_vector<T> & b) {
        sample_vector<long long int> r(a.size());
        #pragma omp simd
        for(size_t i = 0; i < a.size(); i++) r[i] = std::fdim(a[i],b[i]);
        return r;
    }    
    #undef fmax
    template<class T>
    sample_vector<T> fmax(const sample_vector<T> & a, const sample_vector<T> & b) {
        sample_vector<long long int> r(a.size());
        #pragma omp simd
        for(size_t i = 0; i < a.size(); i++) r[i] = std::fmax(a[i],b[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> fmin(const sample_vector<T> & a, const sample_vector<T> & b) {
        sample_vector<long long int> r(a.size());
        #pragma omp simd
        for(size_t i = 0; i < a.size(); i++) r[i] = std::fmin(a[i],b[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> fma(const sample_vector<T> & a, const sample_vector<T> & b, const sample_vector<T> & c) {
        sample_vector<long long int> r(a.size());
        #pragma omp simd
        for(size_t i = 0; i < a.size(); i++) r[i] = std::fma(a[i],b[i],c[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> fabs(const sample_vector<T> & a) {
        sample_vector<T> r(a.size());
        #pragma omp simd
        for(size_t i = 0; i < a.size(); i++) r[i] = std::fabs(a[i]);
        return r;
    }    
#else


template<typename T>
    sample_matrix<T> matmul(sample_matrix<T> & a, sample_matrix<T> & b)
    {
        sample_matrix<T> r = a * b;
        return r;
    }

    
    template<typename T>
    sample_vector<T> sqr(sample_vector<T> & a) {
        sample_vector<T> r(a.size());                
        cppmkl::vsqr(a,r);
        return r;
    }
    template<typename T>
    sample_vector<T> abs(const sample_vector<T> & a) {
        sample_vector<T> r(a.size());
        cppmkl::vabs(a,r);
        return r;
    }
    template<typename T>
    sample_vector<T> inv(const sample_vector<T> & a) {
        sample_vector<T> r(a.size());
        cppmkl::vinv(a,r);
        return r;            
    }
    template<typename T>
    sample_vector<T> sqrt(const sample_vector<T> & a) {
        sample_vector<T> r(a.size());
        cppmkl::vsqrt(a,r);
        return r;
    }
    template<typename T>
    sample_vector<T> rsqrt(const sample_vector<T> & a) {
        sample_vector<T> r(a.size());
        cppmkl::vinvsqrt(a,r);
        return r;
    }
    template<typename T>
    sample_vector<T> cbrt(const sample_vector<T> & a) {
        sample_vector<T> r(a.size());
        cppmkl::vcbrt(a,r);
        return r;
    }
    template<typename T>
    sample_vector<T> rcbrt(const sample_vector<T> & a) {
        sample_vector<T> r(a.size());
        cppmkl::vinvcbrt(a,r);
        return r;
    }
    template<typename T>
    sample_vector<T> pow(const sample_vector<T> & a,const sample_vector<T> & b) {
        sample_vector<T> r(a.size());
        cppmkl::vpow(a,b,r);
        return r;
    }
    template<typename T>
    sample_vector<T> pow2o3(const sample_vector<T> & a) {
        sample_vector<T> r(a.size());
        cppmkl::vpow2o3(a,r);
        return r;
    }
    template<typename T>
    sample_vector<T> pow3o2(const sample_vector<T> & a) {
        sample_vector<T> r(a.size());
        cppmkl::vpow3o2(a,r);
        return r;
    }
    template<typename T>
    sample_vector<T> pow(const sample_vector<T> & a,const T b) {
        sample_vector<T> r(a.size());
        cppmkl::vpowx(a,b,r);
        return r;
    }
    template<typename T>
    sample_vector<T> hypot(const sample_vector<T> & a,const sample_vector<T> & b) {
        sample_vector<T> r(a.size());
        cppmkl::vhypot(a,b,r);
        return r;
    }
    template<typename T>
    sample_vector<T> exp(const sample_vector<T> & a) {
        sample_vector<T> r(a.size());
        cppmkl::vexp(a,r);
        return r;
    }
    template<typename T>
    sample_vector<T> exp2(const sample_vector<T> & a) {
        sample_vector<T> r(a.size());
        cppmkl::vexp2(a,r);
        return r;
    }
    template<typename T>
    sample_vector<T> exp10(const sample_vector<T> & a) {
        sample_vector<T> r(a.size());
        cppmkl::vexp10(a,r);
        return r;
    }
    template<typename T>
    sample_vector<T> expm1(const sample_vector<T> & a) {
        sample_vector<T> r(a.size());
        cppmkl::vexpm1(a,r);
        return r;
    }
    template<typename T>
    sample_vector<T> ln(const sample_vector<T> & a) {
        sample_vector<T> r(a.size());
        cppmkl::vln(a,r);
        return r;
    }
    template<typename T>
    sample_vector<T> log10(const sample_vector<T> & a) {
        sample_vector<T> r(a.size());
        cppmkl::vlog10(a,r);
        return r;
    }
    template<typename T>
    sample_vector<T> log2(const sample_vector<T> & a) {
        sample_vector<T> r(a.size());
        cppmkl::vlog2(a,r);
        return r;
    }
    template<typename T>
    sample_vector<T> logb(const sample_vector<T> & a) {
        sample_vector<T> r(a.size());
        cppmkl::vlogb(a,r);
        return r;
    }
    template<typename T>
    sample_vector<T> log1p(const sample_vector<T> & a) {
        sample_vector<T> r(a.size());
        cppmkl::vlog1p(a,r);
        return r;
    }
    template<typename T>
    sample_vector<T> cos(const sample_vector<T> & a) {
        sample_vector<T> r(a.size());
        cppmkl::vcos(a,r);
        return r;
    }
    template<typename T>
    sample_vector<T> sin(const sample_vector<T> & a) {
        sample_vector<T> r(a.size());
        cppmkl::vsin(a,r);
        return r;            
    }
    template<typename T>
    sample_vector<T> tan(const sample_vector<T> & a) {
        sample_vector<T> r(a.size());
        cppmkl::vtan(a,r);
        return r;
    }
    template<typename T>
    sample_vector<T> cosh(const sample_vector<T> & a) {
        sample_vector<T> r(a.size());
        cppmkl::vcosh(a,r);
        return r;
    }
    template<typename T>
    sample_vector<T> sinh(const sample_vector<T> & a) {
        sample_vector<T> r(a.size());
        cppmkl::vsinh(a,r);
        return r;            
    }
    template<typename T>
    sample_vector<T> tanh(const sample_vector<T> & a) {
        sample_vector<T> r(a.size());
        cppmkl::vtanh(a,r);
        return r;            
    }
    template<typename T>
    sample_vector<T> acos(const sample_vector<T> & a) {
        sample_vector<T> r(a.size());
        cppmkl::vacos(a,r);
        return r;
    }
    template<typename T>
    sample_vector<T> asin(const sample_vector<T> & a) {
        sample_vector<T> r(a.size());
        cppmkl::vasin(a,r);
        return r;            
    }
    template<typename T>
    sample_vector<T> atan(const sample_vector<T> & a) {
        sample_vector<T> r(a.size());
        cppmkl::vatan(a,r);
        return r;
    }
    template<typename T>
    sample_vector<T> atan2(const sample_vector<T> & a,const sample_vector<T> &n) {
        sample_vector<T> r(a.size());
        cppmkl::vatan2(a,n,r);
        return r;
    }
    template<typename T>
    sample_vector<T> acosh(const sample_vector<T> & a) {
        sample_vector<T> r(a.size());
        cppmkl::vacosh(a,r);
        return r;
    }
    template<typename T>
    sample_vector<T> asinh(const sample_vector<T> & a) {
        sample_vector<T> r(a.size());
        cppmkl::vasinh(a,r);
        return r;            
    }
    template<typename T>
    sample_vector<T> atanh(const sample_vector<T> & a) {
        sample_vector<T> r(a.size());
        cppmkl::vatanh(a,r);
        return r;
    }        
    template<typename T>
    void sincos(const sample_vector<T> & a, sample_vector<T> & b, sample_vector<T> & r) {        
        cppmkl::vsincos(a,b,r);        
    }
    template<typename T>
    sample_vector<T> erf(const sample_vector<T> & a) {
        sample_vector<T> r(a.size());
        cppmkl::verf(a,r);
        return r;
    }
    template<typename T>
    sample_vector<T> erfinv(const sample_vector<T> & a) {
        sample_vector<T> r(a.size());
        cppmkl::verfinv(a,r);
        return r;
    }
    template<typename T>
    sample_vector<T> erfc(const sample_vector<T> & a) {
        sample_vector<T> r(a.size());
        cppmkl::verfc(a,r);
        return r;
    }
    template<typename T>
    sample_vector<T> cdfnorm(const sample_vector<T> & a) {
        sample_vector<T> r(a.size());
        cppmkl::vcdfnorm(a,r);
        return r;
    }    
    template<typename T>
    sample_vector<T> cdfnorminv(const sample_vector<T> & a) {
        sample_vector<T> r(a.size());
        cppmkl::vcdfnorm(a,r);
        return r;
    }    
    template<typename T>
    sample_vector<T> floor(const sample_vector<T> & a) {
        sample_vector<T> r(a.size());
        cppmkl::vfloor(a,r);
        return r;
    }    
    template<typename T>
    sample_vector<T> ceil(const sample_vector<T> & a) {
        sample_vector<T> r(a.size());
        cppmkl::vceil(a,r);
        return r;
    }    
    template<typename T>
    sample_vector<T> trunc(const sample_vector<T> & a) {
        sample_vector<T> r(a.size());
        cppmkl::vtrunc(a,r);
        return r;
    }    
    template<typename T>
    sample_vector<T> round(const sample_vector<T> & a) {
        sample_vector<T> r(a.size());
        cppmkl::vround(a,r);
        return r;        
    }    
    template<typename T>
    sample_vector<T> nearbyint(const sample_vector<T> & a) {
        sample_vector<T> r(a.size());
        cppmkl::vnearbyint(a,r);
        return r;
    }    
    template<typename T>
    sample_vector<T> rint(const sample_vector<T> & a) {
        sample_vector<T> r(a.size());
        cppmkl::vrint(a,r);
        return r;
    }    
    template<typename T>
    sample_vector<T> fmod(const sample_vector<T> & a, sample_vector<T> & b) {
        sample_vector<T> r(a.size());
        cppmkl::vmodf(a,b,r);
        return r;
    }    
    
    template<typename T>
    sample_vector<T> CIS(const sample_vector<T> & a) {
        sample_vector<T> r(a.size());
        cppmkl::vCIS(a,r);
        return r;
    }    
    template<typename T>
    sample_vector<T> cospi(const sample_vector<T> & a) {
        sample_vector<T> r(a.size());
        cppmkl::vcospi(a,r);
        return r;
    }    
    template<typename T>
    sample_vector<T> sinpi(const sample_vector<T> & a) {
        sample_vector<T> r(a.size());
        cppmkl::vsinpi(a,r);
        return r;
    }    
    template<typename T>
    sample_vector<T> tanpi(const sample_vector<T> & a) {
        sample_vector<T> r(a.size());
        cppmkl::vtanpi(a,r);
        return r;
    }    
    template<typename T>
    sample_vector<T> acospi(const sample_vector<T> & a) {
        sample_vector<T> r(a.size());
        cppmkl::vacospi(a,r);
        return r;
    }    
    template<typename T>
    sample_vector<T> asinpi(const sample_vector<T> & a) {
        sample_vector<T> r(a.size());
        cppmkl::vasinpi(a,r);
        return r;
    }    
    template<typename T>
    sample_vector<T> atanpi(const sample_vector<T> & a) {
        sample_vector<T> r(a.size());
        cppmkl::vatanpi(a,r);
        return r;
    }    
    template<typename T>
    sample_vector<T> atan2pi(const sample_vector<T> & a, sample_vector<T> & b) {
        sample_vector<T> r(a.size());
        cppmkl::vatan2pi(a,b,r);
        return r;
    }    
    template<typename T>
    sample_vector<T> cosd(const sample_vector<T> & a) {
        sample_vector<T> r(a.size());
        cppmkl::vcosd(a.size(),a.data(),r.data());
        return r;
    }    
    template<typename T>
    sample_vector<T> sind(const sample_vector<T> & a) {
        sample_vector<T> r(a.size());
        cppmkl::vsind(a.size(),a.data(),r.data());
        return r;
    }    
    template<typename T>
    sample_vector<T> tand(const sample_vector<T> & a) {
        sample_vector<T> r(a.size());
        cppmkl::vtand(a.size(),a.data(),r.data());
        return r;
    }       
    template<typename T>
    sample_vector<T> lgamma(const sample_vector<T> & a) {
        sample_vector<T> r(a.size());
        cppmkl::vlgamma(a,r);
        return r;
    }       
    template<typename T>
    sample_vector<T> tgamma(const sample_vector<T> & a) {
        sample_vector<T> r(a.size());
        cppmkl::vtgamma(a,r);
        return r;
    }       
    template<typename T>
    sample_vector<T> expint1(const sample_vector<T> & a) {
        sample_vector<T> r(a.size());
        cppmkl::vexpint1(a,r);
        return r;
    }       
    template<typename T>
    sample_matrix<T> sqr(sample_matrix<T> & a) {
        sample_matrix<T> r(a.rows(),a.cols());
        cppmkl::vsqr(a,r);
        return r;
    }
    template<typename T>
    sample_matrix<T> abs(const sample_matrix<T> & a) {
        sample_matrix<T> r(a.rows(),a.cols());
        cppmkl::vabs(a,r);
        return r;
    }
    template<typename T>
    sample_matrix<T> inv(const sample_matrix<T> & a) {
        sample_matrix<T> r(a.rows(),a.cols());
        cppmkl::vinv(a,r);
        return r;            
    }
    template<typename T>
    sample_matrix<T> sqrt(const sample_matrix<T> & a) {
        sample_matrix<T> r(a.rows(),a.cols());
        cppmkl::vsqrt(a,r);
        return r;
    }
    template<typename T>
    sample_matrix<T> rsqrt(const sample_matrix<T> & a) {
        sample_matrix<T> r(a.rows(),a.cols());
        cppmkl::vinvsqrt(a,r);
        return r;
    }
    template<typename T>
    sample_matrix<T> cbrt(const sample_matrix<T> & a) {
        sample_matrix<T> r(a.rows(),a.cols());
        cppmkl::vcbrt(a,r);
        return r;
    }
    template<typename T>
    sample_matrix<T> rcbrt(const sample_matrix<T> & a) {
        sample_matrix<T> r(a.rows(),a.cols());
        cppmkl::vinvcbrt(a,r);
        return r;
    }
    template<typename T>
    sample_matrix<T> pow(const sample_matrix<T> & a,const sample_matrix<T> & b) {
        sample_matrix<T> r(a.rows(),a.cols());
        cppmkl::vpow(a,b,r);
        return r;
    }
    template<typename T>
    sample_matrix<T> pow2o3(const sample_matrix<T> & a) {
        sample_matrix<T> r(a.rows(),a.cols());
        cppmkl::vpow2o3(a,r);
        return r;
    }
    template<typename T>
    sample_matrix<T> pow3o2(const sample_matrix<T> & a) {
        sample_matrix<T> r(a.rows(),a.cols());
        cppmkl::vpow3o2(a,r);
        return r;
    }
    template<typename T>
    sample_matrix<T> pow(const sample_matrix<T> & a,const T b) {
        sample_matrix<T> r(a.rows(),a.cols());
        cppmkl::vpowx(a,b,r);
        return r;
    }
    template<typename T>
    sample_matrix<T> hypot(const sample_matrix<T> & a,const sample_matrix<T> & b) {
        sample_matrix<T> r(a.rows(),a.cols());
        cppmkl::vhypot(a,b,r);
        return r;
    }
    template<typename T>
    sample_matrix<T> exp(const sample_matrix<T> & a) {
        sample_matrix<T> r(a.rows(),a.cols());
        cppmkl::vexp(a,r);
        return r;
    }
    template<typename T>
    sample_matrix<T> exp2(const sample_matrix<T> & a) {
        sample_matrix<T> r(a.rows(),a.cols());
        cppmkl::vexp2(a,r);
        return r;
    }
    template<typename T>
    sample_matrix<T> exp10(const sample_matrix<T> & a) {
        sample_matrix<T> r(a.rows(),a.cols());
        cppmkl::vexp10(a,r);
        return r;
    }
    template<typename T>
    sample_matrix<T> expm1(const sample_matrix<T> & a) {
        sample_matrix<T> r(a.rows(),a.cols());
        cppmkl::vexpm1(a,r);
        return r;
    }
    template<typename T>
    sample_matrix<T> ln(const sample_matrix<T> & a) {
        sample_matrix<T> r(a.rows(),a.cols());
        cppmkl::vln(a,r);
        return r;
    }
    template<typename T>
    sample_matrix<T> log10(const sample_matrix<T> & a) {
        sample_matrix<T> r(a.rows(),a.cols());
        cppmkl::vlog10(a,r);
        return r;
    }
    template<typename T>
    sample_matrix<T> log2(const sample_matrix<T> & a) {
        sample_matrix<T> r(a.rows(),a.cols());
        cppmkl::vlog2(a,r);        
        return r;
    }
    template<typename T>
    sample_matrix<T> logb(const sample_matrix<T> & a) {
        sample_matrix<T> r(a.rows(),a.cols());
        cppmkl::vlogb(a,r);
        return r;
    }
    template<typename T>
    sample_matrix<T> log1p(const sample_matrix<T> & a) {
        sample_matrix<T> r(a.rows(),a.cols());
        cppmkl::vlog1p(a,r);
        return r;
    }
    template<typename T>
    sample_matrix<T> cos(const sample_matrix<T> & a) {
        sample_matrix<T> r(a.rows(),a.cols());
        cppmkl::vcos(a,r);
        return r;
    }
    template<typename T>
    sample_matrix<T> sin(const sample_matrix<T> & a) {
        sample_matrix<T> r(a.rows(),a.cols());
        cppmkl::vsin(a,r);
        return r;            
    }
    template<typename T>
    sample_matrix<T> tan(const sample_matrix<T> & a) {
        sample_matrix<T> r(a.rows(),a.cols());
        cppmkl::vtan(a,r);
        return r;
    }
    template<typename T>
    sample_matrix<T> cosh(const sample_matrix<T> & a) {
        sample_matrix<T> r(a.rows(),a.cols());
        cppmkl::vcosh(a,r);
        return r;
    }
    template<typename T>
    sample_matrix<T> sinh(const sample_matrix<T> & a) {
        sample_matrix<T> r(a.rows(),a.cols());
        cppmkl::vsinh(a,r);
        return r;            
    }
    template<typename T>
    sample_matrix<T> tanh(const sample_matrix<T> & a) {
        sample_matrix<T> r(a.rows(),a.cols());
        cppmkl::vtanh(a,r);
        return r;            
    }
    template<typename T>
    sample_matrix<T> acos(const sample_matrix<T> & a) {
        sample_matrix<T> r(a.rows(),a.cols());
        cppmkl::vacos(a,r);
        return r;
    }
    template<typename T>
    sample_matrix<T> asin(const sample_matrix<T> & a) {
        sample_matrix<T> r(a.rows(),a.cols());
        cppmkl::vasin(a,r);
        return r;            
    }
    template<typename T>
    sample_matrix<T> atan(const sample_matrix<T> & a) {
        sample_matrix<T> r(a.rows(),a.cols());
        cppmkl::vatan(a,r);
        return r;
    }
    template<typename T>
    sample_matrix<T> atan2(const sample_matrix<T> & a,const sample_matrix<T> &n) {
        sample_matrix<T> r(a.rows(),a.cols());
        cppmkl::vatan2(a,n,r);
        return r;
    }
    template<typename T>
    sample_matrix<T> acosh(const sample_matrix<T> & a) {
        sample_matrix<T> r(a.rows(),a.cols());
        cppmkl::vacosh(a,r);
        return r;
    }
    template<typename T>
    sample_matrix<T> asinh(const sample_matrix<T> & a) {
        sample_matrix<T> r(a.rows(),a.cols());
        cppmkl::vasinh(a,r);
        return r;            
    }
    template<typename T>
    sample_matrix<T> atanh(const sample_matrix<T> & a) {
        sample_matrix<T> r(a.rows(),a.cols());
        cppmkl::vatanh(a,r);
        return r;
    }        
    template<typename T>
    void sincos(const sample_matrix<T> & a, sample_matrix<T> & b, sample_matrix<T> & r) {        
        cppmkl::vsincos(a,b,r);
    }
    template<typename T>
    sample_matrix<T> erf(const sample_matrix<T> & a) {
        sample_matrix<T> r(a.rows(),a.cols());
        cppmkl::verf(a,r);
        return r;
    }
    template<typename T>
    sample_matrix<T> erfinv(const sample_matrix<T> & a) {
        sample_matrix<T> r(a.rows(),a.cols());
        cppmkl::verfinv(a,r);
        return r;
    }
    template<typename T>
    sample_matrix<T> erfc(const sample_matrix<T> & a) {
        sample_matrix<T> r(a.rows(),a.cols());
        cppmkl::verfc(a,r);
        return r;
    }
    template<typename T>
    sample_matrix<T> cdfnorm(const sample_matrix<T> & a) {
        sample_matrix<T> r(a.rows(),a.cols());
        cppmkl::vcdfnorm(a,r);
        return r;
    }    
    template<typename T>
    sample_matrix<T> cdfnorminv(const sample_matrix<T> & a) {
        sample_matrix<T> r(a.rows(),a.cols());
        cppmkl::vcdfnorm(a,r);
        return r;
    }    
    template<typename T>
    sample_matrix<T> floor(const sample_matrix<T> & a) {
        sample_matrix<T> r(a.rows(),a.cols());
        cppmkl::vfloor(a,r);
        return r;
    }    
    template<typename T>
    sample_matrix<T> ceil(const sample_matrix<T> & a) {
        sample_matrix<T> r(a.rows(),a.cols());
        cppmkl::vceil(a,r);
        return r;
    }    
    template<typename T>
    sample_matrix<T> trunc(const sample_matrix<T> & a) {
        sample_matrix<T> r(a.rows(),a.cols());
        cppmkl::vtrunc(a,r);
        return r;
    }    
    template<typename T>
    sample_matrix<T> round(const sample_matrix<T> & a) {
        sample_matrix<T> r(a.rows(),a.cols());
        cppmkl::vround(a,r);
        return r;        
    }    
    template<typename T>
    sample_matrix<T> nearbyint(const sample_matrix<T> & a) {
        sample_matrix<T> r(a.rows(),a.cols());
        cppmkl::vnearbyint(a,r);
        return r;
    }    
    template<typename T>
    sample_matrix<T> rint(const sample_matrix<T> & a) {
        sample_matrix<T> r(a.rows(),a.cols());
        cppmkl::vrint(a,r);
        return r;
    }    
    template<typename T>
    sample_matrix<T> fmod(const sample_matrix<T> & a, sample_matrix<T> & b) {
        sample_matrix<T> r(a.rows(),a.cols());
        cppmkl::vmodf(a,b,r);
        return r;
    }    
    template<typename T>
    sample_matrix<T> mulbyconj(const sample_matrix<T> & a, const sample_matrix<T> & b) {
        sample_matrix<T> r(a.rows(),a.cols());
        cppmkl::vmulbyconj(a,b,r);
        return r;
    }    
    template<typename T>
    sample_matrix<T> conj(const sample_matrix<T> & a) {
        sample_matrix<T> r(a.rows(),a.cols());
        cppmkl::vconj(a,r);
        return r;
    }    
    template<typename T>
    sample_matrix<T> arg(const sample_matrix<T> & a) {
        sample_matrix<T> r(a.rows(),a.cols());
        cppmkl::varg(a,r);
        return r;
    }    
    template<typename T>
    sample_matrix<T> CIS(const sample_matrix<T> & a) {
        sample_matrix<T> r(a.rows(),a.cols());
        cppmkl::vCIS(a,r);
        return r;
    }    
    template<typename T>
    sample_matrix<T> cospi(const sample_matrix<T> & a) {
        sample_matrix<T> r(a.rows(),a.cols());
        cppmkl::vcospi(a,r);
        return r;
    }    
    template<typename T>
    sample_matrix<T> sinpi(const sample_matrix<T> & a) {
        sample_matrix<T> r(a.rows(),a.cols());
        cppmkl::vsinpi(a,r);
        return r;
    }    
    template<typename T>
    sample_matrix<T> tanpi(const sample_matrix<T> & a) {
        sample_matrix<T> r(a.rows(),a.cols());
        cppmkl::vtanpi(a,r);
        return r;
    }    
    template<typename T>
    sample_matrix<T> acospi(const sample_matrix<T> & a) {
        sample_matrix<T> r(a.rows(),a.cols());
        cppmkl::vacospi(a,r);
        return r;
    }    
    template<typename T>
    sample_matrix<T> asinpi(const sample_matrix<T> & a) {
        sample_matrix<T> r(a.rows(),a.cols());
        cppmkl::vasinpi(a,r);
        return r;
    }    
    template<typename T>
    sample_matrix<T> atanpi(const sample_matrix<T> & a) {
        sample_matrix<T> r(a.rows(),a.cols());
        cppmkl::vatanpi(a,r);
        return r;
    }    
    template<typename T>
    sample_matrix<T> atan2pi(const sample_matrix<T> & a, sample_matrix<T> & b) {
        sample_matrix<T> r(a.rows(),a.cols());
        cppmkl::vatan2pi(a,b,r);
        return r;
    }    
    template<typename T>
    sample_matrix<T> cosd(const sample_matrix<T> & a) {
        sample_matrix<T> r(a.rows(),a.cols());
        cppmkl::vcosd(a.size(),a.data(),r.data());
        return r;
    }    
    template<typename T>
    sample_matrix<T> sind(const sample_matrix<T> & a) {
        sample_matrix<T> r(a.rows(),a.cols());
        cppmkl::vsind(a.size(),a.data(),r.data());
        return r;
    }    
    template<typename T>
    sample_matrix<T> tand(const sample_matrix<T> & a) {
        sample_matrix<T> r(a.rows(),a.cols());
        cppmkl::vtand(a.size(),a.data(),r.data());
        return r;
    }       
    template<typename T>
    sample_matrix<T> lgamma(const sample_matrix<T> & a) {
        sample_matrix<T> r(a.rows(),a.cols());
        cppmkl::vlgamma(a,r);
        return r;
    }       
    template<typename T>
    sample_matrix<T> tgamma(const sample_matrix<T> & a) {
        sample_matrix<T> r(a.rows(),a.cols());
        cppmkl::vtgamma(a,r);
        return r;
    }       
    template<typename T>
    sample_matrix<T> expint1(const sample_matrix<T> & a) {
        sample_matrix<T> r(a.rows(),a.cols());
        cppmkl::vexpint1(a,r);
        return r;
    }       

    template<typename T>
    sample_vector<T> copy(const sample_vector<T> & a) {
        sample_vector<T> r(a.size());
        cppmkl::cblas_copy(a.size(),a.data(),1,r.data(),1);
        return r;
    }       

    template<typename T> T sum(const sample_vector<T> & a) {        
        return cppmkl::cblas_asum(a.size(), a.data(),1);        
    }       

    template<typename T>
    sample_vector<T> add(const sample_vector<T> & a, const sample_vector<T> & b) {
        sample_vector<T> r(b);
        cppmkl::cblas_axpy(a.size(),1.0,a.data(),1,r.data(),1);
        return r;
    }       
    template<typename T>
    sample_vector<T> sub(const sample_vector<T> & a, const sample_vector<T> & b) {
        sample_vector<T> r(a);
        cppmkl::cblas_axpy(a.size(),-1.0,b.data(),1,r.data(),1);
        return r;
    }       
    template<typename T>
    T dot(const sample_vector<T> & a, sample_vector<T> & b) {
        return cppmkl::cblas_dot(a,b);
    }       
    template<typename T>
    T nrm2(const sample_vector<T> & a) {
        sample_vector<T> r(a);
        return cppmkl::cblas_nrm2(a);        
    }       
    
    template<typename T>
    void scale(sample_vector<T> & x, T alpha) {
        cppmkl::cblas_scal(x.size(),alpha,x.data(),1);
    }


    template<typename T>
    size_t min_index(const sample_matrix<T> & m) { return cppmkl::cblas_iamin(m.size(),m.data(),1); }
    template<typename T>
    size_t max_index(const sample_matrix<T> & m) { return cppmkl::cblas_iamax(m.size(),m.data(),1); }
    
    template<typename T>
    size_t min_index(const sample_vector<T> & v) { return cppmkl::cblas_iamin(v.size(),v.data(),1); }
    template<typename T>
    size_t max_index(const sample_vector<T> & v) { return cppmkl::cblas_iamax(v.size(),v.data(),1); }

    template<typename T>
    sample_vector<T> linspace(size_t n, T start, T inc=1) {
        sample_vector<T> r(n);
        for(size_t i = 0; i < n; i++) {
            r[i] = start + i*inc;
        }
        return r;
    }
    template<typename T>
    sample_vector<T> linspace(T start, T end, T inc=1) {
        size_t n = (end - start)/inc;
        sample_vector<T> r(n);
        for(size_t i = 0; i < n; i++) {
            r[i] = start + i*inc;
        }
        return r;
    }

    template<typename T>
    sample_vector<T> operator * (T a, const sample_vector<T> & b) {
        sample_vector<T> x(b.size());
        x.fill(a);
        sample_vector<T> r(b.size());
        cppmkl::vmul(x,b,r);            
        return r;
    }
    template<typename T>
    sample_vector<T> operator / (T a, const sample_vector<T> & b) {
        sample_vector<T> x(b.size());
        x.fill(a);
        sample_vector<T> r(b.size());
        cppmkl::vdiv(x,b,r);            
        return r;
    }
    template<typename T>
    sample_vector<T> operator + (T a, const sample_vector<T> & b) {
        sample_vector<T> x(b.size());
        x.fill(a);
        sample_vector<T> r(b.size());
        cppmkl::vadd(x,b,r);            
        return r;
    }
    template<typename T>
    sample_vector<T> operator - (T a, const sample_vector<T> & b) {
        sample_vector<T> x(b.size());
        x.fill(a);
        sample_vector<T> r(b.size());
        cppmkl::vsub(x,b,r);            
        return r;
    }

    template<typename T>
    sample_vector<T> operator * (const sample_vector<T> & a, T b) {
        sample_vector<T> x(a.size());
        x.fill(b);
        sample_vector<T> r(a.size());
        cppmkl::vmul(a,x,r);            
        return r;
    }
    template<typename T>
    sample_vector<T> operator / (const sample_vector<T> & a , T b) {
        sample_vector<T> x(a.size());
        x.fill(b);
        sample_vector<T> r(a.size());
        cppmkl::vdiv(a,x,r);            
        return r;
    }
    template<typename T>
    sample_vector<T> operator + (const sample_vector<T> & a, T b) {
        sample_vector<T> x(a.size());
        x.fill(b);
        sample_vector<T> r(a.size());
        cppmkl::vadd(a,x,r);            
        return r;
    }
    template<typename T>
    sample_vector<T> operator - (const sample_vector<T> & a, T b) {
        sample_vector<T> x(a.size());
        x.fill(b);
        sample_vector<T> r(a.size());
        cppmkl::vsub(a,x,r);            
        return r;
    }
    
    template<typename T>
    sample_vector<T> operator - (const sample_vector<T> & a, const sample_vector<T> & b) {        
        sample_vector<T> r(a.size());
        cppmkl::vsub(a,b,r);            
        return r;
    }
    template<typename T>
    sample_vector<T> operator + (const sample_vector<T> & a, const sample_vector<T> & b) {        
        sample_vector<T> r(a.size());
        cppmkl::vadd(a,b,r);                    
        return r;
    }
    template<typename T>
    sample_vector<T> operator * (const sample_vector<T> & a, const sample_vector<T> & b) {        
        sample_vector<T> r(a.size());
        cppmkl::vmul(a,b,r);            
        return r;
    }
    template<typename T>
    sample_vector<T> operator / (const sample_vector<T> & a, const sample_vector<T> & b) {        
        sample_vector<T> r(a.size());
        cppmkl::vdiv(a,b,r);            
        return r;
    }

    template<typename T>
    sample_matrix<T> operator * (const sample_matrix<T> & a, const sample_matrix<T> & b) {
        assert(a.N == b.M);
        sample_matrix<T> r(a.rows(),b.cols());
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
    sample_matrix<T> operator + (const sample_matrix<T> & a, const sample_matrix<T> & b) {
        assert(a.rows() == b.rows() && a.cols() == b.cols());
        sample_matrix<T> r(a.rows(),a.cols());
        cppmkl::vadd(a,b,r);            
        return r;
    }
    template<typename T>
    sample_matrix<T> operator / (const sample_matrix<T> & a, const sample_matrix<T> & b) {
        assert(a.rows() == b.rows() && a.cols() == b.cols());
        sample_matrix<T> r(a.rows(),a.cols());
        cppmkl::vdiv(a,b,r);            
        return r;
    }        
    template<typename T>
    sample_matrix<T> operator - (const sample_matrix<T> & a, const sample_matrix<T> & b) {
        assert(a.rows() == b.rows() && a.cols() == b.cols());
        sample_matrix<T> r(a.rows(),a.cols());
        cppmkl::vsub(a,b,r);            
        return r;
    }        

    template<typename T>
    sample_matrix<T> operator * (T a, const sample_matrix<T> & b) {
        sample_matrix<T> x(b.M,b.N);
        x.fill(a);
        sample_matrix<T> r(b.M,b.N);
        cppmkl::vmul(x,b,r);
        return r;
    }
    
    template<typename T>
    sample_matrix<T> operator + (T a, const sample_matrix<T> & b) {
        sample_matrix<T> x(b.M,b.N);
        x.fill(a);
        sample_matrix<T> r(b.M,b.N);
        cppmkl::vadd(x,b,r);
        return r;
    }
    template<typename T>
    sample_matrix<T> operator - (T a, const sample_matrix<T> & b) {
        sample_matrix<T> x(b.M,b.N);
        x.fill(a);
        sample_matrix<T> r(b.M,b.N);
        cppmkl::vsub(x,b,r);
        return r;
    }
    template<typename T>
    sample_matrix<T> operator / (T a, const sample_matrix<T> & b) {
        sample_matrix<T> x(b.M,b.N);
        x.fill(a);
        sample_matrix<T> r(b.M,b.N);
        cppmkl::vdiv(x,b,r);
        return r;
    }
    template<typename T>
    sample_matrix<T> operator * (const sample_matrix<T> & a, T b) {
        sample_matrix<T> x(a.M,a.N);
        x.fill(b);
        sample_matrix<T> r(a.M,a.N);
        cppmkl::vmul(a,x,r);
        return r;
    }
    
    template<typename T>
    sample_matrix<T> operator + (const sample_matrix<T> & a, T b) {
        sample_matrix<T> x(a.M,a.N);
        x.fill(b);
        sample_matrix<T> r(a.M,a.N);
        cppmkl::vadd(a,x,r);
        return r;
    }
    template<typename T>
    sample_matrix<T> operator - (const sample_matrix<T> & a, T b) {
        sample_matrix<T> x(a.M,a.N);
        x.fill(b);
        sample_matrix<T> r(a.M,a.N);
        cppmkl::vsub(a,x,r);
        return r;
    }
    template<typename T>
    sample_matrix<T> operator / (const sample_matrix<T> & a, T b) {
        sample_matrix<T> x(a.M,a.N);
        x.fill(b);
        sample_matrix<T> r(a.M,a.N);
        cppmkl::vdiv(a,x,r);
        return r;
    }    
    template<typename T>
    sample_matrix<T> hadamard(const sample_matrix<T> & a, const sample_matrix<T> & b) {
        sample_matrix<T> r(a.rows(),a.cols());
        assert(a.rows() == b.rows() && a.cols() == b.cols());
        cppmkl::vmul(a,b,r);
        return r;
    }      
    template<typename T>
    sample_vector<T> operator *(const sample_vector<T> &a, const sample_matrix<T> &b) {
        sample_vector<T> r(a.size());        
        sgemv(CblasRowMajor,CblasNoTrans,b.rows(),b.cols(),1.0,b.data(),b.cols(),a.data(),1,1.0,r.data(),1);
        return r;
    }      
    template<typename T>
    sample_vector<T> operator *(const sample_matrix<T> &a, const sample_vector<T> &b) {
        sample_vector<T> r(b.size());        
        sgemv(CblasRowMajor,CblasNoTrans,a.rows(),a.cols(),1.0,a.data(),a.cols(),b.data(),1,1.0,r.data(),1);
        return r;
    }      
#endif    
}