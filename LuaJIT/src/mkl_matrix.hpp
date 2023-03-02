#pragma once

#include "Std.h"
#include <ccomplex>
#include <complex>
#include "cppmkl/cppmkl_allocator.h"
#include "cppmkl/cppmkl_vml.h"
#include "cppmkl/cppmkl_cblas.h"

namespace MKL
{
    
    template<typename T> struct Matrix;
    template<typename T>
    struct MatrixView 
    {
        Matrix<T> * matrix;
        size_t row;
        MatrixView(Matrix<T> * m, size_t r) {
            matrix = m;
            row = r;
        }

        T& operator[](size_t i);
        T  __getitem(size_t i);
        void __setitem(size_t i, T x);
    };

    template<typename T>
    struct Matrix
    {
        Vector<T> matrix;
        size_t M;
        size_t N;

        Matrix() { M = N = 0; }
        Matrix(size_t m, size_t n) {
            matrix.resize(m*n);
            M = m;
            N = n;
            assert(M > 0);
            assert(N > 0);
        }
        Matrix(const Matrix<T> & m) {
            matrix = m.matrix;
            M = m.M;
            N = m.N;
        }
        Matrix(T * ptr, size_t m, size_t n)
        {
            matrix.resize(m*n);
            M = m;
            N = n;
            memcpy(matrix.data(),ptr,m*n*sizeof(T));
        }

        Matrix<T>& operator = (const Matrix<T> & m) {
            matrix = m.matrix;
            M = m.M;
            N = m.N;
            return *this;
        }

        size_t size() const { return matrix.size(); }
        size_t rows() const { return M; }
        size_t cols() const { return N; }

        void resize(size_t r, size_t c) {
            M = r;
            N = c;
            matrix.resize(r*c);
        }
        T& operator()(size_t i, size_t j) { return matrix[i*N + j]; }
        
        MatrixView<T> operator[](size_t row) { return MatrixView<T>(&matrix,row); }

        Matrix<T> operator + (const Matrix<T> & m) {
            Matrix<T> r(rows(),cols());
            r.matrix = matrix + m.matrix;
            return r;
        }
        Matrix<T> operator - (const Matrix<T> & m) {
            Matrix<T> r(rows(),cols());
            r.matrix = matrix - m.matrix;
            return r;
        }
        Matrix<T> operator * (Matrix<T> & b) {
            assert(N == b.M);
            Matrix<T> r(rows(),b.cols());
            int m = M;
            int k = N;
            int n = b.N;
            cppmkl::cblas_gemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,m,n,k,1.0,matrix.vector.data(),m,b.matrix.vector.data(),n,0.0,r.matrix.vector.data(),n);
            return r;
        }
                

        Matrix<T> hadamard(Matrix<T> & m) {
            Matrix<T> r(rows(),cols());
            assert(rows() == m.rows() && cols() == m.cols());
            cppmkl::vmul(matrix.vector,m.matrix.vector,r.matrix.vector);
            return r;
        }

        Matrix<T> transpose() {
            Matrix<T> r(cols(),rows());
            for(size_t i = 0; i < rows(); i++)
                for(size_t j = 0; j < cols(); j++)
                    r(j,i) = matrix[i*N + j];
            return r;
        }
        void zero() {
            memset(matrix.vector.data(),0x00,size()*sizeof(T));
        }
        void fill(T x) {
            for(size_t i = 0; i < size(); i++) matrix[i] = x;
        }
        void ones() {
            fill((T)1);
        }
        void random(T min = T(0), T max = T(1)) {
            Std::StdRandomUniform r;
            for(size_t i = 0; i < size(); i++) matrix[i] = r.random(min,max);
        }
        void identity() {
            size_t x = 0;
            zero();
            for(size_t i = 0; i < rows(); i++)
            {
                matrix[i*N + x++] = 0;
            }
        }
        Vector<T>& get_matrix() { return matrix; }

    };

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
    Matrix<T> matmul(Matrix<T> & a, Matrix<T> & b)
    {
        Matrix<T> r = a * b;
        return r;
    }

    template<typename T>
    Matrix<T> sqr(Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vsqr(a.matrix.vector,r.matrix.vector);
        return r;
    }
    template<typename T>
    Matrix<T> abs(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vabs(a.matrix.vector,r.matrix.vector);
        return r;
    }
    template<typename T>
    Matrix<T> inv(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());;
        cppmkl::vinv(a.matrix.vector,r.matrix.vector);
        return r;            
    }
    template<typename T>
    Matrix<T> sqrt(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());;
        cppmkl::vsqrt(a.matrix.vector,r.matrix.vector);
        return r;
    }
    template<typename T>
    Matrix<T> rsqrt(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());;
        cppmkl::vinvsqrt(a.matrix.vector,r.matrix.vector);
        return r;
    }
    template<typename T>
    Matrix<T> cbrt(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());;
        cppmkl::vcbrt(a.matrix.vector,r.matrix.vector);
        return r;
    }
    template<typename T>
    Matrix<T> rcbrt(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());;
        cppmkl::vinvcbrt(a.matrix.vector,r.matrix.vector);
        return r;
    }
    template<typename T>
    Matrix<T> pow(const Matrix<T> & a,const Matrix<T> & b) {
        Matrix<T> r(a.rows(),a.cols());;
        cppmkl::vpow(a.matrix.vector,b.matrix.vector,r.matrix.vector);
        return r;
    }
    template<typename T>
    Matrix<T> pow2o3(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());;
        cppmkl::vpow2o3(a.matrix.vector,r.matrix.vector);
        return r;
    }
    template<typename T>
    Matrix<T> pow3o2(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());;
        cppmkl::vpow3o2(a.matrix.vector,r.matrix.vector);
        return r;
    }
    template<typename T>
    Matrix<T> pow(const Matrix<T> & a,const T b) {
        Matrix<T> r(a.rows(),a.cols());;
        cppmkl::vpowx(a.matrix.vector,b,r.matrix.vector);
        return r;
    }
    template<typename T>
    Matrix<T> hypot(const Matrix<T> & a,const Matrix<T> & b) {
        Matrix<T> r(a.rows(),a.cols());;
        cppmkl::vhypot(a.matrix.vector,b.matrix.vector,r.matrix.vector);
        return r;
    }
    template<typename T>
    Matrix<T> exp(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());;
        cppmkl::vexp(a.matrix.vector,r.matrix.vector);
        return r;
    }
    template<typename T>
    Matrix<T> exp2(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());;
        cppmkl::vexp2(a.matrix.vector,r.matrix.vector);
        return r;
    }
    template<typename T>
    Matrix<T> exp10(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());;
        cppmkl::vexp10(a.matrix.vector,r.matrix.vector);
        return r;
    }
    template<typename T>
    Matrix<T> expm1(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());;
        cppmkl::vexpm1(a.matrix.vector,r.matrix.vector);
        return r;
    }
    template<typename T>
    Matrix<T> ln(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());;
        cppmkl::vln(a.matrix.vector,r.matrix.vector);
        return r;
    }
    template<typename T>
    Matrix<T> log10(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());;
        cppmkl::vlog10(a.matrix.vector,r.matrix.vector);
        return r;
    }
    template<typename T>
    Matrix<T> log2(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());;
        cppmkl::vlog2(a.matrix.vector,r.matrix.vector);
        return r;
    }
    template<typename T>
    Matrix<T> logb(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());;
        cppmkl::vlogb(a.matrix.vector,r.matrix.vector);
        return r;
    }
    template<typename T>
    Matrix<T> log1p(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());;
        cppmkl::vlog1p(a.matrix.vector,r.matrix.vector);
        return r;
    }
    template<typename T>
    Matrix<T> cos(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());;
        cppmkl::vcos(a.matrix.vector,r.matrix.vector);
        return r;
    }
    template<typename T>
    Matrix<T> sin(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());;
        cppmkl::vsin(a.matrix.vector,r.matrix.vector);
        return r;            
    }
    template<typename T>
    Matrix<T> tan(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());;
        cppmkl::vtan(a.matrix.vector,r.matrix.vector);
        return r;
    }
    template<typename T>
    Matrix<T> cosh(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());;
        cppmkl::vcosh(a.matrix.vector,r.matrix.vector);
        return r;
    }
    template<typename T>
    Matrix<T> sinh(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());;
        cppmkl::vsinh(a.matrix.vector,r.matrix.vector);
        return r;            
    }
    template<typename T>
    Matrix<T> tanh(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());;
        cppmkl::vtanh(a.matrix.vector,r.matrix.vector);
        return r;            
    }
    template<typename T>
    Matrix<T> acos(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());;
        cppmkl::vacos(a.matrix.vector,r.matrix.vector);
        return r;
    }
    template<typename T>
    Matrix<T> asin(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());;
        cppmkl::vasin(a.matrix.vector,r.matrix.vector);
        return r;            
    }
    template<typename T>
    Matrix<T> atan(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());;
        cppmkl::vatan(a.matrix.vector,r.matrix.vector);
        return r;
    }
    template<typename T>
    Matrix<T> atan2(const Matrix<T> & a,const Matrix<T> &n) {
        Matrix<T> r(a.rows(),a.cols());;
        cppmkl::vatan2(a.matrix.vector,n.matrix.vector,r.matrix.vector);
        return r;
    }
    template<typename T>
    Matrix<T> acosh(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());;
        cppmkl::vacosh(a.matrix.vector,r.matrix.vector);
        return r;
    }
    template<typename T>
    Matrix<T> asinh(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());;
        cppmkl::vasinh(a.matrix.vector,r.matrix.vector);
        return r;            
    }
    template<typename T>
    Matrix<T> atanh(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());;
        cppmkl::vatanh(a.matrix.vector,r.matrix.vector);
        return r;
    }        
    template<typename T>
    void sincos(const Matrix<T> & a, Matrix<T> & b, Matrix<T> & r) {        
        cppmkl::vsincos(a.matrix.vector,b.matrix.vector,r.matrix.vector);
    }
    template<typename T>
    Matrix<T> erf(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());;
        cppmkl::verf(a.matrix.vector,r.matrix.vector);
        return r;
    }
    template<typename T>
    Matrix<T> erfinv(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());;
        cppmkl::verfinv(a.matrix.vector,r.matrix.vector);
        return r;
    }
    template<typename T>
    Matrix<T> erfc(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());;
        cppmkl::verfc(a.matrix.vector,r.matrix.vector);
        return r;
    }
    template<typename T>
    Matrix<T> cdfnorm(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());;
        cppmkl::vcdfnorm(a.matrix.vector,r.matrix.vector);
        return r;
    }    
    template<typename T>
    Matrix<T> cdfnorminv(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());;
        cppmkl::vcdfnorm(a.matrix.vector,r.matrix.vector);
        return r;
    }    
    template<typename T>
    Matrix<T> floor(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());;
        cppmkl::vfloor(a.matrix.vector,r.matrix.vector);
        return r;
    }    
    template<typename T>
    Matrix<T> ceil(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());;
        cppmkl::vceil(a.matrix.vector,r.matrix.vector);
        return r;
    }    
    template<typename T>
    Matrix<T> trunc(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());;
        cppmkl::vtrunc(a.matrix.vector,r.matrix.vector);
        return r;
    }    
    template<typename T>
    Matrix<T> round(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());;
        cppmkl::vround(a.matrix.vector,r.matrix.vector);
        return r;        
    }    
    template<typename T>
    Matrix<T> nearbyint(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());;
        cppmkl::vnearbyint(a.matrix.vector,r.matrix.vector);
        return r;
    }    
    template<typename T>
    Matrix<T> rint(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());;
        cppmkl::vrint(a.matrix.vector,r.matrix.vector);
        return r;
    }    
    template<typename T>
    Matrix<T> fmod(const Matrix<T> & a, Matrix<T> & b) {
        Matrix<T> r(a.rows(),a.cols());;
        cppmkl::vmodf(a.matrix.vector,b.matrix.vector,r.matrix.vector);
        return r;
    }    
    template<typename T>
    Matrix<T> mulbyconj(const Matrix<T> & a, const Matrix<T> & b) {
        Matrix<T> r(a.rows(),a.cols());;
        cppmkl::vmulbyconj(a.matrix.vector,b.matrix.vector,r.matrix.vector);
        return r;
    }    
    template<typename T>
    Matrix<T> conj(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());;
        cppmkl::vconj(a.matrix.vector,r.matrix.vector);
        return r;
    }    
    template<typename T>
    Matrix<T> arg(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());;
        cppmkl::varg(a.matrix.vector,r.matrix.vector);
        return r;
    }    
    template<typename T>
    Matrix<T> CIS(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());;
        cppmkl::vCIS(a.matrix.vector,r.matrix.vector);
        return r;
    }    
    template<typename T>
    Matrix<T> cospi(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());;
        cppmkl::vcospi(a.matrix.vector,r.matrix.vector);
        return r;
    }    
    template<typename T>
    Matrix<T> sinpi(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());;
        cppmkl::vsinpi(a.matrix.vector,r.matrix.vector);
        return r;
    }    
    template<typename T>
    Matrix<T> tanpi(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());;
        cppmkl::vtanpi(a.matrix.vector,r.matrix.vector);
        return r;
    }    
    template<typename T>
    Matrix<T> acospi(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());;
        cppmkl::vacospi(a.matrix.vector,r.matrix.vector);
        return r;
    }    
    template<typename T>
    Matrix<T> asinpi(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());;
        cppmkl::vasinpi(a.matrix.vector,r.matrix.vector);
        return r;
    }    
    template<typename T>
    Matrix<T> atanpi(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());;
        cppmkl::vatanpi(a.matrix.vector,r.matrix.vector);
        return r;
    }    
    template<typename T>
    Matrix<T> atan2pi(const Matrix<T> & a, Matrix<T> & b) {
        Matrix<T> r(a.rows(),a.cols());;
        cppmkl::vatan2pi(a.matrix.vector,b.matrix.vector,r.matrix.vector);
        return r;
    }    
    template<typename T>
    Matrix<T> cosd(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());;
        cppmkl::vcosd(a.size(),a.matrix.vector.data(),r.matrix.vector.data());
        return r;
    }    
    template<typename T>
    Matrix<T> sind(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());;
        cppmkl::vsind(a.size(),a.matrix.vector.data(),r.matrix.vector.data());
        return r;
    }    
    template<typename T>
    Matrix<T> tand(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());;
        cppmkl::vtand(a.size(),a.matrix.vector.data(),r.matrix.vector.data());
        return r;
    }       
    template<typename T>
    Matrix<T> lgamma(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());;
        cppmkl::vlgamma(a.matrix.vector,r.matrix.vector);
        return r;
    }       
    template<typename T>
    Matrix<T> tgamma(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());;
        cppmkl::vtgamma(a.matrix.vector,r.matrix.vector);
        return r;
    }       
    template<typename T>
    Matrix<T> expint1(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());;
        cppmkl::vexpint1(a.matrix.vector,r.matrix.vector);
        return r;
    }       
}
