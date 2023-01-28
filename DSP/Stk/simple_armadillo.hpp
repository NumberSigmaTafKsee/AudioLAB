#pragma once

#include<armadillo>
#include<vector>
#include<complex>

namespace Armadillo
{

    template<typename T>
    struct Vector : public arma::Col<T> {
        arma::Col<T> vector;

        Vector() = default;
        Vector(size_t n) : arma::Col<T>(n) {}
        Vector(const Vector<T> & v) : arma::Col<T>(v) {}
        Vector(const arma::Col<T> & v) : arma::Col<T>(v) {}
        //Vector(std::initializer_list<T> & l) { vector = arma::Col<T>(l.begin(),l.end()}; }
        ~Vector() = default;

        T    __getitem__(size_t i) { return (*this)(i); }
        void __setitem__(size_t i, T v) { (*this)(i) = v; }


        using arma::Col<T>::operator =;
        using arma::Col<T>::operator -;
        using arma::Col<T>::operator +;
        using arma::Col<T>::operator *;
        using arma::Col<T>::operator /;

        using arma::Col<T>::at;
        using arma::Col<T>::operator ();
        using arma::Col<T>::operator [];

        using arma::Col<T>::size;
        using arma::Col<T>::zeros;
        using arma::Col<T>::ones;
        using arma::Col<T>::randu;
        using arma::Col<T>::random;
        using arma::Col<T>::randn;
        using arma::Col<T>::fill;
        using arma::Col<T>::replace;
        using arma::Col<T>::clamp;
        using arma::Col<T>::set_size;
        using arma::Col<T>::resize;
        using arma::Col<T>::copy_size;
        using arma::Col<T>::reset;
        using arma::Col<T>::data;
        using arma::Col<T>::print;
        using arma::Col<T>::raw_print;
        using arma::Col<T>::brief_print;        

    };

    template<typename T>
    struct RowVector : public arma::Row<T> 
    {    
        RowVector() = default;
        RowVector(size_t n) : arma::Row<T> (n) {}
        RowVector(const RowVector<T> & v) : arma::Row<T>(v) {}
        RowVector(const arma::Row<T> & v) : arma::Row<T> (v) {}
        //RowVector(std::initializer_list<T> & l) { vector = arma::Row<T>(l}; }
        ~RowVector() = default;

        T    __getitem__(size_t i) { return (*this)(i); }
        void __setitem__(size_t i, T v) { (*this)(i) = v; }


        using arma::Row<T>::operator =;
        using arma::Row<T>::operator -;
        using arma::Row<T>::operator +;
        using arma::Row<T>::operator *;
        using arma::Row<T>::operator /;

        using arma::Row<T>::at;
        using arma::Row<T>::operator ();
        using arma::Row<T>::operator [];

        using arma::Row<T>::size;
        using arma::Row<T>::zeros;
        using arma::Row<T>::ones;
        using arma::Row<T>::randu;
        using arma::Row<T>::random;
        using arma::Row<T>::randn;
        using arma::Row<T>::fill;
        using arma::Row<T>::replace;
        using arma::Row<T>::clamp;
        using arma::Row<T>::set_size;
        using arma::Row<T>::resize;
        using arma::Row<T>::copy_size;
        using arma::Row<T>::reset;
        using arma::Row<T>::data;
        using arma::Row<T>::print;
        using arma::Row<T>::raw_print;
        using arma::Row<T>::brief_print;        
    };

    
    template<typename T>
    struct Matrix : arma::Mat<T> {
        

        Matrix() = default;
        Matrix(size_t rows, size_t cols) : arma::Mat<T>(rows,cols) {}
        Matrix(const Matrix<T> & m) : arma::Mat<T>(m) {}
        Matrix(const arma::Mat<T> & m) : arma::Mat<T>(m) {}
        //Matrix(const std::vector<T> & v, size_t rows, size_t cols, bool copy_mem = true, bool strict = false) { matrix = arma::Mat<T>(v.data(),rows,cols,copy_mem,strict); }
        Matrix(const std::initializer_list<T> & l) : arma::Mat<T>(l) {}
        ~Matrix() = default;

                        
        using arma::Mat<T>::operator =;
        //using arma::Mat<T>::operator -;
        //using arma::Mat<T>::operator +;
        //using arma::Mat<T>::operator *;
        //using arma::Mat<T>::operator /;

        using arma::Mat<T>::at;
        using arma::Mat<T>::operator ();
        using arma::Mat<T>::operator [];

        using arma::Mat<T>::size;
        using arma::Mat<T>::zeros;
        using arma::Mat<T>::ones;
        using arma::Mat<T>::randu;
        //using arma::Mat<T>::random;
        using arma::Mat<T>::randn;
        using arma::Mat<T>::fill;
        using arma::Mat<T>::replace;
        using arma::Mat<T>::clamp;
        using arma::Mat<T>::set_size;
        using arma::Mat<T>::resize;
        using arma::Mat<T>::copy_size;
        using arma::Mat<T>::reset;
        //using arma::Mat<T>::data;        

        using arma::Mat<T>::rows;
        //using arma::Mat<T>::Mats;
        using arma::Mat<T>::row;
        using arma::Mat<T>::Mat;
        using arma::Mat<T>::as_row;
        //using arma::Mat<T>::as_Mat;
        using arma::Mat<T>::eye;
        using arma::Mat<T>::reshape;
        //using arma::Mat<T>::Matptr;
        using arma::Mat<T>::t;

        using arma::Mat<T>::min;
        using arma::Mat<T>::max;
        using arma::Mat<T>::index_min;
        using arma::Mat<T>::index_max;
        using arma::Mat<T>::eval;
        using arma::Mat<T>::is_empty;
        using arma::Mat<T>::is_trimatu;
        using arma::Mat<T>::is_trimatl;
        using arma::Mat<T>::is_diagmat;
        using arma::Mat<T>::is_square;
        using arma::Mat<T>::is_symmetric;
        using arma::Mat<T>::is_zero;
        using arma::Mat<T>::is_finite;
        using arma::Mat<T>::has_inf;
        using arma::Mat<T>::has_nan;
        
        using arma::Mat<T>::diag;
        using arma::Mat<T>::insert_rows;
        //using arma::Mat<T>::insert_Mats;
        using arma::Mat<T>::shed_row;
        //using arma::Mat<T>::shed_Mat;
        using arma::Mat<T>::shed_rows;
        //using arma::Mat<T>::shed_Mats;
        using arma::Mat<T>::swap_rows;
        //using arma::Mat<T>::swap_Mats;
        using arma::Mat<T>::swap;
        using arma::Mat<T>::submat;
        using arma::Mat<T>::save;
        using arma::Mat<T>::load;                            
        using arma::Mat<T>::print;
        using arma::Mat<T>::raw_print;
        using arma::Mat<T>::brief_print;        
        
        Matrix<T> cwiseProduct(const  Matrix<T> & b)
        {
            Matrix<T> r = arma::cumprod(*this,b);
            return r;
        }
                
    };

    using arma::linspace;
    using arma::logspace;
    using arma::regspace;
    using arma::randperm;
    using arma::eye;

    using arma::operator -;
    using arma::operator +;
    using arma::operator *;
    using arma::operator /;
    
    template<typename T>
    T zero(const T & x) {
        T r(x);
        r.zeros();
        return r;
    }
    template<typename T>
    T ones(const T & x) {
        T r(x);
        r.ones();
        return r;
    }
    template<typename T>
    T random(const T & x) {
        T r(x);
        r.random();
        return r;
    }
    template<typename T>
    T randu(const T & x) {
        T r(x);
        r.randu();
        return r;
    }
    

    template<typename T>
    Vector<T> random_vector(size_t n, T min=(T)0, T max = (T)0) {
        Vector<T> r(n);
        r.random();
        r = min + (max-min)*r;
        return r;
    }

    using arma::toeplitz;
    using arma::circ_toeplitz;
    
    using arma::accu;
    using arma::arg;
    using arma::sum;
    using arma::affmul;   
    using arma::cdot;
    using arma::approx_equal;
    
    using arma::as_scalar;
    using arma::clamp;
    using arma::cond;
    using arma::cumsum;
    using arma::cumprod;
    using arma::det;
    using arma::diagmat;
    using arma::diagvec;
    using arma::diff;
    using arma::dot;
    using arma::norm_dot;
    using arma::inplace_trans;
    using arma::min;
    using arma::max;
    using arma::normalise;
    using arma::prod;
    using arma::powmat;
    using arma::rank;
    using arma::repelem;
    using arma::repmat;
    using arma::reshape;
    using arma::reverse;
    using arma::shift;
    using arma::shuffle;
    using arma::size;
    using arma::sort;
    using arma::symmatu;
    using arma::trace;
    using arma::trans;
    using arma::trimatu;
    using arma::vectorise;
    using arma::abs;
    using arma::exp;
    using arma::exp2;
    using arma::exp10;
    using arma::expm1;
    using arma::trunc_exp;
    using arma::log;
    using arma::log2;
    using arma::log1p;
    using arma::log10;
    using arma::trunc_log;
    using arma::pow;
    using arma::square;
    using arma::sqrt;
    using arma::floor;
    using arma::ceil;
    using arma::round;
    using arma::trunc;
    using arma::erf;
    using arma::erfc;
    using arma::tgamma;
    using arma::lgamma;
    using arma::sign;
    using arma::cos;
    using arma::sin;
    using arma::tan;
    using arma::acos;
    using arma::asin;
    using arma::atan;
    using arma::cosh;
    using arma::sinh;
    using arma::tanh;
    using arma::acosh;
    using arma::asinh;
    using arma::atanh;
    using arma::atan2;
    using arma::hypot;
    using arma::sinc;
    
    template<typename T>
    void addToEachRow(Matrix<T> & m, const RowVector<T> & v) {
        for(size_t i = 0; i < m.rows(); i++)
            m.matrix.row(i) += v.vector;
    }


    template<typename T>
    void addToEachCol(Matrix<T> & m, const Vector<T> & v) {
        for(size_t i = 0; i < m.cols(); i++)
            m.matrix.col(i) += v.vector;
    }

// complex    
    //using arma::st;
    
    using arma::real;
    using arma::imag;
    using arma::log_det;
    using arma::log_det_sympd;
    using arma::logmat_sympd;
    using arma::logmat;
    using arma::sqrtmat;
    using arma::sqrtmat_sympd;
    using arma::chol;
    using arma::eig_sym;
    using arma::eig_gen;
    using arma::eig_pair;
    using arma::hess;
    using arma::inv;
    using arma::inv_sympd;
    using arma::lu;
    using arma::null;
    using arma::orth;
    using arma::pinv;
    using arma::qr;
    using arma::qr_econ;
    using arma::schur;
    using arma::solve;
    using arma::svd;
    using arma::conv;
    using arma::conv2;
    using arma::fft;
    using arma::ifft;
    using arma::fft2;
    using arma::ifft2;
    using arma::interp1;
    using arma::interp2;
    using arma::polyfit;
    using arma::polyval;
    using arma::mean;
    using arma::median;
    using arma::stddev;
    using arma::var;
    using arma::range;
    using arma::cov;
    using arma::cor;
    using arma::hist;
    using arma::histc;
    using arma::quantile;
    using arma::princomp;    
};
