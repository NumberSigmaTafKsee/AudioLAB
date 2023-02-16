#define EIGEN_USE_MKL_ALL
//#define EIGEN_USE_BLAS
//#define EIGEN_USE_LAPACKE
#include "Eigen/Core"
//#include "unsupported/Eigen/AutoDiff"
#include <cfloat>
#include <random>
#include <algorithm>
#include <chrono>
#include <iostream>
#include <vector>
#include <functional>
#include <boost/math/differentiation/autodiff.hpp>

template<typename T, typename CoeffType>
struct Arrayoid1 : public Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>
{

};

template<typename T, typename CoeffType>
struct Arrayoid2 : public Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>
{

};

template<typename T, typename CoeffType>
struct Vectoroid : public Eigen::Matrix<T,1,Eigen::Dynamic,Eigen::RowMajor>
{

};

// despite T is there it is always a kind of Eigen::Matrix
template<typename T, typename CoeffType>
struct Matroid : public Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>
{
    using M = Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>;
    Matroid() = default;    
    Matroid(size_t i, size_t j) : M(i,j) {}
    Matroid(const Matroid& m) : M(m) {}
    Matroid(const M & m) : M(m) {}

    using M::operator ();
    using M::operator =;
    using M::operator +=;
    using M::operator -=;
    using M::operator *=;        
    using M::operator +;
    using M::operator -;
    using M::operator *;
    using M::operator /;
    //using M::operator <<;
    //using M::operator <;
    //using M::operator <=;


    using M::size;
    using M::resize;
    using M::fill;
    using M::data;        
    using M::setZero;
    using M::setOnes;
    using M::setRandom;            

    using M::setIdentity;
    using M::head;
    using M::tail;
    using M::segment;
    using M::block;
    using M::row;
    using M::col;
    using M::rows;
    using M::cols;        
    using M::leftCols;
    using M::middleCols;
    using M::rightCols;
    using M::topRows;
    using M::middleRows;
    using M::bottomRows;
    using M::topLeftCorner;
    using M::topRightCorner;
    using M::bottomLeftCorner;
    using M::bottomRightCorner;
    //using M::seq;
    //using M::seqN;
    //using M::A;
    //using M::v;        
    using M::adjoint;
    using M::transpose;
    using M::diagonal;
    using M::eval;
    using M::asDiagonal;        
    
    using M::replicate;
    using M::reshaped;        
    using M::select;
    using M::cwiseProduct;
    using M::cwiseQuotient;
    using M::cwiseSqrt;
    using M::cwiseInverse;        
    using M::cwiseMin;        
    using M::cwiseMax; 
    using M::cwiseAbs;
    using M::cwiseAbs2;
    using M::unaryExpr;
    using M::array;

    using M::minCoeff;
    using M::maxCoeff;
    using M::sum;
    using M::colwise;
    using M::rowwise;
    using M::trace;
    using M::all;
    using M::any;

    using M::norm;
    using M::squaredNorm;
    
    
    using M::real;
    using M::imag;

    using M::ldlt;
    using M::llt;
    using M::lu;
    //using M::qr;
    //using M::svd;
    using M::eigenvalues;   
    
    using M::begin;
    using M::end;

    void print() const { 
        for(size_t i = 0; i < rows(); i++)
            for(size_t j = 0; j < cols(); j++)
                std::cout << (*this)(i,j) << std::endl; 
    }
    
    Matroid& operator = (const CoeffType & c) {
        for(size_t i = 0; i < rows(); i++)
            for(size_t j = 0; j < cols(); j++)
            {
                (*this)(i,j).fill(c);
            }
        return *this;
    }
    Matroid& operator += (const CoeffType & c) {
        
        for(size_t i = 0; i < rows(); i++)
            for(size_t j = 0; j < cols(); j++)
            {
        
                (*this)(i,j).array() += c;
            }
        return *this;
    }
    Matroid& operator -= (const CoeffType & c) {
        
        for(size_t i = 0; i < rows(); i++)
            for(size_t j = 0; j < cols(); j++)
            {
        
                (*this)(i,j) -= c;
            }
        return *this;
    }
    Matroid& operator *= (const CoeffType & c) {
        
        for(size_t i = 0; i < rows(); i++)
            for(size_t j = 0; j < cols(); j++)
            {
        
                (*this)(i,j) *= c;
            }
        return *this;
    }
    Matroid& operator /= (const CoeffType & c) {
        
        for(size_t i = 0; i < rows(); i++)
            for(size_t j = 0; j < cols(); j++)
            {        
                (*this)(i,j) /= c;
            }
        return *this;
    }
};

int main()
{
    Matroid<Eigen::MatrixXf,float> m(3,3);
    Eigen::MatrixXf x(3,3);
    x.setIdentity();
    m.fill(x);
    m += 2.0f;
    m.print();
}