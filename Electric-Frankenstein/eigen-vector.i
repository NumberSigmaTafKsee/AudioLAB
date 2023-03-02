%{
#include <Eigen/Eigen>
#include <iostream>
#include <vector>
%}

namespace Eigen
{
    template<typename T>
    struct Vector : public Eigen::Matrix<T,1,Eigen::Dynamic,Eigen::RowMajor>
    {
        Vector() = default;        
        Vector(size_t n) : Eigen::Matrix<T,1,Eigen::Dynamic,Eigen::RowMajor>(n) {}
        Vector(const Vector<T> & v) { *this = v; }
        Vector(const std::vector<T> & d) {
            resize(d.size());
            memcpy(data(),d.data(),d.size()*sizeof(T));
        }

        using Eigen::Matrix<T,1,Eigen::Dynamic,Eigen::RowMajor>::setLinSpaced;

        using Eigen::Matrix<T,1,Eigen::Dynamic,Eigen::RowMajor>::operator [];
        using Eigen::Matrix<T,1,Eigen::Dynamic,Eigen::RowMajor>::operator ();
        using Eigen::Matrix<T,1,Eigen::Dynamic,Eigen::RowMajor>::operator =;
        using Eigen::Matrix<T,1,Eigen::Dynamic,Eigen::RowMajor>::operator +=;
        using Eigen::Matrix<T,1,Eigen::Dynamic,Eigen::RowMajor>::operator -=;
        using Eigen::Matrix<T,1,Eigen::Dynamic,Eigen::RowMajor>::operator *=;        
        using Eigen::Matrix<T,1,Eigen::Dynamic,Eigen::RowMajor>::operator +;
        using Eigen::Matrix<T,1,Eigen::Dynamic,Eigen::RowMajor>::operator -;
        using Eigen::Matrix<T,1,Eigen::Dynamic,Eigen::RowMajor>::operator *;
        using Eigen::Matrix<T,1,Eigen::Dynamic,Eigen::RowMajor>::operator /;
        using Eigen::Matrix<T,1,Eigen::Dynamic,Eigen::RowMajor>::operator <<;
        //using Eigen::Matrix<T,1,Eigen::Dynamic,Eigen::RowMajor>::operator <;
        //using Eigen::Matrix<T,1,Eigen::Dynamic,Eigen::RowMajor>::operator <=;


        using Eigen::Matrix<T,1,Eigen::Dynamic,Eigen::RowMajor>::size;
        using Eigen::Matrix<T,1,Eigen::Dynamic,Eigen::RowMajor>::resize;
        using Eigen::Matrix<T,1,Eigen::Dynamic,Eigen::RowMajor>::fill;
        using Eigen::Matrix<T,1,Eigen::Dynamic,Eigen::RowMajor>::data;
        using Eigen::Matrix<T,1,Eigen::Dynamic,Eigen::RowMajor>::cols;
        
        using Eigen::Matrix<T,1,Eigen::Dynamic,Eigen::RowMajor>::setZero;
        using Eigen::Matrix<T,1,Eigen::Dynamic,Eigen::RowMajor>::setOnes;
        using Eigen::Matrix<T,1,Eigen::Dynamic,Eigen::RowMajor>::setRandom;
        
        using Eigen::Matrix<T,1,Eigen::Dynamic,Eigen::RowMajor>::sum;
        using Eigen::Matrix<T,1,Eigen::Dynamic,Eigen::RowMajor>::dot;
        using Eigen::Matrix<T,1,Eigen::Dynamic,Eigen::RowMajor>::cross;

        using Eigen::Matrix<T,1,Eigen::Dynamic,Eigen::RowMajor>::real;
        using Eigen::Matrix<T,1,Eigen::Dynamic,Eigen::RowMajor>::imag;
        
        using Eigen::Matrix<T,1,Eigen::Dynamic,Eigen::RowMajor>::cwiseProduct;
        using Eigen::Matrix<T,1,Eigen::Dynamic,Eigen::RowMajor>::cwiseQuotient;
        using Eigen::Matrix<T,1,Eigen::Dynamic,Eigen::RowMajor>::cwiseSqrt;
        using Eigen::Matrix<T,1,Eigen::Dynamic,Eigen::RowMajor>::cwiseInverse;        
        using Eigen::Matrix<T,1,Eigen::Dynamic,Eigen::RowMajor>::cwiseMin;        
        using Eigen::Matrix<T,1,Eigen::Dynamic,Eigen::RowMajor>::cwiseMax; 
        using Eigen::Matrix<T,1,Eigen::Dynamic,Eigen::RowMajor>::cwiseAbs;
        using Eigen::Matrix<T,1,Eigen::Dynamic,Eigen::RowMajor>::cwiseAbs2;
        using Eigen::Matrix<T,1,Eigen::Dynamic,Eigen::RowMajor>::unaryExpr;
        using Eigen::Matrix<T,1,Eigen::Dynamic,Eigen::RowMajor>::array;                                         
    };
}