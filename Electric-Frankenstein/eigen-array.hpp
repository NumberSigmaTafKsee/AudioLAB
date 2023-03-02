#pragma once
#include <vector>

    // it has to be special for 1d array    
    template<class T>
    class Array : public Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>
    {
    public:

        Array() = default;        
        Array(size_t n) : Eigen::Array<T,1,Eigen::Dynamic>(n) {}
        Array(const Array<T> & v) { *this = v; }
        Array(const std::vector<T> & d) {
            resize(d.size());
            memcpy(data(),d.data(),d.size()*sizeof(T));
        }
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::operator +;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::operator -;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::operator *;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::operator /;
        
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::operator [];
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::operator ();
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::operator =;

        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::data;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::size;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::cols;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::resize;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::fill;        
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::setZero;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::setOnes;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::setRandom;        
        
        /*                                        
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::operator +=;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::operator -=;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::operator *=;        
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::operator <<;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::operator <;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::operator <=;

        
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::matrix;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::inverse;        
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::pow;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::square;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::cube;        
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::sqrt;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::rsqrt;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::exp;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::log;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::log1p;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::log10;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::max;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::min;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::abs;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::abs2;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::sin;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::cos;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::tan;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::asin;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::acos;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::atan;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::sinh;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::cosh;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::tanh;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::asinh;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::acosh;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::atanh;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::ceil;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::floor;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::round;
        */
        
		T __getitem__(size_t i) { return (*this)[i]; }
		void __setitem__(size_t i, const T & val) { (*this)[i] = val; }

    };
