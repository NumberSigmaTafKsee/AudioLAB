
#pragma once
#include<Eigen/Core>
#include<iostream>
#include<vector>

    template<typename T>
    struct Array : public Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>
    {
        Array() = default;
        Array(size_t n) : Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>(n) {}
        Array(const Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor> &a ) :
            Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>(a) {}
        
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::operator [];
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::operator ();
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::operator =;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::operator +=;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::operator -=;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::operator *=;        
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::operator +;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::operator -;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::operator *;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::operator /;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::operator <<;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::operator <;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::operator <=;

        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::data;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::size;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::cols;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::resize;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::fill;        
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::setZero;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::setOnes;
        using Eigen::Array<T,1,Eigen::Dynamic,Eigen::RowMajor>::setRandom;        
        
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
        
        
        Array<T> clamp(T min = (T)0.0, T max = (T)1.0) {
            Array<T> r(*this);
            r = Array<T>(r.cwiseMin(max).cwiseMax(min));
            return r;
        }

        void print() const { std::cout << *this << std::endl; }
    };

    using ArrayXf = Array<float>;
    using ArrayXd = Array<double>;

    template<typename T>
    Array<T> abs(const Array<T> & array) { return Array<T>(array.abs()); }
    template<typename T>
    Array<T> inverse(const Array<T> & array) { return Array<T>(array.inverse());   }
    template<typename T>
    Array<T> exp(const Array<T> & array) { return Array<T>(array.exp());   }
    template<typename T>
    Array<T> log(const Array<T> & array) { return Array<T>(array.log());   }
    template<typename T>
    Array<T> log1p(const Array<T> & array) { return Array<T>(array.log1p());  }
    template<typename T>
    Array<T> log10(const Array<T> & array) { return Array<T>(array.log10());  }
    template<typename T>
    Array<T> pow(const Array<T> & array,const T& b) { return Array<T>(array.pow(b));  }
    template<typename T>
    Array<T> sqrt(const Array<T> & array) { return Array<T>(array.sqrt());  }
    template<typename T>
    Array<T> rsqrt(const Array<T> & array) { return Array<T>(array.rsqrt());  }
    template<typename T>
    Array<T> square(const Array<T> & array) { return Array<T>(array.square());  }
    template<typename T>
    Array<T> cube(const Array<T> & array) { return Array<T>(array.cube());  }
    template<typename T>
    Array<T> abs2(const Array<T> & array) { return Array<T>(array.abs2());  }
    template<typename T>
    Array<T> sin(const Array<T> & array) { return Array<T>(array.sin());  }
    template<typename T>
    Array<T> cos(const Array<T> & array) { return Array<T>(array.cos());  }
    template<typename T>
    Array<T> tan(const Array<T> & array) { return Array<T>(array.tan()); }
    template<typename T>
    Array<T> asin(const Array<T> & array) { return Array<T>(array.asin());  }
    template<typename T>
    Array<T> acos(const Array<T> & array) { return Array<T>(array.acos());  }
    template<typename T>
    Array<T> atan(const Array<T> & array) { return Array<T>(array.atan());  }
    template<typename T>
    Array<T> sinh(const Array<T> & array) { return Array<T>(array.sinh());  }
    template<typename T>
    Array<T> cosh(const Array<T> & array) { return Array<T>(array.cosh());  }
    template<typename T>
    Array<T> tanh(const Array<T> & array) { return Array<T>(array.tanh());  }    
    template<typename T>
    Array<T> ceil(const Array<T> & array) { return Array<T>(array.ceil());  }
    template<typename T>
    Array<T> floor(const Array<T> & array) { return Array<T>(array.floor());  }
    template<typename T>
    Array<T> round(const Array<T> & array) { return Array<T>(array.round());  }


    template<typename T> Array<T> asinh(const Array<T> & array) { return Array<T>(array.asinh()); }
    template<typename T> Array<T> acosh(const Array<T> & array) { return Array<T>(array.acosh()); }
    template<typename T> Array<T> atanh(const Array<T> & array) { return Array<T>(array.atanh()); }
    template<typename T> Array<T> rint(const Array<T> & array)  { return Array<T>(array.rint());  }

    template<typename T>
    size_t size(Array<T> & array) { return array.size(); }  
    template<typename T>
    void    random(int i, Array<T> & array) { array.random(i); }
    template<typename T>
    void    random(Array<T> & array) { array.random(); }    
    template<typename T>
    void    fill(Array<T> & array, T v) { array.fill(v); }
    template<typename T>
    size_t  cols(Array<T> & array) { return array.cols(); }

    template<typename T>
    Array<T> operator + (T a, Array<T> & b) {
        Array<T> r(b);
        r.array = a + b.array;
        return r;
    }
    template<typename T>
    Array<T> operator - (T a, Array<T> & b) {
        Array<T> r(b);
        r.array = a - b.array;
        return r;
    }
    template<typename T>
    Array<T> operator * (T a, Array<T> & b) {
        Array<T> r(b);
        r.array = a * b.array;
        return r;
    }

    template<typename T>
    Array<T> scalar_array(size_t n, T v) {
        Array<T> r(n);
        r.fill(v);
        return r;
    }

    
    template<typename T>
    void resize(Array<T> & a, size_t n) { a.resize(n); }

    template<typename T>
    void println(Array<T> & a) { a.print(); }



    
    using ArrayXXf = Array<float>;
    using ArrayXXd = Array<double>;


    template<typename T> 
    struct Array2D : public Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>
    {
        Array2D() = default;
        Array2D(size_t n) : Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>(n) {}
        Array2D(const Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> &a ) :
            Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>(a) {}

        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::operator [];
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::operator ();
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::operator =;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::operator +=;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::operator -=;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::operator *=;        
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::operator +;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::operator -;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::operator *;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::operator /;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::operator <<;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::operator <;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::operator <=;

        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::data;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::size;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::resize;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::fill;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::rows;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::cols;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::row;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::col;
        
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::setZero;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::setOnes;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::setRandom;        
        
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::matrix;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::inverse;        
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::pow;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::square;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::cube;        
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::sqrt;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::rsqrt;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::exp;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::log;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::log1p;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::log10;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::max;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::min;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::abs;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::abs2;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::sin;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::cos;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::tan;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::asin;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::acos;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::atan;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::sinh;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::cosh;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::tanh;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::asinh;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::acosh;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::atanh;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::ceil;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::floor;
        using Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::round;
                
        Array2D<T> clamp(T min = (T)0.0, T max = (T)1.0) {
            Array2D<T> r(*this);
            r = Array2D<T>(r.cwiseMin(max).cwiseMax(min));
            return r;
        }

        void print() const { std::cout << *this << std::endl; }
    };

    template<typename T>
    Array2D<T> abs(const Array2D<T> & array) { return Array2D<T>(array.abs()); }
    template<typename T>
    Array2D<T> inverse(const Array2D<T> & array) { return Array2D<T>(array.inverse());   }
    template<typename T>
    Array2D<T> exp(const Array2D<T> & array) { return Array2D<T>(array.exp());   }
    template<typename T>
    Array2D<T> log(const Array2D<T> & array) { return Array2D<T>(array.log());   }
    template<typename T>
    Array2D<T> log1p(const Array2D<T> & array) { return Array2D<T>(array.log1p());  }
    template<typename T>
    Array2D<T> log10(const Array2D<T> & array) { return Array2D<T>(array.log10());  }
    template<typename T>
    Array2D<T> pow(const Array2D<T> & array,const T& b) { return Array2D<T>(array.pow(b));  }
    template<typename T>
    Array2D<T> sqrt(const Array2D<T> & array) { return Array2D<T>(array.sqrt());  }
    template<typename T>
    Array2D<T> rsqrt(const Array2D<T> & array) { return Array2D<T>(array.rsqrt());  }
    template<typename T>
    Array2D<T> square(const Array2D<T> & array) { return Array2D<T>(array.square());  }
    template<typename T>
    Array2D<T> cube(const Array2D<T> & array) { return Array2D<T>(array.cube());  }
    template<typename T>
    Array2D<T> abs2(const Array2D<T> & array) { return Array2D<T>(array.abs2());  }
    template<typename T>
    Array2D<T> sin(const Array2D<T> & array) { return Array2D<T>(array.sin());  }
    template<typename T>
    Array2D<T> cos(const Array2D<T> & array) { return Array2D<T>(array.cos());  }
    template<typename T>
    Array2D<T> tan(const Array2D<T> & array) { return Array2D<T>(array.tan()); }
    template<typename T>
    Array2D<T> asin(const Array2D<T> & array) { return Array2D<T>(array.asin());  }
    template<typename T>
    Array2D<T> acos(const Array2D<T> & array) { return Array2D<T>(array.acos());  }
    template<typename T>
    Array2D<T> atan(const Array2D<T> & array) { return Array2D<T>(array.atan());  }
    template<typename T>
    Array2D<T> sinh(const Array2D<T> & array) { return Array2D<T>(array.sinh());  }
    template<typename T>
    Array2D<T> cosh(const Array2D<T> & array) { return Array2D<T>(array.cosh());  }
    template<typename T>
    Array2D<T> tanh(const Array2D<T> & array) { return Array2D<T>(array.tanh());  }
    template<typename T>
    Array2D<T> ceil(const Array2D<T> & array) { return Array2D<T>(array.ceil());  }
    template<typename T>
    Array2D<T> floor(const Array2D<T> & array) { return Array2D<T>(array.floor());  }
    template<typename T>
    Array2D<T> round(const Array2D<T> & array) { return Array2D<T>(array.round());  }

    template<typename T> Array2D<T> asinh(const Array2D<T> & array) { return Array2D<T>(array.asinh()); }
    template<typename T> Array2D<T> acosh(const Array2D<T> & array) { return Array2D<T>(array.acosh()); }
    template<typename T> Array2D<T> atanh(const Array2D<T> & array) { return Array2D<T>(array.atanh()); }
    template<typename T> Array2D<T> rint(const Array2D<T> & array)  { return Array2D<T>(array.rint());  }

    template<typename T>
    size_t size(const Array2D<T> & array) { return array.size(); }  

    template<typename T>
    void println(Array2D<T> & array) { array.print(); }  

    template<typename T>
    size_t rows(const Array2D<T> & array) { return array.rows(); }

    template<typename T>
    size_t cols(const Array2D<T> & array) { return array.cols(); }

    template<typename T>
    void random(int i,int j,Array2D<T> & array) { array.random(i,j); }

    template<typename T>
    void random(Array2D<T> & array) { array.random(); }    

    template<typename T>
    void fill(Array2D<T> & array,T v) { array.fill(v); }
        
    template<typename T>
    void resize(Array2D<T> & array, size_t i, size_t j) { array.resize(i,j); }
        

    template<typename T>
    Array2D<T> operator + (T a, Array2D<T> & b) {
        Array2D<T> r(b);
        r.array = a + b.array;
        return r;
    }
    template<typename T>
    Array2D<T> operator - (T a, Array2D<T> & b) {
        Array2D<T> r(b);
        r.array = a - b.array;
        return r;
    }
    template<typename T>
    Array2D<T> operator * (T a, Array2D<T> & b) {
        Array2D<T> r(b);
        r.array = a * b.array;
        return r;
    }

    template<typename T>
    Array2D<T> scalar_array2d(size_t n, T v) {
        Array2D<T> r(n);
        r.fill(v);
        return r;
    }


    template<typename T> using RowVector = Eigen::Matrix<T,1,Eigen::Dynamic,Eigen::RowMajor>;
    template<typename T> using ColVector = Eigen::Matrix<T,Eigen::Dynamic,1,Eigen::ColMajor>;

    template<typename T, int X=Eigen::Dynamic,int Y=Eigen::Dynamic, int Z=Eigen::RowMajor>
    struct Matrix : public Eigen::Matrix<T,X,Y,Z>
    {
        Matrix() {}
        Matrix(int x) : Eigen::Matrix<T,X,Y,Z>(x,y) {}
        Matrix(int x, int y) : Eigen::Matrix<T,X,Y,Z>(x,y) {}
        Matrix(int x, int y, const T& val) : Eigen::Matrix<T,X,Y,Z>(x,y) { fill(val);  }
        Matrix(const RowVector<T> & v) : Eigen::Matrix<T,X,Y,Z>(v) {}
        Matrix(const ColVector<T> & v) : Eigen::Matrix<T,X,Y,Z>(v) {}
        Matrix(const Matrix<T>& m) : Eigen::Matrix<T,X,Y,Z>(m) {}        
        Matrix(const Array<T> & a) : Eigen::Matrix<T,X,Y,Z>(a) {}        
        Matrix(const Array2D<T> & a) : Eigen::Matrix<T,X,Y,Z>(a) {}        
        Matrix(const std::vector<T> & data, size_t r, size_t c) {       
            resize(r,c);
            for(size_t i = 0; i < r; i++)
                for(size_t j = 0; j < c; j++)
                    (*this)(i,j) = data[i*c + j];        
        }

        using Eigen::Matrix<T,X,Y,Z>::operator ();
        using Eigen::Matrix<T,X,Y,Z>::operator =;
        using Eigen::Matrix<T,X,Y,Z>::operator +=;
        using Eigen::Matrix<T,X,Y,Z>::operator -=;
        using Eigen::Matrix<T,X,Y,Z>::operator *=;        
        using Eigen::Matrix<T,X,Y,Z>::operator +;
        using Eigen::Matrix<T,X,Y,Z>::operator -;
        using Eigen::Matrix<T,X,Y,Z>::operator *;
        using Eigen::Matrix<T,X,Y,Z>::operator /;
        using Eigen::Matrix<T,X,Y,Z>::operator <<;
        //using Eigen::Matrix<T,X,Y,Z>::operator <;
        //using Eigen::Matrix<T,X,Y,Z>::operator <=;


        using Eigen::Matrix<T,X,Y,Z>::size;
        using Eigen::Matrix<T,X,Y,Z>::resize;
        using Eigen::Matrix<T,X,Y,Z>::fill;
        using Eigen::Matrix<T,X,Y,Z>::data;        
        using Eigen::Matrix<T,X,Y,Z>::setZero;
        using Eigen::Matrix<T,X,Y,Z>::setOnes;
        using Eigen::Matrix<T,X,Y,Z>::setRandom;            

        using Eigen::Matrix<T,X,Y,Z>::setIdentity;
        using Eigen::Matrix<T,X,Y,Z>::head;
        using Eigen::Matrix<T,X,Y,Z>::tail;
        using Eigen::Matrix<T,X,Y,Z>::segment;
        using Eigen::Matrix<T,X,Y,Z>::block;
        using Eigen::Matrix<T,X,Y,Z>::row;
        using Eigen::Matrix<T,X,Y,Z>::col;
        using Eigen::Matrix<T,X,Y,Z>::rows;
        using Eigen::Matrix<T,X,Y,Z>::cols;        
        using Eigen::Matrix<T,X,Y,Z>::leftCols;
        using Eigen::Matrix<T,X,Y,Z>::middleCols;
        using Eigen::Matrix<T,X,Y,Z>::rightCols;
        using Eigen::Matrix<T,X,Y,Z>::topRows;
        using Eigen::Matrix<T,X,Y,Z>::middleRows;
        using Eigen::Matrix<T,X,Y,Z>::bottomRows;
        using Eigen::Matrix<T,X,Y,Z>::topLeftCorner;
        using Eigen::Matrix<T,X,Y,Z>::topRightCorner;
        using Eigen::Matrix<T,X,Y,Z>::bottomLeftCorner;
        using Eigen::Matrix<T,X,Y,Z>::bottomRightCorner;
        //using Eigen::Matrix<T,X,Y,Z>::seq;
        //using Eigen::Matrix<T,X,Y,Z>::seqN;
        //using Eigen::Matrix<T,X,Y,Z>::A;
        //using Eigen::Matrix<T,X,Y,Z>::v;        
        using Eigen::Matrix<T,X,Y,Z>::adjoint;
        using Eigen::Matrix<T,X,Y,Z>::transpose;
        using Eigen::Matrix<T,X,Y,Z>::diagonal;
        using Eigen::Matrix<T,X,Y,Z>::eval;
        using Eigen::Matrix<T,X,Y,Z>::asDiagonal;        
        
        using Eigen::Matrix<T,X,Y,Z>::replicate;
        using Eigen::Matrix<T,X,Y,Z>::reshaped;        
        using Eigen::Matrix<T,X,Y,Z>::select;
        using Eigen::Matrix<T,X,Y,Z>::cwiseProduct;
        using Eigen::Matrix<T,X,Y,Z>::cwiseQuotient;
        using Eigen::Matrix<T,X,Y,Z>::cwiseSqrt;
        using Eigen::Matrix<T,X,Y,Z>::cwiseInverse;        
        using Eigen::Matrix<T,X,Y,Z>::cwiseMin;        
        using Eigen::Matrix<T,X,Y,Z>::cwiseMax; 
        using Eigen::Matrix<T,X,Y,Z>::cwiseAbs;
        using Eigen::Matrix<T,X,Y,Z>::cwiseAbs2;
        using Eigen::Matrix<T,X,Y,Z>::unaryExpr;
        using Eigen::Matrix<T,X,Y,Z>::array;

        using Eigen::Matrix<T,X,Y,Z>::minCoeff;
        using Eigen::Matrix<T,X,Y,Z>::maxCoeff;
        using Eigen::Matrix<T,X,Y,Z>::sum;
        using Eigen::Matrix<T,X,Y,Z>::colwise;
        using Eigen::Matrix<T,X,Y,Z>::rowwise;
        using Eigen::Matrix<T,X,Y,Z>::trace;
        using Eigen::Matrix<T,X,Y,Z>::all;
        using Eigen::Matrix<T,X,Y,Z>::any;

        using Eigen::Matrix<T,X,Y,Z>::norm;
        using Eigen::Matrix<T,X,Y,Z>::squaredNorm;
        
        
        using Eigen::Matrix<T,X,Y,Z>::real;
        using Eigen::Matrix<T,X,Y,Z>::imag;

        using Eigen::Matrix<T,X,Y,Z>::ldlt;
        using Eigen::Matrix<T,X,Y,Z>::llt;
        using Eigen::Matrix<T,X,Y,Z>::lu;
        //using Eigen::Matrix<T,X,Y,Z>::qr;
        //using Eigen::Matrix<T,X,Y,Z>::svd;
        using Eigen::Matrix<T,X,Y,Z>::eigenvalues;   
        
        void print() const { std::cout << *this << std::endl; }
        
    };

    template<typename T> using Vector = Matrix<T,1,Eigen::Dynamic,Eigen::RowMajor>;
    using VectorXf = Vector<float>;
    using VectorXd = Vector<double>;

    
    template<typename T> Matrix<T> abs(const Matrix<T> & matrix) { matrix.array().abs().matrix();}
    template<typename T> Matrix<T> inverse(const Matrix<T> & matrix) { matrix.array().inverse().matrix(); }
    template<typename T> Matrix<T> exp(const Matrix<T> & matrix) { matrix.array().exp().matrix();  }
    template<typename T> Matrix<T> log(const Matrix<T> & matrix) { matrix.array().log().matrix();  }
    template<typename T> Matrix<T> log1p(const Matrix<T> & matrix) { matrix.array().log1p().matrix(); }
    template<typename T> Matrix<T> log10(const Matrix<T> & matrix) { matrix.array().log10().matrix(); }
    template<typename T> Matrix<T> pow(const Matrix<T> & matrix,const T& b) { matrix.array().pow(b).matrix();}
    template<typename T> Matrix<T> sqrt(const Matrix<T> & matrix) { matrix.array().sqrt().matrix();}
    template<typename T> Matrix<T> rsqrt(const Matrix<T> & matrix) { matrix.array().rsqrt().matrix(); }
    template<typename T> Matrix<T> square(const Matrix<T> & matrix) { matrix.array().square().matrix();}
    template<typename T> Matrix<T> cube(const Matrix<T> & matrix) { matrix.array().cube().matrix(); }
    template<typename T> Matrix<T> abs2(const Matrix<T> & matrix) { matrix.array().abs2().matrix(); }
    template<typename T> Matrix<T> sin(const Matrix<T> & matrix) { matrix.array().sin().matrix();}
    template<typename T> Matrix<T> cos(const Matrix<T> & matrix) { matrix.array().cos().matrix(); }
    template<typename T> Matrix<T> tan(const Matrix<T> & matrix) { matrix.array().tan().matrix();}
    template<typename T> Matrix<T> asin(const Matrix<T> & matrix) { matrix.array().asin().matrix(); }
    template<typename T> Matrix<T> acos(const Matrix<T> & matrix) { matrix.array().acos().matrix(); }
    template<typename T> Matrix<T> atan(const Matrix<T> & matrix) { matrix.array().atan().matrix(); }
    template<typename T> Matrix<T> sinh(const Matrix<T> & matrix) { matrix.array().sinh().matrix(); }
    template<typename T> Matrix<T> cosh(const Matrix<T> & matrix) { matrix.array().cosh().matrix(); }
    template<typename T> Matrix<T> tanh(const Matrix<T> & matrix) { matrix.array().tanh().matrix(); }
    template<typename T> Matrix<T> ceil(const Matrix<T> & matrix) { matrix.array().ceil().matrix(); }
    template<typename T> Matrix<T> floor(const Matrix<T> & matrix) { matrix.array().floor().matrix(); }
    template<typename T> Matrix<T> round(const Matrix<T> & matrix) { matrix.array().round().matrix(); }
    template<typename T> Matrix<T> asinh(const Matrix<T> & matrix) { matrix.array().asinh().matrix(); }
    template<typename T> Matrix<T> acosh(const Matrix<T> & matrix) { matrix.array().acosh().matrix(); }
    template<typename T> Matrix<T> atanh(const Matrix<T> & matrix) { matrix.array().atanh().matrix(); }
    template<typename T> Matrix<T> rint(const Matrix<T> & matrix) { matrix.array().rint().matrix(); }


    template<typename T> 
    void random(Matrix<T> & matrix, int x, int y) { matrix.random(x,y); }
    template<typename T> 
    void random(Matrix<T> & matrix) { matrix.random(); }
    template<typename T> 
    void identity(Matrix<T> & matrix,int x, int y) { matrix.identity(x,y); }
    template<typename T> 
    void identity(Matrix<T> & matrix) { matrix.identity(); }
    template<typename T> 
    void zero(Matrix<T> & matrix,int x, int y) { matrix.zero(x,y); }
    template<typename T> 
    void zero(Matrix<T> & matrix) { matrix.zero(); }
    template<typename T> 
    void ones(Matrix<T> & matrix,int x, int y) { matrix.ones(x,y); }
    template<typename T> 
    void ones(Matrix<T> & matrix) { matrix.ones(); }

    template<typename T> 
    T get(Matrix<T> & matrix,int i, int j) { return matrix.get(i,j); }

    template<typename T> 
    void set(Matrix<T> & matrix,int i, int j, T v) { matrix.set(i,j,v); }

    template<typename T> 
    T norm(Matrix<T> & matrix) { return matrix.norm(); }

    template<typename T> 
    T squaredNorm(Matrix<T> & matrix) { return matrix.squaredNorm(); }

    template<typename T> 
    bool all(Matrix<T> & matrix) { return matrix.all(); }

    template<typename T> 
    bool allFinite(Matrix<T> & matrix) { return matrix.allFinite(); }

    template<typename T> 
    bool any(Matrix<T> & matrix) { return matrix.any(); }

    template<typename T> 
    bool count(Matrix<T> & matrix) { return matrix.count(); }

    template<typename T> 
    size_t rows(Matrix<T> & matrix) { return matrix.rows(); }

    template<typename T> 
    size_t cols(Matrix<T> & matrix) { return matrix.cols(); }

    template<typename T> 
    void resize(Matrix<T> & matrix,int x, int y) { matrix.resize(x,y); }

    template<typename T> 
    void normalize(Matrix<T> & matrix) { matrix.normalize(); }

    template<typename T> 
    Matrix<T>  normalized(Matrix<T> & matrix) { return matrix.normalized(); }    


    template<typename T> 
    void fill(Matrix<T> & matrix,T v) { matrix.fill(v); }

    template<typename T> 
    Matrix<T>  eval(Matrix<T> & matrix) { return Matrix<T> (matrix.eval()); }

    template<typename T> 
    bool hasNaN(Matrix<T> & matrix) { return matrix.hasNaN(); }

    template<typename T> 
    size_t innerSize(Matrix<T> & matrix) { return matrix.innerSize(); }

    template<typename T> 
    size_t outerSize(Matrix<T> & matrix) { return matrix.outerSize(); }    

    template<typename T> 
    bool isMuchSmallerThan(Matrix<T> & matrix,const Matrix<T> & n, T v) { return matrix.isMuchSmallerThan(n.matrix,v); }

    template<typename T> 
    bool isOnes(Matrix<T> & matrix) { return matrix.isOnes(); }

    template<typename T> 
    bool isZero(Matrix<T> & matrix) { return matrix.isZero(); }

    template<typename T> 
    Matrix<T>  adjoint(Matrix<T> & matrix)  { return matrix.adjoint(); }

    template<typename T> 
    Matrix<T>  transpose(Matrix<T> & matrix) { return matrix.tranpose(); }

    template<typename T> 
    Matrix<T>  diagonal(Matrix<T> & matrix) { return matrix.diagonal(); }        

    template<typename T> 
    Matrix<T>  reverse(Matrix<T> & matrix) { return matrix.revese(); }    

    template<typename T> 
    Matrix<T>  replicate(Matrix<T> & matrix, size_t i, size_t j) { return matrix.replicate(i,j); }
        

    template<typename T> 
    T sum(Matrix<T> & matrix)    {        
        return matrix.sum();        
    }

    template<typename T> 
    T prod(Matrix<T> & matrix)    {        
        return matrix.prod();        
    }

    template<typename T> 
    T mean(Matrix<T> & matrix)    {        
        return matrix.mean();        
    }

    template<typename T> 
    T minCoeff(Matrix<T> & matrix)    {        
        return matrix.minCoeff();        
    }

    template<typename T> 
    T maxCoeff(Matrix<T> & matrix)    {        
        return matrix.maxCoeff();        
    }    

    template<typename T> 
    T trace(Matrix<T> & matrix)    {        
        return matrix.trace();        
    }

    template<typename T> 
    Matrix<T>  addToEachRow(Matrix<T>  & m, Matrix<T>  & v)    {
        Matrix<T>  r(m);        
        r.matrix = r.matrix.rowwise() + RowVector<T>(v.matrix).vector;
    
    }

    template<typename T> 
    Matrix<T>  cwiseAbs(Matrix<T>  & matrix)    {
        Matrix<T>  r = matrix.cwiseAbs();
    
    }

    template<typename T> 
    Matrix<T>  cwiseAbs2(Matrix<T>  & matrix)    {
        Matrix<T>  r = matrix.cwiseAbs2();
    
    }

    template<typename T> 
    Matrix<T>  cwiseProduct(Matrix<T>  & matrix,const Matrix<T> & q)    {
        Matrix<T>  r = matrix.cwiseProduct(q.matrix); 
    
    }

    template<typename T> 
    Matrix<T>  cwiseQuotient(Matrix<T>  & matrix, const Matrix<T> & q)    {
        Matrix<T>  r = matrix.cwiseQuotient(q.matrix); 
    
    }

    template<typename T> 
    Matrix<T>  cwiseInverse(Matrix<T>  & matrix)    {
        Matrix<T>  r = matrix.cwiseInverse();
    
    }

    template<typename T> 
    Matrix<T>  cwiseSqrt(Matrix<T>  & matrix)    {
        Matrix<T>  r = matrix.cwiseSqrt();
    
    }

    template<typename T> 
    Matrix<T>  cwiseMax(Matrix<T>  & matrix, Matrix<T> & q)    {
        Matrix<T>  r = matrix.cwiseMin(q.matrix);
            
    }

    template<typename T> 
    Matrix<T>  cwiseMin(Matrix<T>  & matrix, Matrix<T> & q)    {
        Matrix<T>  r = matrix.cwiseMin(q.matrix);
    
    }


    template<typename T> 
    Matrix<T>  slice(Matrix<T>  & matrix,int first_r,int first_c, int last_r=-1, int last_c=-1)    {
        return matrix.slice(first_r,first_c,last_r,last_c);
    }

    template<typename T> 
    Matrix<T>  sliceN1(Matrix<T>  & matrix,int first_r,int first_rn, int first_c, int last_c=-1)    {        
        return matrix.sliceN1(first_r,first_rn,first_c,last_c);
    }

    template<typename T> 
    Matrix<T>  sliceN2(Matrix<T>  & matrix,int first_r,int first_c, int first_cn, int last_r=-1)    {                
        return matrix.sliceN2(first_r, first_c, first_cn, last_r);
    }

    template<typename T> 
    Matrix<T>  slicen(Matrix<T>  & matrix,int first_r,int first_rn, int first_c, int first_cn)    {        
        return matrix.slicen(first_r,first_rn,first_c,first_cn);
    }

    template<typename T> 
    Array2D <T> array(Matrix<T>  & matrix) { return matrix.array(); }
    

    using MatrixXf = Matrix<float>;
    using MatrixXd = Matrix<double>;

/*
    template<typename T>
    RowVector<T> get_row_vector(Array<T> & array) { return RowVector<T>(array.matrix()); }

    template<typename T>
    ColVector<T> get_col_vector(Array<T> & array) { return ColVector<T>(array.matrix()); }
*/


