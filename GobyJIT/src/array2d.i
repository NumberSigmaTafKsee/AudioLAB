%module array2d
%{
#include "Eigen"
#include <iostream>
using namespace Eigen;
%}

%inline %{
    template<typename T, int X, int Y>
    struct ArrayView {
        Eigen::Array<T,X,Y> * array;
        size_t row;

        ArrayView(size_t r, Eigen::Array<T,X,Y> *a) {
            array = a;
            row   = r;
        }
        
        T __getitem(size_t c) { return (*array)(row,c); }
        void __setitem(size_t c, T val) { (*array)(row,c) = val; }
    };
%}
namespace Eigen 
{    
    template<class T, int _Rows=Dynamic, int _Cols=Dynamic>
    class Array
    {
    public:

        Array();        
        Array(const Array & a);
        Array(size_t r, size_t c);
        //Array(size_t size, const T & val);
        ~Array();
        
        Array<T,_Rows,_Cols>& operator = (const Array & a);
        
        T& operator()(size_t r, size_t c);
        T& operator[](size_t r, size_t c);

        size_t size() const;
        void   resize(size_t r, size_t c);
        
        void    setRandom();
        void    fill(T v);
        size_t  cols() const;
            
        
        %extend {            
            // lua is 1 index like fortran
            
            ArrayView<T,_Rows,_Cols>      __getitem(size_t i) { return ArrayView<T,_Rows,_Cols>(i,$self); }

            void print()
            {
                std::cout << *$self << std::endl;
            }
        }

            
        Array<T,_Rows,_Cols> operator + (const Array& b);        
        Array<T,_Rows,_Cols> operator - (const Array& b);
        Array<T,_Rows,_Cols> operator * (const Array& b);
        Array<T,_Rows,_Cols> operator / (const Array& b); 
         
    };
}

%inline %{

namespace Ops {
    template<typename T>
    Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> abs(const Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> & array) { return array.abs(); }
    template<typename T>
    Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> inverse(const Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> & array) { return array.inverse();   }
    template<typename T>
    Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> exp(const Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> & array) { return array.exp();   }
    template<typename T>
    Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> log(const Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> & array) { return array.log();   }
    template<typename T>
    Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> log1p(const Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> & array) { return array.log1p();  }
    template<typename T>
    Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> log10(const Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> & array) { return array.log10();  }
    template<typename T>
    Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> pow(const Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> & array,const T& b) { return array.pow(b);  }
    template<typename T>
    Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> sqrt(const Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> & array) { return array.sqrt();  }
    template<typename T>
    Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> rsqrt(const Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> & array) { return array.rsqrt();  }
    template<typename T>
    Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> square(const Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> & array) { return array.square();  }
    template<typename T>
    Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> cube(const Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> & array) { return array.cube();  }
    template<typename T>
    Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> abs2(const Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> & array) { return array.abs2();  }
    template<typename T>
    Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> sin(const Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> & array) { return array.sin();  }
    template<typename T>
    Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> cos(const Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> & array) { return array.cos();  }
    template<typename T>
    Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> tan(const Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> & array) { return array.tan(); }
    template<typename T>
    Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> asin(const Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> & array) { return array.asin();  }
    template<typename T>
    Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> acos(const Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> & array) { return array.acos();  }
    template<typename T>
    Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> atan(const Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> & array) { return array.atan();  }
    template<typename T>
    Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> sinh(const Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> & array) { return array.sinh();  }
    template<typename T>
    Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> cosh(const Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> & array) { return array.cosh();  }
    template<typename T>
    Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> tanh(const Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> & array) { return array.tanh();  }    
    template<typename T>
    Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> ceil(const Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> & array) { return array.ceil();  }
    template<typename T>
    Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> floor(const Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> & array) { return array.floor();  }
    template<typename T>
    Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> round(const Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> & array) { return array.round();  }


    template<typename T> Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> asinh(const Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> & array) { return array.asinh(); }
    template<typename T> Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> acosh(const Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> & array) { return array.acosh(); }
    template<typename T> Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> atanh(const Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> & array) { return array.atanh(); }
    template<typename T> Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> rint(const Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> & array)  { return array.rint();  }

    template<typename T>
    size_t size(Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> & array) { return array.size(); }  
    //template<typename T>
    //void    random(int i, int j, Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> & array) { array.resize(i,j); array.setRandom(); }
    template<typename T>
    void    random(Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> & array) { array.setRandom(); }    
    template<typename T>
    void    fill(Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> & array, T v) { array.fill(v); }
    template<typename T>
    size_t  cols(Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> & array) { return array.cols(); }

    template<typename T>
    Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> operator + (T a, Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> & b) {
        Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> r(b);
        r = a + b;
        return r;
    }
    template<typename T>
    Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> operator - (T a, Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> & b) {
        Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> r(b);
        r = a - b;
        return r;
    }
    template<typename T>
    Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> operator * (T a, Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> & b) {
        Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> r(b);
        r = a * b;
        return r;
    }

    template<typename T>
    Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> scalar_array(size_t n, T v) {
        Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> r(n);
        r.fill(v);
        return r;
    }
    
    template<typename T>
    void resize(Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> & a, size_t r, size_t c) { a.resize(r,c); }
    
}
%}

%template(flt_array) Eigen::Array<float,Eigen::Dynamic,Eigen::Dynamic>;
%template(absf) Ops::abs<float>;
%template(abs2f) Ops::abs2<float>;
%template(inversef) Ops::inverse<float>;
%template(expf) Ops::exp<float>;
%template(logf) Ops::log<float>;
%template(log1pf) Ops::log1p<float>;
%template(log10f) Ops::log10<float>;
%template(powf) Ops::pow<float>;
%template(sqrtf) Ops::sqrt<float>;
%template(rsqrtf) Ops::rsqrt<float>;
%template(square) Ops::square<float>;
%template(cube) Ops::cube<float>;
%template(sinf) Ops::sin<float>;
%template(cosf) Ops::cos<float>;
%template(tanf) Ops::tan<float>;
%template(asinf) Ops::asin<float>;
%template(acosf) Ops::acos<float>;
%template(atanf) Ops::atan<float>;
%template(sinhf) Ops::sinh<float>;
%template(coshf) Ops::cosh<float>;
%template(tanhf) Ops::tanh<float>;
%template(asinhf) Ops::asinh<float>;
%template(acoshf) Ops::acosh<float>;
%template(atanhf) Ops::atanh<float>;
%template(floorf) Ops::floor<float>;
%template(ceilf) Ops::ceil<float>;
%template(roundf) Ops::round<float>;
%template(rintf) Ops::rint<float>;
%template(sizef) Ops::size<float>;
%template(randomf) Ops::random<float>;
%template(fillf) Ops::fill<float>;
%template(colsf) Ops::cols<float>;
%template(resizef) Ops::resize<float>;

