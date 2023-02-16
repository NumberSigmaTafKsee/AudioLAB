%{
#include "Eigen/Core"
#include <iostream>
%}

%include "std_complex.i"
%include "std_vector.i"

namespace Eigen 
{    

    // this is the base
    // it will have to specialized for the types and classes
    // array    
    // complex
    // integer

    template<class T, int _Rows=Dynamic, int _Cols=Dynamic, int _Z = RowMajor>
    class Array
    {
    public:

        Array();        
        Array(const Array & a);
        Array(size_t size);
        //Array(size_t size, const T & val);
        ~Array();
        
        Array<T,_Rows,_Cols,_Z>& operator = (const Array & a);
        
        T& operator()(size_t i, size_t j);
        T& operator[](size_t i);

        T*      data();
        size_t  size() const;
        void    resize(size_t i, size_t j);        
        void    setRandom();
        void    setZero();
        void    setOnes();
        void    fill(T v);
        size_t  rows() const;
        size_t  cols() const;
            
        
        %extend {            
            // lua is 1 index like fortran
            void    __setitem__(size_t i, T v) { (*$self)[i] = v; }
            T       __getitem__(size_t i) { return (*$self)[i]; }

            void println()
            {
                std::cout << *$self << std::endl;
            }
        }

        Array<T,_Rows,_Cols,_Z> operator - ();    

        Array<T,_Rows,_Cols,_Z> operator + (const Array<T,_Rows,_Cols,_Z>& b);        
        Array<T,_Rows,_Cols,_Z> operator - (const Array<T,_Rows,_Cols,_Z>& b);
        Array<T,_Rows,_Cols,_Z> operator * (const Array<T,_Rows,_Cols,_Z>& b);
        Array<T,_Rows,_Cols,_Z> operator / (const Array<T,_Rows,_Cols,_Z>& b); 

        Array<T,_Rows,_Cols,_Z> operator + (const T& b);        
        Array<T,_Rows,_Cols,_Z> operator - (const T& b);
        Array<T,_Rows,_Cols,_Z> operator * (const T& b);
        Array<T,_Rows,_Cols,_Z> operator / (const T& b); 

        Array<T,_Rows,_Cols,_Z>& operator = (const Array<T,_Rows,_Cols,_Z> & a);
    
        
        void Random(int i,int j);        
        
    };
}

    template<typename T,int X, int Y, int Z>  Eigen::Array<T,X,Y,Z> abs2( Eigen::Array<T,X,Y,Z> & m) { return m.abs2(); }
    template<typename T,int X, int Y, int Z>  Eigen::Array<T,X,Y,Z> inverse( Eigen::Array<T,X,Y,Z> & m) { return m.inverse(); }
    template<typename T,int X, int Y, int Z>  Eigen::Array<T,X,Y,Z> exp( Eigen::Array<T,X,Y,Z> & m) { return m.exp(); }
    template<typename T,int X, int Y, int Z>  Eigen::Array<T,X,Y,Z> log( Eigen::Array<T,X,Y,Z> & m) { return m.log(); }
    template<typename T,int X, int Y, int Z>  Eigen::Array<T,X,Y,Z> log1p( Eigen::Array<T,X,Y,Z> & m) { return m.log1p(); }
    template<typename T,int X, int Y, int Z>  Eigen::Array<T,X,Y,Z> log10( Eigen::Array<T,X,Y,Z> & m) { return m.log10(); }
    template<typename T,int X, int Y, int Z>  Eigen::Array<T,X,Y,Z> pow( Eigen::Array<T,X,Y,Z> & m,  Eigen::Array<T,X,Y,Z> & p) { return m.pow(p); }
    template<typename T,int X, int Y, int Z>  Eigen::Array<T,X,Y,Z> pow( Eigen::Array<T,X,Y,Z> & m,  T p) { return m.pow(p); }
    template<typename T,int X, int Y, int Z>  Eigen::Array<T,X,Y,Z> sqrt( Eigen::Array<T,X,Y,Z> & m) { return m.sqrt(); }
    template<typename T,int X, int Y, int Z>  Eigen::Array<T,X,Y,Z> rsqrt( Eigen::Array<T,X,Y,Z> & m) { return m.rsqrt(); }
    template<typename T,int X, int Y, int Z>  Eigen::Array<T,X,Y,Z> square( Eigen::Array<T,X,Y,Z> & m) { return m.square(); }
    template<typename T,int X, int Y, int Z>  Eigen::Array<T,X,Y,Z> sin( Eigen::Array<T,X,Y,Z> & m) { return m.sin(); }
    template<typename T,int X, int Y, int Z>  Eigen::Array<T,X,Y,Z> cos( Eigen::Array<T,X,Y,Z> & m) { return m.cos(); }
    template<typename T,int X, int Y, int Z>  Eigen::Array<T,X,Y,Z> tan( Eigen::Array<T,X,Y,Z> & m) { return m.tan(); }
    template<typename T,int X, int Y, int Z>  Eigen::Array<T,X,Y,Z> asin( Eigen::Array<T,X,Y,Z> & m) { return m.asin(); }
    template<typename T,int X, int Y, int Z>  Eigen::Array<T,X,Y,Z> acos( Eigen::Array<T,X,Y,Z> & m) { return m.acos(); }
    template<typename T,int X, int Y, int Z>  Eigen::Array<T,X,Y,Z> atan( Eigen::Array<T,X,Y,Z> & m) { return m.atan(); }
    template<typename T,int X, int Y, int Z>  Eigen::Array<T,X,Y,Z> sinh( Eigen::Array<T,X,Y,Z> & m) { return m.sinh(); }
    template<typename T,int X, int Y, int Z>  Eigen::Array<T,X,Y,Z> cosh( Eigen::Array<T,X,Y,Z> & m) { return m.cosh(); }
    template<typename T,int X, int Y, int Z>  Eigen::Array<T,X,Y,Z> tanh( Eigen::Array<T,X,Y,Z> & m) { return m.tanh(); }
    template<typename T,int X, int Y, int Z>  Eigen::Array<T,X,Y,Z> ceil( Eigen::Array<T,X,Y,Z> & m) { return m.ceil(); }
    template<typename T,int X, int Y, int Z>  Eigen::Array<T,X,Y,Z> floor( Eigen::Array<T,X,Y,Z> & m) { return m.floor(); }
    template<typename T,int X, int Y, int Z>  Eigen::Array<T,X,Y,Z> round( Eigen::Array<T,X,Y,Z> & m) { return m.round(); }
    
    template<typename T,int X, int Y, int Z>
    size_t size(Eigen::Array<T,X,Y,Z> & array) { return array.size(); }  
    template<typename T,int X, int Y, int Z>
    void    random(int i, Eigen::Array<T,X,Y,Z> & array) { array.resize(i); array.setRandom(i); }
    template<typename T,int X, int Y, int Z>
    void    random(Eigen::Array<T,X,Y,Z> & array) { array.setRandom(); }    
    template<typename T,int X, int Y, int Z>
    void    fill(Eigen::Array<T,X,Y,Z> & array, T v) { array.fill(v); }
    template<typename T,int X, int Y, int Z>
    size_t  cols(Eigen::Array<T,X,Y,Z> & array) { return array.cols(); }

    template<typename T,int X, int Y, int Z>
    Eigen::Array<T,X,Y,Z> operator + (T a, Eigen::Array<T,X,Y,Z> & b) {
        Eigen::Array<T,X,Y,Z> r(b);
        r = a + b;
        return r;
    }
    template<typename T,int X, int Y, int Z>
    Eigen::Array<T,X,Y,Z> operator - (T a, Eigen::Array<T,X,Y,Z> & b) {
        Eigen::Array<T,X,Y,Z> r(b);
        r = a - b;
        return r;
    }
    template<typename T,int X, int Y, int Z>
    Eigen::Array<T,X,Y,Z> operator * (T a, Eigen::Array<T,X,Y,Z> & b) {
        Eigen::Array<T,X,Y,Z> r(b);
        r = a * b;
        return r;
    }

    template<typename T,int X, int Y, int Z>
    Eigen::Array<T,X,Y,Z> scalar_array(size_t n, T v) {
        Eigen::Array<T,X,Y,Z> r(n);
        r.fill(v);
        return r;
    }

    template<typename T,int X, int Y, int Z>
    void resize(Eigen::Array<T,X,Y,Z> & a, size_t n) { a.resize(n); }    



