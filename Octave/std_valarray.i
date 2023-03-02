// very basic resizeable array with slice
%{
#include <valarray>
#include <algorithm>
%}

namespace std 
{
    
    class gslice
    {
        gslice();
        gslice(size_t,const valarray<size_t>&,const valarray<size_t>&);
        gslice(const gslice &);
        
        size_t start() const;
        size_t size() const;
        size_t stride() const;
    };

    
    class slice
    {
        slice();
        slice(size_t,const valarray<size_t>&,const valarray<size_t>&);
        slice(const slice &);
        
        size_t start() const;
        size_t size() const;
        size_t stride() const;
    };

    template <class T> class gslice_array {
    public:
        typedef T value_type;
        void operator=   (const valarray<T>&) const;
        void operator*=  (const valarray<T>&) const;
        void operator/=  (const valarray<T>&) const;
        void operator%=  (const valarray<T>&) const;
        void operator+=  (const valarray<T>&) const;
        void operator-=  (const valarray<T>&) const;
        void operator^=  (const valarray<T>&) const;
        void operator&=  (const valarray<T>&) const;
        void operator|=  (const valarray<T>&) const;
        void operator<<= (const valarray<T>&) const;
        void operator>>= (const valarray<T>&) const;
        void operator=(const T&) const;
        
        gslice_array(const gslice_array&);
        const gslice_array& operator= (const gslice_array&) const;
        ~gslice_array();
    };
    template <class T> class indirect_array {
    public:
        typedef T value_type;
        void operator=   (const valarray<T>&) const;
        void operator*=  (const valarray<T>&) const;
        void operator/=  (const valarray<T>&) const;
        void operator%=  (const valarray<T>&) const;
        void operator+=  (const valarray<T>&) const;
        void operator-=  (const valarray<T>&) const;
        void operator^=  (const valarray<T>&) const;
        void operator&=  (const valarray<T>&) const;
        void operator|=  (const valarray<T>&) const;
        void operator<<= (const valarray<T>&) const;
        void operator>>= (const valarray<T>&) const;
        void operator=(const T&);

        indirect_array(const gslice_array&);
        const indirect_array& operator= (const indirect_array&) const;
        ~indirect_array();
    };

        
    template <class T> class mask_array {
    public:
        typedef T value_type;
        void operator=   (const valarray<T>&) const;
        void operator*=  (const valarray<T>&) const;
        void operator/=  (const valarray<T>&) const;
        void operator%=  (const valarray<T>&) const;
        void operator+=  (const valarray<T>&) const;
        void operator-=  (const valarray<T>&) const;
        void operator^=  (const valarray<T>&) const;
        void operator&=  (const valarray<T>&) const;
        void operator|=  (const valarray<T>&) const;
        void operator<<= (const valarray<T>&) const;
        void operator>>= (const valarray<T>&) const;
        void operator=(const T&) const;

        mask_array(const mask_array&);
        const mask_array& operator= (const mask_array&) const;
        ~mask_array();
    };

    template <class T> class slice_array {
    public:
        typedef T value_type;
        void operator=   (const valarray<T>&) const;
        void operator*=  (const valarray<T>&) const;
        void operator/=  (const valarray<T>&) const;
        void operator%=  (const valarray<T>&) const;
        void operator+=  (const valarray<T>&) const;
        void operator-=  (const valarray<T>&) const;
        void operator^=  (const valarray<T>&) const;
        void operator&=  (const valarray<T>&) const;
        void operator|=  (const valarray<T>&) const;
        void operator<<= (const valarray<T>&) const;
        void operator>>= (const valarray<T>&) const;
        void operator=(const T&) const;

        slice_array(const slice_array&);
        const slice_array& operator= (const slice_array&) const;
        ~slice_array();
    };

        

    template<typename T>
    class valarray {
    public:
        valarray();
        valarray(size_t count);
        valarray(const T & val, size_t count);
        valarray(const T* vals, size_t count);
        valarray(const valarray& other);
        valarray(const std::slice_array<T> &sa);
        valarray(const std::gslice_array<T> &ga);
        valarray(const std::mask_array<T> & ma);
        valarray(const std::indirect_array<T> & ia);
    
        valarray<T>& operator =  (const valarray<T> & other);

        %extend {

            // lua starts at 1 like fortran
            T       __getitem__(size_t i) { return (*$self)[i-1]; }
            void    __setitem__(size_t i, const T& v) { (*$self)[i-1] = v; }

            valarray       __getitem__(const gslice& i) { return (*$self)[i]; }
            void    __setitem__(const gslice& i, const T& v) { (*$self)[i] = v; }

            valarray       __getitem__(const slice& i) { return (*$self)[i]; }
            void    __setitem__(const slice& i, const T& v) { (*$self)[i] = v; }

            valarray<T> __add__(const valarray<T> & b) { return *$self + b; }
            valarray<T> __sub__(const valarray<T> & b) { return *$self - b; }
            valarray<T> __mul__(const valarray<T> & b) { return *$self * b; }      
            valarray<T> __div__(const valarray<T> & b) { return *$self / b; }
            valarray<T> __unm__(const valarray<T> & b) { return -*$self; }
            
            valarray<T> __pow__(const valarray<T> & b) { return std::pow(*$self,b); }
            valarray<T> __pow__(const T & b) { return std::pow(*$self,b); }

            //bool __eq__(const valarray<T> & b) { return *$self == b; }
            //bool __le__(const valarray<T> & b) { return *$self <= b; }
            //bool __lt__(const valarray<T> & b) { return *$self < b; }

        }

        size_t size() const;
        void resize(size_t count, T value = T() );
        T sum() const;
        T min() const;
        T max() const;

        valarray<T> shift(int count) const;
        valarray<T> cshift(int count) const;
        valarray<T> apply(T func(T) ) const;
        valarray<T> apply(T func(const T&) ) const;
    };

    template<typename T> swap( valarray<T> & lhs, valarray<T> & rhs);
    template<typename T> sort( valarray<T> & a);
    template<typename T> copy( valarray<T> & dst, valarray<T> & src);
    template<typename T> random_shuffle(valarray<T> & o);

    template<typename T> valarray<T> abs(const valarray<T> & va);
    template<typename T> valarray<T> exp(const valarray<T> & va);
    template<typename T> valarray<T> log(const valarray<T> & va);        
    template<typename T> valarray<T> log10(const valarray<T> & va);

    template<typename T> valarray<T> sqrt(const valarray<T> & va);
    
    template<typename T> valarray<T> sin(const valarray<T> & va);
    template<typename T> valarray<T> cos(const valarray<T> & va);
    template<typename T> valarray<T> tan(const valarray<T> & va);

    template<typename T> valarray<T> asin(const valarray<T> & va);
    template<typename T> valarray<T> acos(const valarray<T> & va);
    template<typename T> valarray<T> atan(const valarray<T> & va);

    template<typename T> valarray<T> sinh(const valarray<T> & va);
    template<typename T> valarray<T> cosh(const valarray<T> & va);
    template<typename T> valarray<T> tanh(const valarray<T> & va);

    template<typename T> valarray<T> asinh(const valarray<T> & va);
    template<typename T> valarray<T> acosh(const valarray<T> & va);
    template<typename T> valarray<T> atanh(const valarray<T> & va);

    template<typename T> valarray<T> pow(const valarray<T> & va, const valarray<T> & p);
    template<typename T> valarray<T> pow(const valarray<T> & va, const T& p);
    template<typename T> valarray<T> pow(const T & va, const valarray<T> & p);

    template<typename T> valarray<T> atan2(const valarray<T> & va, const valarray<T> & p);
    template<typename T> valarray<T> atan2(const valarray<T> & va, const T & p);
    template<typename T> valarray<T> atan2(const T & va, const valarray<T> & p);

}