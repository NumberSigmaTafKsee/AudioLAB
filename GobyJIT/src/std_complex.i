%{
#include <complex>
%}

namespace std 
{
    template<typename T>
    class complex {
    public:        
        complex(const T& re = T(), const T& im = T());        
        complex(const complex<T>& other);

        T real() const;
        T imag() const;
        void real(T re);
        void imag(T im);

        complex<T>& operator = (const complex<T> &b) { *$self = b; return *this; }

        %extend {
            complex<T> __add__(const complex<T> &b) { return *$self + b; }
            complex<T> __sub__(const complex<T> &b) { return *$self - b; }
            complex<T> __div__(const complex<T> &b) { return *$self / b; }        
            complex<T> __mul__(const complex<T> &b) { return *$self * b; }
            complex<T> __pow__(const complex<T> &b) { return std::pow(*$self,b); }
            complex<T> __pow__(const T &b) { return std::pow(*$self,b); }
            bool       __eq__(const complex<T> &b) { return *$self == b; }
            //bool       __lt__(const complex<T> &b) { return *$self < b; }
            //bool       __le__(const complex<T> &b) { return *$self <= b; }
        }

    };

    template<typename T> T real(const complex<T> & z);
    template<typename T> T imag(const complex<T> & z);

    template<typename T> T abs(const complex<T> & z);
    template<typename T> T arg(const complex<T> & z);
    template<typename T> T norm(const complex<T> & z);
    template<typename T> complex<T> proj(const complex<T> & z);
    template<typename T> complex<T> polar(const T& r, const T& theta=T());

    template<typename T> complex<T> exp(const complex<T> & z);
    template<typename T> complex<T> log(const complex<T> & z);
    template<typename T> complex<T> log10(const complex<T> & z);
    template<typename T> complex<T> pow(const complex<T>& x, const complex<T> & z);
    //template<typename T> complex<T> pow(const complex<T>& x, const T& z);
    //template<typename T> complex<T> pow(const T& x, const complex<T> & z);

    template<typename T> complex<T> sqrt(const complex<T> & z);
    
    template<typename T> complex<T> sin(const complex<T> & z);
    template<typename T> complex<T> cos(const complex<T> & z);
    template<typename T> complex<T> tan(const complex<T> & z);

    template<typename T> complex<T> asin(const complex<T> & z);
    template<typename T> complex<T> acos(const complex<T> & z);
    template<typename T> complex<T> atan(const complex<T> & z);

    template<typename T> complex<T> sinh(const complex<T> & z);
    template<typename T> complex<T> cosh(const complex<T> & z);
    template<typename T> complex<T> tanh(const complex<T> & z);

    template<typename T> complex<T> asinh(const complex<T> & z);
    template<typename T> complex<T> acosh(const complex<T> & z);
    template<typename T> complex<T> atanh(const complex<T> & z);
};


%template(crealf) std::real<float>;
%template(creald) std::real<double>;
%template(cimagf) std::imag<float>;
%template(cimagd) std::imag<double>;
%template(cabsf) std::abs<float>;
%template(cabsd) std::abs<double>;
%template(cargf) std::arg<float>;
%template(cargd) std::arg<double>;
%template(cnormf) std::norm<float>;
%template(cnormd) std::norm<double>;
%template(cprojf) std::proj<float>;
%template(cprojf) std::proj<double>;
%template(cpolarf) std::polar<float>;
%template(cpolard) std::polar<double>;
%template(cexpf) std::exp<float>;
%template(cexpd) std::exp<double>;
%template(clogf) std::log<float>;
%template(clogd) std::log<double>;
%template(clog10f) std::log10<float>;
%template(clog10d) std::log10<double>;
%template(cpowf) std::pow<float>;
%template(cpowd) std::pow<double>;
%template(csqrtf) std::sqrt<float>;
%template(csqrtd) std::sqrt<double>;
%template(csinf) std::sin<float>;
%template(csind) std::sin<double>;
%template(ccosf) std::cos<float>;
%template(ccosd) std::cos<double>;
%template(ctanf) std::tan<float>;
%template(ctand) std::tan<double>;
%template(casinf) std::asin<float>;
%template(casind) std::asin<double>;
%template(cacosf) std::acos<float>;
%template(cacosd) std::acos<double>;
%template(catanf) std::atan<float>;
%template(catand) std::atan<double>;
%template(csinhf) std::sinh<float>;
%template(csinhd) std::sinh<double>;
%template(ccoshf) std::cosh<float>;
%template(ccoshd) std::cosh<double>;
%template(ctanhf) std::tanh<float>;
%template(ctanhd) std::tanh<double>;
%template(casinhf) std::asinh<float>;
%template(casinhd) std::asinh<double>;
%template(cacoshf) std::acosh<float>;
%template(cacoshd) std::acosh<double>;
%template(catanhf) std::atanh<float>;
%template(catanhd) std::atanh<double>;
