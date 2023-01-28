#pragma once

namespace kfr3
{
    template<typename T>
    struct complex : public kfr::complex<T>
    {
        constexpr complex()  = default;
        constexpr complex(T r)  : kfr::complex<T>(r,0) {}
        constexpr complex(T r, T i)  : kfr::complex<T>(r,i) {}
        constexpr complex(const complex&)  = default;
        constexpr complex(const std::complex<T> & c) : kfr::complex<T>(c.real(),c.imag()) {}
        
        using kfr::complex<T>::real;
        using kfr::complex<T>::imag;
        using kfr::complex<T>::operator =;       

        operator std::complex<T>() {
            return std::complex<T>(real(),imag());
        } 
    };

    template<typename T> complex<T> csin(const complex<T> & v) { return kfr::csin(v); }
    template<typename T> complex<T> ccos(const complex<T> & v) { return kfr::ccos(v); }
    // it should probably be static-cast
    template<typename T> kfr::univector<T> ctan(const kfr::univector<complex<T>> & v) { return std::tan((std::complex<T>)(v)); }

    template<typename T> complex<T> csinh(const complex<T> & v) { return kfr::csinh(v); }
    template<typename T> complex<T> ccosh(const complex<T> & v) { return kfr::ccosh(v); }
    template<typename T> kfr::univector<T> ctanh(const kfr::univector<complex<T>> & v) { return std::tanh((std::complex<T>)(v)); }

    template<typename T> complex<T> cabssqr(const complex<T> & v) { return kfr::cabssqr(v); }
    template<typename T> complex<T> cabs(const complex<T> & v) { return kfr::cabs(v); }
    template<typename T> complex<T> carg(const complex<T> & v) { return kfr::carg(v); }
    
    template<typename T> complex<T> clog(const complex<T> & v) { return kfr::clog(v); }
    template<typename T> complex<T> clog2(const complex<T> & v) { return kfr::clog2(v); }
    template<typename T> complex<T> clog10(const complex<T> & v) { return kfr::clog10(v); }

    template<typename T> complex<T> cexp(const complex<T> & v) { return kfr::cexp(v); }
    template<typename T> complex<T> cexp2(const complex<T> & v) { return kfr::cexp2(v); }
    template<typename T> complex<T> cexp10(const complex<T> & v) { return kfr::cexp10(v); }

    template<typename T> complex<T> polar(const complex<T> & v) { return kfr::polar(v); }
    template<typename T> complex<T> cartesian(const complex<T> & v) { return kfr::cartesian(v); }
    //template<typename T> kfr::univector<T> cabsdup(const kfr::univector<complex<T>> & v) { return kfr::cabsdup(v); }

    template<typename T> complex<T> csqrt(const complex<T> & v) { return kfr::csqrt(v); }
    template<typename T> complex<T> csqr(const complex<T> & v) { return kfr::csqr(v); }
    template<typename T> complex<T> pow(const complex<T> & a, const complex<T> &b)
    {
        return kfr::pow(a,b);
    }

    template<typename T>
    struct sample_vector : public kfr::univector<T>
    {
        using base = kfr::univector<T>;
        using kfr::univector<T>::size;
        using kfr::univector<T>::resize;
        using kfr::univector<T>::data;
        using kfr::univector<T>::operator [];
        using kfr::univector<T>::operator =;
        
    
        sample_vector() = default;
        sample_vector(size_t n) : base(n) {}
        
        operator std::vector<T>() {
            std::vector<T> r(size());
            memcpy(r.data(),data(),size()*sizeof(T));
        }
        T& operator()(size_t i) { return (*this)[i]; }
                
    };

    template<typename X> X rol(X x, X y) { return kfr::rol(x,y); }
    template<typename X> X ror(X x, X y) { return kfr::ror(x,y); }
    template<typename X> X shl(X x, X y) { return kfr::shl(x,y); }
    template<typename X> X shr(X x, X y) { return kfr::rol(x,y); }

    template<typename T> sample_vector<T> bitwiseand(const sample_vector<T> & a, const sample_vector<T> & b) { sample_vector<T> r; r = kfr::bitwiseand(a,b); return r; }
    template<typename T> sample_vector<T> bitwiseandnot(const sample_vector<T> & a,const sample_vector<T> & b) { sample_vector<T> r; r = kfr::bitwiseandnot(a,b); return r; }
    template<typename T> sample_vector<T> bitwisenot(const sample_vector<T> & a) { sample_vector<T> r; r = kfr::bitwisenot(a); return r; }
    template<typename T> sample_vector<T> bitwiseor(const sample_vector<T> & a,const sample_vector<T> & b) { sample_vector<T> r; r = kfr::bitwiseor(a,b); return r; }
    template<typename T> sample_vector<T> bitwisexor(const sample_vector<T> & a,const sample_vector<T> & b) { sample_vector<T> r; r = kfr::bitwisexor(a,b); return r; }

    template<typename T>    
    sample_vector<T> linspace(T start, T stop, size_t size, bool endpoint=false,bool trunc=false)
    {   
        sample_vector<T> r; 
        r = kfr::linspace(start,stop,size,endpoint,trunc); 
        return r; 
    }

    template<typename T>    
    sample_vector<T> pad(const sample_vector<T> & in, const T & fill_value = T(0))
    { 
        sample_vector<T> r; 
        r = kfr::padded(in,fill_value); 
        return r; 
    }

    template<typename T>    
    sample_vector<T> slice(const sample_vector<T> & v, size_t start, size_t end=kfr::max_size_t)
    {   
        sample_vector<T> r;        
        r = v.slice(start,end);
        return r;
    }

    template<typename T>    
    sample_vector<T> truncate(const sample_vector<T> & v, size_t size)
    {   
        sample_vector<T> r; 
        r = v.truncate();
        return r;
    }

    template<typename T>    
    sample_vector<T> reverse(const sample_vector<T> & v)
    {   
        sample_vector<T> r;         
        r = kfr::reverse(v);
        return r;
    }


    template<typename T>    
    T& ringbuf_read(sample_vector<T> &v,size_t & cursor, T& value) { v.ringbuf_read(cursor,value); return value; }

    template<typename T>    
    void ringbuf_write(sample_vector<T> &v, size_t & cursor, T& value) { v.ringbuf_write(cursor,value); }
    
    template<typename T> sample_vector<T> abs(const sample_vector<T>& v) { return kfr::abs(v); }
    template<typename T> sample_vector<T> add(const sample_vector<T> & a,const sample_vector<T> & b) { sample_vector<T> r; r = kfr::add(a,b); return r; }
    template<typename T> sample_vector<T> add(const sample_vector<T> & a,const T & b) { sample_vector<T> r; r = kfr::add(a,b); return r; }
    template<typename T> sample_vector<T> absmax(const sample_vector<T> & a,const sample_vector<T> & b) { sample_vector<T> r; r = kfr::absmax(a,b); return r; }
    template<typename T> sample_vector<T> absmax(const sample_vector<T> & a,const T & b) { sample_vector<T> r; r = kfr::absmax(a,b); return r; }
    template<typename T> sample_vector<T> absmin(const sample_vector<T> & a,const sample_vector<T> & b) { sample_vector<T> r; r = kfr::absmin(a,b); return r; }
    template<typename T> sample_vector<T> absmin(const sample_vector<T> & a,const T & b) { sample_vector<T> r; r = kfr::absmin(a,b); return r; }    
    template<typename T> sample_vector<T> clamp(const sample_vector<T> & a,const sample_vector<T> & lo, const sample_vector<T> &hi) { sample_vector<T> r; r = kfr::clamp(a,lo,hi); return r; }
    template<typename T> sample_vector<T> clamp(const sample_vector<T> & a,const T& lo, const T &hi) { sample_vector<T> r; r = kfr::clamp(a,lo,hi); return r; }
    template<typename T> sample_vector<T> cube(const sample_vector<T> & a) { sample_vector<T> r; r = kfr::cub(a); return r; }
    template<typename T> sample_vector<T> div(const sample_vector<T> & a,const sample_vector<T> & b) { sample_vector<T> r; r = kfr::div(a,b); return r; }    
    template<typename T> sample_vector<T> fmadd(const sample_vector<T> & a,const sample_vector<T> & y, const sample_vector<T> & z) { sample_vector<T> r; r = kfr::fmadd(a,y,z); return r; }
    template<typename T> sample_vector<T> fmsub(const sample_vector<T> & a,const sample_vector<T> & y, const sample_vector<T> & z) { sample_vector<T> r; r = kfr::fmsub(a,y,z); return r; }    
    template<typename T> sample_vector<T> max(const sample_vector<T> & a,const sample_vector<T> & b) { sample_vector<T> r; r = kfr::max(a,b); return r; }
    template<typename T> sample_vector<T> max(const sample_vector<T> & a, const T & b) { sample_vector<T> r; r = kfr::max(a,b); return r; }
    template<typename T> sample_vector<T> min(const sample_vector<T> & a, const sample_vector<T> & b) { sample_vector<T> r; r = kfr::min(a,b); return r; }
    template<typename T> sample_vector<T> min(const sample_vector<T> & a, const T & b) { sample_vector<T> r; r = kfr::min(a,b); return r; }
    template<typename T> sample_vector<T> mix(const sample_vector<T> & a, const T& c, const sample_vector<T> & y) { sample_vector<T> r; r = kfr::mix(c,a,y); return r; }
    template<typename T> sample_vector<T> mixs(const sample_vector<T> & a, const T& c, const sample_vector<T> & y) { sample_vector<T> r; r = kfr::mixs(c,a,y); return r; }
    template<typename T> sample_vector<T> mul(const sample_vector<T> & a,const sample_vector<T> & b) { sample_vector<T> r; r = kfr::mul(a,b); return r; }
    template<typename T> sample_vector<T> mul(const sample_vector<T> & a, const T & b) { sample_vector<T> r; r = kfr::mul(a,b); return r; }
    template<typename T> sample_vector<T> neg(const sample_vector<T> & a) { sample_vector<T> r; r = kfr::neg(a); return r; }        
    template<typename T> sample_vector<T> sqr(const sample_vector<T> & a) { sample_vector<T> r; r = kfr::sqr(a); return r; }
    template<typename T> sample_vector<T> sqrt(const sample_vector<T> & a) { sample_vector<T> r; r = kfr::sqrt(a); return r; }
    template<typename T> sample_vector<T> exp(const sample_vector<T> & a) { sample_vector<T> r; r = kfr::exp(a); return r; }
    template<typename T> sample_vector<T> exp10(const sample_vector<T> & a) { sample_vector<T> r; r = kfr::exp10(a); return r; }
    template<typename T> sample_vector<T> exp2(const sample_vector<T> & a) { sample_vector<T> r; r = kfr::exp2(a); return r; }
    template<typename T> sample_vector<T> exp_fmadd(const sample_vector<T> & a,const sample_vector<T> & y, const sample_vector<T> & z) { sample_vector<T> r; r = kfr::exp_fmadd(a,y,z); return r; }
    template<typename T> sample_vector<T> log(const sample_vector<T> & a) { sample_vector<T> r; r = kfr::log(a); return r; }
    template<typename T> sample_vector<T> log10(const sample_vector<T> & a) { sample_vector<T> r; r = kfr::log10(a); return r; }
    template<typename T> sample_vector<T> log2(const sample_vector<T> & a) { sample_vector<T> r; r = kfr::log2(a); return r; }
    template<typename T> sample_vector<T> log_fmadd(const sample_vector<T> & a,const sample_vector<T> & y, const sample_vector<T> & z) { sample_vector<T> r; r = kfr::log_fmadd(a,y,z); return r; }
    template<typename T> sample_vector<T> logb(const sample_vector<T> & a) { sample_vector<T> r; r = kfr::logb(a); return r; }
    template<typename T> sample_vector<T> logm(const sample_vector<T> & a,const sample_vector<T> & b) { sample_vector<T> r; r = kfr::logm(a,b); return r; }
    template<typename T> sample_vector<T> logn(const sample_vector<T> & a,const sample_vector<T> & b) { sample_vector<T> r; r = kfr::logn(a,b); return r; }
    template<typename T> sample_vector<T> pow(const sample_vector<T> & a,const T & y) { sample_vector<T> r; r = kfr::pow(a,y); return r; }
    template<typename T> sample_vector<T> pow(const sample_vector<T> & a,const sample_vector<T> & y) { sample_vector<T> r; r = kfr::pow(a,y); return r; }
    template<typename T> sample_vector<T> root(const sample_vector<T> & a,const sample_vector<T> & y) { sample_vector<T> r; r = kfr::root(a,y); return r; }
    template<typename T> sample_vector<T> floor(const sample_vector<T> & a) { sample_vector<T> r; r = kfr::floor(a); return r; }        
    template<typename T> sample_vector<T> acos(const sample_vector<T> & a) { sample_vector<T> r; r = kfr::acos(a); return r; }
    template<typename T> sample_vector<T> asin(const sample_vector<T> & a) { sample_vector<T> r; r = kfr::asin(a); return r; }
    template<typename T> sample_vector<T> atan(const sample_vector<T> & a) { sample_vector<T> r; r = kfr::atan(a); return r; }
    template<typename T> sample_vector<T> atan2(const sample_vector<T> & a,const T & b) { sample_vector<T> r; r = kfr::atan2(a,b); return r; }
    template<typename T> sample_vector<T> atan2(const sample_vector<T> & a,const sample_vector<T> & b) { sample_vector<T> r; r = kfr::atan2(a,b); return r; }
    template<typename T> sample_vector<T> atan2deg(const sample_vector<T> & a,const T & b) { sample_vector<T> r; r = kfr::atan2deg(a,b); return r; }
    template<typename T> sample_vector<T> atan2deg(const sample_vector<T> & a,const sample_vector<T> & b) { sample_vector<T> r; r = kfr::atan2deg(a,b); return r; }
    template<typename T> sample_vector<T> atandeg(const sample_vector<T> & a) { sample_vector<T> r; r = kfr::atandeg(a); return r; }
    template<typename T> sample_vector<T> cos(const sample_vector<T> & a) { sample_vector<T> r; r = kfr::cos(a); return r; }
    template<typename T> sample_vector<T> sin(const sample_vector<T> & a) { sample_vector<T> r; r = kfr::sin(a); return r; }    
    template<typename T> sample_vector<T> cosdeg(const sample_vector<T> & a) { sample_vector<T> r; r = kfr::cosdeg(a); return r; }        
    template<typename T> sample_vector<T> sindeg(const sample_vector<T> & a) { sample_vector<T> r; r = kfr::cosdeg(a); return r; }    
    template<typename T> sample_vector<T> sinc(const sample_vector<T> & a) { sample_vector<T> r; r = kfr::sinc(a); return r; }
    template<typename T> sample_vector<T> tan(const sample_vector<T> & a) { sample_vector<T> r; r = kfr::tan(a); return r; }        
    template<typename T> sample_vector<T> cosh(const sample_vector<T> & a) { sample_vector<T> r; r = kfr::cosh(a); return r; }
    template<typename T> sample_vector<T> coth(const sample_vector<T> & a) { sample_vector<T> r; r = kfr::coth(a); return r; }    
    template<typename T> sample_vector<T> sinh(const sample_vector<T> & a) { sample_vector<T> r; r = kfr::sinh(a); return r; }    
    template<typename T> sample_vector<T> tanh(const sample_vector<T> & a) { sample_vector<T> r; r = kfr::tanh(a); return r; }
    template<typename T> sample_vector<T> gamma(const sample_vector<T> & a) { sample_vector<T> r; r = kfr::gamma(a); return r; }

    template<typename T> T absmaxof(const sample_vector<T> & a) { return kfr::absmaxof(a); }
    template<typename T> T absminof(const sample_vector<T> & a) { return kfr::absminof(a); }
    template<typename T> T dot(const sample_vector<T> & a,const sample_vector<T> & b) { return kfr::dotproduct(a,b); }
    template<typename T> T maxof(const sample_vector<T> & a) { return kfr::maxof(a); }
    template<typename T> T minof(const sample_vector<T> & a) { return kfr::minof(a); }
    template<typename T> T mean(const sample_vector<T> & a) { return kfr::mean(a); }
    template<typename T> T product(const sample_vector<T> & a) { return kfr::product(a); }
    template<typename T> T rms(const sample_vector<T> & a) { return kfr::rms(a); }
    template<typename T> T sum(const sample_vector<T> & a) { return kfr::sum(a); }
    template<typename T> T sumsqr(const sample_vector<T> & a) { return kfr::sumsqr(a); }

    // doesn't compile
    //template<typename T>    
    //sample_vector<T> div(const sample_vector<T> & a,const T b) { sample_vector<T> r; r = kfr::div<T>(a,b); return r; }

    template<typename T>    
    sample_vector<T> ipow(const sample_vector<T> & v, int base) { sample_vector<T> r; r = kfr::ipow(v,base); return r; }

    template<typename T>    
    T kcos2x(const T s, const T c) {return kfr::cos2x<T>(s,c); }

    template<typename T>    
    T kcos3x(const T & s, const T & c) {return kfr::cos3x(s,c); }

    template<typename T>    
    T ksin2x(const T & s, const T & c) {return kfr::sin2x(s,c); }

    template<typename T>    
    T ksin3x(const T & s, const T & c) {return kfr::sin3x(s,c); }

    template<typename T>    
    sample_vector<T> cossin(const sample_vector<T> & v) { sample_vector<T> r; r = kfr::cossin(v); return r; }

    template<typename T>    
    sample_vector<T> sincos(const sample_vector<T> & v) 
    { 
        sample_vector<T> r; 
        for(size_t i = 0; i < v.size(); i++) r[i] = kfr::sincos(v[i]);
        return r; }

    template<typename T>    
    T kreciprocal(const T & v) { return kfr::reciprocal(v); }

    template<typename T>    
    T rem(const T v,const T b) { return kfr::rem(v,b); }    

    template<typename T>    
    T satadd(const T v,const T y) { return kfr::satadd(v,y); }

    template<typename T>    
    T satsub(const T v,const T  y) { return kfr::satsub(v,y); }

    //? dont know how to make these work yet.
    template<typename T>    
    sample_vector<T> fastcos(const sample_vector<T> & v) { 
        sample_vector<T> r; 
        for(size_t i = 0; i < v.size(); i++) r[i] = kfr::fastcos(v[i]);
        return r; }

    template<typename T>    
    sample_vector<T> fastcosdeg(const sample_vector<T> & v) { sample_vector<T> r; r = kfr::fastcosdeg(v); return r; }

    template<typename T>    
    sample_vector<T> fastsin(const sample_vector<T> & v) { sample_vector<T> r; r = kfr::fastsin(v); return r; }

    template<typename T>    
    sample_vector<T> fastsindeg(const sample_vector<T> & v) { sample_vector<T> r; r = kfr::fastsindeg(v); return r; }        

    template<typename T>    
    sample_vector<T> coshsinh(const sample_vector<T> & v) { sample_vector<T> r; r = kfr::coshsinh(v); return r; }

    template<typename T>    
    sample_vector<T> sinhcosh(const sample_vector<T> & v) { sample_vector<T> r; r = kfr::sinhcosh(v); return r; }

    template<typename T>    
    sample_vector<T> cossindeg(const sample_vector<T> & v) { sample_vector<T> r; r = kfr::cossindeg(v); return r; }    

    template<typename T>    
    sample_vector<T> sincosdeg(const sample_vector<T> & v) { sample_vector<T> r; r = kfr::sincosdeg(v); return r; }    

    // I dont understand the kfr random at all yet
    template<typename T>    
    sample_vector<T> random(size_t s) 
    {
        std::default_random_engine generator;
        std::uniform_real_distribution<T> distrib(0.0,1.0);
        sample_vector<T> r(s);    
        for(size_t i = 0; i < s; i++)
            r[i] = distrib(generator);
        return r;
    }   

    template<typename T>    
    sample_vector<T> random(size_t s, T min, T max) 
    {
        std::default_random_engine generator;
        std::uniform_real_distribution<T> distrib(min,max);
        sample_vector<T> r(s);    
        for(size_t i = 0; i < s; i++)
            r[i] = distrib(generator);
        return r;
    }  

    #if PYTHON_IS_INSTALLED
    template<typename T>
    void plot_save(const sample_vector<T> & v, const std::string& name="", const std::string& options="") {
            kfr::plot_save(name,v,options);
        }

    template<typename T>    
    void plot_show(const sample_vector<T> & v, const std::string& name="", const std::string&  options="") {
        kfr::plot_show(name,v,options);
    }
    #endif

    // ?
    //template<typename T> sample_vector<T> make_univec(const T * data, size_t s) { return sample_vector<T>(kfr::make_univector<T>(data,s));  }    
    template<typename T>
    struct sample_matrix : public sample_vector<T>
    {        
        size_t M,N;

        sample_matrix() = default;
        sample_matrix(size_t m, size_t n) : sample_vector<T>(m*n) { M = m; N = n; }
        
        using base = sample_vector<T>;
        using kfr::univector<T>::size;
        using kfr::univector<T>::resize;
        using kfr::univector<T>::data;
        using kfr::univector<T>::operator [];

        sample_matrix<T>& operator = (const sample_matrix<T>& v) {
            v = v.v;
            M = v.M;
            N = v.N;
            return *this;
        }
        void resize(size_t i, size_t j) {
            this->resize(i*j);
            M = i;
            N = j;
        }
        
        T& operator()(size_t i, size_t j) { return (*this)[i*N + j]; }        
        T operator()(size_t i, size_t j) const { return (*this)[i*N + j]; }       
    };

    template<typename T>
    sample_vector<T> make_univec(kfr::univector<T> & r) {        
        sample_vector<T> x(r);        
        return x;
    }

    template<typename T>
    sample_vector<T> make_window_hann(size_t s) {
        return make_univec(kfr::univector<T>(kfr::window_hann<T>(s)));
    }
    template<typename T>
    sample_vector<T> make_window_hamming(size_t s) {
        return make_univec(kfr::univector<T>(kfr::window_hamming<T>(s)));
    }
    template<typename T>
    sample_vector<T> make_window_blackman(size_t s, const T alpha=T(0.16), kfr::window_symmetry symmetry = kfr::window_symmetry::symmetric) {
        return make_univec(kfr::univector<T>(kfr::window_blackman<T>(s,alpha,symmetry)));
    }
    template<typename T>
    sample_vector<T> make_window_blackman_harris(size_t s, kfr::window_symmetry symmetry = kfr::window_symmetry::symmetric) {
        return make_univec(kfr::univector<T>(kfr::window_blackman_harris<T>(s,symmetry)));
    }
    template<typename T>
    sample_vector<T> make_window_gaussian(size_t s, const T alpha=T(0.25)) {
        return make_univec(kfr::univector<T>(kfr::window_gaussian<T>(s,alpha)));
    }
    template<typename T>
    sample_vector<T> make_window_triangular(size_t s) {
        return make_univec(kfr::univector<T>(kfr::window_triangular<T>(s)));
    }
    template<typename T>
    sample_vector<T> make_window_bartlett(size_t s) {
        return make_univec(kfr::univector<T>(kfr::window_bartlett<T>(s)));
    }
    template<typename T>
    sample_vector<T> make_window_cosine(size_t s) {
        return make_univec(kfr::univector<T>(kfr::window_cosine<T>(s)));
    }
    template<typename T>
    sample_vector<T> make_window_bartlett_hann(size_t s) {
        return make_univec(kfr::univector<T>(kfr::window_bartlett_hann<T>(s)));
    }
    template<typename T>
    sample_vector<T> make_window_bohman(size_t s) {
        return make_univec(kfr::univector<T>(kfr::window_bohman<T>(s)));
    }
    template<typename T>
    sample_vector<T> make_window_lanczos(size_t s) {
        return make_univec(kfr::univector<T>(kfr::window_lanczos<T>(s)));
    }
    template<typename T>
    sample_vector<T> make_window_flattop(size_t s) {
        return make_univec(kfr::univector<T>(kfr::window_flattop<T>(s)));
    }
    template<typename T>
    sample_vector<T> make_window_rectangular(size_t s) {
        return make_univec(kfr::univector<T>(kfr::window_rectangular<T>(s)));
    }
    template<typename T>
    sample_vector<T> make_window_kaiser(size_t s, const T beta = T(0.5)) {
        return make_univec(kfr::univector<T>(kfr::window_kaiser<T>(s,beta)));
    }

    template<typename T>
    T energy_to_loudness(T energy) {
        return kfr::energy_to_loudness(energy);
    }
    template<typename T>
    T loudness_to_energy(T loudness) {
        return kfr::loudness_to_energy(loudness);
    }    
    template<typename T> T normalize_frequency(T f, T sample_rate) {
        return f/sample_rate;
    }
    
    template<typename T>
    T amp_to_dB(const T & in) {
        return kfr::amp_to_dB(in);
    }        

}