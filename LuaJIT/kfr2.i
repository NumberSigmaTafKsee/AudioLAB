%module kfr2f
%{

typedef float SampleType;

#include <cmath>
#include <vector>
#include <complex>
#include <iostream>
#include <random> 

#include <kfr/kfr.h>
#include <kfr/dft.hpp>
#include <kfr/io.hpp>
#include <kfr/math.hpp>
#include "kfrcore.hpp"

#include "KfrDSP/KfrDsp.hpp"
#include "IIRFilters.hpp"

using namespace KfrDSP1;
using namespace Filters;
%}


%include "stdint.i"
%include "std_complex.i"
%include "std_vector.i"
%include "std_math.i"
%include "std_list.i"
%include "std_map.i"

%include "kfr_biquad.hpp"
%include "kfr_bessel.hpp"
%include "kfr_butterworth.hpp"
%include "kfr_chebyshev1.hpp"
%include "kfr_chebyshev2.hpp"
%include "kfr_convolve.hpp"
%include "kfr_dct.hpp"
%include "kfr_dft.hpp"
%include "kfr_fileio.hpp"
%include "kfr_fir.hpp"
%include "kfr_resample.hpp"
%include "kfr_window.hpp"
%include "kfrcore.hpp"




namespace kfr {
    using b8   = bool;
    using f32  = float;
    using f64  = double;
    using i8   = int8_t;
    using i16  = int16_t;
    using i32  = int32_t;
    using i64  = int64_t;
    using u8   = uint8_t;
    using u16  = uint16_t;
    using u32  = uint32_t;
    using u64  = uint64_t;
    using umax = uint64_t;
    using imax = int64_t;
    using fmax = double;
    using f80  = long double;
    using fbase = SampleType;
}

/*
%inline %{
    typedef int8_t   i8;
    typedef uint8_t   u8;
    typedef int16_t  i16;
    typedef uint16_t  u16;
    typedef int32_t  i32;
    typedef uint32_t u32;
    typedef signed long   i64;
    typedef unsigned long  u64;
    typedef float    f32;
    typedef double   f64;
%};
*/
namespace kfr {

    enum class audio_sample_type
    {
        unknown,
        i8,
        i16,
        i24,
        i32,
        i64,
        f32,
        f64,
        first_float = f32
    };

    enum class sample_rate_conversion_quality : int
    {
        draft   = 4,
        low     = 6,
        normal  = 8,
        high    = 10,
        perfect = 12,
    };

    enum class biquad_type
    {
        lowpass,
        highpass,
        bandpass,
        bandstop,
        peak,
        notch,
        lowshelf,
        highshelf
    };

    enum class Speaker : int
    {
        None          = -1,
        Mono          = 0,
        M             = static_cast<int>(Mono),
        Left          = 1,
        L             = static_cast<int>(Left),
        Right         = 2,
        R             = static_cast<int>(Right),
        Center        = 3,
        C             = static_cast<int>(Center),
        Lfe           = 4,
        Ls            = 5,
        LeftSurround  = static_cast<int>(Ls),
        Rs            = 6,
        RightSurround = static_cast<int>(Rs),
        Lc            = 7,
        Rc            = 8,
        S             = 9,
        Cs            = static_cast<int>(S),
        Sl            = 10,
        Sr            = 11,
        Tm            = 12,
        Tfl           = 13,
        Tfc           = 14,
        Tfr           = 15,
        Trl           = 16,
        Trc           = 17,
        Trr           = 18,
        Lfe2          = 19
    };

    enum class SpeakerArrangement : int
    {
        None           = -1,
        Mono           = 0,
        Stereo         = 1,
        StereoSurround = 2,
        StereoCenter   = 3,
        StereoSide     = 4,
        StereoCLfe     = 5,
        Cine30         = 6,
        Music30        = 7,
        Cine31         = 8,
        Music31        = 9,
        Cine40         = 10,
        Music40        = 11,
        Cine41         = 12,
        Music41        = 13,
        Arr50          = 14,
        Arr51          = 15,
        Cine60         = 16,
        Music60        = 17,
        Cine61         = 18,
        Music61        = 19,
        Cine70         = 20,
        Music70        = 21,
        Cine71         = 22,
        Music71        = 23,
        Cine80         = 24,
        Music80        = 25,
        Cine81         = 26,
        Music81        = 27,
        Arr102         = 28
    };

    /// @brief Seek origin
    enum class seek_origin : int
    {
        current = SEEK_CUR, ///< From the current position
        begin   = SEEK_SET, ///< From the beginning
        end     = SEEK_END, ///< From the end
    };

    struct audio_format
    {
        size_t channels        = 2;
        audio_sample_type type = audio_sample_type::i16;
        kfr::fmax samplerate        = 44100;
        bool use_w64           = false;
    };

    struct audio_format_and_length : audio_format
    {        
        constexpr audio_format_and_length();
        constexpr audio_format_and_length(const audio_format& fmt);

        imax length = 0; // in samples
    };

    constexpr size_t audio_sample_sizeof(audio_sample_type type);
    constexpr size_t audio_sample_bit_depth(audio_sample_type type);

    
    /*
    void CMT_ARCH_NAME::deinterleave(SampleType* out[], const SampleType *in, size_t channels, size_t size);
    void CMT_ARCH_NAME::interleave(SampleType* out, const SampleType* in[], size_t channels, size_t size);
    void CMT_ARCH_NAME::convert(SampleType* out, const SampleTYpe* in, size_t size);
    */
   
    struct fraction
    {
        fraction(i64 num = 0, i64 den = 1);
        void normalize();
        
        i64 numerator;
        i64 denominator;

        fraction operator+() const;
        fraction operator-() const;

        //explicit operator bool() const;
        //explicit operator double() const;
        //explicit operator float() const;
        //explicit operator kfr::signed long long() const;
    };    

    template <typename T>
    struct complex
    {        
        constexpr complex()  = default;
        constexpr complex(T re)  : re(re), im(0) {}
        constexpr complex(T re, T im)  : re(re), im(im) {}
        constexpr complex(const complex&)  = default;
        
        constexpr complex& operator=(const complex&)  = default;
        constexpr complex& operator=(complex&&)  = default;
        constexpr const T& real() const  { return re; }
        constexpr const T& imag() const  { return im; }
        constexpr void real(T value)  { re = value; }
        constexpr void imag(T value)  { im = value; }    
    };
}

%inline %{
    namespace Ops 
    {
        /*
        template<typename T> kfr::univector<T> csin(const kfr::univector<kfr::complex<T>> & v) { return kfr::csin(v); }
        template<typename T> kfr::univector<T> ccos(const kfr::univector<kfr::complex<T>> & v) { return kfr::ccos(v); }
        //template<typename T> kfr::univector<T> ctan(const kfr::univector<kfr::complex<T>> & v) { return kfr::ctan(v); }

        template<typename T> kfr::univector<T> csinh(const kfr::univector<kfr::complex<T>> & v) { return kfr::csinh(v); }
        template<typename T> kfr::univector<T> ccosh(const kfr::univector<kfr::complex<T>> & v) { return kfr::ccosh(v); }
        //template<typename T> kfr::univector<T> ctanh(const kfr::univector<kfr::complex<T>> & v) { return kfr::ctanh(v); }

        template<typename T> kfr::univector<T> cabssqr(const kfr::univector<kfr::complex<T>> & v) { return kfr::cabssqr(v); }
        template<typename T> kfr::univector<T> cabs(const kfr::univector<kfr::complex<T>> & v) { return kfr::cabs(v); }
        template<typename T> kfr::univector<T> carg(const kfr::univector<kfr::complex<T>> & v) { return kfr::carg(v); }
        
        template<typename T> kfr::univector<T> clog(const kfr::univector<kfr::complex<T>> & v) { return kfr::clog(v); }
        template<typename T> kfr::univector<T> clog2(const kfr::univector<kfr::complex<T>> & v) { return kfr::clog2(v); }
        template<typename T> kfr::univector<T> clog10(const kfr::univector<kfr::complex<T>> & v) { return kfr::clog10(v); }

        template<typename T> kfr::univector<T> cexp(const kfr::univector<kfr::complex<T>> & v) { return kfr::cexp(v); }
        template<typename T> kfr::univector<T> cexp2(const kfr::univector<kfr::complex<T>> & v) { return kfr::cexp2(v); }
        template<typename T> kfr::univector<T> cexp10(const kfr::univector<kfr::complex<T>> & v) { return kfr::cexp10(v); }

        template<typename T> kfr::univector<T> polar(const kfr::univector<kfr::complex<T>> & v) { return kfr::polar(v); }
        template<typename T> kfr::univector<T> cartesian(const kfr::univector<kfr::complex<T>> & v) { return kfr::cartesian(v); }
        //template<typename T> kfr::univector<T> cabsdup(const kfr::univector<kfr::complex<T>> & v) { return kfr::cabsdup(v); }

        template<typename T> kfr::univector<T> csqrt(const kfr::univector<kfr::complex<T>> & v) { return kfr::csqrt(v); }
        template<typename T> kfr::univector<T> csqr(const kfr::univector<kfr::complex<T>> & v) { return kfr::csqr(v); }
        */
        template<typename T> kfr::complex<T> csin(const kfr::complex<T> & v) { return kfr::csin(v); }
        template<typename T> kfr::complex<T> ccos(const kfr::complex<T> & v) { return kfr::ccos(v); }
        //template<typename T> kfr::univector<T> ctan(const kfr::univector<kfr::complex<T>> & v) { return kfr::ctan(v); }

        template<typename T> kfr::complex<T> csinh(const kfr::complex<T> & v) { return kfr::csinh(v); }
        template<typename T> kfr::complex<T> ccosh(const kfr::complex<T> & v) { return kfr::ccosh(v); }
        //template<typename T> kfr::univector<T> ctanh(const kfr::univector<kfr::complex<T>> & v) { return kfr::ctanh(v); }

        template<typename T> kfr::complex<T> cabssqr(const kfr::complex<T> & v) { return kfr::cabssqr(v); }
        template<typename T> kfr::complex<T> cabs(const kfr::complex<T> & v) { return kfr::cabs(v); }
        template<typename T> kfr::complex<T> carg(const kfr::complex<T> & v) { return kfr::carg(v); }
        
        template<typename T> kfr::complex<T> clog(const kfr::complex<T> & v) { return kfr::clog(v); }
        template<typename T> kfr::complex<T> clog2(const kfr::complex<T> & v) { return kfr::clog2(v); }
        template<typename T> kfr::complex<T> clog10(const kfr::complex<T> & v) { return kfr::clog10(v); }

        template<typename T> kfr::complex<T> cexp(const kfr::complex<T> & v) { return kfr::cexp(v); }
        template<typename T> kfr::complex<T> cexp2(const kfr::complex<T> & v) { return kfr::cexp2(v); }
        template<typename T> kfr::complex<T> cexp10(const kfr::complex<T> & v) { return kfr::cexp10(v); }

        template<typename T> kfr::complex<T> polar(const kfr::complex<T> & v) { return kfr::polar(v); }
        template<typename T> kfr::complex<T> cartesian(const kfr::complex<T> & v) { return kfr::cartesian(v); }
        //template<typename T> kfr::univector<T> cabsdup(const kfr::univector<kfr::complex<T>> & v) { return kfr::cabsdup(v); }

        template<typename T> kfr::complex<T> csqrt(const kfr::complex<T> & v) { return kfr::csqrt(v); }
        template<typename T> kfr::complex<T> csqr(const kfr::complex<T> & v) { return kfr::csqr(v); }
    }
%}

%template(vectorf32)  std::vector<f32>;
%template(vectorf64)  std::vector<f64>;
%template(vectori8)   std::vector<i8>;
%template(vectorui8)  std::vector<u8>;
%template(vectori16)  std::vector<i16>;
%template(vectorui16) std::vector<u16>;
%template(vectori32)  std::vector<i32>;
%template(vectorui32) std::vector<u32>;
%template(vectori64)  std::vector<i64>;
%template(vectorui64)  std::vector<u64>;

%template(cvector32) std::vector<kfr::complex<float>>;
%template(cvector64) std::vector<kfr::complex<double>>;

%template(complex32) kfr::complex<float>;
%template(complex64) kfr::complex<double>;


namespace kfr
{    
    // single channel vector (SampleVector)
    template<typename T> 
    struct univector
    {        
        univector() {}
        univector(size_t s);
        univector(const univector<T> & u);        
      
        size_t size() const;
        void resize(size_t s);
      
        %extend {
            // lua is 1 based like fortran
            T       __getitem__(size_t i) { assert(i > 0) ; return (*$self)[i-1]; }
            void    __setitem__(size_t i, const T & val) { assert(i > 0);(*$self)[i-1] = val; }

            univector<T> __add__(const univector& b) { return *$self + b; }
            univector<T> __sub__(const univector& b) { return *$self - b; }
            univector<T> __mul__(const univector& b) { return *$self * b; }
            univector<T> __div__(const univector& b) { return *$self / b; }
            univector<T> __unm__() { return -*$self; }
            //univector<SampleType> __pow__(const SampleType& b) { return pow(*$self,b); }
            //univector<double> __pow__(const double& b) { return pow(*$self,b); }
            //bool         __eq__(const univector& b) { return (bool)*$self == b; }
            //bool         __lt__(const univector& b) { return *$self < b; }
            //bool         __le__(const univector& b) { return *$self <= b; }

            void fill(const T& val) { for(size_t i = 0; i < $self->size(); i++) (*$self)[i] = val; }
            void print() const { kfr::println(*$self); }
        }
                
        T& at(size_t pos);
        T& front();
        T& back();
        T* data();
            
        univector<T>& operator = (const univector<T>& u);
        
    };    

    // single channel vector (SampleVector)
    template<typename T> 
    struct univector2d
    {        
        univector2d();
        univector2d(size_t s);        

        size_t size() const;
        void resize(size_t s);
        
        %extend {
            // lua is 1 based like fortran
            univector<T>  __getitem__(size_t i) { assert(i > 0) ; return (*$self)[i-1]; }
            void __setitem__(size_t i, const univector<T> & val) { assert(i > 0);(*$self)[i-1] = val; }
        }
            
        univector2d<T>& operator = (const univector2d<T>& u);
        
    };    
}


%inline %{

    namespace Ops 
    {
        kfr::univector2d<SampleType> deinterleave(const kfr::univector<SampleType> & v) {
            kfr::univector2d<SampleType> r(2);
            r[0].resize(v.size()/2);
            r[1].resize(v.size()/2);
            for(size_t i = 0; i < v.size()/2; i++)
            {
                r[0][i] = v[i*2];
                r[1][i] = v[i*2+1];
            }
            return r;
        }
        kfr::univector<SampleType> interleave(const kfr::univector2d<SampleType> & v) {
            kfr::univector<SampleType> r(v[0].size()*2);
            for(size_t i = 0; i < v[0].size(); i++)
            {
                r[2*i]   = v[0][i];
                r[2*i+1] = v[1][i];
            }
            return r;
        }
        kfr::univector<SampleType> to_univector(const std::vector<SampleType> & v) {
            kfr::univector<SampleType> r(v.size());
            std::copy(v.begin(),v.end(),r.begin());
            return r;
        }

        std::vector<SampleType> to_vector(const kfr::univector<SampleType> & v) {
            std::vector<SampleType> r(v.size());
            std::copy(v.begin(),v.end(),r.begin());
            return r;
        }

        template<typename X> X rol(X x, X y) { return kfr::rol(x,y); }
        template<typename X> X ror(X x, X y) { return kfr::ror(x,y); }
        template<typename X> X shl(X x, X y) { return kfr::shl(x,y); }
        template<typename X> X shr(X x, X y) { return kfr::rol(x,y); }

        template<typename T> kfr::univector<T> bitwiseand(const kfr::univector<T> & a, const kfr::univector<T> & b) { kfr::univector<T> r; r = kfr::bitwiseand(a,b); return r; }
        template<typename T> kfr::univector<T> bitwiseandnot(const kfr::univector<T> & a,const kfr::univector<T> & b) { kfr::univector<T> r; r = kfr::bitwiseandnot(a,b); return r; }
        template<typename T> kfr::univector<T> bitwisenot(const kfr::univector<T> & a) { kfr::univector<T> r; r = kfr::bitwisenot(a); return r; }
        template<typename T> kfr::univector<T> bitwiseor(const kfr::univector<T> & a,const kfr::univector<T> & b) { kfr::univector<T> r; r = kfr::bitwiseor(a,b); return r; }
        template<typename T> kfr::univector<T> bitwisexor(const kfr::univector<T> & a,const kfr::univector<T> & b) { kfr::univector<T> r; r = kfr::bitwisexor(a,b); return r; }

        template<typename T>    
        kfr::univector<T> linspace(T start, T stop, size_t size, bool endpoint=false,bool trunc=false)
        {   
            kfr::univector<T> r; 
            r = kfr::linspace(start,stop,size,endpoint,trunc); 
            return r; 
        }

        template<typename T>    
        kfr::univector<T> pad(const kfr::univector<T> & in, const T & fill_value = T(0))
        { 
            kfr::univector<T> r; 
            r = kfr::padded(in,fill_value); 
            return r; 
        }

        template<typename T>    
        kfr::univector<T> slice(const kfr::univector<T> & v, size_t start, size_t end=kfr::max_size_t)
        {   
            kfr::univector<T> r;        
            r = v.slice(start,end);
            return r;
        }

        template<typename T>    
        kfr::univector<T> truncate(const kfr::univector<T> & v, size_t size)
        {   
            kfr::univector<T> r; 
            r = v.truncate();
            return r;
        }

        template<typename T>    
        kfr::univector<T> reverse(const kfr::univector<T> & v)
        {   
            kfr::univector<T> r;         
            r = kfr::reverse(v);
            return r;
        }


        template<typename T>    
        T& ringbuf_read(kfr::univector<T> &v,size_t & cursor, T& value) { v.ringbuf_read(cursor,value); return value; }

        template<typename T>    
        void ringbuf_write(kfr::univector<T> &v, size_t & cursor, T& value) { v.ringbuf_write(cursor,value); }
        
        template<typename T> kfr::univector<T> abs(const kfr::univector<T>& v) { return kfr::abs(v); }
        template<typename T> kfr::univector<T> add(const kfr::univector<T> & a,const kfr::univector<T> & b) { kfr::univector<T> r; r = kfr::add(a,b); return r; }
        template<typename T> kfr::univector<T> add(const kfr::univector<T> & a,const T & b) { kfr::univector<T> r; r = kfr::add(a,b); return r; }
        template<typename T> kfr::univector<T> absmax(const kfr::univector<T> & a,const kfr::univector<T> & b) { kfr::univector<T> r; r = kfr::absmax(a,b); return r; }
        template<typename T> kfr::univector<T> absmax(const kfr::univector<T> & a,const T & b) { kfr::univector<T> r; r = kfr::absmax(a,b); return r; }
        template<typename T> kfr::univector<T> absmin(const kfr::univector<T> & a,const kfr::univector<T> & b) { kfr::univector<T> r; r = kfr::absmin(a,b); return r; }
        template<typename T> kfr::univector<T> absmin(const kfr::univector<T> & a,const T & b) { kfr::univector<T> r; r = kfr::absmin(a,b); return r; }    
        template<typename T> kfr::univector<T> clamp(const kfr::univector<T> & a,const kfr::univector<T> & lo, const kfr::univector<T> &hi) { kfr::univector<T> r; r = kfr::clamp(a,lo,hi); return r; }
        template<typename T> kfr::univector<T> clamp(const kfr::univector<T> & a,const T& lo, const T &hi) { kfr::univector<T> r; r = kfr::clamp(a,lo,hi); return r; }
        template<typename T> kfr::univector<T> cube(const kfr::univector<T> & a) { kfr::univector<T> r; r = kfr::cub(a); return r; }
        template<typename T> kfr::univector<T> div(const kfr::univector<T> & a,const kfr::univector<T> & b) { kfr::univector<T> r; r = kfr::div(a,b); return r; }    
        template<typename T> kfr::univector<T> fmadd(const kfr::univector<T> & a,const kfr::univector<T> & y, const kfr::univector<T> & z) { kfr::univector<T> r; r = kfr::fmadd(a,y,z); return r; }
        template<typename T> kfr::univector<T> fmsub(const kfr::univector<T> & a,const kfr::univector<T> & y, const kfr::univector<T> & z) { kfr::univector<T> r; r = kfr::fmsub(a,y,z); return r; }    
        template<typename T> kfr::univector<T> max(const kfr::univector<T> & a,const kfr::univector<T> & b) { kfr::univector<T> r; r = kfr::max(a,b); return r; }
        template<typename T> kfr::univector<T> max(const kfr::univector<T> & a, const T & b) { kfr::univector<T> r; r = kfr::max(a,b); return r; }
        template<typename T> kfr::univector<T> min(const kfr::univector<T> & a, const kfr::univector<T> & b) { kfr::univector<T> r; r = kfr::min(a,b); return r; }
        template<typename T> kfr::univector<T> min(const kfr::univector<T> & a, const T & b) { kfr::univector<T> r; r = kfr::min(a,b); return r; }
        template<typename T> kfr::univector<T> mix(const kfr::univector<T> & a, const T& c, const kfr::univector<T> & y) { kfr::univector<T> r; r = kfr::mix(c,a,y); return r; }
        template<typename T> kfr::univector<T> mixs(const kfr::univector<T> & a, const T& c, const kfr::univector<T> & y) { kfr::univector<T> r; r = kfr::mixs(c,a,y); return r; }
        template<typename T> kfr::univector<T> mul(const kfr::univector<T> & a,const kfr::univector<T> & b) { kfr::univector<T> r; r = kfr::mul(a,b); return r; }
        template<typename T> kfr::univector<T> mul(const kfr::univector<T> & a, const T & b) { kfr::univector<T> r; r = kfr::mul(a,b); return r; }
        template<typename T> kfr::univector<T> neg(const kfr::univector<T> & a) { kfr::univector<T> r; r = kfr::neg(a); return r; }        
        template<typename T> kfr::univector<T> sqr(const kfr::univector<T> & a) { kfr::univector<T> r; r = kfr::sqr(a); return r; }
        template<typename T> kfr::univector<T> sqrt(const kfr::univector<T> & a) { kfr::univector<T> r; r = kfr::sqrt(a); return r; }
        template<typename T> kfr::univector<T> exp(const kfr::univector<T> & a) { kfr::univector<T> r; r = kfr::exp(a); return r; }
        template<typename T> kfr::univector<T> exp10(const kfr::univector<T> & a) { kfr::univector<T> r; r = kfr::exp10(a); return r; }
        template<typename T> kfr::univector<T> exp2(const kfr::univector<T> & a) { kfr::univector<T> r; r = kfr::exp2(a); return r; }
        template<typename T> kfr::univector<T> exp_fmadd(const kfr::univector<T> & a,const kfr::univector<T> & y, const kfr::univector<T> & z) { kfr::univector<T> r; r = kfr::exp_fmadd(a,y,z); return r; }
        template<typename T> kfr::univector<T> log(const kfr::univector<T> & a) { kfr::univector<T> r; r = kfr::log(a); return r; }
        template<typename T> kfr::univector<T> log10(const kfr::univector<T> & a) { kfr::univector<T> r; r = kfr::log10(a); return r; }
        template<typename T> kfr::univector<T> log2(const kfr::univector<T> & a) { kfr::univector<T> r; r = kfr::log2(a); return r; }
        template<typename T> kfr::univector<T> log_fmadd(const kfr::univector<T> & a,const kfr::univector<T> & y, const kfr::univector<T> & z) { kfr::univector<T> r; r = kfr::log_fmadd(a,y,z); return r; }
        template<typename T> kfr::univector<T> logb(const kfr::univector<T> & a) { kfr::univector<T> r; r = kfr::logb(a); return r; }
        template<typename T> kfr::univector<T> logm(const kfr::univector<T> & a,const kfr::univector<T> & b) { kfr::univector<T> r; r = kfr::logm(a,b); return r; }
        template<typename T> kfr::univector<T> logn(const kfr::univector<T> & a,const kfr::univector<T> & b) { kfr::univector<T> r; r = kfr::logn(a,b); return r; }
        template<typename T> kfr::univector<T> pow(const kfr::univector<T> & a,const T & y) { kfr::univector<T> r; r = kfr::pow(a,y); return r; }
        template<typename T> kfr::univector<T> pow(const kfr::univector<T> & a,const kfr::univector<T> & y) { kfr::univector<T> r; r = kfr::pow(a,y); return r; }
        template<typename T> kfr::univector<T> root(const kfr::univector<T> & a,const kfr::univector<T> & y) { kfr::univector<T> r; r = kfr::root(a,y); return r; }
        template<typename T> kfr::univector<T> floor(const kfr::univector<T> & a) { kfr::univector<T> r; r = kfr::floor(a); return r; }        
        template<typename T> kfr::univector<T> acos(const kfr::univector<T> & a) { kfr::univector<T> r; r = kfr::acos(a); return r; }
        template<typename T> kfr::univector<T> asin(const kfr::univector<T> & a) { kfr::univector<T> r; r = kfr::asin(a); return r; }
        template<typename T> kfr::univector<T> atan(const kfr::univector<T> & a) { kfr::univector<T> r; r = kfr::atan(a); return r; }
        template<typename T> kfr::univector<T> atan2(const kfr::univector<T> & a,const T & b) { kfr::univector<T> r; r = kfr::atan2(a,b); return r; }
        template<typename T> kfr::univector<T> atan2(const kfr::univector<T> & a,const kfr::univector<T> & b) { kfr::univector<T> r; r = kfr::atan2(a,b); return r; }
        template<typename T> kfr::univector<T> atan2deg(const kfr::univector<T> & a,const T & b) { kfr::univector<T> r; r = kfr::atan2deg(a,b); return r; }
        template<typename T> kfr::univector<T> atan2deg(const kfr::univector<T> & a,const kfr::univector<T> & b) { kfr::univector<T> r; r = kfr::atan2deg(a,b); return r; }
        template<typename T> kfr::univector<T> atandeg(const kfr::univector<T> & a) { kfr::univector<T> r; r = kfr::atandeg(a); return r; }
        template<typename T> kfr::univector<T> cos(const kfr::univector<T> & a) { kfr::univector<T> r; r = kfr::cos(a); return r; }
        template<typename T> kfr::univector<T> sin(const kfr::univector<T> & a) { kfr::univector<T> r; r = kfr::sin(a); return r; }    
        template<typename T> kfr::univector<T> cosdeg(const kfr::univector<T> & a) { kfr::univector<T> r; r = kfr::cosdeg(a); return r; }        
        template<typename T> kfr::univector<T> sindeg(const kfr::univector<T> & a) { kfr::univector<T> r; r = kfr::cosdeg(a); return r; }    
        template<typename T> kfr::univector<T> sinc(const kfr::univector<T> & a) { kfr::univector<T> r; r = kfr::sinc(a); return r; }
        template<typename T> kfr::univector<T> tan(const kfr::univector<T> & a) { kfr::univector<T> r; r = kfr::tan(a); return r; }        
        template<typename T> kfr::univector<T> cosh(const kfr::univector<T> & a) { kfr::univector<T> r; r = kfr::cosh(a); return r; }
        template<typename T> kfr::univector<T> coth(const kfr::univector<T> & a) { kfr::univector<T> r; r = kfr::coth(a); return r; }    
        template<typename T> kfr::univector<T> sinh(const kfr::univector<T> & a) { kfr::univector<T> r; r = kfr::sinh(a); return r; }    
        template<typename T> kfr::univector<T> tanh(const kfr::univector<T> & a) { kfr::univector<T> r; r = kfr::tanh(a); return r; }
        template<typename T> kfr::univector<T> gamma(const kfr::univector<T> & a) { kfr::univector<T> r; r = kfr::gamma(a); return r; }

        template<typename T> T absmaxof(const kfr::univector<T> & a) { return kfr::absmaxof(a); }
        template<typename T> T absminof(const kfr::univector<T> & a) { return kfr::absminof(a); }
        template<typename T> T dot(const kfr::univector<T> & a,const kfr::univector<T> & b) { return kfr::dotproduct(a,b); }
        template<typename T> T maxof(const kfr::univector<T> & a) { return kfr::maxof(a); }
        template<typename T> T minof(const kfr::univector<T> & a) { return kfr::minof(a); }
        template<typename T> T mean(const kfr::univector<T> & a) { return kfr::mean(a); }
        template<typename T> T product(const kfr::univector<T> & a) { return kfr::product(a); }
        template<typename T> T rms(const kfr::univector<T> & a) { return kfr::rms(a); }
        template<typename T> T sum(const kfr::univector<T> & a) { return kfr::sum(a); }
        template<typename T> T sumsqr(const kfr::univector<T> & a) { return kfr::sumsqr(a); }

        // doesn't compile
        //template<typename T>    
        //kfr::univector<T> div(const kfr::univector<T> & a,const T b) { kfr::univector<T> r; r = kfr::div<T>(a,b); return r; }

        template<typename T>    
        kfr::univector<T> ipow(const kfr::univector<T> & v, int base) { kfr::univector<T> r; r = kfr::ipow(v,base); return r; }

        template<typename T>    
        T kcos2x(const T s, const T c) {return kfr::cos2x<SampleType>(s,c); }

        template<typename T>    
        T kcos3x(const T & s, const T & c) {return kfr::cos3x(s,c); }

        template<typename T>    
        T ksin2x(const T & s, const T & c) {return kfr::sin2x(s,c); }

        template<typename T>    
        T ksin3x(const T & s, const T & c) {return kfr::sin3x(s,c); }

        template<typename T>    
        kfr::univector<T> cossin(const kfr::univector<T> & v) { kfr::univector<T> r; r = kfr::cossin(v); return r; }

        template<typename T>    
        kfr::univector<T> sincos(const kfr::univector<T> & v) 
        { 
            kfr::univector<T> r; 
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
        kfr::univector<T> fastcos(const kfr::univector<T> & v) { 
            kfr::univector<T> r; 
            for(size_t i = 0; i < v.size(); i++) r[i] = kfr::fastcos(v[i]);
            return r; }

        template<typename T>    
        kfr::univector<T> fastcosdeg(const kfr::univector<T> & v) { kfr::univector<T> r; r = kfr::fastcosdeg(v); return r; }

        template<typename T>    
        kfr::univector<T> fastsin(const kfr::univector<T> & v) { kfr::univector<T> r; r = kfr::fastsin(v); return r; }

        template<typename T>    
        kfr::univector<T> fastsindeg(const kfr::univector<T> & v) { kfr::univector<T> r; r = kfr::fastsindeg(v); return r; }        

        template<typename T>    
        kfr::univector<T> coshsinh(const kfr::univector<T> & v) { kfr::univector<T> r; r = kfr::coshsinh(v); return r; }

        template<typename T>    
        kfr::univector<T> sinhcosh(const kfr::univector<T> & v) { kfr::univector<T> r; r = kfr::sinhcosh(v); return r; }

        template<typename T>    
        kfr::univector<T> cossindeg(const kfr::univector<T> & v) { kfr::univector<T> r; r = kfr::cossindeg(v); return r; }    

        template<typename T>    
        kfr::univector<T> sincosdeg(const kfr::univector<T> & v) { kfr::univector<T> r; r = kfr::sincosdeg(v); return r; }    

        // I dont understand the kfr random at all yet
        template<typename T>    
        kfr::univector<T> random(size_t s) 
        {
            std::default_random_engine generator;
            std::uniform_real_distribution<T> distrib(0.0,1.0);
            kfr::univector<T> r(s);    
            for(size_t i = 0; i < s; i++)
                r[i] = distrib(generator);
            return r;
        }   

        template<typename T>    
        kfr::univector<T> random(size_t s, T min, T max) 
        {
            std::default_random_engine generator;
            std::uniform_real_distribution<T> distrib(min,max);
            kfr::univector<T> r(s);    
            for(size_t i = 0; i < s; i++)
                r[i] = distrib(generator);
            return r;
        }  

        template<typename T>
        void plot_save(const kfr::univector<T> & v, const std::string& name="", const std::string& options="") {
                kfr::plot_save(name,v,options);
            }

        template<typename T>    
        void plot_show(const kfr::univector<T> & v, const std::string& name="", const std::string&  options="") {
            kfr::plot_show(name,v,options);
        }
    
        // ?
        //template<typename T> kfr::univector<T> make_univec(const T * data, size_t s) { return kfr::univector<T>(kfr::make_univector<T>(data,s));  }    
    }
%}

%template(csin)  Ops::csin<SampleType>;
%template(ccos)  Ops::ccos<SampleType>;
%template(csinh) Ops::csinh<SampleType>;
%template(ccosh) Ops::ccosh<SampleType>;
%template(cabssqr) Ops::cabssqr<SampleType>;
%template(cabs) Ops::cabs<SampleType>;
%template(carg) Ops::carg<SampleType>;
%template(clog) Ops::clog<SampleType>;
%template(clog10) Ops::clog10<SampleType>;
%template(clog2) Ops::clog2<SampleType>;
%template(cexp) Ops::cexp<SampleType>;
%template(cexp2) Ops::cexp2<SampleType>;
%template(cexp10) Ops::cexp10<SampleType>;
%template(cpolar) Ops::polar<SampleType>;
%template(ccartesian) Ops::cartesian<SampleType>;
%template(csqrt) Ops::csqrt<SampleType>;
%template(csqr) Ops::csqr<SampleType>;

%template(abs) Ops::abs<SampleType>;
%template(add) Ops::add<SampleType>;
%template(absmax) Ops::absmax<SampleType>;
%template(absmin) Ops::absmin<SampleType>;
%template(clamp) Ops::clamp<SampleType>;
%template(cube) Ops::cube<SampleType>;
%template(div)  Ops::div<SampleType>;
%template(fmadd) Ops::fmadd<SampleType>;
%template(fmsub) Ops::fmsub<SampleType>;
%template(max) Ops::max<SampleType>;
%template(min) Ops::min<SampleType>;
%template(mix) Ops::mix<SampleType>;
%template(mixs) Ops::mixs<SampleType>;
%template(mul) Ops::mul<SampleType>;
%template(neg) Ops::neg<SampleType>;
%template(sqr) Ops::sqr<SampleType>;
%template(sqrt) Ops::sqrt<SampleType>;
%template(exp) Ops::exp<SampleType>;
%template(exp10) Ops::exp10<SampleType>;
%template(exp2) Ops::exp2<SampleType>;
%template(exp_fmadd) Ops::exp_fmadd<SampleType>;
%template(log) Ops::log<SampleType>;
%template(log10) Ops::log10<SampleType>;
%template(log2) Ops::log2<SampleType>;
%template(log_fmadd) Ops::log_fmadd<SampleType>;
%template(logb) Ops::logb<SampleType>;
%template(logm) Ops::logm<SampleType>;
%template(logn) Ops::logn<SampleType>;
%template(pow) Ops::pow<SampleType>;
%template(root) Ops::root<SampleType>;
%template(floor) Ops::floor<SampleType>;
%template(acos) Ops::acos<SampleType>;
%template(asin) Ops::asin<SampleType>;
%template(atan) Ops::atan<SampleType>;
%template(atan2) Ops::atan2<SampleType>;
%template(cos) Ops::cos<SampleType>;
%template(sin) Ops::sin<SampleType>;
%template(tan) Ops::tan<SampleType>;
%template(cosh) Ops::cosh<SampleType>;
%template(coth) Ops::coth<SampleType>;
%template(sinh) Ops::sinh<SampleType>;
%template(tanh) Ops::tanh<SampleType>;
%template(atan2deg) Ops::atan2deg<SampleType>;
%template(cosdeg) Ops::cosdeg<SampleType>;
%template(sindeg) Ops::sindeg<SampleType>;
%template(sinc) Ops::sinc<SampleType>;
%template(gamma) Ops::gamma<SampleType>;
%template(absmaxo) Ops::absmaxof<SampleType>;
%template(dot) Ops::dot<SampleType>;
%template(maxo) Ops::maxof<SampleType>;
%template(mino) Ops::minof<SampleType>;
%template(mean) Ops::mean<SampleType>;
%template(prdocut) Ops::product<SampleType>;
%template(rms) Ops::rms<SampleType>;
%template(sum) Ops::sum<SampleType>;
%template(sumsqr) Ops::sumsqr<SampleType>;
%template(ipow) Ops::ipow<SampleType>;
%template(cos2x) Ops::kcos2x<SampleType>;
%template(sin2x) Ops::ksin2x<SampleType>;
%template(cos3x) Ops::kcos3x<SampleType>;
%template(sin3x) Ops::ksin3x<SampleType>;
%template(reciprocal) Ops::kreciprocal<SampleType>;

%template(linspace) Ops::linspace<SampleType>;
%template(pad)      Ops::pad<SampleType>;
%template(slice)    Ops::slice<SampleType>;
%template(truncate) Ops::truncate<SampleType>;
%template(reverse)  Ops::reverse<SampleType>;
%template(ringbuf_read)  Ops::ringbuf_read<SampleType>;
%template(ringbuf_write) Ops::ringbuf_write<SampleType>;
%template (random)       Ops::random<SampleType>;

// Plot
%template (plot_save) Ops::plot_save<SampleType>;
%template (plot_show) Ops::plot_show<SampleType>;

/*
%template(acosh) acosh<SampleType>;
%template(asinh) asinh<SampleType>;
%template(atanh) atanh<SampleType>;
%template(cbrt) cbrt<SampleType>;
%template(ceil) ceil<SampleType>;
%template(copysign) copysign<SampleType>;
%template(er) erf<SampleType>;
%template(erfc) erfc<SampleType>;
%template(expm1) expm1<SampleType>;
%template(fdim) fdim<SampleType>;
%template(fma) fma<SampleType>;
%template(fmax) fmax<SampleType>;
%template(fmin) fmin<SampleType>;
%template(fmod) fmod<SampleType>;
%template(fpclassify) fpclassify<SampleType>;
%template(hypot) hypot<SampleType>;
%template(ilogb) ilogb<SampleType>;
%template(isfinite) isfinite<SampleType>;
%template(isgreater) isgreater<SampleType>;
%template(isgreaterequal) isgreaterequal<SampleType>;
%template(isin) isinf<SampleType>;
%template(isless) isless<SampleType>;
%template(islessequal) islessequal<SampleType>;
%template(isnan) isnan<SampleType>;
%template(isnormal) isnormal<SampleType>;
%template(isunordered) isunordered<SampleType>;
%template(ldexp) ldexp<SampleType>;
%template(lgamma) lgamma<SampleType>;
%template(llrint) llrint<SampleType>;
%template(llround) llround<SampleType>;
%template(log1p) log1p<SampleType>;
%template(lrint) lrint<SampleType>;
%template(lround) lround<SampleType>;
%template(nan) nan<SampleType>;
%template(nanf) nanf<SampleType>;
%template(nanl) nanl<SampleType>;
%template(nearbyint) nearbyint<SampleType>;
%template(nextafter) nextafter<SampleType>;
%template(nexttoward) nexttoward<SampleType>;
%template(remainder) remainder<SampleType>;
%template(rint) rint<SampleType>;
%template(round) round<SampleType>;
%template(scalbln) scalbln<SampleType>;
%template(scalbn) scalbn<SampleType>;
%template(square) square<SampleType>;
%template(tgamma) tgamma<SampleType>;
%template(trunc) trunc<SampleType>;
*/


/*
%template (fastcos)     Ops::fastcos<SampleType>;
%template (fastcosdeg)  Ops::fastcosdeg<f32>;
%template (fastsin)     Ops::fastsin<f32>;
%template (fastsindeg)  Ops::fastsindeg<f32>;
%template (coshsinh)    Ops::coshsinh<f32>;
%template (sinhcosh)    Ops::sinhcosh<f32>;
%template (cossindeg)   Ops::cossindeg<f32>;
%template (sincosdeg)   Ops::fastsindeg<f32>;
*/
/*
// these work but only for integer
%template(sataddi64) Ops::satadd<i64>;
%template(roli32) Ops::rol<i32>;
%template(ro4i32) Ops::ror<i32>;
%template(shli32) Ops::shl<i32>;
%template(shri32) Ops::shr<i32>;
%template(remi32)    Ops::rem<i32>;
%template(sataddi32) Ops::satadd<i32>;
%template(satsubi32) Ops::satsub<i32>;
%template(bitwiseandi32) Ops::bitwiseand<i32>;
%template(bitwiseori32) Ops::bitwiseor<i32>;
%template(bitwisexori32) Ops::bitwisexor<i32>;
%template(bitwiseandnoti32) Ops::bitwiseandnot<i32>;
%template(bitwisenoti32) Ops::bitwisenot<i32>;
*/

%inline %{

    template<typename T>
    struct SampleVector
    {
        std::vector<kfr::univector<T>> samples;
        size_t                         channels;

        SampleVector(size_t channels) {
            samples.resize(channels);
            this->channels = channels;
        }
        
        T& operator()(size_t ch, size_t i) { return samples[ch][i]; }
        
        size_t num_channels() const { return channels; }
        size_t size() const { return samples[0].size(); }
        
        kfr::univector<T> get_channel(size_t channel) { return samples[channel]; }
        void set_channel(size_t channel, kfr::univector<T> & v) { samples[channel] = v; }

        kfr::univector<T> __getitem(size_t i ) { return samples[i]; }
        void __setitem(size_t i, kfr::univector<T> & v) { samples[i] = v; }

    };

    /*
    template<typename T>
    SampleVector<T> deinterleave(size_t channels, kfr::univector<T> & v) {
        SampleVector<T> r(channels);        
        for(size_t i = 0; i < channels; i++) {
            r.samples[i].resize(v.size()/channels)
            for(size_t j = i; j < v.size(); j += channels)
                r[channels][j] = v[j];
        }
        return r;
    }
    template<typename T>
    void interleave(SampleVector<T> & samples, kfr::univector<T> & out) {
        out.resize(samples.channels * samples[0].size());        
        for(size_t i = 0; i < samples.channels; i++)            
            for(size_t j = i; j < samples[i].size(); i+=samples.channels)
                out[j*channels + i] = samples[i][j];
        }
    */
    template<typename T>
    void copy(kfr::univector<T> & dst, std::vector<T> & src) {
        std::copy(src.begin(),src.end(),dst.begin());
    }
    template<typename T>
    void copy(std::vector<T> & dst, kfr::univector<T> & src) {
        std::copy(src.begin(),src.end(),dst.begin());
    }

%}



%inline %{

        template <typename T> T f_note_to_hertz(const T& input) {
            return kfr::note_to_hertz<T>(input);
        }    
        template <typename T> T f_hertz_to_note(const T& input) {
            return kfr::hertz_to_note<T>(input);
        }    
        template <typename T> T f_amp_to_dB(const T& input) {
            return kfr::amp_to_dB<T>(input);
        }    
        template <typename T> T f_dB_to_amp(const T& input) {
            return kfr::dB_to_amp<T>(input);
        }    
        template <typename T> T f_power_to_dB(const T& input) {
            return kfr::power_to_dB<T>(input);
        }    
        template <typename T> T f_dB_to_power(const T& input) {
            return kfr::dB_to_power<T>(input);
        }    
        
        /*
        template<typename T> kfr::complex<T> goertzal(kfr::complex<T> & result, T  omega) {
            kfr::complex<T> r(result);
            kfr::goertzal(r,omega);
            return r;
        }
        */
        

        template <typename T> T waveshaper_hardclip(T & input, double clip_level) 
        {            
            return kfr::waveshaper_hardclip(input,clip_level);
        }
        template <typename T> kfr::univector<T> waveshaper_hardclip(kfr::univector<T> & input, double clip_level) 
        {            
            kfr::univector r(input.size());
            for(size_t i = 0; i < input.size(); i++)
                r[i] = kfr::waveshaper_hardclip(input[i],clip_level);
            return r;
        }

        template <typename T> T waveshaper_tanh(T & input, double sat) 
        {            
            return kfr::waveshaper_tanh(input,sat);
        }
        template <typename T> kfr::univector<T> waveshaper_tanh(kfr::univector<T> & input, double sat) 
        {            
            kfr::univector r(input.size());
            for(size_t i = 0; i < input.size(); i++)
                r[i] = kfr::waveshaper_tanh(input[i],sat);
            return r;
        }

        template <typename T> T waveshaper_saturate_I(T & input, double sat) 
        {            
            return kfr::waveshaper_saturate_I(input,sat);
        }
        template <typename T> kfr::univector<T> waveshaper_saturate_I(kfr::univector<T> & input, double sat) 
        {            
            kfr::univector r(input.size());
            for(size_t i = 0; i < input.size(); i++)
                r[i] = kfr::waveshaper_saturate_I(input[i],sat);
            return r;
        }

        template <typename T> T waveshaper_saturate_II(T & input, double sat) 
        {            
            return kfr::waveshaper_saturate_II(input, sat);
        }
        template <typename T> kfr::univector<T> waveshaper_saturate_II(kfr::univector<T> & input, double sat) 
        {            
            kfr::univector r(input.size());
            for(size_t i = 0; i < input.size(); i++)
                r[i] = kfr::waveshaper_saturate_II(input[i],sat);
            return r;
        }
        /*
        template <typename T> T waveshaper_poly(T & input) 
        {            
            return kfr::waveshaper_poly(input);
        }
        template <typename T> kfr::univector<T> waveshaper_poly(kfr::univector<T> & input) 
        {            
            kfr::univector r(input.size());
            for(size_t i = 0; i < input.size(); i++)
                r[i] = kfr::waveshaper_poly(input[i]);
            return r;
        }
        */

%}

// don't really want any of these from DSP:: but the stuff in Kfr can be awkward to use directly in the script languages
%template(waveshaper_hardclip) waveshaper_hardclip<SampleType>;
%template(waveshaper_tanh) waveshaper_tanh<SampleType>;
%template(waveshaper_saturateI) waveshaper_saturate_I<SampleType>;
%template(waveshaper_saturateII) waveshaper_saturate_II<SampleType>;
//%template(waveshaper_poly) waveshaper_poly<SampleType>;


%template(note_to_hertz) f_note_to_hertz<SampleType>;
%template(hertz_to_note) f_hertz_to_note<SampleType>;
%template(amp_to_dB) f_amp_to_dB<SampleType>;
%template(dB_to_amp) f_dB_to_amp<SampleType>;
%template(power_to_dB) f_power_to_dB<SampleType>;
%template(dB_to_power) f_dB_to_power<SampleType>;

%template(dcremove) DSP::dcremove<SampleType>;

%template(window_hann) DSP::make_window_hann<SampleType>;
%template(window_hamming) DSP::make_window_hamming<SampleType>;
%template(window_blackman) DSP::make_window_blackman<SampleType>;
%template(window_blackman_harris) DSP::make_window_blackman_harris<SampleType>;
%template(window_gaussian) DSP::make_window_gaussian<SampleType>;
%template(window_triangular) DSP::make_window_triangular<SampleType>;
%template(window_bartlett) DSP::make_window_bartlett<SampleType>;
%template(window_cosine) DSP::make_window_cosine<SampleType>;
%template(window_bartlett_hann) DSP::make_window_bartlett_hann<SampleType>;
%template(window_bohman) DSP::make_window_bohman<SampleType>;
%template(window_lanczos) DSP::make_window_lanczos<SampleType>;
%template(window_flattop) DSP::make_window_flattop<SampleType>;
%template(window_kaiser) DSP::make_window_kaiser<SampleType>;

%template(window_hann_ptr) DSP::make_window_hann_ptr<SampleType>;
%template(window_hamming_ptr) DSP::make_window_hamming_ptr<SampleType>;
%template(window_blackman_ptr) DSP::make_window_blackman_ptr<SampleType>;
%template(window_blackman_harris_ptr) DSP::make_window_blackman_harris_ptr<SampleType>;
%template(window_gaussian_ptr) DSP::make_window_gaussian_ptr<SampleType>;
%template(window_triangular_ptr) DSP::make_window_triangular_ptr<SampleType>;
%template(window_bartlett_ptr) DSP::make_window_bartlett_ptr<SampleType>;
%template(window_cosine_ptr) DSP::make_window_cosine_ptr<SampleType>;
%template(window_bartlett_hann_ptr) DSP::make_window_bartlett_hann_ptr<SampleType>;
%template(window_bohman_ptr) DSP::make_window_bohman_ptr<SampleType>;
%template(window_lanczos_ptr) DSP::make_window_lanczos_ptr<SampleType>;
%template(window_flattop_ptr) DSP::make_window_flattop_ptr<SampleType>;
%template(window_kaiser_ptr) DSP::make_window_kaiser_ptr<SampleType>;

%template(energy_to_loudness) DSP::energy_to_loudness<SampleType>;
%template(loudness_to_energy) DSP::loudness_to_energy<SampleType>;


%template (sinewave) DSP::sinewave<SampleType>;
%template (squarewave) DSP::squarewave<SampleType>;
%template (trianglewave) DSP::trianglewave<SampleType>;
%template (sawtoothwave) DSP::sawtoothwave<SampleType>;

%template (generate_sine) DSP::generate_sin<SampleType>;
%template (generate_linear) DSP::generate_linear<SampleType>;
%template (generate_exp) DSP::generate_exp<SampleType>;
%template (generate_exp2) DSP::generate_exp2<SampleType>;
%template (generate_cossin) DSP::generate_cossin<SampleType>;


%template (resample) DSP::resample<SampleType>;
%template (convert_sample) DSP::convert_sample<SampleType>;
//%template (amp_to_dB) DSP::amp_to_dB<SampleType>;

//%template(interleave) DSP::do_interleave<SampleType>;
//%template(deinterleave) DSP::do_deinterleave<SampleType>; 

%inline %{
    // needs to return channels and samplerate
    template<typename T>
    kfr::univector<T> wav_load(const char * filename, kfr::audio_format& fmt) {
        DSP::WavReader<T> r(filename);        
        kfr::univector<T> v(r.size());
        r.read(v);
        fmt = r.format();
        return v;
    }
    template<typename T>
    void wav_write(kfr::univector<T> & v, const char * filename, size_t channels, double sample_rate, bool use_w64=false) {
        DSP::WavWriter<T> w(filename,kfr::audio_format{channels,kfr::audio_sample_type::f32,sample_rate,use_w64});
        w.write(v);        
    }
    template<typename T>
    kfr::univector<T> mp3_load(const char * filename) {
        DSP::MP3Reader<T> r(filename);
        kfr::univector<T> v(r.size());
        r.read(v);
        return v;
    }
    template<typename T>
    kfr::univector<T> flac_load(const char * filename) {
        DSP::FlacReader<T> r(filename);
        kfr::univector<T> v(r.size());
        r.read(v);
        return v;
    }
%}
%template(float_wavreader) DSP::WavReader<SampleType>;
%template(float_wavwriter) DSP::WavWriter<SampleType>;
%template(float_mp3_reader) DSP::MP3Reader<SampleType>;
%template(float_flac_reader) DSP::FlacReader<SampleType>;

%template(wav_load) wav_load<SampleType>;
%template(wav_save) wav_write<SampleType>;
%template(mp3_load) mp3_load<SampleType>;
%template(flac_load) flac_load<SampleType>;

//////////////////////////////////////////////////////////////////////////
// Biquad
//////////////////////////////////////////////////////////////////////////
namespace kfr
{
    template <typename T>
    struct biquad_params
    {        
        constexpr biquad_params(const biquad_params<T>& bq);
        constexpr static bool is_pod;
        
        constexpr biquad_params();
        constexpr biquad_params(T a0, T a1, T a2, T b0, T b1, T b2);

        T a0;
        T a1;
        T a2;
        T b0;
        T b1;
        T b2;
        biquad_params<T> normalized_a0() const;
        biquad_params<T> normalized_b0() const;
        biquad_params<T> normalized_all() const;
    };

    template <typename T>
        struct zpk
        {
            univector<complex<T>> z;
            univector<complex<T>> p;
            T k;
        };

    namespace CMT_ARCH_NAME
    {

        template <typename T, size_t maxfiltercount = 64>
        class biquad_filter
        {
        public:       
            biquad_filter(const biquad_params<T>* bq, size_t count);     
            biquad_filter(const std::vector<biquad_params<T>>& bq);
            void apply(kfr::univector<T> & vector);
            void apply(kfr::univector<T> & dst, const kfr::univector<T> & src);
            void apply(T * dst, const T* src, size_t n);
            void apply(T * buffer, size_t n);
            void reset();                    
        };

        template <typename T = fbase> biquad_params<T> biquad_allpass(T frequency, T Q);
        template <typename T = fbase> biquad_params<T> biquad_lowpass(T frequency, T Q);
        template <typename T = fbase> biquad_params<T> biquad_highpass(T frequency, T Q);
        template <typename T = fbase> biquad_params<T> biquad_bandpass(T frequency, T Q);
        template <typename T = fbase> biquad_params<T> biquad_notch(T frequency, T Q);
        template <typename T = fbase> biquad_params<T> biquad_peak(T frequency, T Q, T gain);
        template <typename T = fbase> biquad_params<T> biquad_lowshelf(T frequency, T gain);
        template <typename T = fbase> biquad_params<T> biquad_highshelf(T frequency, T gain);       
        
        template <typename T> zpk<T> chebyshev1(int N, T rp);
        template <typename T> zpk<T> chebyshev2(int N, T rs);
        template <typename T> zpk<T> butterworth(int N);
        template <typename T> zpk<T> bessel(int N);

        template <typename T> zpk<T> iir_lowpass(const zpk<T>& filter, T frequency, T fs = 2.0);
        template <typename T> zpk<T> iir_highpass(const zpk<T>& filter, T frequency, T fs = 2.0);
        template <typename T> zpk<T> iir_bandpass(const zpk<T>& filter, T lowfreq, T highfreq, T fs = 2.0);
        template <typename T> zpk<T> iir_bandstop(const zpk<T>& filter, T lowfreq, T highfreq, T fs = 2.0);
 
        template <typename T> std::vector<biquad_params<T>> to_sos(const zpk<T>& filter);                
    }
    template <typename T> expression_pointer<T> make_kfilter(int samplerate);
}

%template(Biquad)           DSP::Biquad<SampleType>;
%template(biquadparams)     DSP::BiQuadParams<SampleType>;
%template(biquad)    DSP::biquad<SampleType>;
%template(notch_params)     DSP::notch_params<SampleType>;
%template(lowpass_params)   DSP::lowpass_params<SampleType>;
%template(highpass_params)  DSP::highpass_params<SampleType>;
%template(peak_params)      DSP::peak_params<SampleType>;
%template(lowshelf_params)  DSP::lowshelf_params<SampleType>;
%template(highshelf_params) DSP::highshelf_params<SampleType>;
%template(bandpass_params)  DSP::bandpass_params<SampleType>;
%template(notch_filter)     DSP::NotchFilter<SampleType>;
%template(lowpass_filter)   DSP::LowPassFilter<SampleType>;
%template(highpass_filter)  DSP::HighPassFilter<SampleType>;
%template(bandpass_filter)  DSP::BandPassFilter<SampleType>;
%template(peak_filter) DSP::PeakFilter<SampleType>;
%template(lowshelf_filter)  DSP::LowShelfFilter<SampleType>;
%template(highshelf_filter) DSP::HighShelfFilter<SampleType>;

%template(normalize_frequency) DSP::normalize_frequency<SampleType>;

%template(lowpassfilter)   DSP::lowpassfilter<SampleType>;
%template(highpassfilter)  DSP::highpassfilter<SampleType>;
%template(bandpassfilter)  DSP::bandpassfilter<SampleType>;
%template(peakfilter)      DSP::peakfilter<SampleType>;
%template(lowshelffilter)  DSP::lowshelffilter<SampleType>;
%template(highshelffilter) DSP::highshelffilter<SampleType>;
%template(notchfilter)     DSP::notchfilter<SampleType>;

%template (bessel_filter) DSP::BesselFilter<SampleType>;
%template (bessel_lowpass_filter) DSP::BesselLowPassFilter<SampleType>;
%template (bessel_highpass_filter) DSP::BesselHighPassFilter<SampleType>;
%template (bessel_bandpass_filter) DSP::BesselBandPassFilter<SampleType>;
%template (bessel_bandstop_filter) DSP::BesselBandStopFilter<SampleType>;

%template (butterworth_filter) DSP::ButterworthFilter<SampleType>;
%template (butterworth_lowpass_filter) DSP::ButterworthLowPassFilter<SampleType>;
%template (butterworth_highpass_filter) DSP::ButterworthHighPassFilter<SampleType>;
%template (butterworth_bandpass_filter) DSP::ButterworthBandPassFilter<SampleType>;
%template (butterworth_bandstop_filter) DSP::ButterworthBandStopFilter<SampleType>;

%template (chevyshev1_filter)         DSP::Chebyshev1Filter<SampleType>;
%template (chevyshev1_lowpass_filter)  DSP::Chebyshev1LowPassFilter<SampleType>;
%template (chevyshev1_highpass_filter) DSP::Chebyshev1HighPassFilter<SampleType>;
%template (chevyshev1_bandpass_filter) DSP::Chebyshev1BandPassFilter<SampleType>;
%template (chevyshev1_bandptop_filter) DSP::Chebyshev1BandStopFilter<SampleType>;

%template (chevyshev2_filter) DSP::Chebyshev2Filter<SampleType>;
%template (chevyshev2_lowpass_filter) DSP::Chebyshev2LowPassFilter<SampleType>;
%template (chevyshev2_highpass_filter) DSP::Chebyshev2HighPassFilter<SampleType>;
%template (chevyshev2_bandpass_filter) DSP::Chebyshev2BandPassFilter<SampleType>;
%template (chevyshev2_bandstop_filter) DSP::Chebyshev2BandStopFilter<SampleType>;

%template(bessel_lowpass)  DSP::bessel_lowpass<SampleType>;
%template(bessel_highpass) DSP::bessel_highpass<SampleType>;
%template(bessel_bandpass) DSP::bessel_bandpass<SampleType>;
%template(bessel_bandstop) DSP::bessel_bandstop<SampleType>;

%template(butterworth_lowpass) DSP::butterworth_lowpass<SampleType>;
%template(butterworth_highpass) DSP::butterworth_highpass<SampleType>;
%template(butterworth_bandpass) DSP::butterworth_bandpass<SampleType>;
%template(butterworth_bandstop) DSP::butterworth_bandstop<SampleType>;

%template(chebyshev1_lowpass)  DSP::chebyshev1_lowpass<SampleType>;
%template(chebyshev1_highpass) DSP::chebyshev1_highpass<SampleType>;
%template(chebyshev1_bandpass) DSP::chebyshev1_bandpass<SampleType>;
%template(chebyshev1_bandstop) DSP::chebyshev1_bandstop<SampleType>;

%template(chebyshev2_lowpass)  DSP::chebyshev2_lowpass<SampleType>;
%template(chebyshev2_highpass) DSP::chebyshev2_highpass<SampleType>;
%template(chebyshev2_bandpass) DSP::chebyshev2_bandpass<SampleType>;
%template(chebyshev2_bandstop) DSP::chebyshev2_bandstop<SampleType>;

%template(biquad_params) kfr::biquad_params<SampleType>;
%template(biquad_filter) kfr::CMT_ARCH_NAME::biquad_filter<SampleType>;
%template(biquad_allpass) kfr::CMT_ARCH_NAME::biquad_allpass<SampleType>;
%template(biquad_lowpass) kfr::CMT_ARCH_NAME::biquad_lowpass<SampleType>;
%template(biquad_highpass) kfr::CMT_ARCH_NAME::biquad_highpass<SampleType>;
%template(biquad_bandpass) kfr::CMT_ARCH_NAME::biquad_bandpass<SampleType>;
%template(biquad_notch) kfr::CMT_ARCH_NAME::biquad_notch<SampleType>;
%template(biquad_peak) kfr::CMT_ARCH_NAME::biquad_peak<SampleType>;
%template(biquad_lowshelf) kfr::CMT_ARCH_NAME::biquad_lowshelf<SampleType>;
%template(biquad_highshelf) kfr::CMT_ARCH_NAME::biquad_highshelf<SampleType>;

%template(zpk)              kfr::zpk<SampleType>;
%template(chebyshev1)       kfr::CMT_ARCH_NAME::chebyshev1<SampleType>;
%template(chebyshev2)       kfr::CMT_ARCH_NAME::chebyshev2<SampleType>;
%template(butterworth)      kfr::CMT_ARCH_NAME::butterworth<SampleType>;
%template(bessel)           kfr::CMT_ARCH_NAME::bessel<SampleType>;
%template(iir_lowpass)      kfr::CMT_ARCH_NAME::iir_lowpass<SampleType>;
%template(iir_highpass)      kfr::CMT_ARCH_NAME::iir_highpass<SampleType>;
%template(iir_bandpass)      kfr::CMT_ARCH_NAME::iir_bandpass<SampleType>;
%template(iir_bandstop)      kfr::CMT_ARCH_NAME::iir_bandstop<SampleType>;
%template(to_sos)            kfr::CMT_ARCH_NAME::to_sos<SampleType>;
//%template(make_kfilter)      kfr::make_kfilter<SampleType>;

%template (biquad_params_vector) std::vector<kfr::biquad_params<SampleType>>;


namespace kfr::CMT_ARCH_NAME
{
        
    template <typename T> univector<T> convolve(const univector<T>& src1, const univector<T>& src2);
    template <typename T> univector<T> correlate(const univector<T>& src1, const univector<T>& src2);        
    template <typename T> univector<T> autocorrelate(const univector<T>& src);

}
namespace kfr {
    template <typename T> class convolve_filter
    {
    public:
        
        explicit convolve_filter(size_t size, size_t block_size = 1024);
        explicit convolve_filter(const univector<T>& data, size_t block_size = 1024);
        
        void set_data(const univector<T>& data);
        void reset() final;
        /// Apply filter to multiples of returned block size for optimal processing efficiency.
        size_t input_block_size() const { return block_size; }    

        void apply(univector<T> & in);
        void apply(univector<T> & in, univector<T> & out);
        
        void apply(T* in,size_t);
        void apply(T* dst, T* src,size_t size);

    };

    template <typename T> kfr::filter<T>* make_convolve_filter(const univector<T>& taps, size_t block_size);    
    
}

//////////////////////////////////////////////////////////////////////////
// Convolution
//////////////////////////////////////////////////////////////////////////

%template (convolve_filter) kfr::convolve_filter<SampleType>;
%template (complex_convolve_filter) kfr::convolve_filter<kfr::complex<SampleType>>;
%template (make_convolve_filter) kfr::make_convolve_filter<SampleType>;

%template (conv)  kfr::CMT_ARCH_NAME::convolve<SampleType>;
%template (complex_conv)  kfr::CMT_ARCH_NAME::convolve<kfr::complex<SampleType>>;
%template (acorr) kfr::CMT_ARCH_NAME::autocorrelate<SampleType>;
%template (complex_acorr) kfr::CMT_ARCH_NAME::autocorrelate<kfr::complex<SampleType>>;
%template (xcorr) kfr::CMT_ARCH_NAME::correlate<SampleType>;
%template (complex_xcorr) kfr::CMT_ARCH_NAME::correlate<kfr::complex<SampleType>>;




namespace kfr 
{
    
    enum class dft_type
    {
        both,
        direct,
        inverse
    };

    enum class dft_order
    {
        normal,
        internal, // possibly bit/digit-reversed, implementation-defined, faster to compute
    };

    enum class dft_pack_format
    {
        Perm, // {X[0].r, X[N].r}, ... {X[i].r, X[i].i}, ... {X[N-1].r, X[N-1].i}
        CCs // {X[0].r, 0}, ... {X[i].r, X[i].i}, ... {X[N-1].r, X[N-1].i},  {X[N].r, 0}
    };    

    template <typename T>
    struct dft_plan
    {
        size_t size;
        size_t temp_size;

        explicit dft_plan(size_t size, dft_order order = dft_order::normal);
        ~dft_plan() {}

        void dump() const;
        void execute(complex<T>* out, const complex<T>* in, u8* temp,bool inverse = false) const;        
        
        template <bool inverse> void execute(complex<T>* out, const complex<T>* in, u8* temp, cbool_t<inverse> inv) const;

        
        void execute(univector<complex<T>>& out, const univector<complex<T>>& in,univector<u8>& temp, bool inverse = false) const;                
        //void execute(univector<complex<T>>& out, const univector<complex<T>>& in,univector<u8>& temp, cbool_t<kfr::dft_type::inverse> inv) const;

        
        size_t data_size;
        
    };

    template <typename T>
    struct dft_plan_real
    {
        size_t size;
    
        explicit dft_plan_real(size_t size, dft_pack_format fmt = dft_pack_format::CCs);
                                
        void execute(complex<T>* out, const T* in, u8* temp) const;
        void execute(T* out, const complex<T>* in, u8* temp) const;
        void execute(univector<complex<T>>& out, const univector<T>& in,univector<u8>& temp) const;        
        void execute(univector<T>& out, const univector<complex<T>>& in,univector<u8>& temp) const;
                
    };

    /// @brief DCT type 2 (unscaled)
    template <typename T>
    struct dct_plan
    {
        dct_plan(size_t size);
        void execute(T* out, const T* in, u8* temp, bool inverse = false) const;
        
        void execute(univector<T>& out, const univector<T>& in,univector<u8>& temp, bool inverse = false) const;
    };

    
}

%template(dft_plan) kfr::dft_plan<SampleType>;
%template(dft_plan_real) kfr::dft_plan_real<SampleType>;
%template(dct_plan) kfr::dct_plan<SampleType>;

%template(dft) DSP::run_dft<SampleType>;
%template(dft_real) DSP::run_realdft<SampleType>;
%template(idft) DSP::run_idft<SampleType>;
%template(idft_real) DSP::run_irealdft<SampleType>;

%template (DCTPlan) DSP::DCTPlan<SampleType>;
%template (DFTPlan) DSP::DFTPlan<SampleType>;
%template (DFTRealPlan) DSP::DFTRealPlan<SampleType>;


namespace kfr::CMT_ARCH_NAME
{
    template <typename T, size_t Size>
    using fir_taps = univector<T, Size>;

    template <typename T, typename U = T>
    class fir_filter 
    {
    public:
        fir_filter(const univector<T>& taps);
        void set_taps(const univector<T>& taps);
        void reset() final;            
        void apply(kfr::univector<T> & vector);
        void apply(const kfr::univector<T> & src, kfr::univector<T> & dst);           
    };

    template<typename T> void fir_lowpass(univector<T>& taps, T cutoff,
                                const expression_pointer<T>& window, bool normalize = true);
    template<typename T> void fir_highpass(univector<T>& taps, T cutoff, const expression_pointer<T>& window, bool normalize = true);
    template<typename T> void fir_bandpass(univector<T>& taps, T frequency1, T frequency2,
                                    const expression_pointer<T>& window, bool normalize = true);
    template<typename T> void fir_bandstop(univector<T>& taps, T frequency1, T frequency2,
                                    const expression_pointer<T>& window, bool normalize = true);
}

namespace kfr {
    template <typename U, typename T> fir_filter<U>* make_fir_filter(const univector<T>& taps);
}

%template (fir_filter) DSP::FIRFilter<SampleType>;
%template (fir_bandpass_filter) DSP::FIRBandpassFilter<SampleType>;
%template (fir_lowpass_filter) DSP::FIRLowpassFilter<SampleType>;
%template (fir_highpass_filter) DSP::FIRHighpassFilter<SampleType>;
%template (fir_bandstop_filter) DSP::FIRBandstopFilter<SampleType>;
%template (fir_lowpass)  DSP::fir_lowpass<SampleType>;
%template (fir_highpass) DSP::fir_highpass<SampleType>;
%template (fir_bandpass) DSP::fir_bandpass<SampleType>;
%template (fir_bandstop) DSP::fir_bandstop<SampleType>;

namespace kfr::CMT_ARCH_NAME 
{
    using resample_quality = sample_rate_conversion_quality;

    template <typename T>
    struct samplerate_converter
    {    
        using itype = i64;
        using ftype = subtype<T>;
        static size_t filter_order(sample_rate_conversion_quality quality);
        static ftype sidelobe_attenuation(sample_rate_conversion_quality quality);
        static ftype transition_width(sample_rate_conversion_quality quality);
        static ftype window_param(sample_rate_conversion_quality quality);        
        samplerate_converter(sample_rate_conversion_quality quality, itype interpolation_factor, itype decimation_factor, ftype scale = ftype(1), ftype cutoff = 0.5f);
        itype input_position_to_intermediate(itype in_pos) const;
        itype output_position_to_intermediate(itype out_pos) const;
        itype input_position_to_output(itype in_pos) const;
        itype output_position_to_input(itype out_pos) const;
        itype output_size_for_input(itype input_size) const;
        itype input_size_for_output(itype output_size) const;
        size_t skip(size_t output_size, univector<T> input);        
        size_t process(kfr::univector<T> output, univector<T> input);
        double get_fractional_delay() const;
        size_t get_delay() const;

        ftype kaiser_beta;
        itype depth;
        itype taps;
        size_t order;
        itype interpolation_factor;
        itype decimation_factor;
        //univector<T> filter;
        //univector<T> delay;
        itype input_position;
        itype output_position;
    };
    // Deprecated in 0.9.2        

    template<typename T>
    samplerate_converter<T> resampler(sample_rate_conversion_quality quality,
                                                size_t interpolation_factor, size_t decimation_factor,
                                                T scale, T cutoff);
}


%template(samplerate_converter) kfr::CMT_ARCH_NAME::samplerate_converter<SampleType>;
%template(resampler) kfr::CMT_ARCH_NAME::resampler<SampleType>;


////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////
namespace kfr {
    template <typename T>
    struct audio_reader : public abstract_reader<T>
    {
        /// @brief Reads interleaved audio
        using abstract_reader<T>::read;

        univector2d<T> read_channels();
        univector2d<T> read_channels(size_t size);
        /// @brief Returns audio format description
        virtual const audio_format_and_length& format() const = 0;
    };

    template <typename T>
    struct audio_writer : public abstract_writer<T>
    {
        /// @brief Writes interleaved audio
        using abstract_writer<T>::write;

        template <univector_tag Tag1, univector_tag Tag2>
        size_t write_channels(const univector2d<T, Tag1, Tag2>& data);
        /// @brief Returns audio format description
        virtual const audio_format_and_length& format() const = 0;

        /// @brief Finishes writing and closes underlying writer
        virtual void close() = 0;
    };
    
    template <typename T>
    struct audio_writer_wav : audio_writer<T>
    {
        /// @brief Constructs WAV writer using target writer and format
        audio_writer_wav(std::shared_ptr<abstract_writer<>>&& writer, const audio_format& fmt);
        ~audio_writer_wav() override;
        using audio_writer<T>::write;

        /// @brief Write data to underlying binary writer
        /// data is PCM samples in interleaved format
        /// size is the number of samples (PCM frames * channels)
        size_t write(const T* data, size_t size) override;
        void close() override;
        const audio_format_and_length& format() const override;
        imax tell() const override;
        bool seek(imax, seek_origin) override;
    };

    /// @brief WAV format reader
    template <typename T>
    struct audio_reader_wav : audio_reader<T>
    {
        using audio_reader<T>::read;

        /// @brief Constructs WAV reader
        audio_reader_wav(std::shared_ptr<abstract_reader<>>&& reader);
        ~audio_reader_wav() override;
        const audio_format_and_length& format() const override;
        size_t read(T* data, size_t size) override;
        imax tell() const override;
        bool seek(imax offset, seek_origin origin) override;
    };

    /// @brief FLAC format reader
    template <typename T>
    struct audio_reader_flac : audio_reader<T>
    {
        /// @brief Constructs FLAC reader
        audio_reader_flac(std::shared_ptr<abstract_reader<>>&& reader);
        ~audio_reader_flac() override;
        const audio_format_and_length& format() const override;
        size_t read(T* data, size_t size) override;
        imax tell() const override;
        bool seek(imax offset, seek_origin origin) override;
    };

    template <typename T>
    struct audio_reader_mp3 : audio_reader<T>
    {
        /// @brief Constructs MP3 reader
        audio_reader_mp3(std::shared_ptr<abstract_reader<>>&& reader);
        ~audio_reader_mp3() override;

        drmp3_config config{ 0, 0 };

        /// @brief Returns audio format description
        const audio_format_and_length& format() const override;
        size_t read(T* data, size_t size) override;
        imax tell() const override;
        bool seek(imax offset, seek_origin origin) override;
    };
} // namespace kfr



%include "KfrDSP/KfrAudio.hpp"
%include "KfrDSP/KfrCombFilter.hpp"
%include "KfrDSP/KfrFIRCombFilter.hpp"
%include "KfrDSP/KfrIIRCombFilter.hpp"
%include "KfrDSP/KfrMultiTapCombFilter.hpp"
%include "KfrDSP/KfrMultiTapFirCombFilter.hpp"
%include "KfrDSP/KfrMultitapIIRCombFilter.hpp"

%include "KfrDSP/KfrDelayLine.hpp"
%include "KfrDSP/KfrMultitapDelayLine.hpp"
%include "KfrDSP/KfrNoise.hpp"
%include "KfrDSP/KfrUtils.hpp"

%include "KfrDSP/KfrBiquads.hpp"
%include "KfrDSP/KfrRBJ.hpp"
%include "KfrDSP/KfrZolzer.hpp"

%include "KfrDSP/KfrIIRBesselFilter.hpp"
%include "KfrDSP/KfrIIRButterworthFilter.hpp"
%include "KfrDSP/KfrIIRChebyshev1Filter.hpp"
%include "KfrDSP/KfrIIRChebyshev2Filter.hpp"
%include "KfrDSP/KfrIIR.hpp"

%include "KfrDSP/KfrFIRFilter.hpp"
%include "KfrDSP/KfrFIRBandpassFilter.hpp"
%include "KfrDSP/KfrFIRBandstopFilter.hpp"
%include "KfrDSP/KfrFIRLowpassFilter.hpp"
%include "KfrDSP/KfrFIRHighpassFilter.hpp"
%include "KfrDSP/KfrFIR.hpp"
%include "IIRFilters.hpp"

%include "KfrDSP/KfrWindows.hpp"
%include "KfrDSP/KfrConvolution.hpp"
%include "KfrDSP/KfrDFT.hpp"
%include "KfrDSP/KfrFunctions.hpp"
%include "KfrDSP/KfrResample.hpp"

%template(SampleVector) kfr::univector<SampleType>;
//%template(SampleMatrix)  kfr::univector2d<SampleType>;
%template(ComplexVector) kfr::univector<kfr::complex<SampleType>>;
//%template(ComplexMatrix) kfr::univector2d<kfr::complex<SampleType>>;

%template(ConvolutionFilter) KfrDSP1::ConvolutionFilter<SampleType>;
%template(StereoConvolutionFilter) KfrDSP1::StereoConvolutionFilter<SampleType>;

%template(CombFilter) KfrDSP1::CombFilter<SampleType>;
%template(FirCombFilter) KfrDSP1::FIRCombFilter<SampleType>;
%template(IirCombFilter) KfrDSP1::IIRCombFilter<SampleType>;

%template(MultitapCombFilter)    KfrDSP1::MultiTapCombFilter<SampleType>;
%template(MultitapFirCombFilter) KfrDSP1::MultiTapFIRCombFilter<SampleType>;
%template(MultitapIirCombFilter) KfrDSP1::MultiTapIIRCombFilter<SampleType>;

%template(DelayLine) KfrDSP1::DelayLine<SampleType>;
%template(MultitapDelayLine) KfrDSP1::MultiTapDelayLine<SampleType>;

%template(FirFilter) KfrDSP1::FIRFilter<SampleType>;
%template(FirLowpassFilter) KfrDSP1::FIRLowpassFilter<SampleType>;
%template(FirHighpassFilter) KfrDSP1::FIRHighpassFilter<SampleType>;
%template(FirBandpassFilter) KfrDSP1::FIRBandpassFilter<SampleType>;
%template(FirBandstopFilter) KfrDSP1::FIRBandstopFilter<SampleType>;


%inline %{
float sampleRate=44100.0;
float invSampleRate=1.0/sampleRate;
Std::RandomMersenne noise;
%}