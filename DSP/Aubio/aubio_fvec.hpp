#pragma once

#define SQR(x) ((x)*(x))

/** Window types */
typedef enum
{
  win_ones,
  win_rectangle,
  win_hamming,
  win_hanning,
  win_hanningz,
  win_blackman,
  win_blackman_harris,
  win_gaussian,
  win_welch,
  win_parzen,
  win_default = win_hanningz,
} window_type;


template<typename T>
struct FVec : public std::vector<T>
{
    using base = std::vector<T>;
    using base::operator [];        
    using base::size;
    using base::resize;
    using base::max_size;
    using base::capacity;
    using base::empty;
    using base::reserve;
    using base::shrink_to_fit;
    using base::at;
    using base::front;
    using base::back;
    using base::data;
    using base::assign;
    using base::push_back;
    using base::pop_back;
    using base::insert;
    using base::erase;
    using base::swap;
    using base::clear;
    using base::emplace;
    using base::emplace_back;

    FVec(size_t i) : std::vector<T>(i) {}

    void set_sample(size_t pos, T v) {
        (*this)[i] = v;
    }
    T&   get_sample(size_t pos) {
        return (*this)[i];
    }

    void fill(T v) {
        #pragma omp simd
        for(size_t i = 0; i < size(); i++) (*this)[i] = v;
    }
    void zero() {
        fill((T)0);
    }
    void ones() {
        fill((T)1);
    }

    void weight(const FVec &weight) {
        size_t length = std::min(size(), weight.size());
    #if defined(HAVE_INTEL_IPP)
        ippsMul(data(), weight.data(), data(), (int)size();
    #elif defined(HAVE_ACCELERATE)
        vDSP_vmul( data(), 1, weight.data(), 1, data(), 1, size() );
    #else
        size_t j;
        #pragma omp simd
        for (j = 0; j < length; j++) {
            (*this)[j] *= weight[j];
        }
    #endif 
    }
    void print() {
        std::cout << "VECTOR[" << size() << "]";        
        for(size_t i = 0; i < size(); i++)
            std::cout << (*this)[i] << ",";
        std::cout << std::endl;
    }

    FVec<T>& operator +=(const FVec<T> & b) {
        #pragma omp simd
        for(size_t i = 0; i < a.size(); i++) (*this) += b[i];
        return *this;
    }   
    
    FVec<T>& operator -=(const FVec<T> & b) {
        #pragma omp simd
        for(size_t i = 0; i < a.size(); i++) (*this) -= b[i];
        return *this;
    }   
    FVec<T>& operator *=(const FVec<T> & b) {
        #pragma omp simd
        for(size_t i = 0; i < a.size(); i++) (*this) *= b[i];
        return *this;
    }   
    FVec<T>& operator /=(const FVec<T> & b) {
        #pragma omp simd
        for(size_t i = 0; i < a.size(); i++) (*this) /= b[i];
        return *this;
    }   

    
    FVec<T>& operator += (const T& b) {
        #pragma omp simd
        for(size_t i = 0; i < a.size(); i++) (*this) += b;
        return *this;
    }   
    FVec<T>& operator -= (const T& b) {
        #pragma omp simd
        for(size_t i = 0; i < a.size(); i++) (*this) -= b;
        return *this;
    }   
    FVec<T>& operator *= (const T& b) {
        #pragma omp simd
        for(size_t i = 0; i < a.size(); i++) (*this) *= b;
        return *this;
    }   
    FVec<T>& operator /= (const T& b) {
        #pragma omp simd
        for(size_t i = 0; i < a.size(); i++) (*this) /= b;
        return *this;
    }   

};

template<typename T>
std::ostream& operator << (std::ostream & o, const FVec<T> & m )
{
    for(size_t i = 0; i < m.size(); i++)
    {        
        o << m(i) << ",";        
    }
    o << std::endl;
    return o;
}

template<typename T>
bool operator == (const FVec<T> & a, const FVec<T> & b) {
    return std::equal(a.begin(),a.end(),b.end());
}

template<class T>
FVec<T> operator +(const FVec<T> & a, const FVec<T> & b) {
    FVec<T> r(a.size());
    #pragma omp simd
    for(size_t i = 0; i < a.size(); i++) r[i] = a[i] + b[i];
    return r;
}   
template<class T>
FVec<T> operator -(const FVec<T> & a, const FVec<T> & b) {
    FVec<T> r(a.size());
    #pragma omp simd
    for(size_t i = 0; i < a.size(); i++) r[i] = a[i] - b[i];
    return r;
}   
template<class T>
FVec<T> operator *(const FVec<T> & a, const FVec<T> & b) {
    FVec<T> r(a.size());
    #pragma omp simd
    for(size_t i = 0; i < a.size(); i++) r[i] = a[i] * b[i];
    return r;
}   
template<class T>
FVec<T> operator /(const FVec<T> & a, const FVec<T> & b) {
    FVec<T> r(a.size());
    #pragma omp simd
    for(size_t i = 0; i < a.size(); i++) r[i] = a[i] / b[i];
    return r;
}   
template<class T>
FVec<T> operator %(const FVec<T> & a, const FVec<T> & b) {
    FVec<T> r(a.size());
    #pragma omp simd
    for(size_t i = 0; i < a.size(); i++) r[i] = std::fmod(a[i],b[i]);
    return r;
}   
template<class T>
FVec<T> operator +(const FVec<T> & a, const T& b) {
    FVec<T> r(a.size());
    #pragma omp simd
    for(size_t i = 0; i < a.size(); i++) r[i] = a[i] + b;
    return r;
}   
template<class T>
FVec<T> operator -(const FVec<T> & a, const T& b) {
    FVec<T> r(a.size());
    #pragma omp simd
    for(size_t i = 0; i < a.size(); i++) r[i] = a[i] - b;
    return r;
}   
template<class T>
FVec<T> operator *(const FVec<T> & a, const T& b) {
    FVec<T> r(a.size());
    #pragma omp simd
    for(size_t i = 0; i < a.size(); i++) r[i] = a[i] * b;
    return r;
}   
template<class T>
FVec<T> operator / (const FVec<T> & a, const T& b) {
    FVec<T> r(a.size());
    #pragma omp simd
    for(size_t i = 0; i < a.size(); i++) r[i] = a[i] / b;
    return r;
}   
template<class T>
FVec<T> operator %(const FVec<T> & a, const T& b) {
    FVec<T> r(a.size());
    #pragma omp simd
    for(size_t i = 0; i < a.size(); i++) r[i] = std::fmod(a[i],b);
    return r;
}   
template<class T>
FVec<T> operator +(const T & a, const FVec<T>& b) {
    FVec<T> r(a.size());
    #pragma omp simd
    for(size_t i = 0; i < a.size(); i++) r[i] = a + b[i];
    return r;
}   
template<class T>
FVec<T> operator -(const T & a, const FVec<T>& b) {
    FVec<T> r(a.size());
    #pragma omp simd
    for(size_t i = 0; i < a.size(); i++) r[i] = a - b[i];
    return r;
}   
template<class T>
FVec<T> operator *(const T & a, const FVec<T>& b) {
    FVec<T> r(a.size());
    #pragma omp simd
    for(size_t i = 0; i < a.size(); i++) r[i] = a * b[i];
    return r;
}   
template<class T>
FVec<T> operator /(const T & a, const FVec<T>& b) {
    FVec<T> r(a.size());
    #pragma omp simd
    for(size_t i = 0; i < a.size(); i++) r[i] = a / b[i];
    return r;
}   
template<typename T>
void copy_fvec(FVec<T> & dst, const FVec<T> & src) {
    std::copy(src.begin(),src.end(),dst.begin());
}

template<typename T>
void copy_vector(FVec<T> & dst, size_t n, const T * src) {
    std::copy(&src[0],&src[n-1],dst.begin());
}
template<typename T>
FVec<T> slice_fvec(size_t start, size_t end, const FVec<T> & src) {
    FVec<T> r(end-start);
    std::copy(src.begin()+start,src.begin()+end,r.begin());
    return r;
}

template<typename T>
FVec<T> normalize(const FVec<T> & a) {
    FVec<T> r(a);        
    auto max = std::max_element(r.begin(),r.end());
    if(*max > 0) for(size_t i = 0; i < r.size(); i++) r[i] /= *max;
    return r;
}
template<class T>
FVec<T> cos(const FVec<T> & v) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::cos(v[i]);
    return r;
}
template<class T>
FVec<T> sin(const FVec<T> & v) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::sin(v[i]);
    return r;
}    
template<class T>
FVec<T> tan(const FVec<T> & v) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::tan(v[i]);
    return r;
}

template<class T>
FVec<T> acos(const FVec<T> & v) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::acos(v[i]);
    return r;
}
template<class T>
FVec<T> asin(const FVec<T> & v) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::asin(v[i]);
    return r;
}    
template<class T>
FVec<T> atan(const FVec<T> & v) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::atan(v[i]);
    return r;
}    
template<class T>
FVec<T> atan2(const FVec<T> & v, const T value) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::atan2(v[i], value);
    return r;
}    
template<class T>
FVec<T> atan2(const FVec<T> & v, const FVec<T> value) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::atan2(v[i], value[i]);
    return r;
}    
template<class T>
FVec<T> cosh(const FVec<T> & v) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::cosh(v[i]);
    return r;
}
template<class T>
FVec<T> sinh(const FVec<T> & v) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::sinh(v[i]);
    return r;
}    
template<class T>
FVec<T> tanh(const FVec<T> & v) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::tanh(v[i]);
    return r;
}

template<class T>
FVec<T> acosh(const FVec<T> & v) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::acosh(v[i]);
    return r;
}
template<class T>
FVec<T> asinh(const FVec<T> & v) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::asinh(v[i]);
    return r;
}    
template<class T>
FVec<T> atanh(const FVec<T> & v) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::atanh(v[i]);
    return r;
}    

template<class T>
FVec<T> exp(const FVec<T> & v) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::exp(v[i]);
    return r;
}    
template<class T>
FVec<T> log(const FVec<T> & v) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::log(v[i]);
    return r;
}    
template<class T>
FVec<T> log10(const FVec<T> & v) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::log10(v[i]);
    return r;
}    
template<class T>
FVec<T> exp2(const FVec<T> & v) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::exp2(v[i]);
    return r;
}    
template<class T>
FVec<T> expm1(const FVec<T> & v) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::expm1(v[i]);
    return r;
}    
template<class T>
FVec<T> ilogb(const FVec<T> & v) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::ilogb(v[i]);
    return r;
}    
template<class T>
FVec<T> log2(const FVec<T> & v) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::log2(v[i]);
    return r;
}    
template<class T>
FVec<T> log1p(const FVec<T> & v) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::log1p(v[i]);
    return r;
}    
template<class T>
FVec<T> logb(const FVec<T> & v) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::logb(v[i]);
    return r;
}    
template<class T>
FVec<T> scalbn(const FVec<T> & v, const FVec<int> & x) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::scalbn(v[i],x[i]);
    return r;
}    
template<class T>
FVec<T> scalbn(const FVec<T> & v, const int x) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::scalbn(v[i],x);
    return r;
}    
template<class T>
FVec<T> scalbln(const FVec<T> & v, const FVec<long int> & x) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::scalbln(v[i],x[i]);
    return r;
}    
template<class T>
FVec<T> scalbln(const FVec<T> & v, const long int x) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::scalbln(v[i],x);
    return r;
}    
template<class T>
FVec<T> pow(const FVec<T> & v, const FVec<T> & x) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::pow(v[i],x[i]);
    return r;
}    
template<class T>
FVec<T> pow(const FVec<T> & v, const T x) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::pow(v[i],x);
    return r;
}    
template<class T>
FVec<T> pow(const T x, const FVec<T> & v) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::pow(x,v[i]);
    return r;
}    
template<class T>
FVec<T> sqrt(const FVec<T> & v) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::sqrt(v[i]);
    return r;
}    
template<class T>
FVec<T> cbrt(const FVec<T> & v) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::cbrt(v[i]);
    return r;
}    
template<class T>
FVec<T> hypot(const FVec<T> & v, const FVec<T> & x) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::hypot(v[i],x[i]);
    return r;
}    
template<class T>
FVec<T> hypot(const FVec<T> & v, const T x) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::hypot(v[i],x);
    return r;
}    
template<class T>
FVec<T> hypot(const T x, const FVec<T> & v) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::hypot(x,v[i]);
    return r;
}    
template<class T>
FVec<T> erf(const FVec<T> & v) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::erf(v[i]);
    return r;
}    
template<class T>
FVec<T> erfc(const FVec<T> & v) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::erfc(v[i]);
    return r;
}    
template<class T>
FVec<T> tgamma(const FVec<T> & v) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::tgamma(v[i]);
    return r;
}    
template<class T>
FVec<T> lgamma(const FVec<T> & v) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::lgamma(v[i]);
    return r;
}    
template<class T>
FVec<T> ceil(const FVec<T> & v) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::ceil(v[i]);
    return r;
}    
template<class T>
FVec<T> floor(const FVec<T> & v) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::floor(v[i]);
    return r;
}    
template<class T>
FVec<T> fmod(const FVec<T> & v, const FVec<T> & x) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::fmod(v[i],x[i]);
    return r;
}    
template<class T>
FVec<T> fmod(const FVec<T> & v, const T x) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::fmod(v[i],x);
    return r;
}    
template<class T>
FVec<T> fmod(const T x, const FVec<T> & v) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::fmod(x,v[i]);
    return r;
}    
template<class T>
FVec<T> trunc(const FVec<T> & v) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::trunc(v[i]);
    return r;
}    
template<class T>
FVec<T> round(const FVec<T> & v) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::round(v[i]);
    return r;
}    
template<class T>
FVec<long int> lround(const FVec<T> & v) {
    FVec<long int> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::lround(v[i]);
    return r;
}    
template<class T>
FVec<long long int> llround(const FVec<T> & v) {
    FVec<long long int> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::llround(v[i]);
    return r;
}    
template<class T>
FVec<T> nearbyint(const FVec<T> & v) {
    FVec<long long int> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::nearbyint(v[i]);
    return r;
}    
template<class T>
FVec<T> remainder(const FVec<T> & a, const FVec<T> & b) {
    FVec<long long int> r(a.size());
    #pragma omp simd
    for(size_t i = 0; i < a.size(); i++) r[i] = std::remainder(a[i],b[i]);
    return r;
}    
template<class T>
FVec<T> copysign(const FVec<T> & a, const FVec<T> & b) {
    FVec<long long int> r(a.size());
    #pragma omp simd
    for(size_t i = 0; i < a.size(); i++) r[i] = std::copysign(a[i],b[i]);
    return r;
}    
template<class T>
FVec<T> fdim(const FVec<T> & a, const FVec<T> & b) {
    FVec<long long int> r(a.size());
    #pragma omp simd
    for(size_t i = 0; i < a.size(); i++) r[i] = std::fdim(a[i],b[i]);
    return r;
}    
#undef fmax
template<class T>
FVec<T> fmax(const FVec<T> & a, const FVec<T> & b) {
    FVec<long long int> r(a.size());
    #pragma omp simd
    for(size_t i = 0; i < a.size(); i++) r[i] = std::fmax(a[i],b[i]);
    return r;
}    
template<class T>
FVec<T> fmin(const FVec<T> & a, const FVec<T> & b) {
    FVec<long long int> r(a.size());
    #pragma omp simd
    for(size_t i = 0; i < a.size(); i++) r[i] = std::fmin(a[i],b[i]);
    return r;
}    
template<class T>
FVec<T> fma(const FVec<T> & a, const FVec<T> & b, const FVec<T> & c) {
    FVec<long long int> r(a.size());
    #pragma omp simd
    for(size_t i = 0; i < a.size(); i++) r[i] = std::fma(a[i],b[i],c[i]);
    return r;
}    
template<class T>
FVec<T> fabs(const FVec<T> & a) {
    FVec<T> r(a.size());
    #pragma omp simd
    for(size_t i = 0; i < a.size(); i++) r[i] = std::fabs(a[i]);
    return r;
}    

template<class T>
FVec<T> cos(const FVec<T> & v) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::cos(v[i]);
    return r;
}
template<class T>
FVec<T> sin(const FVec<T> & v) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::sin(v[i]);
    return r;
}    
template<class T>
FVec<T> tan(const FVec<T> & v) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::tan(v[i]);
    return r;
}

template<class T>
FVec<T> acos(const FVec<T> & v) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::acos(v[i]);
    return r;
}
template<class T>
FVec<T> asin(const FVec<T> & v) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::asin(v[i]);
    return r;
}    
template<class T>
FVec<T> atan(const FVec<T> & v) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::atan(v[i]);
    return r;
}    
template<class T>
FVec<T> atan2(const FVec<T> & v, const T value) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::atan2(v[i], value);
    return r;
}    
template<class T>
FVec<T> atan2(const FVec<T> & v, const FVec<T> value) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::atan2(v[i], value[i]);
    return r;
}    
template<class T>
FVec<T> cosh(const FVec<T> & v) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::cosh(v[i]);
    return r;
}
template<class T>
FVec<T> sinh(const FVec<T> & v) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::sinh(v[i]);
    return r;
}    
template<class T>
FVec<T> tanh(const FVec<T> & v) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::tanh(v[i]);
    return r;
}

template<class T>
FVec<T> acosh(const FVec<T> & v) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::acosh(v[i]);
    return r;
}
template<class T>
FVec<T> asinh(const FVec<T> & v) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::asinh(v[i]);
    return r;
}    
template<class T>
FVec<T> atanh(const FVec<T> & v) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::atanh(v[i]);
    return r;
}    

template<class T>
FVec<T> exp(const FVec<T> & v) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::exp(v[i]);
    return r;
}    
template<class T>
FVec<T> log(const FVec<T> & v) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::log(v[i]);
    return r;
}    
template<class T>
FVec<T> log10(const FVec<T> & v) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::log10(v[i]);
    return r;
}    
template<class T>
FVec<T> exp2(const FVec<T> & v) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::exp2(v[i]);
    return r;
}    
template<class T>
FVec<T> expm1(const FVec<T> & v) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::expm1(v[i]);
    return r;
}    
template<class T>
FVec<T> ilogb(const FVec<T> & v) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::ilogb(v[i]);
    return r;
}    
template<class T>
FVec<T> log2(const FVec<T> & v) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::log2(v[i]);
    return r;
}    
template<class T>
FVec<T> log1p(const FVec<T> & v) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::log1p(v[i]);
    return r;
}    
template<class T>
FVec<T> logb(const FVec<T> & v) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::logb(v[i]);
    return r;
}    
template<class T>
FVec<T> scalbn(const FVec<T> & v, const FVec<int> & x) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::scalbn(v[i],x[i]);
    return r;
}    
template<class T>
FVec<T> scalbn(const FVec<T> & v, const int x) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::scalbn(v[i],x);
    return r;
}    
template<class T>
FVec<T> scalbln(const FVec<T> & v, const FVec<long int> & x) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::scalbln(v[i],x[i]);
    return r;
}    
template<class T>
FVec<T> scalbln(const FVec<T> & v, const long int x) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::scalbln(v[i],x);
    return r;
}    
template<class T>
FVec<T> pow(const FVec<T> & v, const FVec<T> & x) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::pow(v[i],x[i]);
    return r;
}    
template<class T>
FVec<T> pow(const FVec<T> & v, const T x) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::pow(v[i],x);
    return r;
}    
template<class T>
FVec<T> pow(const T x, const FVec<T> & v) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::pow(x,v[i]);
    return r;
}    
template<class T>
FVec<T> sqrt(const FVec<T> & v) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::sqrt(v[i]);
    return r;
}    
template<class T>
FVec<T> cbrt(const FVec<T> & v) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::cbrt(v[i]);
    return r;
}    
template<class T>
FVec<T> hypot(const FVec<T> & v, const FVec<T> & x) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::hypot(v[i],x[i]);
    return r;
}    
template<class T>
FVec<T> hypot(const FVec<T> & v, const T x) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::hypot(v[i],x);
    return r;
}    
template<class T>
FVec<T> hypot(const T x, const FVec<T> & v) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::hypot(x,v[i]);
    return r;
}    
template<class T>
FVec<T> erf(const FVec<T> & v) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::erf(v[i]);
    return r;
}    
template<class T>
FVec<T> erfc(const FVec<T> & v) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::erfc(v[i]);
    return r;
}    
template<class T>
FVec<T> tgamma(const FVec<T> & v) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::tgamma(v[i]);
    return r;
}    
template<class T>
FVec<T> lgamma(const FVec<T> & v) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::lgamma(v[i]);
    return r;
}    
template<class T>
FVec<T> ceil(const FVec<T> & v) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::ceil(v[i]);
    return r;
}    
template<class T>
FVec<T> floor(const FVec<T> & v) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::floor(v[i]);
    return r;
}    
template<class T>
FVec<T> fmod(const FVec<T> & v, const FVec<T> & x) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::fmod(v[i],x[i]);
    return r;
}    
template<class T>
FVec<T> fmod(const FVec<T> & v, const T x) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::fmod(v[i],x);
    return r;
}    
template<class T>
FVec<T> fmod(const T x, const FVec<T> & v) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::fmod(x,v[i]);
    return r;
}    
template<class T>
FVec<T> trunc(const FVec<T> & v) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::trunc(v[i]);
    return r;
}    
template<class T>
FVec<T> round(const FVec<T> & v) {
    FVec<T> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::round(v[i]);
    return r;
}    
template<class T>
FVec<long int> lround(const FVec<T> & v) {
    FVec<long int> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::lround(v[i]);
    return r;
}    
template<class T>
FVec<long long int> llround(const FVec<T> & v) {
    FVec<long long int> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::llround(v[i]);
    return r;
}    
template<class T>
FVec<T> nearbyint(const FVec<T> & v) {
    FVec<long long int> r(v.size());
    #pragma omp simd
    for(size_t i = 0; i < v.size(); i++) r[i] = std::nearbyint(v[i]);
    return r;
}    
template<class T>
FVec<T> remainder(const FVec<T> & a, const FVec<T> & b) {
    FVec<long long int> r(a.size());
    #pragma omp simd
    for(size_t i = 0; i < a.size(); i++) r[i] = std::remainder(a[i],b[i]);
    return r;
}    
template<class T>
FVec<T> copysign(const FVec<T> & a, const FVec<T> & b) {
    FVec<long long int> r(a.size());
    #pragma omp simd
    for(size_t i = 0; i < a.size(); i++) r[i] = std::copysign(a[i],b[i]);
    return r;
}    
template<class T>
FVec<T> fdim(const FVec<T> & a, const FVec<T> & b) {
    FVec<long long int> r(a.size());
    #pragma omp simd
    for(size_t i = 0; i < a.size(); i++) r[i] = std::fdim(a[i],b[i]);
    return r;
}    
#undef fmax
template<class T>
FVec<T> fmax(const FVec<T> & a, const FVec<T> & b) {
    FVec<long long int> r(a.size());
    #pragma omp simd
    for(size_t i = 0; i < a.size(); i++) r[i] = std::fmax(a[i],b[i]);
    return r;
}    
template<class T>
FVec<T> fmin(const FVec<T> & a, const FVec<T> & b) {
    FVec<long long int> r(a.size());
    #pragma omp simd
    for(size_t i = 0; i < a.size(); i++) r[i] = std::fmin(a[i],b[i]);
    return r;
}    
template<class T>
FVec<T> fma(const FVec<T> & a, const FVec<T> & b, const FVec<T> & c) {
    FVec<long long int> r(a.size());
    #pragma omp simd
    for(size_t i = 0; i < a.size(); i++) r[i] = std::fma(a[i],b[i],c[i]);
    return r;
}    
template<class T>
FVec<T> fabs(const FVec<T> & a) {
    FVec<T> r(a.size());
    #pragma omp simd
    for(size_t i = 0; i < a.size(); i++) r[i] = std::fabs(a[i]);
    return r;
}    

template<typename T>
T unwrap2pi (T phase)
{
    /* mod(phase+pi,-2pi)+pi */
    return phase + TWO_PI * (1. + std::floor(-(phase + PI) / TWO_PI));
}


template<typename T>
T mean (const FVec<T> & s)
{
    T_t tmp = 0.0;
    #if defined(HAVE_INTEL_IPP)
        ippsMean(s.data(), (int)s.size(), &tmp);
        return tmp;
    #elif defined(HAVE_ACCELERATE)
        vDSP_meanv(s.data(), 1, &tmp, s.size());
        return tmp;
    #else
        size_t j;
        #pragma omp simd
        for (j = 0; j < s.size(); j++) {
            tmp += s[j];
        }
        return tmp / (T)(s.size());
    #endif
}

template<typename T>
T sum (const FVec<T>& s)
{
    T_t tmp = 0.0;
    #if defined(HAVE_INTEL_IPP)
        ippsSum(s.data(), (int)s.size(), &tmp);
    #elif defined(HAVE_ACCELERATE)
        vDSP_sve(s.data(), 1, &tmp, s.size());
    #else
        fsize_t j;
        #pragma omp simd
        for (j = 0; j < s.size(); j++) {
            tmp += s[j];
        }
    #endif
    return tmp;
}

template<typename T>
T max (const FVec<T>& s)
{
    #if defined(HAVE_INTEL_IPP)
        T tmp = 0.;
        ippsMax( s.data(), (int)s.size(), &tmp);
    #elif defined(HAVE_ACCELERATE)
        T tmp = 0.;
        vDSP_maxv( s.data(), 1, &tmp, s.size() );
    #else
        size_t j;
        T tmp = s.data()[0];
        for (j = 1; j < s.size(); j++) {
            tmp = (tmp > s[j]) ? tmp : s[j];
        }
    #endif
    return tmp;
}

template<typename T>
T min (const FVec<T>& s)
{
    #if defined(HAVE_INTEL_IPP)
        T tmp = 0.;
        ippsMin(s.data(), (int)s.size(), &tmp);
    #elif defined(HAVE_ACCELERATE)
        T tmp = 0.;
        vDSP_minv(s.data(), 1, &tmp, s.size());
    #else
        size_t j;
        T tmp = s[0];
        #pragma omp simd
        for (j = 1; j < s.size(); j++) {
            tmp = (tmp < s[j]) ? tmp : s[j];
        }
    #endif
    return tmp;
}

template<typename T>
size_t min_elem (const FVec<T>& s)
{
#ifndef HAVE_ACCELERATE
  size_t j, pos = 0.;
  T tmp = s[0];
  #pragma omp simd
  for (j = 0; j < s.size(); j++) {
    pos = (tmp < s[j]) ? pos : j;
    tmp = (tmp < s[j]) ? tmp : s[j];
  }
#else
  T tmp = 0.;
  vDSP_Length pos = 0;
  vDSP_minvi(s, 1, &tmp, &pos, s.size());
#endif
  return )size_t)pos;
}

template<typename T>
size_t max_elem (const FVec<T>& s)
{
#ifndef HAVE_ACCELERATE
  size_t j, pos = 0;
  T tmp = 0.0;
  #pragma omp simd
  for (j = 0; j < s.size(); j++) {
    pos = (tmp > s[j]) ? pos : j;
    tmp = (tmp > s[j]) ? tmp : s[j];
  }
#else
  T_t tmp = 0.;
  vDSP_Length pos = 0;
  vDSP_maxvi(s.data(), 1, &tmp, &pos, s.size());
#endif
  return (fsize_t)pos;
}

template<typename T>
void shift (const FVec<T>& s)
{
  size_t half = s.size() / 2, start = half, j;
  // if length is odd, middle element is moved to the end
  if (2 * half < s.size()) start ++;
#ifndef HAVE_BLAS
  for (j = 0; j < half; j++) {
    std::swap(s[j], s[j + start]);
  }
#else
  cblas_swap(half, s.data(), 1, s.data() + start, 1);
#endif
  if (start != half) {
    for (j = 0; j < half; j++) {
      std::swap(s[j + start - 1], s[j + start]);
    }
  }
}

template<typename T>
void ishift (const FVec<T> & s)
{
  size_t half = s.size() / 2, start = half, j;
  // if length is odd, middle element is moved to the beginning
  if (2 * half < s.size()) start ++;
#ifndef HAVE_BLAS
    #pragma omp simd
    for (j = 0; j < half; j++) {
        std::swap(s[j], s[j + start]);
    }
#else
    cblas_swap(half, s->data, 1, s->data + start, 1);
#endif
    #pragma omp simd
    if (start != half) {
        for (j = 0; j < half; j++) {
            std::swap(s[half], s[j]);
        }
    }
}

template<typename T>
void clamp(const FVec<T>& in, T absmax) {
  size_t i;
  absmax = fabs(absmax);
  #pragma omp simd  
  for (i = 0; i < in->length; i++) in[i] = std::clamp(in[i],-absmax,absmax);  
}


template<typename T>
T level_lin (const FVec<T>& f)
{
    T energy = 0.;
    #ifndef HAVE_BLAS
        size_t j;
        #pragma omp simd
        for (j = 0; j < f.size(); j++) {
            energy += SQR (f[j]);
        }
    #else
        energy = cblas_dot(f.size(), f.data(), 1, f.data(), 1);
    #endif
    return energy / f.size();
}

template<typename T>
T local_hfc (const FVec<T>& v)
{
  T hfc = 0.;
  size_t j;
  #pragma omp simd
  for (j = 0; j < v->length; j++) {
    hfc += (j + 1) * v[j];
  }
  return hfc;
}

template<typename T>
void min_removal (const FVec<T>& v)
{
  T v_min = min(v);
  v += -v_min;  
}

template<typename T>
T alpha_norm (const FVec<T>& o, T alpha)
{
  size_t j;
  T tmp = 0.;
  #pragma omp simd
  for (j = 0; j < o.size(); j++) {
    tmp += std::pow(std::fabs(o[j]), alpha);
  }
  return std::pow(tmp / o.size(), 1. / alpha);
}

template<typename T>
void alpha_normalise (const FVec<T>& o, T alpha)
{
  size_t j;
  T norm = alpha_norm (o, alpha);
  o /= norm;  
}


template<typename T>
T median (const FVec<T>& input) {
  size_t n = input->length;
  T * arr = input.data();
  size_t low, high ;
  size_t median;
  size_t middle, ll, hh;

  low = 0 ; high = n-1 ; median = (low + high) / 2;
  #pragma omp simd
  for (;;) {
    if (high <= low) /* One element only */
      return arr[median] ;

    if (high == low + 1) {  /* Two elements only */
      if (arr[low] > arr[high])
        std::swap(arr[low], arr[high]) ;
      return arr[median] ;
    }

    /* Find median of low, middle and high items; swap into position low */
    middle = (low + high) / 2;
    if (arr[middle] > arr[high])    std::swap(arr[middle], arr[high]);
    if (arr[low]    > arr[high])    std::swap(arr[low],    arr[high]);
    if (arr[middle] > arr[low])     std::swap(arr[middle], arr[low]) ;

    /* Swap low item (now in position middle) into position (low+1) */
    std::swap(arr[middle], arr[low+1]) ;

    /* Nibble from each end towards middle, swapping items when stuck */
    ll = low + 1;
    hh = high;
    for (;;) {
      do ll++; while (arr[low] > arr[ll]) ;
      do hh--; while (arr[hh]  > arr[low]) ;

      if (hh < ll)
        break;

      std::swap(arr[ll], arr[hh]) ;
    }

    /* Swap middle item (in position low) back into correct position */
    std::swap(arr[low], arr[hh]) ;

    /* Re-set active partition */
    if (hh <= median)
      low = ll;
    if (hh >= median)
      high = hh - 1;
  }
}


template<typename T>
T moving_thres (const FVec<T>& vec, const FVec<T>& tmpvec,
    fsize_t post, fsize_t pre, fsize_t pos)
{
  fsize_t k;
  T *medar = (T *) tmpvec;
  fsize_t win_length = post + pre + 1;
  fsize_t length = vec.size();
  /* post part of the buffer does not exist */
  if (pos < post + 1) {
    for (k = 0; k < post + 1 - pos; k++)
      medar[k] = 0.;            /* 0-padding at the beginning */
    #pragma omp simd
    for (k = post + 1 - pos; k < win_length; k++)
      medar[k] = vec[k + pos - post];
    /* the buffer is fully defined */
  } else if (pos + pre < length) {
    #pragma omp simd
    for (k = 0; k < win_length; k++)
      medar[k] = vec[k + pos - post];
    /* pre part of the buffer does not exist */
  } else {
    #pragma omp simd
    for (k = 0; k < length - pos + post; k++)
      medar[k] = vec[k + pos - post];
    #pragma omp simd
    for (k = length - pos + post; k < win_length; k++)
      medar[k] = 0.;            /* 0-padding at the end */
  }
  return median (tmpvec);
}

template<typename T>
void fvec_adapt_thres(const FVec<T>& vec, const FVec<T>& tmp,
    size_t post, size_t pre) {
  size_t length = vec.size(), j;  
  for (j=0;j<length;j++) {
    vec[j] -= moving_thres(vec, tmp, post, pre, j);
  }
}

template<typename T>
T quadratic_peak_pos (const FVec<T>& x, size_t pos) {
  T s0, s1, s2; 
  size_t x0, x2;
  T half = .5, two = 2.;
  if (pos == 0 || pos == x.size() - 1) return pos;
  x0 = (pos < 1) ? pos : pos - 1;
  x2 = (pos + 1 < x.size()) ? pos + 1 : pos;
  if (x0 == pos) return (x[pos] <= x[x2]) ? pos : x2;
  if (x2 == pos) return (x[pos] <= x[x0]) ? pos : x0;
  s0 = x[x0];
  s1 = x[pos];
  s2 = x[x2];
  return pos + half * (s0 - s2 ) / (s0 - two * s1 + s2);
}


template<typename T>
T quadratic_peak_mag (const FVec<T>& *x, T pos) {
  T x0, x1, x2;
  size_t index = (size_t)(pos - .5) + 1;
  if (pos >= x.size() || pos < 0.) return 0.;
  if ((T)index == pos) return x[index];
  x0 = x[index - 1];
  x1 = x[index];
  x2 = x[index + 1];
  return x1 - .25 * (x0 - x2) * (pos - index);
}

template<typename T>
size_t peakpick(const FVec<T>& onset, size_t pos) {
  size_t tmp=0;
  tmp = (onset[pos] > onset[pos-1]
      &&  onset[pos] > onset[pos+1]
      &&  onset[pos] > 0.);
  return tmp;
}

template<typename T>
T quadfrac (T s0, T s1, T s2, T pf)
{
  T tmp =
      s0 + (pf / 2.) * (pf * (s0 - 2. * s1 + s2) - 3. * s0 + 4. * s1 - s2);
  return tmp;
}

template<typename T>
T freqtomidi (T freq)
{
  T midi;
  if (freq < 2. || freq > 100000.) return 0.; // avoid nans and infs
  /* log(freq/A-2)/log(2) */
  midi = freq / 6.875;
  midi = LOG (midi) / 0.6931471805599453;
  midi *= 12;
  midi -= 3;
  return midi;
}

template<typename T>
T miditofreq (T midi)
{
  T freq;
  if (midi > 140.) return 0.; // avoid infs
  freq = (midi + 3.) / 12.;
  freq = EXP (freq * 0.6931471805599453);
  freq *= 6.875;
  return freq;
}

template<typename T>
T bintofreq (T bin, T samplerate, T fftsize)
{
  T freq = samplerate / fftsize;
  return freq * MAX(bin, 0);
}

template<typename T>
T bintomidi (T bin, T samplerate, T fftsize)
{
  T midi = bintofreq (bin, samplerate, fftsize);
  return freqtomidi (midi);
}

template<typename T>
T freqtobin (T freq, T samplerate, T fftsize)
{
  T bin = fftsize / samplerate;
  return MAX(freq, 0) * bin;
}

template<typename T>
T miditobin (T midi, T samplerate, T fftsize)
{
  T freq = miditofreq (midi);
  return freqtobin (freq, samplerate, fftsize);
}


size_t is_power_of_two (size_t a)
{
  if ((a & (a - 1)) == 0) {
    return 1;
  } else {
    return 0;
  }
}


size_t next_power_of_two (size_t a)
{
  size_t i = 1;
  while (i < a) i <<= 1;
  return i;
}

size_t power_of_two_order (size_t a)
{
  int order = 0;
  int temp = next_power_of_two(a);
  while (temp >>= 1) {
    ++order;
  }
  return order;
}

template<typename T>
T db_spl (const FVec<T>& o)
{
  return 10. * LOG10 (level_lin (o));
}

template<typename T>
size_t silence_detection (const FVec<T>& * o, T threshold)
{
  return (db_spl (o) < threshold);
}

template<typename T>
T level_detection (const FVec<T>& * o, T threshold)
{
  T db_spl = db_spl (o);
  if (db_spl < threshold) {
    return 1.;
  } else {
    return db_spl;
  }
}

template<typename T>
T zero_crossing_rate (FVec<T>& * input)
{
  size_t j;
  size_t zcr = 0;
  for (j = 1; j < input->length; j++) {
    // previous was strictly negative
    if (input->data[j - 1] < 0.) {
      // current is positive or null
      if (input->data[j] >= 0.) {
        zcr += 1;
      }
      // previous was positive or null
    } else {
      // current is strictly negative
      if (input->data[j] < 0.) {
        zcr += 1;
      }
    }
  }
  return zcr / (T) input->length;
}

template<typename T>
void autocorr (const FVec<T>& * input, FVec<T>& * output)
{
  size_t i, j, length = input->length;
  T *data, *acf;
  T tmp = 0;
  data = input->data;
  acf = output->data;
  for (i = 0; i < length; i++) {
    tmp = 0.;
    for (j = i; j < length; j++) {
      tmp += data[j - i] * data[j];
    }
    acf[i] = tmp / (T) (length - i);
  }
}

template<typename T>
FVec<T> create_window (size_t n, window_type wintype) {
    FVec<T> win(n);
    T * w = win.data();
    size_t i, size =n;
  
    switch(wintype) {
    case win_ones:
        win.ones();      
        break;
    case win_rectangle:
        win.fill(.5);
        break;
    case win_hamming:
        #pragma omp simd
        for (i=0;i<size;i++)
            w[i] = 0.54 - 0.46 * std::cos(TWO_PI * i / (size));
        break;
    case win_hanning:
        #pragma omp simd
        for (i=0;i<size;i++)
            w[i] = 0.5 - (0.5 * std::cos(TWO_PI * i / (size)));
        break;
    case win_hanningz:
        #pragma omp simd
        for (i=0;i<size;i++)
            w[i] = 0.5 * (1.0 - std::(TWO_PI * i / (size)));
        break;
    case win_blackman:
        #pragma omp simd
        for (i=0;i<size;i++)
            w[i] = 0.42
            - 0.50 * std::cos(    TWO_PI*i/(size-1.0))
            + 0.08 * std::cos(2.0*TWO_PI*i/(size-1.0));
        break;
    case win_blackman_harris:
    #pragma omp simd
        for (i=0;i<size;i++)
            w[i] = 0.35875
            - 0.48829 * std::cos(    TWO_PI*i/(size-1.0))
            + 0.14128 * std::cos(2.0*TWO_PI*i/(size-1.0))
            - 0.01168 * std::cos(3.0*TWO_PI*i/(size-1.0));
        break;
    case win_gaussian:
        {
        T a, b, c = 0.5;
        size_t n;
        #pragma omp simd
        for (n = 0; n < size; n++)
        {
            a = (n-c*(size-1))/(SQR(c)*(size-1));
            b = -c*SQR(a);
            w[n] = std::exp(b);
        }
        }
        break;
    case win_welch:
        #pragma omp simd
        for (i=0;i<size;i++)
            w[i] = 1.0 - SQR((2.*i-size)/(size+1.0));
        break;
    case win_parzen:
        #pragma omp simd
        for (i=0;i<size;i++)
            w[i] = 1.0 - std::fabs((2.f*i-size)/(size+1.0f));
        break;
    default:
        break;
    }  
}

T
hztomel (T freq)
{
  const T lin_space = 3./200.;
  const T split_hz = 1000.;
  const T split_mel = split_hz * lin_space;
  const T log_space = 27./LOG(6400/1000.);
  if (freq < 0) {
    AUBIO_WRN("hztomel: input frequency should be >= 0\n");
    return 0;
  }
  if (freq < split_hz)
  {
    return freq * lin_space;
  } else {
    return split_mel + log_space * LOG (freq / split_hz);
  }

}

template<typename T>
T meltohz (T mel)
{
  const T lin_space = 200./3.;
  const T split_hz = 1000.;
  const T split_mel = split_hz / lin_space;
  const T logSpacing = std::pow(6400/1000., 1/27.);
  if (mel < 0) {
    std::cerr << "meltohz: input mel should be >= 0\n";
    return 0;
  }
  if (mel < split_mel) {
    return lin_space * mel;
  } else {
    return split_hz * std::pow(logSpacing, mel - split_mel);
  }
}

template<typename T>
T hztomel_htk (T freq)
{
  const T split_hz = 700.;
  const T log_space = 1127.;
  if (freq < 0) {
    std::cerr << "hztomel_htk: input frequency should be >= 0\n";
    return 0;
  }
  return log_space * std::log(1 + freq / split_hz);
}

template<typename T>
T meltohz_htk (T mel)
{
  const T split_hz = 700.;
  const T log_space = 1./1127.;
  if (mel < 0) {
    std::cerr << "meltohz_htk: input frequency should be >= 0\n";
    return 0;
  }
  return split_hz * (std::exp( mel * log_space) - 1.);
}
