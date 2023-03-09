#pragma once

namespace AudioDSP
{
    template<typename T>
    struct sample_vector : public std::vector<T,Allocator::aligned_allocator<T,64>>
    {
        using base = std::vector<T,Allocator::aligned_allocator<T,64>>;
        sample_vector() = default;
        sample_vector(size_t i) : base(i) {}

        using base::operator [];        
        using base::begin;
        using base::end;
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

        void copy(size_t n, T * p)
        {
            resize(n);
            std::memcpy(data(),p,n*sizeof(T));
        }
        void fill(const T v) {
            std::fill(begin(),end(),v);
        }
        
        void print() {
            std::cout << "VECTOR[" << size() << "]";
            for(size_t i = 0; i < size(); i++)
                std::cout << (*this)[i] << ",";
            std::cout << std::endl;
        }

        #ifdef SWIG
        %extend {
            T __getitem__(size_t i) { return (*$self)[i-1]; }
            void __setitem(size_t i, const T v) { (*$self)[i-1] = v; }
                                    
            sample_vector<T> __add__(const sample_vector<T> & v)
            {
                assert($self->size() == v.size());
                sample_vector<T> r(v.size());
                #pragma omp simd
                for(size_t i = 0; i < v.size(); i++) r[i] = (*$self)[i] + v[i];
                return r;
            }
        
            sample_vector<T> __sub__(const sample_vector<T> & v)
            {
                assert($self->size() == v.size());
                sample_vector<T> r(v.size());
                #pragma omp simd
                for(size_t i = 0; i < v.size(); i++) r[i] = (*$self)[i] - v[i];
                return r;
            }
            sample_vector<T> __mul__(const sample_vector<T> & v)
            {
                assert($self->size() == v.size());
                sample_vector<T> r(v.size());
                #pragma omp simd
                for(size_t i = 0; i < v.size(); i++) r[i] = (*$self)[i] * v[i];
                return r;
            }
            sample_vector<T> __div__(const sample_vector<T> & v)
            {
                assert($self->size() == v.size());
                sample_vector<T> r(v.size());
                #pragma omp simd
                for(size_t i = 0; i < v.size(); i++) r[i] = (*$self)[i] / v[i];
                return r;
            }
            sample_vector<T> __pow__(const sample_vector<T> & v)
            {
                assert($self->size() == v.size());
                sample_vector<T> r(v.size());
                #pragma omp simd
                for(size_t i = 0; i < v.size(); i++) r[i] = std::pow((*$self)[i],v[i]);
                return r;
            }
            sample_vector<T> __add__(const T x)
            {
                assert($self->size() == $self->size());
                sample_vector<T> r($self->size());
                #pragma omp simd
                for(size_t i = 0; i < $self->size(); i++) r[i] = (*$self)[i] + x;
                return r;
            }
            sample_vector<T> __sub__(const T x)
            {
                assert($self->size() == $self->size());
                sample_vector<T> r($self->size());
                #pragma omp simd
                for(size_t i = 0; i < $self->size(); i++) r[i] = (*$self)[i] - x;
                return r;
            }
            sample_vector<T> __mul__(const T x)
            {
                assert($self->size() == $self->size());
                sample_vector<T> r($self->size());
                #pragma omp simd
                for(size_t i = 0; i < $self->size(); i++) r[i] = (*$self)[i] * x;
                return r;
            }
            sample_vector<T> __div__(const T x)
            {
                assert($self->size() == $self->size());
                sample_vector<T> r($self->size());
                #pragma omp simd
                for(size_t i = 0; i < $self->size(); i++) r[i] = (*$self)[i] / x;
                return r;
            }
            sample_vector<T> __pow__(const T x)
            {
                assert($self->size() == $self->size());
                sample_vector<T> r($self->size());
                #pragma omp simd
                for(size_t i = 0; i < $self->size(); i++) r[i] = std::pow((*$self)[i],x);
                return r;
            }
            sample_vector<T> __neg__()
            {
                assert($self->size() == $self->size());
                sample_vector<T> r($self->size());
                #pragma omp simd
                for(size_t i = 0; i < $self->size(); i++) r[i] = -(*$self)[i];
                return r;
            }
            //bool __eq__(const sample_vector<T> & v);
        }
        #endif       
    };

    template<class T>
    sample_vector<T> abs(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::fabs(v[i]);
        return r;
    }
    template<class T>
    sample_vector<T> cube(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = v[i]*v[i]*v[i];
        return r;
    }
    template<class T>
    sample_vector<T> sqr(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = v[i]*v[i];
        return r;
    }
    template<class T>
    sample_vector<T> cos(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::cos(v[i]);
        return r;
    }
    template<class T>
    sample_vector<T> sin(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::sin(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> tan(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::tan(v[i]);
        return r;
    }

    template<class T>
    sample_vector<T> acos(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::acos(v[i]);
        return r;
    }
    template<class T>
    sample_vector<T> asin(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::asin(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> atan(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::atan(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> atan2(const sample_vector<T> & v, const T value) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::atan2(v[i], value);
        return r;
    }    
    template<class T>
    sample_vector<T> atan2(const sample_vector<T> & v, const sample_vector<T> value) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::atan2(v[i], value[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> cosh(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::cosh(v[i]);
        return r;
    }
    template<class T>
    sample_vector<T> sinh(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::sinh(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> tanh(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::tanh(v[i]);
        return r;
    }

    template<class T>
    sample_vector<T> acosh(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::acosh(v[i]);
        return r;
    }
    template<class T>
    sample_vector<T> asinh(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::asinh(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> atanh(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::atanh(v[i]);
        return r;
    }    

    template<class T>
    sample_vector<T> exp(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::exp(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> log(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::log(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> log10(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::log10(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> exp2(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::exp2(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> expm1(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::expm1(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> ilogb(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::ilogb(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> log2(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::log2(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> log1p(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::log1p(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> logb(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::logb(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> scalbn(const sample_vector<T> & v, const sample_vector<int> & x) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::scalbn(v[i],x[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> scalbn(const sample_vector<T> & v, const int x) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::scalbn(v[i],x);
        return r;
    }    
    template<class T>
    sample_vector<T> scalbln(const sample_vector<T> & v, const sample_vector<long int> & x) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::scalbln(v[i],x[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> scalbln(const sample_vector<T> & v, const long int x) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::scalbln(v[i],x);
        return r;
    }    
    template<class T>
    sample_vector<T> pow(const sample_vector<T> & v, const sample_vector<T> & x) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::pow(v[i],x[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> pow(const sample_vector<T> & v, const T x) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::pow(v[i],x);
        return r;
    }    
    template<class T>
    sample_vector<T> pow(const T x, const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::pow(x,v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> sqrt(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::sqrt(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> cbrt(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::cbrt(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> hypot(const sample_vector<T> & v, const sample_vector<T> & x) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::hypot(v[i],x[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> hypot(const sample_vector<T> & v, const T x) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::hypot(v[i],x);
        return r;
    }    
    template<class T>
    sample_vector<T> hypot(const T x, const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::hypot(x,v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> erf(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::erf(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> erfc(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::erfc(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> tgamma(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::tgamma(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> lgamma(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::lgamma(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> ceil(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::ceil(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> floor(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::floor(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> fmod(const sample_vector<T> & v, const sample_vector<T> & x) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::fmod(v[i],x[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> fmod(const sample_vector<T> & v, const T x) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::fmod(v[i],x);
        return r;
    }    
    template<class T>
    sample_vector<T> fmod(const T x, const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::fmod(x,v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> trunc(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::trunc(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> round(const sample_vector<T> & v) {
        sample_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::round(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<long int> lround(const sample_vector<T> & v) {
        sample_vector<long int> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::lround(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<long long int> llround(const sample_vector<T> & v) {
        sample_vector<long long int> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::llround(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> nearbyint(const sample_vector<T> & v) {
        sample_vector<long long int> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::nearbyint(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> remainder(const sample_vector<T> & a, const sample_vector<T> & b) {
        sample_vector<T> r(a.size());
        #pragma omp simd
        for(size_t i = 0; i < a.size(); i++) r[i] = std::remainder(a[i],b[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> copysign(const sample_vector<T> & a, const sample_vector<T> & b) {
        sample_vector<T> r(a.size());
        #pragma omp simd
        for(size_t i = 0; i < a.size(); i++) r[i] = std::copysign(a[i],b[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> fdim(const sample_vector<T> & a, const sample_vector<T> & b) {
        sample_vector<T> r(a.size());
        #pragma omp simd
        for(size_t i = 0; i < a.size(); i++) r[i] = std::fdim(a[i],b[i]);
        return r;
    }    
    #undef fmax
    template<class T>
    sample_vector<T> fmax(const sample_vector<T> & a, const sample_vector<T> & b) {
        sample_vector<T> r(a.size());
        #pragma omp simd
        for(size_t i = 0; i < a.size(); i++) r[i] = std::fmax(a[i],b[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> fmin(const sample_vector<T> & a, const sample_vector<T> & b) {
        sample_vector<T> r(a.size());
        #pragma omp simd
        for(size_t i = 0; i < a.size(); i++) r[i] = std::fmin(a[i],b[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> fma(const sample_vector<T> & a, const sample_vector<T> & b, const sample_vector<T> & c) {
        sample_vector<T> r(a.size());
        #pragma omp simd
        for(size_t i = 0; i < a.size(); i++) r[i] = std::fma(a[i],b[i],c[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> fabs(const sample_vector<T> & a) {
        sample_vector<T> r(a.size());
        #pragma omp simd
        for(size_t i = 0; i < a.size(); i++) r[i] = std::fabs(a[i]);
        return r;
    }    

    template<class T>
    sample_vector<T> operator +(const sample_vector<T> & a, const sample_vector<T> & b) {
        sample_vector<T> r(a.size());
        #pragma omp simd
        for(size_t i = 0; i < a.size(); i++) r[i] = a[i] + b[i];
        return r;
    }   
    template<class T>
    sample_vector<T> operator -(const sample_vector<T> & a, const sample_vector<T> & b) {
        sample_vector<T> r(a.size());
        #pragma omp simd
        for(size_t i = 0; i < a.size(); i++) r[i] = a[i] - b[i];
        return r;
    }   
    template<class T>
    sample_vector<T> operator *(const sample_vector<T> & a, const sample_vector<T> & b) {
        sample_vector<T> r(a.size());
        #pragma omp simd
        for(size_t i = 0; i < a.size(); i++) r[i] = a[i] * b[i];
        return r;
    }   
    template<class T>
    sample_vector<T> operator /(const sample_vector<T> & a, const sample_vector<T> & b) {
        sample_vector<T> r(a.size());
        #pragma omp simd
        for(size_t i = 0; i < a.size(); i++) r[i] = a[i] / b[i];
        return r;
    }   
    template<class T>
    sample_vector<T> operator %(const sample_vector<T> & a, const sample_vector<T> & b) {
        sample_vector<T> r(a.size());
        #pragma omp simd
        for(size_t i = 0; i < a.size(); i++) r[i] = std::fmod(a[i],b[i]);
        return r;
    }   
    template<class T>
    sample_vector<T> operator +(const sample_vector<T> & a, const T& b) {
        sample_vector<T> r(a.size());
        #pragma omp simd
        for(size_t i = 0; i < a.size(); i++) r[i] = a[i] + b;
        return r;
    }   
    template<class T>
    sample_vector<T> operator -(const sample_vector<T> & a, const T& b) {
        sample_vector<T> r(a.size());
        #pragma omp simd
        for(size_t i = 0; i < a.size(); i++) r[i] = a[i] - b;
        return r;
    }   
    template<class T>
    sample_vector<T> operator *(const sample_vector<T> & a, const T& b) {
        sample_vector<T> r(a.size());
        #pragma omp simd
        for(size_t i = 0; i < a.size(); i++) r[i] = a[i] * b;
        return r;
    }   
    template<class T>
    sample_vector<T> operator / (const sample_vector<T> & a, const T& b) {
        sample_vector<T> r(a.size());
        #pragma omp simd
        for(size_t i = 0; i < a.size(); i++) r[i] = a[i] / b;
        return r;
    }   
    template<class T>
    sample_vector<T> operator %(const sample_vector<T> & a, const T& b) {
        sample_vector<T> r(a.size());
        #pragma omp simd
        for(size_t i = 0; i < a.size(); i++) r[i] = std::fmod(a[i],b);
        return r;
    }   
    template<class T>
    sample_vector<T> operator +(const T & a, const sample_vector<T>& b) {
        sample_vector<T> r(a.size());
        #pragma omp simd
        for(size_t i = 0; i < a.size(); i++) r[i] = a + b[i];
        return r;
    }   
    template<class T>
    sample_vector<T> operator -(const T & a, const sample_vector<T>& b) {
        sample_vector<T> r(a.size());
        #pragma omp simd
        for(size_t i = 0; i < a.size(); i++) r[i] = a - b[i];
        return r;
    }   
    template<class T>
    sample_vector<T> operator *(const T & a, const sample_vector<T>& b) {
        sample_vector<T> r(a.size());
        #pragma omp simd
        for(size_t i = 0; i < a.size(); i++) r[i] = a * b[i];
        return r;
    }   
    template<class T>
    sample_vector<T> operator /(const T & a, const sample_vector<T>& b) {
        sample_vector<T> r(a.size());
        #pragma omp simd
        for(size_t i = 0; i < a.size(); i++) r[i] = a / b[i];
        return r;
    }   
    template<class T>
    T sum(const sample_vector<T>& v ) {
        T r = 0;
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r += v[i];
        return r;
    }   
    template<class T>
    T product(const sample_vector<T>& v ) {
        T r = 1;
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r *= v[i];
        return r;
    }   

    
}