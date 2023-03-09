#pragma once

namespace AudioDSP
{
    
    template<typename T>
    struct complex_vector : public std::vector<std::complex<T>,Allocator::aligned_allocator<std::complex<T>,64>>
    {
        using base = std::vector<std::complex<T>,Allocator::aligned_allocator<std::complex<T>,64>>;
        complex_vector() = default;
        complex_vector(size_t i) : base(i) {}


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
        
        void fill(const std::complex<T>& v) {
            std::fill(begin(),end(),v);
        }
        void print() {
            std::cout << "VECTOR[" << size() << "]";
            for(size_t i = 0; i < size(); i++)
                std::cout << (*this)[i] << ",";
            std::cout << std::endl;
        }
        void copy(size_t n, std::complex<T> * p)
        {
            resize(n);
            std::memcpy(data(),p,n*sizeof(std::complex<T>));
        }
        #ifdef SWIG
        %extend {
            std::complex<T> __getitem__(size_t i ) { return (*$self)[i]; }
            void __setitem__(size_t i, const std::complex<T> & v) { (*$self)[i] = v; }
            std::complex<T>* data() { return $self->data(); }

            complex_vector<T> __add__(const complex_vector<T> & v)
            {
                assert($self->size() == v.size());
                complex_vector<T> r(v.size());
                #pragma omp simd
                for(size_t i = 0; i < v.size(); i++) r[i] = (*$self)[i] + v[i];
                return r;
            }
        
            complex_vector<T> __sub__(const complex_vector<T> & v)
            {
                assert($self->size() == v.size());
                complex_vector<T> r(v.size());
                #pragma omp simd
                for(size_t i = 0; i < v.size(); i++) r[i] = (*$self)[i] - v[i];
                return r;
            }
            complex_vector<T> __mul__(const complex_vector<T> & v)
            {
                assert($self->size() == v.size());
                complex_vector<T> r(v.size());
                #pragma omp simd
                for(size_t i = 0; i < v.size(); i++) r[i] = (*$self)[i] * v[i];
                return r;
            }
            complex_vector<T> __div__(const complex_vector<T> & v)
            {
                assert($self->size() == v.size());
                complex_vector<T> r(v.size());
                #pragma omp simd
                for(size_t i = 0; i < v.size(); i++) r[i] = (*$self)[i] / v[i];
                return r;
            }
            complex_vector<T> __pow__(const complex_vector<T> & v)
            {
                assert($self->size() == v.size());
                complex_vector<T> r(v.size());
                #pragma omp simd
                for(size_t i = 0; i < v.size(); i++) r[i] = std::pow((*$self)[i],v[i]);
                return r;
            }
            complex_vector<T> __add__(const T x)
            {
                assert($self->size() == $self->size());
                complex_vector<T> r($self->size());
                #pragma omp simd
                for(size_t i = 0; i < $self->size(); i++) r[i] = (*$self)[i] + x;
                return r;
            }
            complex_vector<T> __sub__(const T x)
            {
                assert($self->size() == $self->size());
                complex_vector<T> r($self->size());
                #pragma omp simd
                for(size_t i = 0; i < $self->size(); i++) r[i] = (*$self)[i] - x;
                return r;
            }
            complex_vector<T> __mul__(const T x)
            {
                assert($self->size() == $self->size());
                complex_vector<T> r($self->size());
                #pragma omp simd
                for(size_t i = 0; i < $self->size(); i++) r[i] = (*$self)[i] * x;
                return r;
            }
            complex_vector<T> __div__(const T x)
            {
                assert($self->size() == $self->size());
                complex_vector<T> r($self->size());
                #pragma omp simd
                for(size_t i = 0; i < $self->size(); i++) r[i] = (*$self)[i] / x;
                return r;
            }
            complex_vector<T> __pow__(const T x)
            {
                assert($self->size() == $self->size());
                complex_vector<T> r($self->size());
                #pragma omp simd
                for(size_t i = 0; i < $self->size(); i++) r[i] = std::pow((*$self)[i],x);
                return r;
            }
            complex_vector<T> __neg__()
            {
                assert($self->size() == $self->size());
                complex_vector<T> r($self->size());
                #pragma omp simd
                for(size_t i = 0; i < $self->size(); i++) r[i] = -(*$self)[i];
                return r;
            }
            //bool __eq__(const complex_vector<T> & v);
        }
        #endif
    };

    template<class T>
    complex_vector<T> operator +(const complex_vector<T> & a, const complex_vector<T> & b) {
        complex_vector<T> r(a.size());
        #pragma omp simd
        for(size_t i = 0; i < a.size(); i++) r[i] = a[i] + b[i];
        return r;
    }   
    template<class T>
    complex_vector<T> operator -(const complex_vector<T> & a, const complex_vector<T> & b) {
        complex_vector<T> r(a.size());
        #pragma omp simd
        for(size_t i = 0; i < a.size(); i++) r[i] = a[i] - b[i];
        return r;
    }   
    template<class T>
    complex_vector<T> operator *(const complex_vector<T> & a, const complex_vector<T> & b) {
        complex_vector<T> r(a.size());
        #pragma omp simd
        for(size_t i = 0; i < a.size(); i++) r[i] = a[i] * b[i];
        return r;
    }   
    template<class T>
    complex_vector<T> operator /(const complex_vector<T> & a, const complex_vector<T> & b) {
        complex_vector<T> r(a.size());
        #pragma omp simd
        for(size_t i = 0; i < a.size(); i++) r[i] = a[i] / b[i];
        return r;
    }   
    template<class T>
    complex_vector<T> operator %(const complex_vector<T> & a, const complex_vector<T> & b) {
        complex_vector<T> r(a.size());
        #pragma omp simd
        for(size_t i = 0; i < a.size(); i++) r[i] = std::fmod(a[i],b[i]);
        return r;
    }   
    template<class T>
    complex_vector<T> operator +(const complex_vector<T> & a, const std::complex<T>& b) {
        complex_vector<T> r(a.size());
        #pragma omp simd
        for(size_t i = 0; i < a.size(); i++) r[i] = a[i] + b;
        return r;
    }   
    template<class T>
    complex_vector<T> operator -(const complex_vector<T> & a, const std::complex<T>& b) {
        complex_vector<T> r(a.size());
        #pragma omp simd
        for(size_t i = 0; i < a.size(); i++) r[i] = a[i] - b;
        return r;
    }   
    template<class T>
    complex_vector<T> operator *(const complex_vector<T> & a, const std::complex<T>& b) {
        complex_vector<T> r(a.size());
        #pragma omp simd
        for(size_t i = 0; i < a.size(); i++) r[i] = a[i] * b;
        return r;
    }   
    template<class T>
    complex_vector<T> operator / (const complex_vector<T> & a, const std::complex<T>& b) {
        complex_vector<T> r(a.size());
        #pragma omp simd
        for(size_t i = 0; i < a.size(); i++) r[i] = a[i] / b;
        return r;
    }   
    template<class T>
    complex_vector<T> operator %(const complex_vector<T> & a, const std::complex<T>& b) {
        complex_vector<T> r(a.size());
        #pragma omp simd
        for(size_t i = 0; i < a.size(); i++) r[i] = std::fmod(a[i],b);
        return r;
    }   
    template<class T>
    complex_vector<T> operator +(const std::complex<T> & a, const complex_vector<T>& b) {
        complex_vector<T> r(a.size());
        #pragma omp simd
        for(size_t i = 0; i < a.size(); i++) r[i] = a + b[i];
        return r;
    }   
    template<class T>
    complex_vector<T> operator -(const std::complex<T> & a, const complex_vector<T>& b) {
        complex_vector<T> r(a.size());
        #pragma omp simd
        for(size_t i = 0; i < a.size(); i++) r[i] = a - b[i];
        return r;
    }   
    template<class T>
    complex_vector<T> operator *(const std::complex<T> & a, const complex_vector<T>& b) {
        complex_vector<T> r(a.size());
        #pragma omp simd
        for(size_t i = 0; i < a.size(); i++) r[i] = a * b[i];
        return r;
    }   
    template<class T>
    complex_vector<T> operator /(const std::complex<T> & a, const complex_vector<T>& b) {
        complex_vector<T> r(a.size());
        #pragma omp simd
        for(size_t i = 0; i < a.size(); i++) r[i] = a / b[i];
        return r;
    }   
    template<class T>
    std::complex<T> sum(const complex_vector<T>& v ) {
        std::complex<T> r(0,0);
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r += v[i];
        return r;
    }   
    template<class T>
    std::complex<T> product(const complex_vector<T>& v ) {
        std::complex<T> r(1,0);
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r *= v[i];
        return r;
    }   
    template<class T>
    complex_vector<T> pow(const complex_vector<T> & v, const complex_vector<T> & x) {
        complex_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::pow(v[i],x[i]);
        return r;
    }    
    template<class T>
    complex_vector<T> pow(const complex_vector<T> & v, const std::complex<T> x) {
        complex_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::pow(v[i],x);
        return r;
    }    
    template<class T>
    complex_vector<T> pow(const std::complex<T> x, const complex_vector<T> & v) {
        complex_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::pow(x,v[i]);
        return r;
    }    
    template<class T>
    complex_vector<T> sqrt(const complex_vector<T> & v) {
        complex_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::sqrt(v[i]);
        return r;
    }    
      template<class T>
    complex_vector<T> abs(const complex_vector<T> & v) {
        complex_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::fabs(v[i]);
        return r;
    }
    template<class T>
    complex_vector<T> cube(const complex_vector<T> & v) {
        complex_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = v[i]*v[i]*v[i];
        return r;
    }
    template<class T>
    complex_vector<T> sqr(const complex_vector<T> & v) {
        complex_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = v[i]*v[i];
        return r;
    }
    template<class T>
    complex_vector<T> cos(const complex_vector<T> & v) {
        complex_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::cos(v[i]);
        return r;
    }
    template<class T>
    complex_vector<T> sin(const complex_vector<T> & v) {
        complex_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::sin(v[i]);
        return r;
    }    
    template<class T>
    complex_vector<T> tan(const complex_vector<T> & v) {
        complex_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::tan(v[i]);
        return r;
    }

    template<class T>
    complex_vector<T> acos(const complex_vector<T> & v) {
        complex_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::acos(v[i]);
        return r;
    }
    template<class T>
    complex_vector<T> asin(const complex_vector<T> & v) {
        complex_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::asin(v[i]);
        return r;
    }    
    template<class T>
    complex_vector<T> atan(const complex_vector<T> & v) {
        complex_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::atan(v[i]);
        return r;
    }        
    template<class T>
    complex_vector<T> cosh(const complex_vector<T> & v) {
        complex_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::cosh(v[i]);
        return r;
    }
    template<class T>
    complex_vector<T> sinh(const complex_vector<T> & v) {
        complex_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::sinh(v[i]);
        return r;
    }    
    template<class T>
    complex_vector<T> tanh(const complex_vector<T> & v) {
        complex_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::tanh(v[i]);
        return r;
    }

    template<class T>
    complex_vector<T> acosh(const complex_vector<T> & v) {
        complex_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::acosh(v[i]);
        return r;
    }
    template<class T>
    complex_vector<T> asinh(const complex_vector<T> & v) {
        complex_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::asinh(v[i]);
        return r;
    }    
    template<class T>
    complex_vector<T> atanh(const complex_vector<T> & v) {
        complex_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::atanh(v[i]);
        return r;
    }    

    template<class T>
    complex_vector<T> exp(const complex_vector<T> & v) {
        complex_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::exp(v[i]);
        return r;
    }    
    template<class T>
    complex_vector<T> log(const complex_vector<T> & v) {
        complex_vector<T> r(v.size());
        #pragma omp simd
        for(size_t i = 0; i < v.size(); i++) r[i] = std::log(v[i]);
        return r;
    }    
    template<class T>
    sample_vector<T> real(const complex_vector<T> & v) {
        sample_vector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = v[i].real();
        return r;
    }
    template<class T>
    sample_vector<T> imag(const complex_vector<T> & v) {
        sample_vector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = v[i].imag();
        return r;
    }
    template<class T>
    sample_vector<T> arg(const complex_vector<T> & v) {
        sample_vector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::arg(v[i]);
        return r;
    }
    template<class T>
    sample_vector<T> norm(const complex_vector<T> & v) {
        sample_vector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::norm(v[i]());
        return r;
    }
    template<class T>
    complex_vector<T> conj(const complex_vector<T> & v) {
        complex_vector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::conj(v[i]);
        return r;
    }
    template<class T>
    complex_vector<T> proj(const complex_vector<T> & v) {
        complex_vector<T> r(v.size());
        for(size_t i = 0; i < v.size(); i++) r[i] = std::proj(v[i]);
        return r;
    }
    template<class T>
    complex_vector<T> polar(const sample_vector<T> & a, const sample_vector<T> & b) {
        complex_vector<T> r(a.size());
        for(size_t i = 0; i < a.size(); i++) r[i] = std::polar(a[i],b[i]);
        return r;
    }
}