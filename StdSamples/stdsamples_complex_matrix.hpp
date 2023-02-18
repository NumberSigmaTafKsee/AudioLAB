#pragma once

namespace AudioDSP
{
    template<typename T>
    struct complex_matrix : public std::vector<std::complex<T>,Allocator::aligned_allocator<std::complex<T>,64>>
    {
        using base = std::vector<std::complex<T>,Allocator::aligned_allocator<std::complex<T>,64>>;
        size_t M,N;
        complex_matrix() = default;
        complex_matrix(size_t i, size_t j) : base(i) {
            M = i;
            N = j;
            resize(M*N);
        }

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

        std::complex<T>& operator()(size_t i, size_t j) { return (*this)[i*N+j]; };
        std::complex<T>  operator()(size_t i, size_t j) const { return (*this)[i*N+j];};

        size_t rows() const { return M; }
        size_t cols() const { return N; }
        void   resize(size_t i, size_t j) {
            M = i;
            N = j;
            resize(i*j);
        }
        void copy(size_t m, size_t n, std::complex<T> * p)
        {
            resize(m,n);
            std::memcpy(data(),p,m*n*sizeof(std::complex<T>));
        }
        #ifdef SWIG
        %extend {
            complex_matrix<T> __add__(const complex_matrix<T> & v)
            {
                assert($self->size() == v.size());
                complex_matrix<T> r($self->rows(),$self->cols());
                #pragma omp simd
                for(size_t i = 0; i < v.size(); i++) r[i] = (*$self)[i] + v[i];
                return r;
            }
        
            complex_matrix<T> __sub__(const complex_matrix<T> & v)
            {
                assert($self->size() == v.size());
                complex_matrix<T> r($self->rows(),$self->cols());
                #pragma omp simd
                for(size_t i = 0; i < v.size(); i++) r[i] = (*$self)[i] - v[i];
                return r;
            }
            complex_matrix<T> __mul__(const complex_matrix<T> & v)
            {
                assert($self->size() == v.size());
                complex_matrix<T> r($self->rows(),$self->cols());
                #pragma omp parallel for simd                
                for(size_t i = 0; i < $self->rows(); i++)
                {                    
                    for(size_t j = 0; j < v.cols(); j++)
                    {
                        std::complex<T> sum(0,0);
                        for(size_t k = 0; k < $self->cols(); k++)
                        {
                            sum += (*$self)(i,k) * v(i,j);
                        }
                        r(i,j) = sum;
                    }
                }
                return r;
            }
            complex_matrix<T> __div__(const complex_matrix<T> & v)
            {
                assert($self->size() == v.size());
                complex_matrix<T> r($self->rows(),$self->cols());
                #pragma omp simd
                for(size_t i = 0; i < v.size(); i++) r[i] = (*$self)[i] / v[i];
                return r;
            }
            complex_matrix<T> __pow__(const complex_matrix<T> & v)
            {
                assert($self->size() == v.size());
                complex_matrix<T> r($self->rows(),$self->cols());
                #pragma omp simd
                for(size_t i = 0; i < v.size(); i++) r[i] = std::pow((*$self)[i],v[i]);
                return r;
            }
            complex_matrix<T> __add__(const std::complex<T>& x)
            {
                assert($self->size() == $self->size());
                complex_matrix<T> r($self->rows(),$self->cols());
                #pragma omp simd
                for(size_t i = 0; i < $self->size(); i++) r[i] = (*$self)[i] + x;
                return r;
            }
            complex_matrix<T> __sub__(const std::complex<T>& x)
            {
                assert($self->size() == $self->size());
                complex_matrix<T> r($self->rows(),$self->cols());
                #pragma omp simd
                for(size_t i = 0; i < $self->size(); i++) r[i] = (*$self)[i] - x;
                return r;
            }
            complex_matrix<T> __mul__(const std::complex<T>& x)
            {
                assert($self->size() == $self->size());
                complex_matrix<T> r($self->rows(),$self->cols());
                #pragma omp simd
                for(size_t i = 0; i < $self->size(); i++) r[i] = (*$self)[i] * x;
                return r;
            }
            complex_matrix<T> __div__(const std::complex<T>& x)
            {
                assert($self->size() == $self->size());
                complex_matrix<T> r($self->rows(),$self->cols());
                #pragma omp simd
                for(size_t i = 0; i < $self->size(); i++) r[i] = (*$self)[i] / x;
                return r;
            }
            complex_matrix<T> __pow__(const std::complex<T>& x)
            {
                assert($self->size() == $self->size());
                complex_matrix<T> r($self->rows(),$self->cols());
                #pragma omp simd
                for(size_t i = 0; i < $self->size(); i++) r[i] = std::pow((*$self)[i],x);
                return r;
            }
            complex_matrix<T> __neg__()
            {
                assert($self->size() == $self->size());
                complex_matrix<T> r($self->rows(),$self->cols());
                #pragma omp simd
                for(size_t i = 0; i < $self->size(); i++) r[i] = -(*$self)[i];
                return r;
            }
            //bool __eq__(const complex_matrix<T> & v);
        }
        #endif

        void fill(const std::complex<T>& v) {
            std::fill(begin(),end(),v);
        }
        void print() {
            std::cout << "Matrix(" << M << "," << N << ")=\n";
            for(size_t i = 0; i < M; i++)
            {
                for(size_t j = 0; j < N; j++)
                    std::cout << (*this)(i,j) << ",";
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }
    };

    template<typename T>
    complex_matrix<T> matmul(const complex_matrix<T> & a, const complex_matrix<T> & b)
    {
        assert(a.cols() == b.rows());
        complex_matrix<T> r(a.rows(),b.cols());
        #pragma omp parallel for simd                
        for(size_t i = 0; i < a.rows(); i++)
        {                    
            for(size_t j = 0; j < b.cols(); j++)
            {
                std::complex<T> sum = 0;
                for(size_t k = 0; k < a.cols(); k++)
                {
                    sum += a(i,k) * b(i,j);
                }
                r(i,j) = sum;
            }
        }
        return r;
    }
    template<class T>
    complex_matrix<T> operator +(const complex_matrix<T> & a, const complex_matrix<T> & b) {
        complex_matrix<T> r(a.rows(),a.cols());
        #pragma omp simd
        for(size_t i = 0; i < a.rows(); i++) 
        for(size_t j = 0; j < a.cols(); j++)
            r(i,j) = a(i,j) + b(i,j);
        return r;
    }   
    template<class T>
    complex_matrix<T> operator -(const complex_matrix<T> & a, const complex_matrix<T> & b) {
        complex_matrix<T> r(a.rows(),a.cols());
        #pragma omp simd
        for(size_t i = 0; i < a.rows(); i++) 
        for(size_t j = 0; j < a.cols(); j++)
            r(i,j) = a(i,j) - b(i,j);
        return r;
    }   
    template<class T>
    complex_matrix<T> hadamard(const complex_matrix<T> & a, const complex_matrix<T> & b) {
        complex_matrix<T> r(a.rows(),a.cols());
        #pragma omp simd
        for(size_t i = 0; i < a.rows(); i++) 
        for(size_t j = 0; j < a.cols(); j++)
            r(i,j) = a(i,j) * b(i,j);
        return r;
    }   
    template<class T>
    complex_matrix<T> operator *(const complex_matrix<T> & a, const complex_matrix<T> & b) {
        return matmul(a,b);
    }
    template<class T>
    complex_matrix<T> operator /(const complex_matrix<T> & a, const complex_matrix<T> & b) {
        complex_matrix<T> r(a.rows(),a.cols());
        #pragma omp simd
        for(size_t i = 0; i < a.rows(); i++) 
        for(size_t j = 0; j < a.cols(); j++)
            r(i,j) = a(i,j) / b(i,j);
        return r;
    }   
    template<class T>
    complex_matrix<T> operator %(const complex_matrix<T> & a, const complex_matrix<T> & b) {
        complex_matrix<T> r(a.rows(),a.cols());
        #pragma omp simd
        for(size_t i = 0; i < a.rows(); i++) 
        for(size_t j = 0; j < a.cols(); j++)
            r(i,j) = std::fmod(a(i,j),b(i,j));
        return r;
    }   
    template<class T>
    complex_matrix<T> operator +(const complex_matrix<T> & a, const std::complex<T>& b) {
        complex_matrix<T> r(a.rows(),a.cols());
        #pragma omp simd
        for(size_t i = 0; i < a.rows(); i++) 
        for(size_t j = 0; j < a.cols(); j++)
            r(i,j) = a(i,j) + b;
        return r;
    }   
    template<class T>
    complex_matrix<T> operator -(const complex_matrix<T> & a, const std::complex<T>& b) {
        complex_matrix<T> r(a.rows(),a.cols());
        #pragma omp simd
        for(size_t i = 0; i < a.rows(); i++) 
        for(size_t j = 0; j < a.cols(); j++)
            r(i,j) = a(i,j) - b;
        return r;
    }   
    template<class T>
    complex_matrix<T> operator *(const complex_matrix<T> & a, const std::complex<T>& b) {
        complex_matrix<T> r(a.rows(),a.cols());
        #pragma omp simd
        for(size_t i = 0; i < a.rows(); i++) 
        for(size_t j = 0; j < a.cols(); j++)
            r(i,j) = a(i,j) * b;
        return r;
    }   
    template<class T>
    complex_matrix<T> operator / (const complex_matrix<T> & a, const std::complex<T>& b) {
        complex_matrix<T> r(a.rows(),a.cols());
        #pragma omp simd
        for(size_t i = 0; i < a.rows(); i++) 
        for(size_t j = 0; j < a.cols(); j++)
            r(i,j) = a(i,j) / b;
        return r;
    }   
    template<class T>
    complex_matrix<T> operator %(const complex_matrix<T> & a, const std::complex<T>& b) {
        complex_matrix<T> r(a.rows(),a.cols());
        #pragma omp simd
        for(size_t i = 0; i < a.rows(); i++) 
        for(size_t j = 0; j < a.cols(); j++)
            r(i,j) = std::fmod(a(i,j),b);
        return r;
    }   
    template<class T>
    complex_matrix<T> operator +(const std::complex<T> & a, const complex_matrix<T>& b) {
        complex_matrix<T> r(a.rows(),a.cols());
        #pragma omp simd
        for(size_t i = 0; i < a.rows(); i++) 
        for(size_t j = 0; j < a.cols(); j++)        
            r(i,j) = a + b(i,j);
        return r;
    }   
    template<class T>
    complex_matrix<T> operator -(const std::complex<T> & a, const complex_matrix<T>& b) {
        complex_matrix<T> r(a.rows(),a.cols());
        #pragma omp simd
        for(size_t i = 0; i < a.rows(); i++) 
        for(size_t j = 0; j < a.cols(); j++)
            r(i,j) = a - b(i,j);
        return r;
    }   
    template<class T>
    complex_matrix<T> operator *(const std::complex<T> & a, const complex_matrix<T>& b) {
        complex_matrix<T> r(a.rows(),a.cols());
        #pragma omp simd
        for(size_t i = 0; i < a.rows(); i++) 
        for(size_t j = 0; j < a.cols(); j++)
            r(i,j) = a * b(i,j);
        return r;
    }   
    template<class T>
    complex_matrix<T> operator /(const std::complex<T> & a, const complex_matrix<T>& b) {
        complex_matrix<T> r(a.rows(),a.cols());
        #pragma omp simd
        for(size_t i = 0; i < a.rows(); i++) 
        for(size_t j = 0; j < a.cols(); j++)
            r(i,j) = a / b(i,j);
        return r;
    }   
        template<class T>
    complex_matrix<T> pow(const complex_matrix<T> & v, const complex_matrix<T> & x) {
        complex_matrix<T> r(v.rows(),v.cols());
        #pragma omp simd
        for(size_t i = 0; i < v.rows(); i++) 
        for(size_t j = 0; j < v.cols(); j++)
            r(i,j) = std::pow(v(i,j),x(i,j));
        return r;
    }    
    template<class T>
    complex_matrix<T> pow(const complex_matrix<T> & v, const T x) {
        complex_matrix<T> r(v.rows(),v.cols());
        #pragma omp simd
        for(size_t i = 0; i < v.rows(); i++) 
        for(size_t j = 0; j < v.cols(); j++)
            r(i,j) = std::pow(v(i,j),x);
        return r;
    }    
    template<class T>
    complex_matrix<T> pow(const T x, const complex_matrix<T> & v) {
        complex_matrix<T> r(v.rows(),v.cols());
        #pragma omp simd
        for(size_t i = 0; i < v.rows(); i++) 
        for(size_t j = 0; j < v.cols(); j++)
            r(i,j) = std::pow(x,v(i,j));
        return r;
    }    
    template<class T>
    complex_matrix<T> sqrt(const complex_matrix<T> & v) {
        complex_matrix<T> r(v.rows(),v.cols());
        #pragma omp simd
        for(size_t i = 0; i < v.rows(); i++) 
        for(size_t j = 0; j < v.cols(); j++)
            r(i,j) = std::sqrt(v(i,j));
        return r;
    }    
    template<class T>
    complex_matrix<T> abs(const complex_matrix<T> & v) {
        complex_matrix<T> r(v.rows(),v.cols());
        #pragma omp simd
        for(size_t i = 0; i < v.rows(); i++) 
        for(size_t j = 0; j < v.cols(); j++)
            r(i,j) = std::fabs(v(i,j));
        return r;
    }
    template<class T>
    complex_matrix<T> cube(const complex_matrix<T> & v) {
        complex_matrix<T> r(v.rows(),v.cols());
        #pragma omp simd
        for(size_t i = 0; i < v.rows(); i++) 
        for(size_t j = 0; j < v.cols(); j++)
            r(i,j) = v(i,j)*v(i,j)*v(i,j);
        return r;
    }
    template<class T>
    complex_matrix<T> sqr(const complex_matrix<T> & v) {
        complex_matrix<T> r(v.rows(),v.cols());
        #pragma omp simd
        for(size_t i = 0; i < v.rows(); i++) 
        for(size_t j = 0; j < v.cols(); j++)
            r(i,j) = v(i,j)*v(i,j);
        return r;
    }
    template<class T>
    complex_matrix<T> cos(const complex_matrix<T> & v) {
        complex_matrix<T> r(v.rows(),v.cols());
        #pragma omp simd
        for(size_t i = 0; i < v.rows(); i++) 
        for(size_t j = 0; j < v.cols(); j++)
            r(i,j) = std::cos(v(i,j));
        return r;
    }
    template<class T>
    complex_matrix<T> sin(const complex_matrix<T> & v) {
        complex_matrix<T> r(v.rows(),v.cols());
        #pragma omp simd
        for(size_t i = 0; i < v.rows(); i++) 
        for(size_t j = 0; j < v.cols(); j++)
            r(i,j) = std::sin(v(i,j));
        return r;
    }    
    template<class T>
    complex_matrix<T> tan(const complex_matrix<T> & v) {
        complex_matrix<T> r(v.rows(),v.cols());
        #pragma omp simd
        for(size_t i = 0; i < v.rows(); i++) 
        for(size_t j = 0; j < v.cols(); j++)
            r(i,j) = std::tan(v(i,j));
        return r;
    }

    template<class T>
    complex_matrix<T> acos(const complex_matrix<T> & v) {
        complex_matrix<T> r(v.rows(),v.cols());
        #pragma omp simd
        for(size_t i = 0; i < v.rows(); i++) 
        for(size_t j = 0; j < v.cols(); j++)
            r(i,j) = std::acos(v(i,j));
        return r;
    }
    template<class T>
    complex_matrix<T> asin(const complex_matrix<T> & v) {
        complex_matrix<T> r(v.rows(),v.cols());
        #pragma omp simd
        for(size_t i = 0; i < v.rows(); i++) 
        for(size_t j = 0; j < v.cols(); j++)
            r(i,j) = std::asin(v(i,j));
        return r;
    }    
    template<class T>
    complex_matrix<T> atan(const complex_matrix<T> & v) {
        complex_matrix<T> r(v.rows(),v.cols());
        #pragma omp simd
        for(size_t i = 0; i < v.rows(); i++) 
        for(size_t j = 0; j < v.cols(); j++)
            r(i,j) = std::atan(v(i,j));
        return r;
    }        
    template<class T>
    complex_matrix<T> cosh(const complex_matrix<T> & v) {
        complex_matrix<T> r(v.rows(),v.cols());
        #pragma omp simd
        for(size_t i = 0; i < v.rows(); i++) 
        for(size_t j = 0; j < v.cols(); j++)
            r(i,j) = std::cosh(v(i,j));
        return r;
    }
    template<class T>
    complex_matrix<T> sinh(const complex_matrix<T> & v) {
        complex_matrix<T> r(v.rows(),v.cols());
        #pragma omp simd
        for(size_t i = 0; i < v.rows(); i++) 
        for(size_t j = 0; j < v.cols(); j++)
            r(i,j) = std::sinh(v(i,j));
        return r;
    }    
    template<class T>
    complex_matrix<T> tanh(const complex_matrix<T> & v) {
        complex_matrix<T> r(v.rows(),v.cols());
        #pragma omp simd
        for(size_t i = 0; i < v.rows(); i++) 
        for(size_t j = 0; j < v.cols(); j++)
            r(i,j) = std::tanh(v(i,j));
        return r;
    }

    template<class T>
    complex_matrix<T> acosh(const complex_matrix<T> & v) {
        complex_matrix<T> r(v.rows(),v.cols());
        #pragma omp simd
        for(size_t i = 0; i < v.rows(); i++) 
        for(size_t j = 0; j < v.cols(); j++)        
            r(i,j) = std::acosh(v(i,j));
        return r;
    }
    template<class T>
    complex_matrix<T> asinh(const complex_matrix<T> & v) {
        complex_matrix<T> r(v.rows(),v.cols());
        #pragma omp simd
        for(size_t i = 0; i < v.rows(); i++) 
        for(size_t j = 0; j < v.cols(); j++)
            r(i,j) = std::asinh(v(i,j));
        return r;
    }    
    template<class T>
    complex_matrix<T> atanh(const complex_matrix<T> & v) {
        complex_matrix<T> r(v.rows(),v.cols());
        #pragma omp simd
        for(size_t i = 0; i < v.rows(); i++) 
        for(size_t j = 0; j < v.cols(); j++)
            r(i,j) = std::atanh(v(i,j));
        return r;
    }    

    template<class T>
    complex_matrix<T> exp(const complex_matrix<T> & v) {
        complex_matrix<T> r(v.rows(),v.cols());
        #pragma omp simd
        for(size_t i = 0; i < v.rows(); i++) 
        for(size_t j = 0; j < v.cols(); j++)
            r(i,j) = std::exp(v(i,j));
        return r;
    }    
    template<class T>
    complex_matrix<T> log(const complex_matrix<T> & v) {
        complex_matrix<T> r(v.rows(),v.cols());
        #pragma omp simd
        for(size_t i = 0; i < v.rows(); i++) 
        for(size_t j = 0; j < v.cols(); j++)        
            r(i,j) = std::log(v(i,j));
        return r;
    }    
    template<class T>
    complex_matrix<T> log10(const complex_matrix<T> & v) {
        complex_matrix<T> r(v.rows(),v.cols());
        #pragma omp simd
        for(size_t i = 0; i < v.rows(); i++) 
        for(size_t j = 0; j < v.cols(); j++)
            r(i,j) = std::log10(v(i,j));
        return r;
    }    
    template<class T>
    std::complex<T> sum(const complex_matrix<T> & v)
    {
        std::complex<T> r(0,0);
        for(size_t i = 0; i < v.rows(); i++) 
        for(size_t j = 0; j < v.cols(); j++)
            r += v(i,j);
        return r;
    }
    template<class T>
    std::complex<T> product(const complex_matrix<T> & v)
    {
        std::complex<T> r(1,0);
        for(size_t i = 0; i < v.rows(); i++) 
        for(size_t j = 0; j < v.cols(); j++)
            r += v(i,j);
        return r;
    }
        template<class T>
    sample_matrix<T> real(const complex_matrix<T> & v) {
        sample_matrix<T> r(v.rows(),v.cols());
        for(size_t i = 0; i < v.rows(); i++) 
        for(size_t j = 0; i < v.cols(); j++) 
            r(i,j) = v(i,j).real();
        return r;
    }
    template<class T>
    sample_matrix<T> imag(const complex_matrix<T> & v) {
        sample_matrix<T> r(v.rows(),v.cols());
        for(size_t i = 0; i < v.rows(); i++) 
        for(size_t j = 0; i < v.cols(); j++) 
            r(i,j) = v(i,j).imag();
        return r;
    }
    template<class T>
    sample_matrix<T> arg(const complex_matrix<T> & v) {
        sample_matrix<T> r(v.rows(),v.cols());
        for(size_t i = 0; i < v.rows(); i++) 
        for(size_t j = 0; i < v.cols(); j++) 
            r(i,j) = std::arg(v(i,j));
        return r;
    }
    template<class T>
    sample_matrix<T> norm(const complex_matrix<T> & v) {
        sample_matrix<T> r(v.rows(),v.cols());
        for(size_t i = 0; i < v.rows(); i++) 
        for(size_t j = 0; i < v.cols(); j++) 
            r(i,j) = std::norm(v(i,j)());
        return r;
    }
    template<class T>
    complex_matrix<T> conj(const complex_matrix<T> & v) {
        complex_matrix<T> r(v.rows(),v.cols());
        for(size_t i = 0; i < v.rows(); i++) 
        for(size_t j = 0; i < v.cols(); j++) 
            r(i,j) = std::conj(v(i,j));
        return r;
    }
    template<class T>
    complex_matrix<T> proj(const complex_matrix<T> & v) {
        complex_matrix<T> r(v.rows(),v.cols());
        for(size_t i = 0; i < v.rows(); i++) 
        for(size_t j = 0; i < v.cols(); j++) 
            r(i,j) = std::proj(v(i,j));
        return r;
    }
    template<class T>
    complex_matrix<T> polar(const sample_matrix<T> & a, const sample_matrix<T> & b) {
        complex_matrix<T> r(a.rows(),a.cols());
        for(size_t i = 0; i < a.rows(); i++) 
        for(size_t j = 0; i < a.cols(); j++) 
            r(i,j) = std::polar(a[i],b[i]);
        return r;
    }
}