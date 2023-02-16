#pragma once

namespace AudioDSP
{
    template<typename T>
    struct sample_matrix_view
    {
        using base = std::vector<T,Allocator::aligned_allocator<T,64>>;
        base * matrix;
        size_t row;

        sample_matrix_view(base * p, size_t r) : matrix(p),row(r) {}

        T& operator[](size_t j) { return (*matrix)[row + j]; }

        #ifdef SWIG
        %extend {
            T __getitem__(size_t j) { return (*$self)[j]; }
            void __setitem(size_t j, const T v) { (*$self)[j] = v; }
        }
        #endif
    };

    template<typename T>
    struct sample_matrix : public std::vector<T,Allocator::aligned_allocator<T,64>>
    {
        using base = std::vector<T,Allocator::aligned_allocator<T,64>>;
        size_t M,N;

        sample_matrix() = default;
        sample_matrix(size_t i, size_t j) : base(i) {
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

        T& operator()(size_t i, size_t j) { return (*this)[i*N+j]; }
        T  operator()(size_t i, size_t j) const { return (*this)[i*N+j]; }

        #ifdef SWIG
        %extend {
            sample_matrix_view<T> __getitem__(size_t i) {
                sample_matrix_view<T> view($self,i*$self->N);
                return view;
            }
            sample_matrix<T> __add__(const sample_matrix<T> & v)
            {
                assert($self->size() == v.size());
                sample_matrix<T> r($self->rows(),$self->cols());
                #pragma omp simd
                for(size_t i = 0; i < v.size(); i++) r[i] = (*$self)[i] + v[i];
                return r;
            }
        
            sample_matrix<T> __sub__(const sample_matrix<T> & v)
            {
                assert($self->size() == v.size());
                sample_matrix<T> r($self->rows(),$self->cols());
                #pragma omp simd
                for(size_t i = 0; i < v.size(); i++) r[i] = (*$self)[i] - v[i];
                return r;
            }
            sample_matrix<T> __mul__(const sample_matrix<T> & v)
            {
                assert($self->cols() == v.rows());
                sample_matrix<T> r($self->rows(),v.cols());
                #pragma omp parallel for simd                
                for(size_t i = 0; i < $self->rows(); i++)
                {                    
                    for(size_t j = 0; j < v.cols(); j++)
                    {
                        T sum = 0;
                        for(size_t k = 0; k < $self->cols(); k++)
                        {
                            sum += (*$self)(i,k) * v(i,j);
                        }
                        r(i,j) = sum;
                    }
                }
                return r;
            }
            sample_matrix<T> __div__(const sample_matrix<T> & v)
            {
                assert($self->size() == v.size());
                sample_matrix<T> r($self->rows(),$self->cols());
                #pragma omp simd
                for(size_t i = 0; i < v.size(); i++) r[i] = (*$self)[i] / v[i];
                return r;
            }
            sample_matrix<T> __pow__(const sample_matrix<T> & v)
            {
                assert($self->size() == v.size());
                sample_matrix<T> r($self->rows(),$self->cols());
                #pragma omp simd
                for(size_t i = 0; i < v.size(); i++) r[i] = std::pow((*$self)[i],v[i]);
                return r;
            }
            sample_matrix<T> __add__(const T& x)
            {
                assert($self->size() == $self->size());
                sample_matrix<T> r($self->rows(),$self->cols());
                #pragma omp simd
                for(size_t i = 0; i < $self->size(); i++) r[i] = (*$self)[i] + x;
                return r;
            }
            sample_matrix<T> __sub__(const T& x)
            {
                assert($self->size() == $self->size());
                sample_matrix<T> r($self->rows(),$self->cols());
                #pragma omp simd
                for(size_t i = 0; i < $self->size(); i++) r[i] = (*$self)[i] - x;
                return r;
            }
            sample_matrix<T> __mul__(const T& x)
            {
                assert($self->size() == $self->size());
                sample_matrix<T> r($self->rows(),$self->cols());
                #pragma omp simd
                for(size_t i = 0; i < $self->size(); i++) r[i] = (*$self)[i] * x;
                return r;
            }
            sample_matrix<T> __div__(const T& x)
            {
                assert($self->size() == $self->size());
                sample_matrix<T> r($self->rows(),$self->cols());
                #pragma omp simd
                for(size_t i = 0; i < $self->size(); i++) r[i] = (*$self)[i] / x;
                return r;
            }
            sample_matrix<T> __pow__(const T& x)
            {
                assert($self->size() == $self->size());
                sample_matrix<T> r($self->rows(),$self->cols());
                #pragma omp simd
                for(size_t i = 0; i < $self->size(); i++) r[i] = std::pow((*$self)[i],x);
                return r;
            }
            sample_matrix<T> __neg__()
            {
                assert($self->size() == $self->size());
                sample_matrix<T> r($self->rows(),$self->cols());
                #pragma omp simd
                for(size_t i = 0; i < $self->size(); i++) r[i] = -(*$self)[i];
                return r;
            }
        }
        #endif

        size_t rows() const { return M; }
        size_t cols() const { return N; }
        void   resize(size_t i, size_t j) {
            M = i;
            N = j;
            resize(i*j);
        }
        void fill(const T v) {
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
        void copy(size_t m, size_t n, T * p)
        {
            resize(m,n);
            std::memcpy(data(),p,m*n*sizeof(T));
        }        
    };

    template<class T>
    sample_matrix<T> abs(const sample_matrix<T> & v) {
        sample_matrix<T> r(v.rows(),v.cols());
        #pragma omp simd
        for(size_t i = 0; i < v.rows(); i++) 
        for(size_t j = 0; j < v.cols(); j++)
            r(i,j) = std::fabs(v(i,j));
        return r;
    }
    template<class T>
    sample_matrix<T> cube(const sample_matrix<T> & v) {
        sample_matrix<T> r(v.rows(),v.cols());
        #pragma omp simd
        for(size_t i = 0; i < v.rows(); i++) 
        for(size_t j = 0; j < v.cols(); j++)
            r(i,j) = v(i,j)*v(i,j)*v(i,j);
        return r;
    }
    template<class T>
    sample_matrix<T> sqr(const sample_matrix<T> & v) {
        sample_matrix<T> r(v.rows(),v.cols());
        #pragma omp simd
        for(size_t i = 0; i < v.rows(); i++) 
        for(size_t j = 0; j < v.cols(); j++)
        r(i,j) = v(i,j)*v(i,j);
        return r;
    }
    template<class T>
    sample_matrix<T> cos(const sample_matrix<T> & v) {
        sample_matrix<T> r(v.rows(),v.cols());
        #pragma omp simd
        for(size_t i = 0; i < v.rows(); i++) 
        for(size_t j = 0; j < v.cols(); j++)
            r(i,j) = std::cos(v(i,j));
        return r;
    }
    template<class T>
    sample_matrix<T> sin(const sample_matrix<T> & v) {
        sample_matrix<T> r(v.rows(),v.cols());
        #pragma omp simd
        for(size_t i = 0; i < v.rows(); i++) 
        for(size_t j = 0; j < v.cols(); j++)
            r(i,j) = std::sin(v(i,j));
        return r;
    }    
    template<class T>
    sample_matrix<T> tan(const sample_matrix<T> & v) {
        sample_matrix<T> r(v.rows(),v.cols());
        #pragma omp simd
        for(size_t i = 0; i < v.rows(); i++) 
        for(size_t j = 0; j < v.cols(); j++)
            r(i,j) = std::tan(v(i,j));
        return r;
    }

    template<class T>
    sample_matrix<T> acos(const sample_matrix<T> & v) {
        sample_matrix<T> r(v.rows(),v.cols());
        #pragma omp simd
        for(size_t i = 0; i < v.rows(); i++) 
        for(size_t j = 0; j < v.cols(); j++)
            r(i,j) = std::acos(v(i,j));
        return r;
    }
    template<class T>
    sample_matrix<T> asin(const sample_matrix<T> & v) {
        sample_matrix<T> r(v.rows(),v.cols());
        #pragma omp simd
        for(size_t i = 0; i < v.rows(); i++) 
        for(size_t j = 0; j < v.cols(); j++)
            r(i,j) = std::asin(v(i,j));
        return r;
    }    
    template<class T>
    sample_matrix<T> atan(const sample_matrix<T> & v) {
        sample_matrix<T> r(v.rows(),v.cols());
        #pragma omp simd
        for(size_t i = 0; i < v.rows(); i++) 
        for(size_t j = 0; j < v.cols(); j++)
            r(i,j) = std::atan(v(i,j));
        return r;
    }    
    template<class T>
    sample_matrix<T> atan2(const sample_matrix<T> & v, const T value) {
        sample_matrix<T> r(v.rows(),v.cols());
        #pragma omp simd
        for(size_t i = 0; i < v.rows(); i++) 
        for(size_t j = 0; j < v.cols(); j++)
            r(i,j) = std::atan2(v(i,j), value);
        return r;
    }    
    template<class T>
    sample_matrix<T> atan2(const sample_matrix<T> & v, const sample_matrix<T> value) {
        sample_matrix<T> r(v.rows(),v.cols());
        #pragma omp simd
        for(size_t i = 0; i < v.rows(); i++) 
        for(size_t j = 0; j < v.cols(); j++)        
            r(i,j) = std::atan2(v(i,j), value(i,j));
        return r;
    }    
    template<class T>
    sample_matrix<T> cosh(const sample_matrix<T> & v) {
        sample_matrix<T> r(v.rows(),v.cols());
        #pragma omp simd
        for(size_t i = 0; i < v.rows(); i++) 
        for(size_t j = 0; j < v.cols(); j++)
            r(i,j) = std::cosh(v(i,j));
        return r;
    }
    template<class T>
    sample_matrix<T> sinh(const sample_matrix<T> & v) {
        sample_matrix<T> r(v.rows(),v.cols());
        #pragma omp simd
        for(size_t i = 0; i < v.rows(); i++) 
        for(size_t j = 0; j < v.cols(); j++)
            r(i,j) = std::sinh(v(i,j));
        return r;
    }    
    template<class T>
    sample_matrix<T> tanh(const sample_matrix<T> & v) {
        sample_matrix<T> r(v.rows(),v.cols());
        #pragma omp simd
        for(size_t i = 0; i < v.rows(); i++) 
        for(size_t j = 0; j < v.cols(); j++)
            r(i,j) = std::tanh(v(i,j));
        return r;
    }

    template<class T>
    sample_matrix<T> acosh(const sample_matrix<T> & v) {
        sample_matrix<T> r(v.rows(),v.cols());
        #pragma omp simd
        for(size_t i = 0; i < v.rows(); i++) 
        for(size_t j = 0; j < v.cols(); j++)        
            r(i,j) = std::acosh(v(i,j));
        return r;
    }
    template<class T>
    sample_matrix<T> asinh(const sample_matrix<T> & v) {
        sample_matrix<T> r(v.rows(),v.cols());
        #pragma omp simd
        for(size_t i = 0; i < v.rows(); i++) 
        for(size_t j = 0; j < v.cols(); j++)
            r(i,j) = std::asinh(v(i,j));
        return r;
    }    
    template<class T>
    sample_matrix<T> atanh(const sample_matrix<T> & v) {
        sample_matrix<T> r(v.rows(),v.cols());
        #pragma omp simd
        for(size_t i = 0; i < v.rows(); i++) 
        for(size_t j = 0; j < v.cols(); j++)
            r(i,j) = std::atanh(v(i,j));
        return r;
    }    

    template<class T>
    sample_matrix<T> exp(const sample_matrix<T> & v) {
        sample_matrix<T> r(v.rows(),v.cols());
        #pragma omp simd
        for(size_t i = 0; i < v.rows(); i++) 
        for(size_t j = 0; j < v.cols(); j++)
            r(i,j) = std::exp(v(i,j));
        return r;
    }    
    template<class T>
    sample_matrix<T> log(const sample_matrix<T> & v) {
        sample_matrix<T> r(v.rows(),v.cols());
        #pragma omp simd
        for(size_t i = 0; i < v.rows(); i++) 
        for(size_t j = 0; j < v.cols(); j++)        
            r(i,j) = std::log(v(i,j));
        return r;
    }    
    template<class T>
    sample_matrix<T> log10(const sample_matrix<T> & v) {
        sample_matrix<T> r(v.rows(),v.cols());
        #pragma omp simd
        for(size_t i = 0; i < v.rows(); i++) 
        for(size_t j = 0; j < v.cols(); j++)
            r(i,j) = std::log10(v(i,j));
        return r;
    }    
    template<class T>
    sample_matrix<T> exp2(const sample_matrix<T> & v) {
        sample_matrix<T> r(v.rows(),v.cols());
        #pragma omp simd
        for(size_t i = 0; i < v.rows(); i++) 
        for(size_t j = 0; j < v.cols(); j++)
            r(i,j) = std::exp2(v(i,j));
        return r;
    }    
    template<class T>
    sample_matrix<T> expm1(const sample_matrix<T> & v) {
        sample_matrix<T> r(v.rows(),v.cols());
        #pragma omp simd
        for(size_t i = 0; i < v.rows(); i++) 
        for(size_t j = 0; j < v.cols(); j++)
            r(i,j) = std::expm1(v(i,j));
        return r;
    }    
    template<class T>
    sample_matrix<T> ilogb(const sample_matrix<T> & v) {
        sample_matrix<T> r(v.rows(),v.cols());
        #pragma omp simd
        for(size_t i = 0; i < v.rows(); i++) 
        for(size_t j = 0; j < v.cols(); j++)
            r(i,j) = std::ilogb(v(i,j));
        return r;
    }    
    template<class T>
    sample_matrix<T> log2(const sample_matrix<T> & v) {
        sample_matrix<T> r(v.rows(),v.cols());
        #pragma omp simd
        for(size_t i = 0; i < v.rows(); i++) 
        for(size_t j = 0; j < v.cols(); j++)        
            r(i,j) = std::log2(v(i,j));
        return r;
    }    
    template<class T>
    sample_matrix<T> log1p(const sample_matrix<T> & v) {
        sample_matrix<T> r(v.rows(),v.cols());
        #pragma omp simd
        for(size_t i = 0; i < v.rows(); i++) 
        for(size_t j = 0; j < v.cols(); j++)
            r(i,j) = std::log1p(v(i,j));
        return r;
    }    
    template<class T>
    sample_matrix<T> logb(const sample_matrix<T> & v) {
        sample_matrix<T> r(v.rows(),v.cols());
        #pragma omp simd
        for(size_t i = 0; i < v.rows(); i++) 
        for(size_t j = 0; j < v.cols(); j++)
            r(i,j) = std::logb(v(i,j));
        return r;
    }    
    template<class T>
    sample_matrix<T> scalbn(const sample_matrix<T> & v, const sample_matrix<int> & x) {
        sample_matrix<T> r(v.rows(),v.cols());
        #pragma omp simd
        for(size_t i = 0; i < v.rows(); i++) 
        for(size_t j = 0; j < v.cols(); j++)
            r(i,j) = std::scalbn(v(i,j),x(i,j));
        return r;
    }    
    template<class T>
    sample_matrix<T> scalbn(const sample_matrix<T> & v, const int x) {
        sample_matrix<T> r(v.rows(),v.cols());
        #pragma omp simd
        for(size_t i = 0; i < v.rows(); i++) 
        for(size_t j = 0; j < v.cols(); j++)
            r(i,j) = std::scalbn(v(i,j),x);
        return r;
    }    
    template<class T>
    sample_matrix<T> scalbln(const sample_matrix<T> & v, const sample_matrix<long int> & x) {
        sample_matrix<T> r(v.rows(),v.cols());
        #pragma omp simd
        for(size_t i = 0; i < v.rows(); i++) 
        for(size_t j = 0; j < v.cols(); j++)
            r(i,j) = std::scalbln(v(i,j),x(i,j));
        return r;
    }    
    template<class T>
    sample_matrix<T> scalbln(const sample_matrix<T> & v, const long int x) {
        sample_matrix<T> r(v.rows(),v.cols());
        #pragma omp simd
        for(size_t i = 0; i < v.rows(); i++) 
        for(size_t j = 0; j < v.cols(); j++)
            r(i,j) = std::scalbln(v(i,j),x);
        return r;
    }    
    template<class T>
    sample_matrix<T> pow(const sample_matrix<T> & v, const sample_matrix<T> & x) {
        sample_matrix<T> r(v.rows(),v.cols());
        #pragma omp simd
        for(size_t i = 0; i < v.rows(); i++) 
        for(size_t j = 0; j < v.cols(); j++)
            r(i,j) = std::pow(v(i,j),x(i,j));
        return r;
    }    
    template<class T>
    sample_matrix<T> pow(const sample_matrix<T> & v, const T x) {
        sample_matrix<T> r(v.rows(),v.cols());
        #pragma omp simd
        for(size_t i = 0; i < v.rows(); i++) 
        for(size_t j = 0; j < v.cols(); j++)
            r(i,j) = std::pow(v(i,j),x);
        return r;
    }    
    template<class T>
    sample_matrix<T> pow(const T x, const sample_matrix<T> & v) {
        sample_matrix<T> r(v.rows(),v.cols());
        #pragma omp simd
        for(size_t i = 0; i < v.rows(); i++) 
        for(size_t j = 0; j < v.cols(); j++)
            r(i,j) = std::pow(x,v(i,j));
        return r;
    }    
    template<class T>
    sample_matrix<T> sqrt(const sample_matrix<T> & v) {
        sample_matrix<T> r(v.rows(),v.cols());
        #pragma omp simd
        for(size_t i = 0; i < v.rows(); i++) 
        for(size_t j = 0; j < v.cols(); j++)
            r(i,j) = std::sqrt(v(i,j));
        return r;
    }    
    template<class T>
    sample_matrix<T> cbrt(const sample_matrix<T> & v) {
        sample_matrix<T> r(v.rows(),v.cols());
        #pragma omp simd
        for(size_t i = 0; i < v.rows(); i++) 
        for(size_t j = 0; j < v.cols(); j++)
            r(i,j) = std::cbrt(v(i,j));
        return r;
    }    
    template<class T>
    sample_matrix<T> hypot(const sample_matrix<T> & v, const sample_matrix<T> & x) {
        sample_matrix<T> r(v.rows(),v.cols());
        #pragma omp simd
        for(size_t i = 0; i < v.rows(); i++) 
        for(size_t j = 0; j < v.cols(); j++)
            r(i,j) = std::hypot(v(i,j),x(i,j));
        return r;
    }    
    template<class T>
    sample_matrix<T> hypot(const sample_matrix<T> & v, const T x) {
        sample_matrix<T> r(v.rows(),v.cols());
        #pragma omp simd
        for(size_t i = 0; i < v.rows(); i++) 
        for(size_t j = 0; j < v.cols(); j++)
            r(i,j) = std::hypot(v(i,j),x);
        return r;
    }    
    template<class T>
    sample_matrix<T> hypot(const T x, const sample_matrix<T> & v) {
        sample_matrix<T> r(v.rows(),v.cols());
        #pragma omp simd
        for(size_t i = 0; i < v.rows(); i++) 
        for(size_t j = 0; j < v.cols(); j++)
            r(i,j) = std::hypot(x,v(i,j));
        return r;
    }    
    template<class T>
    sample_matrix<T> erf(const sample_matrix<T> & v) {
        sample_matrix<T> r(v.rows(),v.cols());
        #pragma omp simd
        for(size_t i = 0; i < v.rows(); i++) 
        for(size_t j = 0; j < v.cols(); j++)
            r(i,j) = std::erf(v(i,j));
        return r;
    }    
    template<class T>
    sample_matrix<T> erfc(const sample_matrix<T> & v) {
        sample_matrix<T> r(v.rows(),v.cols());
        #pragma omp simd
        for(size_t i = 0; i < v.rows(); i++) 
        for(size_t j = 0; j < v.cols(); j++)
            r(i,j) = std::erfc(v(i,j));
        return r;
    }    
    template<class T>
    sample_matrix<T> tgamma(const sample_matrix<T> & v) {
        sample_matrix<T> r(v.rows(),v.cols());
        #pragma omp simd
        for(size_t i = 0; i < v.rows(); i++) 
        for(size_t j = 0; j < v.cols(); j++)
            r(i,j) = std::tgamma(v(i,j));
        return r;
    }    
    template<class T>
    sample_matrix<T> lgamma(const sample_matrix<T> & v) {
        sample_matrix<T> r(v.rows(),v.cols());
        #pragma omp simd
        for(size_t i = 0; i < v.rows(); i++) 
        for(size_t j = 0; j < v.cols(); j++)
            r(i,j) = std::lgamma(v(i,j));
        return r;
    }    
    template<class T>
    sample_matrix<T> ceil(const sample_matrix<T> & v) {
        sample_matrix<T> r(v.rows(),v.cols());
        #pragma omp simd
        for(size_t i = 0; i < v.rows(); i++) 
        for(size_t j = 0; j < v.cols(); j++)
            r(i,j) = std::ceil(v(i,j));
        return r;
    }    
    template<class T>
    sample_matrix<T> floor(const sample_matrix<T> & v) {
        sample_matrix<T> r(v.rows(),v.cols());
        #pragma omp simd
        for(size_t i = 0; i < v.rows(); i++) 
        for(size_t j = 0; j < v.cols(); j++)
            r(i,j) = std::floor(v(i,j));
        return r;
    }    
    template<class T>
    sample_matrix<T> fmod(const sample_matrix<T> & v, const sample_matrix<T> & x) {
        sample_matrix<T> r(v.rows(),v.cols());
        #pragma omp simd
        for(size_t i = 0; i < v.rows(); i++) 
        for(size_t j = 0; j < v.cols(); j++)
            r(i,j) = std::fmod(v(i,j),x(i,j));
        return r;
    }    
    template<class T>
    sample_matrix<T> fmod(const sample_matrix<T> & v, const T x) {
        sample_matrix<T> r(v.rows(),v.cols());
        #pragma omp simd
        for(size_t i = 0; i < v.rows(); i++) 
        for(size_t j = 0; j < v.cols(); j++)
            r(i,j) = std::fmod(v(i,j),x);
        return r;
    }    
    template<class T>
    sample_matrix<T> fmod(const T x, const sample_matrix<T> & v) {
        sample_matrix<T> r(v.rows(),v.cols());
        #pragma omp simd
        for(size_t i = 0; i < v.rows(); i++) 
        for(size_t j = 0; j < v.cols(); j++)
            r(i,j) = std::fmod(x,v(i,j));
        return r;
    }    
    template<class T>
    sample_matrix<T> trunc(const sample_matrix<T> & v) {
        sample_matrix<T> r(v.rows(),v.cols());
        #pragma omp simd
        for(size_t i = 0; i < v.rows(); i++) 
        for(size_t j = 0; j < v.cols(); j++)
            r(i,j) = std::trunc(v(i,j));
        return r;
    }    
    template<class T>
    sample_matrix<T> round(const sample_matrix<T> & v) {
        sample_matrix<T> r(v.rows(),v.cols());
        #pragma omp simd
        for(size_t i = 0; i < v.rows(); i++) 
        for(size_t j = 0; j < v.cols(); j++)
            r(i,j) = std::round(v(i,j));
        return r;
    }    
    template<class T>
    sample_matrix<long int> lround(const sample_matrix<T> & v) {
        sample_matrix<long int> r(v.rows(),v.cols());
        #pragma omp simd
        for(size_t i = 0; i < v.rows(); i++) 
        for(size_t j = 0; j < v.cols(); j++)
            r(i,j) = std::lround(v(i,j));
        return r;
    }    
    template<class T>
    sample_matrix<long long int> llround(const sample_matrix<T> & v) {
        sample_matrix<long long int> r(v.rows(),v.cols());
        #pragma omp simd
        for(size_t i = 0; i < v.rows(); i++) 
        for(size_t j = 0; j < v.cols(); j++)
            r(i,j) = std::llround(v(i,j));
        return r;
    }    
    template<class T>
    sample_matrix<T> nearbyint(const sample_matrix<T> & v) {
        sample_matrix<long long int> r(v.rows(),v.cols());
        #pragma omp simd
        for(size_t i = 0; i < v.rows(); i++) 
        for(size_t j = 0; j < v.cols(); j++)
            r(i,j) = std::nearbyint(v(i,j));
        return r;
    }    
    template<class T>
    sample_matrix<T> remainder(const sample_matrix<T> & a, const sample_matrix<T> & b) {
        sample_matrix<T> r(a.rows(),a.cols());
        #pragma omp simd
        for(size_t i = 0; i < a.rows(); i++) 
        for(size_t j = 0; j < a.cols(); j++)
            r(i,j) = std::remainder(a(i,j),b(i,j));
        return r;
    }    
    template<class T>
    sample_matrix<T> copysign(const sample_matrix<T> & a, const sample_matrix<T> & b) {
        sample_matrix<T> r(a.rows(),a.cols());
        #pragma omp simd
        for(size_t i = 0; i < a.rows(); i++) 
        for(size_t j = 0; j < a.cols(); j++)
            r(i,j) = std::copysign(a(i,j),b(i,j));
        return r;
    }    
    template<class T>
    sample_matrix<T> fdim(const sample_matrix<T> & a, const sample_matrix<T> & b) {
        sample_matrix<T> r(a.rows(),a.cols());
        #pragma omp simd
        for(size_t i = 0; i < a.rows(); i++) 
        for(size_t j = 0; j < a.cols(); j++)
            r(i,j) = std::fdim(a(i,j),b(i,j));
        return r;
    }    
    #undef fmax
    template<class T>
    sample_matrix<T> fmax(const sample_matrix<T> & a, const sample_matrix<T> & b) {
        sample_matrix<T> r(a.rows(),a.cols());
        #pragma omp simd
        for(size_t i = 0; i < a.rows(); i++) 
        for(size_t j = 0; j < a.cols(); j++)
            r(i,j) = std::fmax(a(i,j),b(i,j));
        return r;
    }    
    template<class T>
    sample_matrix<T> fmin(const sample_matrix<T> & a, const sample_matrix<T> & b) {
        sample_matrix<T> r(a.rows(),a.cols());
        #pragma omp simd
        for(size_t i = 0; i < a.rows(); i++) 
        for(size_t j = 0; j < a.cols(); j++)
            r(i,j) = std::fmin(a(i,j),b(i,j));
        return r;
    }    
    template<class T>
    sample_matrix<T> fma(const sample_matrix<T> & a, const sample_matrix<T> & b, const sample_matrix<T> & c) {
        sample_matrix<T> r(a.rows(),a.cols());
        #pragma omp simd
        for(size_t i = 0; i < a.rows(); i++) 
        for(size_t j = 0; j < a.cols(); j++)
            r(i,j) = std::fma(a(i,j),b(i,j),c(i,j));
        return r;
    }    
    template<class T>
    sample_matrix<T> fabs(const sample_matrix<T> & a) {
        sample_matrix<T> r(a.rows(),a.cols());
        #pragma omp simd
        for(size_t i = 0; i < a.rows(); i++) 
        for(size_t j = 0; j < a.cols(); j++)
            r(i,j) = std::fabs(a(i,j));
        return r;
    }    

    template<typename T>
    sample_matrix<T> matmul(const sample_matrix<T> & a, const sample_matrix<T> & b)
    {
        assert(a.cols() == b.rows());
        sample_matrix<T> r(a.rows(),b.cols());
        #pragma omp parallel for simd                
        for(size_t i = 0; i < a.rows(); i++)
        {                    
            for(size_t j = 0; j < b.cols(); j++)
            {
                T sum = 0;
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
    sample_matrix<T> operator +(const sample_matrix<T> & a, const sample_matrix<T> & b) {
        sample_matrix<T> r(a.rows(),a.cols());
        #pragma omp simd
        for(size_t i = 0; i < a.rows(); i++) 
        for(size_t j = 0; j < a.cols(); j++)
            r(i,j) = a(i,j) + b(i,j);
        return r;
    }   
    template<class T>
    sample_matrix<T> operator -(const sample_matrix<T> & a, const sample_matrix<T> & b) {
        sample_matrix<T> r(a.rows(),a.cols());
        #pragma omp simd
        for(size_t i = 0; i < a.rows(); i++) 
        for(size_t j = 0; j < a.cols(); j++)
            r(i,j) = a(i,j) - b(i,j);
        return r;
    }   
    template<class T>
    sample_matrix<T> hadamard(const sample_matrix<T> & a, const sample_matrix<T> & b) {
        sample_matrix<T> r(a.rows(),a.cols());        
        #pragma omp simd
        for(size_t i = 0; i < a.rows(); i++) 
        for(size_t j = 0; j < a.cols(); j++)
            r(i,j) = a(i,j) * b(i,j);
        return r;
    }   
    template<class T>
    sample_matrix<T> operator *(const sample_matrix<T> & a, const sample_matrix<T> & b) {
        return matmul(a,b);
    }   
    template<class T>
    sample_matrix<T> operator /(const sample_matrix<T> & a, const sample_matrix<T> & b) {
        sample_matrix<T> r(a.rows(),a.cols());
        #pragma omp simd
        for(size_t i = 0; i < a.rows(); i++) 
        for(size_t j = 0; j < a.cols(); j++)
            r(i,j) = a(i,j) / b(i,j);
        return r;
    }   
    template<class T>
    sample_matrix<T> operator %(const sample_matrix<T> & a, const sample_matrix<T> & b) {
        sample_matrix<T> r(a.rows(),a.cols());
        #pragma omp simd
        for(size_t i = 0; i < a.rows(); i++) 
        for(size_t j = 0; j < a.cols(); j++)
            r(i,j) = std::fmod(a(i,j),b(i,j));
        return r;
    }   
    template<class T>
    sample_matrix<T> operator +(const sample_matrix<T> & a, const T& b) {
        sample_matrix<T> r(a.rows(),a.cols());
        #pragma omp simd
        for(size_t i = 0; i < a.rows(); i++) 
        for(size_t j = 0; j < a.cols(); j++)
            r(i,j) = a(i,j) + b;
        return r;
    }   
    template<class T>
    sample_matrix<T> operator -(const sample_matrix<T> & a, const T& b) {
        sample_matrix<T> r(a.rows(),a.cols());
        #pragma omp simd
        for(size_t i = 0; i < a.rows(); i++) 
        for(size_t j = 0; j < a.cols(); j++)
            r(i,j) = a(i,j) - b;
        return r;
    }   
    template<class T>
    sample_matrix<T> operator *(const sample_matrix<T> & a, const T& b) {
        sample_matrix<T> r(a.rows(),a.cols());
        #pragma omp simd
        for(size_t i = 0; i < a.rows(); i++) 
        for(size_t j = 0; j < a.cols(); j++)
            r(i,j) = a(i,j) * b;
        return r;
    }   
    template<class T>
    sample_matrix<T> operator / (const sample_matrix<T> & a, const T& b) {
        sample_matrix<T> r(a.rows(),a.cols());
        #pragma omp simd
        for(size_t i = 0; i < a.rows(); i++) 
        for(size_t j = 0; j < a.cols(); j++)
            r(i,j) = a(i,j) / b;
        return r;
    }   
    template<class T>
    sample_matrix<T> operator %(const sample_matrix<T> & a, const T& b) {
        sample_matrix<T> r(a.rows(),a.cols());
        #pragma omp simd
        for(size_t i = 0; i < a.rows(); i++) 
        for(size_t j = 0; j < a.cols(); j++)
            r(i,j) = std::fmod(a(i,j),b);
        return r;
    }   
    template<class T>
    sample_matrix<T> operator +(const T & a, const sample_matrix<T>& b) {
        sample_matrix<T> r(a.rows(),a.cols());
        #pragma omp simd
        for(size_t i = 0; i < a.rows(); i++) 
        for(size_t j = 0; j < a.cols(); j++)
            r(i,j) = a + b(i,j);
        return r;
    }   
    template<class T>
    sample_matrix<T> operator -(const T & a, const sample_matrix<T>& b) {
        sample_matrix<T> r(a.rows(),a.cols());
        #pragma omp simd
        for(size_t i = 0; i < a.rows(); i++) 
        for(size_t j = 0; j < a.cols(); j++)
            r(i,j) = a - b(i,j);
        return r;
    }   
    template<class T>
    sample_matrix<T> operator *(const T & a, const sample_matrix<T>& b) {
        sample_matrix<T> r(a.rows(),a.cols());
        #pragma omp simd
        for(size_t i = 0; i < a.rows(); i++) 
        for(size_t j = 0; j < a.cols(); j++)
            r(i,j) = a * b(i,j);
        return r;
    }   
    template<class T>
    sample_matrix<T> operator /(const T & a, const sample_matrix<T>& b) {
        sample_matrix<T> r(a.rows(),a.cols());
        #pragma omp simd
        for(size_t i = 0; i < a.rows(); i++) 
        for(size_t j = 0; j < a.cols(); j++)
            r(i,j) = a / b(i,j);
        return r;
    }   
    template<class T>
    T sum(const sample_matrix<T> & v)
    {
        T r = 0;
        for(size_t i = 0; i < v.rows(); i++) 
        for(size_t j = 0; j < v.cols(); j++)
            r += v(i,j);
        return r;
    }
    template<class T>
    T product(const sample_matrix<T> & v)
    {
        T r = 1;
        for(size_t i = 0; i < v.rows(); i++) 
        for(size_t j = 0; j < v.cols(); j++)
            r *= v(i,j);
        return r;
    }
}