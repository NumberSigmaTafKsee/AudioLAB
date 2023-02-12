#pragma once

namespace Casino::MKL
{
    template<typename T> struct ComplexMatrix;
    template<typename T>
    struct ComplexMatrixView 
    {
        ComplexMatrix<T> * matrix;
        size_t row;
        ComplexMatrixView(ComplexMatrix<T> * m, size_t r) {
            matrix = m;
            row = r;
        }

        std::complex<T>& operator[](size_t i);
        std::complex<T>  __getitem__(size_t i);
        void __setitem__(size_t i, std::complex<T> x);
    };

    template<typename T>
    struct ComplexMatrix : public vector_base<std::complex<T>>
    {                
        using vecbase = vector_base<std::complex<T>>;
        using vecbase::size;
        using vecbase::resize;
        using vecbase::data;
        using vecbase::push_back;
        using vecbase::pop_back;
        using vecbase::front;
        using vecbase::back;
        using vecbase::at;
        using vecbase::operator [];
        using vecbase::operator =;

        size_t M;
        size_t N;

        ComplexMatrix() = default;
        ComplexMatrix(size_t i,size_t j) : vecbase(i*j),M(i),N(j) {}        
        ComplexMatrix(const ComplexMatrix<T> & v) : vecbase(v),M(v.M),N(v.N) {}

        size_t rows() const { return M; }
        size_t cols() const { return N; }

        ComplexMatrixView<T> __getitem__(size_t row) { return ComplexMatrixView<T>(this,row); }

        void resize(size_t i, size_t j)
        {
            M = i;
            N = j;
            resize(M*N);
        }
        std::complex<T>& operator()(size_t i, size_t j) {
            return (*this)[i*N + j];
        }
        void fill(const std::complex<T>& c) {
            for(size_t i = 0; i < rows(); i++) 
            for(size_t j = 0; j < cols(); j++)             
                (*this)(i,j) = c;
        }
        void zero() { fill(0.0); }
        void ones() { fill(1.0); }
        
        ComplexMatrix<T> matmul(const ComplexMatrix<T> & b)
        {
            assert(N == b.M);
            ComplexMatrix<T> r(rows(),b.cols());
            r.zero();
                            
            int m = rows();
            int n = b.cols();
            int k = cols();       
            
            cppmkl::cblas_gemm(CblasRowMajor,
                        CblasNoTrans,
                        CblasNoTrans,
                        m,
                        n,
                        k,
                        std::complex<T>(1.0,0.0),
                        this->data(),
                        k,
                        b.data(),
                        n,
                        std::complex<T>(0.0,0.0),
                        r.data(),
                        n);
            return r;
        }


        ComplexMatrix<T>& operator +=  (const ComplexMatrix<T> & v) { 
            cppmkl::vadd(*this,v,*this);
            return *this;
        }
        ComplexMatrix<T>& operator -=  (const ComplexMatrix<T> & v) { 
            cppmkl::vsub(*this,v,*this);
            return *this;
        }
        ComplexMatrix<T>& operator *=  (const ComplexMatrix<T> & v) { 
            *this = matmul(v);
            return *this;
        }
        ComplexMatrix<T>& operator /=  (const ComplexMatrix<T> & v) { 
            cppmkl::vdiv(*this,v,*this);
            return *this;        
	    }

        ComplexMatrix<T>& operator +=  (const std::complex<T> & x) { 
            ComplexMatrix<T> v(rows(),cols());
            v.fill(x);
            cppmkl::vadd(*this,v,*this);
            return *this;
        }
        ComplexMatrix<T>& operator -=  (const std::complex<T> & x) { 
            ComplexMatrix<T> v(rows(),cols());;
            v.fill(x);
            cppmkl::vsub(*this,v,*this);
            return *this;
        }
        ComplexMatrix<T>& operator *=  (const std::complex<T> & x) { 
            ComplexMatrix<T> v(rows(),cols());;
            v.fill(x);
            cppmkl::vmul(*this,v,*this);
            return *this;
        }
        ComplexMatrix<T>& operator /=  (const std::complex<T> & x) { 
            ComplexMatrix<T> v(rows(),cols());;
            v.fill(x);
            cppmkl::vdiv(*this,v,*this);
            return *this;        
	    }

        ComplexMatrix<T> operator - () {
            ComplexMatrix<T> r(*this);
            for(size_t i = 0; i < rows(); i++) 
            for(size_t j = 0; j < cols(); j++)             
                r(i,j) = std::conj(r(i,j));
            return r;
        }
        ComplexMatrix<T>& operator = (const std::complex<T>& v)
        {
            fill(v);
            return *this;
        }

        Matrix<T> real() {
            Matrix<T> r(rows(),cols());
            for(size_t i = 0; i < size(); i++) r[i] = (*this)[i].real();
            return r;
        }
        Matrix<T> imag() {
            Matrix<T> r(rows(),cols());
            for(size_t i = 0; i < size(); i++) r[i] = (*this)[i].imag();
            return r;
        }
        void real(const Matrix<T> & r) {
            for(size_t i = 0; i < size(); i++) (*this)[i].real(r[i]);
        }
        void imag(const Matrix<T> & r) {
            for(size_t i = 0; i < size(); i++) (*this)[i].imag(r[i]);
        }

        Matrix<T> abs() {
            Matrix<T> r(rows(),cols());
            for(size_t i = 0; i < size(); i++) r[i] = std::abs((*this)[i]);
            return r;
        }
        Matrix<T> arg() {
            Matrix<T> r(rows(),cols());
            for(size_t i = 0; i < size(); i++) r[i] = std::abs((*this)[i]);
            return r;
        }
        ComplexMatrix<T> conj() {
            ComplexMatrix<T> r(rows(),cols());
            for(size_t i = 0; i < size(); i++) r[i] = std::conj((*this)[i]);
            return r;
        }
        ComplexMatrix<T> proj() {
            ComplexMatrix<T> r(rows(),cols());
            for(size_t i = 0; i < size(); i++) r[i] = std::proj((*this)[i]);
            return r;
        }
        Matrix<T> norm() {
            Matrix<T> r(rows(),cols());
            for(size_t i = 0; i < size(); i++) r[i] = std::norm((*this)[i]);
            return r;
        }

        #ifdef SWIG
        %extend {
            size_t size() const { return $self->size(); }
            void resize(size_t i) { $self->resize(i); }
            std::complex<T>* data() { return $self->data(); }
            void push_back(const std::complex<T> v) { $self->push_back(v); }
            void pop_back() { $self->pop_back(); }
            std::complex<T> front() { return $self->front(); }
            std::complex<T> back() { return $self->back(); }
            std::complex<T> at(size_t i) { return $self->at(i); }
            
            ComplexMatrixView<T> __getitem__(size_t i) { return ComplexMatrixView<T>(*$self,i); }
            std::complex<T>* data() { return $self->data(); }
            ComplexMatrix<T> __add__(const ComplexMatrix<T> & m) {
                ComplexMatrix<T> r(*$self);
                r += m;
                return r;
            }
            ComplexMatrix<T> __sub__(const ComplexMatrix<T> & m) {
                ComplexMatrix<T> r(*$self);
                r -= m;
                return r;
            }
            ComplexMatrix<T> __mul__(const ComplexMatrix<T> & m) {
                ComplexMatrix<T> r;
                r = $self->matmul(m);
                return r;
            }
            ComplexMatrix<T> __div__(const ComplexMatrix<T> & m) {
                ComplexMatrix<T> r(*$self);
                r /= m;
                return r;
            }
            ComplexMatrix<T> __add__(const std::complex<T> & m) {
                ComplexMatrix<T> r(*$self);
                r += m;
                return r;
            }
            ComplexMatrix<T> __sub__(const std::complex<T> & m) {
                ComplexMatrix<T> r(*$self);
                r -= m;
                return r;
            }
            ComplexMatrix<T> __mul__(const std::complex<T> & m) {
                ComplexMatrix<T> r;
                r *= m;
                return r;
            }
            ComplexMatrix<T> __div__(const std::complex<T> & m) {
                ComplexMatrix<T> r(*$self);
                r /= m;
                return r;
            }
        }
        #endif
        void print() {
            std::cout << "Matrix[" << size() << "]=";
            for(size_t i = 0; i < size()-1; i++) std::cout << (*this)[i] << ",";
            std::cout << (*this)[size()-1] << std::endl;
        }
    };

    
    template<typename T>
    std::complex<T>& ComplexMatrixView<T>::operator[](size_t i) {
        i--;
        return (*matrix)[(row-1)*matrix->N + i];
    }

    template<typename T>    
    std::complex<T>  ComplexMatrixView<T>::__getitem__(size_t i) {
        i--;
        return (*matrix)[(row-1)*matrix->N + i];
    }

    template<typename T>    
    void ComplexMatrixView<T>::__setitem__(size_t i, std::complex<T> v)
    {
        i--;
        (*matrix)[(row-1)*matrix->N + i] = v;
    }

/////////////////////////////////////
// Matrix
/////////////////////////////////////
    template<typename T>
    Matrix<T> operator * (const Matrix<T> & a, const Matrix<T> & b) {
        assert(a.N == b.M);
        Matrix<T> r(a.rows(),b.cols());
        r.zero();
                        
        int m = a.rows();
        int n = b.cols();
        int k = a.cols();       
        
        cppmkl::cblas_gemm(CblasRowMajor,
                    CblasNoTrans,
                    CblasNoTrans,
                    m,
                    n,
                    k,
                    T(1.0),
                    a.data(),
                    k,
                    b.data(),
                    n,
                    T(0.0),
                    r.data(),
                    n);       
        return r;
    } 
    
    template<typename T>
    Matrix<T> operator + (const Matrix<T> & a, const Matrix<T> & b) {
        assert(a.rows() == b.rows() && a.cols() == b.cols());
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vadd(a,b,r);            
        return r;
    }
    template<typename T>
    Matrix<T> operator / (const Matrix<T> & a, const Matrix<T> & b) {
        assert(a.rows() == b.rows() && a.cols() == b.cols());
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vdiv(a,b,r);            
        return r;
    }        
    template<typename T>
    Matrix<T> operator - (const Matrix<T> & a, const Matrix<T> & b) {
        assert(a.rows() == b.rows() && a.cols() == b.cols());
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vsub(a,b,r);            
        return r;
    }        

    template<typename T>
    Matrix<T> operator * (T a, const Matrix<T> & b) {
        Matrix<T> x(b.M,b.N);
        x.fill(a);
        Matrix<T> r(b.M,b.N);
        cppmkl::vmul(x,b,r);
        return r;
    }
    
    template<typename T>
    Matrix<T> operator + (T a, const Matrix<T> & b) {
        Matrix<T> x(b.M,b.N);
        x.fill(a);
        Matrix<T> r(b.M,b.N);
        cppmkl::vadd(x,b,r);
        return r;
    }
    template<typename T>
    Matrix<T> operator - (T a, const Matrix<T> & b) {
        Matrix<T> x(b.M,b.N);
        x.fill(a);
        Matrix<T> r(b.M,b.N);
        cppmkl::vsub(x,b,r);
        return r;
    }
    template<typename T>
    Matrix<T> operator / (T a, const Matrix<T> & b) {
        Matrix<T> x(b.M,b.N);
        x.fill(a);
        Matrix<T> r(b.M,b.N);
        cppmkl::vdiv(x,b,r);
        return r;
    }
    template<typename T>
    Matrix<T> operator * (const Matrix<T> & a, T b) {
        Matrix<T> x(a.M,a.N);
        x.fill(b);
        Matrix<T> r(a.M,a.N);
        cppmkl::vmul(a,x,r);
        return r;
    }
    
    template<typename T>
    Matrix<T> operator + (const Matrix<T> & a, T b) {
        Matrix<T> x(a.M,a.N);
        x.fill(b);
        Matrix<T> r(a.M,a.N);
        cppmkl::vadd(a,x,r);
        return r;
    }
    template<typename T>
    Matrix<T> operator - (const Matrix<T> & a, T b) {
        Matrix<T> x(a.M,a.N);
        x.fill(b);
        Matrix<T> r(a.M,a.N);
        cppmkl::vsub(a,x,r);
        return r;
    }
    template<typename T>
    Matrix<T> operator / (const Matrix<T> & a, T b) {
        Matrix<T> x(a.M,a.N);
        x.fill(b);
        Matrix<T> r(a.M,a.N);
        cppmkl::vdiv(a,x,r);
        return r;
    }    
    template<typename T>
    Matrix<T> hadamard(const Matrix<T> & a, const Matrix<T> & b) {
        Matrix<T> r(a.rows(),a.cols());
        assert(a.rows() == b.rows() && a.cols() == b.cols());
        cppmkl::vmul(a,b,r);
        return r;
    }      
    template<typename T>
    Vector<T> operator *(const Vector<T> &a, const Matrix<T> &b) {
        Vector<T> r(a.size());        
        sgemv(CblasRowMajor,CblasNoTrans,b.rows(),b.cols(),1.0,b.data(),b.cols(),a.data(),1,1.0,r.data(),1);
        return r;
    }      
    template<typename T>
    Vector<T> operator *(const Matrix<T> &a, const Vector<T> &b) {
        Vector<T> r(b.size());        
        sgemv(CblasRowMajor,CblasNoTrans,a.rows(),a.cols(),1.0,a.data(),a.cols(),b.data(),1,1.0,r.data(),1);
        return r;
    }      

    template<typename T>
    Matrix<T> copy(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::cblas_copy(a.size(),a.data(),1,r.data(),1);
        return r;
    }       

    template<typename T> T sum(const Matrix<T> & a) {        
        return cppmkl::cblas_asum(a.size(), a.data(),1);                
    }       

    template<typename T>
    Matrix<T> add(const Matrix<T> & a, const Matrix<T> & b) {
        Matrix<T> r(b);
        cppmkl::cblas_axpy(a.size(),1.0,a.data(),1,r.data(),1);
        return r;
    }       
    template<typename T>
    Matrix<T> sub(const Matrix<T> & a, const Matrix<T> & b) {
        Matrix<T> r(a);
        cppmkl::cblas_axpy(a.size(),-1.0,b.data(),1,r.data(),1);
        return r;
    }       
    template<typename T>
    T dot(const Matrix<T> & a, Matrix<T> & b) {
        return cppmkl::cblas_dot(a,b);
    }       
    template<typename T>
    T nrm2(const Matrix<T> & a) {
        Matrix<T> r(a);
        return cppmkl::cblas_nrm2(a);        
    }       
    
    template<typename T>
    void scale(Matrix<T> & x, T alpha) {
        cppmkl::cblas_scal(x.size(),alpha,x.data(),1);
    }


}