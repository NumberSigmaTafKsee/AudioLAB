#pragma once

namespace Casino::MKL
{
    template<typename T> struct Matrix;
    template<typename T>
    struct MatrixView 
    {
        Matrix<T> * matrix;
        size_t row;
        MatrixView(Matrix<T> * m, size_t r) {
            matrix = m;
            row = r;
        }

        T& operator[](size_t i);
        T  __getitem__(size_t i);
        void __setitem__(size_t i, T x);
    };
    
    template<typename T>
    struct Matrix  : public vector_base<T>
    {                
        using vector_base<T>::size;
        using vector_base<T>::resize;
        using vector_base<T>::data;
        using vector_base<T>::push_back;
        using vector_base<T>::pop_back;
        using vector_base<T>::front;
        using vector_base<T>::back;
        using vector_base<T>::at;
        using vector_base<T>::operator [];
        using vector_base<T>::operator =;

        size_t M;
        size_t N;

        Matrix() { M = N = 0; }
        Matrix(size_t m, size_t n) {
            resize(m,n);
            assert(M > 0);
            assert(N > 0);
        }
        Matrix(const Matrix<T> & m) {
            
            resize(m.M,m.N);
            memcpy(data(),m.data(),size()*sizeof(T));            
        }
        Matrix(T * ptr, size_t m, size_t n)
        {
            resize(m,n);            
            memcpy(data(),ptr,m*n*sizeof(T));
        }
        Matrix<T>& operator = (const T v) {
            fill(v);
            return *this;
        }

        Matrix<T>& operator = (const Matrix<T> & m) {            
            resize(m.M,m.N);
            memcpy(data(),m.data(),size()*sizeof(T));            
            return *this;
        }
        Matrix<T>& operator = (const Vector<T> & m) {            
            resize(1,m.size());
            memcpy(data(),m.data(),size()*sizeof(T));            
            return *this;
        }

        T& operator()(size_t i) { return (*this)[i]; }
        T  operator()(size_t i) const { return (*this)[i]; }

        size_t rows() const { return M; }
        size_t cols() const { return N; }
        
        
        MatrixView<T> __getitem__(size_t row) { return MatrixView<T>(this,row); }

        T get(size_t i,size_t j) { return (*this)(i,j); }
        void set(size_t i, size_t j, T v) { (*this)(i,j) = v; }

        Matrix<T> row(size_t m) { 
            Matrix<T> r(1,cols());
            for(size_t i = 0; i < cols(); i++) r(0,i) = (*this)(m,i);
            return r;
        }
        Matrix<T> col(size_t c) { 
            Matrix<T> r(cols(),1);
            for(size_t i = 0; i < rows(); i++) r(i,0) = (*this)(i,c);
            return r;
        }
        Vector<T> row_vector(size_t m) { 
            Vector<T> r(cols());
            for(size_t i = 0; i < cols(); i++) r[i] = (*this)(m,i);
            return r;
        }
        Vector<T> col_vector(size_t c) { 
            Vector<T> r(cols());
            for(size_t i = 0; i < rows(); i++) r[i] = (*this)(i,c);
            return r;
        }
        void row(size_t m, const Vector<T> & v)
        {            
            for(size_t i = 0; i < cols(); i++) (*this)(m,i) = v[i];
        }
        void col(size_t n, const Vector<T> & v)
        {
            for(size_t i = 0; i < rows(); i++) (*this)(i,n) = v[i];
        }
        
        void resize(size_t r, size_t c) {
            M = r;
            N = c;
            resize(r*c);            
        }

        T& operator()(size_t i, size_t j) {             
            return (*this)[i*N + j]; }

        T  operator()(size_t i, size_t j) const {             
            return (*this)[i*N + j]; }
        
        
        std::ostream& operator << (std::ostream & o )
        {
            for(size_t i = 0; i < rows(); i++)
            {
                for(size_t j = 0; j < cols(); j++)
                    std::cout << (*this)(i,j) << ",";
                std::cout << std::endl;
            }
            return o;
        }

        Matrix<T> operator - () {
            Matrix<T> r(*this);
            return T(-1.0)*r;
        }
        Matrix<T> addToEachRow(Vector<T> & v) {
            Matrix<T> r(*this);
            for(size_t i = 0; i < M; i++)
            {
                for(size_t j = 0; j < N; j++)
                {
                    r(i,j) += v[j];
                }
            }
            return r;
        }
        Matrix<T> addToEachRow(Matrix<T> & v) {
            Matrix<T> r(*this);
            for(size_t i = 0; i < M; i++)
            {
                for(size_t j = 0; j < N; j++)
                {
                    r(i,j) += v(0,j);
                }
            }
            return r;
        }
        Matrix<T> eval() {
            Matrix<T> r(*this);
            return r;
        }
        void printrowscols() const {
            std::cout << "rows=" << rows() << " cols=" << cols() << std::endl;
        }
        
        Matrix<T> matmul(const Matrix<T> & b)
        {            
            assert(N == b.M);
            Matrix<T> r(rows(),b.cols());
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
                        T(1.0),
                        this->data(),
                        k,
                        b.data(),
                        n,
                        T(0.0),
                        r.data(),
                        n);
            return r;
        }
        
        Matrix<T>& operator += (const T b)
        {
            Matrix<T> x(rows(),cols());
            x.fill(b);
            cppmkl::vadd(*this,x,*this);
            return *this;
        }
        Matrix<T>& operator -= (const T b)
        {
            Matrix<T> x(rows(),cols());
            x.fill(b);
            cppmkl::vsub(*this,x,*this);
            return *this;
        }
        Matrix<T>& operator *= (const T b)
        {
            Matrix<T> x(rows(),cols());
            x.fill(b);
            cppmkl::vmul(*this,x,*this);
            return *this;
        }
        Matrix<T>& operator /= (const T b)
        {
            Matrix<T> x(rows(),cols());
            x.fill(b);
            cppmkl::vdiv(*this,x,*this);
            return *this;
        }
        
        Matrix<T>& operator += (const Matrix<T> b)
        {
            assert(rows() == b.rows() && cols() == b.cols());            
            cppmkl::vadd(*this,b,*this);
            return *this;
        }
        Matrix<T>& operator -= (const Matrix<T> b)
        {
            assert(rows() == b.rows() && cols() == b.cols());
            cppmkl::vsub(*this,b,*this);
            return *this;
        }
        Matrix<T>& operator *= (const Matrix<T> b)
        {
            assert(rows() == b.rows() && cols() == b.cols());
            cppmkl::vmul(*this,b,*this);
            return *this;
        }
        Matrix<T>& operator /= (const Matrix<T> b)
        {
            assert(rows() == b.rows() && cols() == b.cols());            
            cppmkl::vdiv(*this,b,*this);
            return *this;
        }
        
        Matrix<T> operator + (const T b)
        {
            Matrix<T> x(rows(),cols()),y(rows(),cols());
            x.fill(b);
            cppmkl::vadd(*this,x,y);
            return y;
        }
        Matrix<T> operator - (const T b)
        {
            Matrix<T> x(rows(),cols()),y(rows(),cols());
            x.fill(b);
            cppmkl::vsub(*this,x,y);
            return y;
        }
        Matrix<T> operator * (const T b)
        {
            Matrix<T> x(rows(),cols()),y(rows(),cols());
            x.fill(b);
            cppmkl::vmul(*this,x,y);
            return y;
        }
        Matrix<T> operator / (const T b)
        {
            Matrix<T> x(rows(),cols()),y(rows(),cols());
            x.fill(b);
            cppmkl::vdiv(*this,x,y);
            return y;
        }
        
        Matrix<T> operator + (const Matrix<T> b)
        {
            Matrix<T> r(rows(),cols());
            assert(rows() == b.rows() && cols() == b.cols());            
            cppmkl::vadd(*this,b,r);
            return r;
        }
        Matrix<T> operator - (const Matrix<T> b)
        {
            Matrix<T> r(rows(),cols());
            assert(rows() == b.rows() && cols() == b.cols());            
            cppmkl::vsub(*this,b,r);
            return r;
        }
        Matrix<T> operator * (const Matrix<T> b)
        {
            return matmul(b);
        }
        Matrix<T> operator / (const Matrix<T> b)
        {
            Matrix<T> r(rows(),cols());
            assert(rows() == b.rows() && cols() == b.cols());            
            cppmkl::vdiv(*this,b,r);
            return r;
        }

        Matrix<T> hadamard(const Matrix<T> & m) {
            Matrix<T> r(rows(),cols());
            assert(rows() == m.rows() && cols() == m.cols());
            cppmkl::vmul(*this,m,r);
            return r;
        }
        
        Matrix<T> transpose() {
            Matrix<T> r(cols(),rows());
            for(size_t i = 0; i < rows(); i++)
                for(size_t j = 0; j < cols(); j++)
                    r(j,i) =(*this)(i,j);
            return r;
        }
        Matrix<T> t() {
            return transpose();
        }
        void transposeInto(Matrix<T> &m) {
            m = transpose();
        }
        void zero() {
            memset(data(),0x00,size()*sizeof(T));
        }
        void zeros() { zero(); }

        void fill(T x) {
            for(size_t i = 0; i < size(); i++) (*this)[i] = x;
        }
        void ones() {
            fill((T)1);
        }
        void minus_ones() {
            fill((T)-1.0);
        }

        T& cell(size_t i, size_t j) { return (*this)(i,j); }
        void swap_cells(size_t i1, size_t j1, size_t i2, size_t j2)
        {
            T x = (*this)(i1,j1);
            (*this)(i1,j1) = (*this)(i2,j2);
            (*this)(i2,j2) = x;
        }
            
        
        void random(T min = T(0), T max = T(1)) {
            Default r;
            for(size_t i = 0; i < size(); i++) (*this)[i] = r.random(min,max);
        }
        void identity() {
            size_t x = 0;
            zero();
            for(size_t i = 0; i < rows(); i++)
            {
                (*this)[i*N + x++] = 1;
            }
        }

        
        void print() const {
            std::cout << "Matrix[" << M << "," << N << "]=";
            for(size_t i = 0; i < rows(); i++) 
            {
                for(size_t j = 0; j < cols(); j++) std::cout << (*this)(i,j) << ",";
                std::cout << std::endl;
            }
        }

        T sum() {        
            return cppmkl::cblas_asum(size(), data(),1);                
        }
        T prod() {
            T p = (T)1.0;
            for(size_t i = 0; i < size(); i++) p *= (*this)[i];
            return p;
        }
        T mean() {
            return sum()/(T)size();
        }
        T geometric_mean() {
            T r = prod();
            return std::pow(r,(T)1.0/(T)size());
        }
        T harmonic_mean() {
            T r = 1.0/sum();
            return size()*std::pow(r,(T)-1.0);
        }
        T stddev() {
            T u = sum();
            T r = 0;
            for(size_t i = 0; i < size(); i++)
                r += std::pow((*this)[i]-u,2.0);
            return r/(T)size();
        }
           
        // eigen compatibility
        void setZero() { zero(); }
        void setOnes() { ones(); }
        void setRandom() { random(); }
        //Matrix<T> matrix() { return eval(); }        
        void setIdentity() { identity(); }

        Matrix<T> cwiseMin(T min) {
            Matrix<T> r(rows(),cols());
            for(size_t i = 0; i < rows(); i++)
            for(size_t j = 0; j < cols(); j++)
                if((*this)(i,j) > min) r(i,j) = min;
                else r(i,j) = (*this)(i,j);
            return r;
        }
        Matrix<T> cwiseMax(T max) {
            Matrix<T> r(rows(),cols());
            for(size_t i = 0; i < rows(); i++)
            for(size_t j = 0; j < cols(); j++)
                if((*this)(i,j) < max) r(i,j) = max;
                else r(i,j) = (*this)(i,j);
            return r;
        }
        Matrix<T> cwiseProd(Matrix<T> & x) {
            Matrix<T> r(*this);
            r *= x;
            return r;
        }
        Matrix<T> cwiseAdd(Matrix<T> & x) {
            Matrix<T> r(*this);
            r += x;
            return r;
        }
        Matrix<T> cwiseSub(Matrix<T> & x) {
            Matrix<T> r(*this);
            r -= x;
            return r;
        }
        Matrix<T> cwiseDiv(Matrix<T> & x) {
            Matrix<T> r(*this);
            r /= x;
            return r;
        }


        #ifdef SWIG
        %extend {
            MatrixView<T> __getitem__(size_t i) { return MatrixView<T>(*$self,i); }
            T* data() { return $self->data(); }
            Matrix<T> __add__(const Matrix<T> & m) {
                Matrix<T> r(*$self);
                r += m;
                return r;
            }
            Matrix<T> __add__(const T m) {
                Matrix<T> r(*$self);
                r += m;
                return r;
            }
            Matrix<T> __sub__(const Matrix<T> & m) {
                Matrix<T> r(*$self);
                r -= m;
                return r;
            }
            Matrix<T> __sub__(const T m) {
                Matrix<T> r(*$self);
                r -= m;
                return r;
            }
            Matrix<T> __mul__(const Matrix<T> & m) {
                return $self->matmul(m);
            }
            Matrix<T> __mul__(const T m) {
                Matrix<T> r(*$self);
                r *= m;
                return r;
            }
            Matrix<T> __div__(const Matrix<T> & m) {
                Matrix<T> r(*$self);
                r /= m;
                return r;
            }
            Matrix<T> __mul__(const T m) {
                Matrix<T> r(*$self);
                r /= m;
                return r;
            }
        }
        #endif
         

    };
    
    template<typename T>
    T& MatrixView<T>::operator[](size_t i) {
        i--;
        return (*matrix)[(row-1)*matrix->N + i];
    }

    template<typename T>    
    T  MatrixView<T>::__getitem__(size_t i) {
        i--;
        return (*matrix)[(row-1)*matrix->N + i];
    }

    template<typename T>    
    void MatrixView<T>::__setitem__(size_t i, T v)
    {
        i--;
        (*matrix)[(row-1)*matrix->N + i] = v;
    }

/////////////////////////////////////
// Matrix
/////////////////////////////////////    
    template<typename T>
    Matrix<T> sqr(Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vsqr(a,r);
        return r;
    }
    template<typename T>
    Matrix<T> abs(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vabs(a,r);
        return r;
    }
    template<typename T>
    Matrix<T> inv(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vinv(a,r);
        return r;            
    }
    template<typename T>
    Matrix<T> sqrt(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vsqrt(a,r);
        return r;
    }
    template<typename T>
    Matrix<T> rsqrt(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vinvsqrt(a,r);
        return r;
    }
    template<typename T>
    Matrix<T> cbrt(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vcbrt(a,r);
        return r;
    }
    template<typename T>
    Matrix<T> rcbrt(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vinvcbrt(a,r);
        return r;
    }
    template<typename T>
    Matrix<T> pow(const Matrix<T> & a,const Matrix<T> & b) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vpow(a,b,r);
        return r;
    }
    template<typename T>
    Matrix<T> pow2o3(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vpow2o3(a,r);
        return r;
    }
    template<typename T>
    Matrix<T> pow3o2(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vpow3o2(a,r);
        return r;
    }
    template<typename T>
    Matrix<T> pow(const Matrix<T> & a,const T b) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vpowx(a,b,r);
        return r;
    }
    template<typename T>
    Matrix<T> hypot(const Matrix<T> & a,const Matrix<T> & b) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vhypot(a,b,r);
        return r;
    }
    template<typename T>
    Matrix<T> exp(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vexp(a,r);
        return r;
    }
    template<typename T>
    Matrix<T> exp2(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vexp2(a,r);
        return r;
    }
    template<typename T>
    Matrix<T> exp10(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vexp10(a,r);
        return r;
    }
    template<typename T>
    Matrix<T> expm1(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vexpm1(a,r);
        return r;
    }
    template<typename T>
    Matrix<T> ln(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vln(a,r);
        return r;
    }
    template<typename T>
    Matrix<T> log(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vln(a,r);
        return r;
    }
    template<typename T>
    Matrix<T> log10(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vlog10(a,r);
        return r;
    }
    template<typename T>
    Matrix<T> log2(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vlog2(a,r);        
        return r;
    }
    template<typename T>
    Matrix<T> logb(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vlogb(a,r);
        return r;
    }
    template<typename T>
    Matrix<T> log1p(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vlog1p(a,r);
        return r;
    }
    template<typename T>
    Matrix<T> cos(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vcos(a,r);
        return r;
    }
    template<typename T>
    Matrix<T> sin(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vsin(a,r);
        return r;            
    }
    template<typename T>
    Matrix<T> tan(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vtan(a,r);
        return r;
    }
    template<typename T>
    Matrix<T> cosh(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vcosh(a,r);
        return r;
    }
    template<typename T>
    Matrix<T> sinh(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vsinh(a,r);
        return r;            
    }
    template<typename T>
    Matrix<T> tanh(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vtanh(a,r);
        return r;            
    }
    template<typename T>
    Matrix<T> acos(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vacos(a,r);
        return r;
    }
    template<typename T>
    Matrix<T> asin(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vasin(a,r);
        return r;            
    }
    template<typename T>
    Matrix<T> atan(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vatan(a,r);
        return r;
    }
    template<typename T>
    Matrix<T> atan2(const Matrix<T> & a,const Matrix<T> &n) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vatan2(a,n,r);
        return r;
    }
    template<typename T>
    Matrix<T> acosh(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vacosh(a,r);
        return r;
    }
    template<typename T>
    Matrix<T> asinh(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vasinh(a,r);
        return r;            
    }
    template<typename T>
    Matrix<T> atanh(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vatanh(a,r);
        return r;
    }        
    template<typename T>
    void sincos(const Matrix<T> & a, Matrix<T> & b, Matrix<T> & r) {        
        cppmkl::vsincos(a,b,r);
    }
    template<typename T>
    Matrix<T> erf(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::verf(a,r);
        return r;
    }
    template<typename T>
    Matrix<T> erfinv(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::verfinv(a,r);
        return r;
    }
    template<typename T>
    Matrix<T> erfc(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::verfc(a,r);
        return r;
    }
    template<typename T>
    Matrix<T> cdfnorm(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vcdfnorm(a,r);
        return r;
    }    
    template<typename T>
    Matrix<T> cdfnorminv(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vcdfnorminv(a,r);
        return r;
    }    
    template<typename T>
    Matrix<T> floor(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vfloor(a,r);
        return r;
    }    
    template<typename T>
    Matrix<T> ceil(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vceil(a,r);
        return r;
    }    
    template<typename T>
    Matrix<T> trunc(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vtrunc(a,r);
        return r;
    }    
    template<typename T>
    Matrix<T> round(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vround(a,r);
        return r;        
    }    
    template<typename T>
    Matrix<T> nearbyint(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vnearbyint(a,r);
        return r;
    }    
    template<typename T>
    Matrix<T> rint(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vrint(a,r);
        return r;
    }    
    template<typename T>
    Matrix<T> fmod(const Matrix<T> & a, Matrix<T> & b) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vmodf(a,b,r);
        return r;
    }    
    template<typename T>
    Matrix<T> mulbyconj(const Matrix<T> & a, const Matrix<T> & b) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vmulbyconj(a,b,r);
        return r;
    }    
    template<typename T>
    Matrix<T> conj(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vconj(a,r);
        return r;
    }    
    template<typename T>
    Matrix<T> arg(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::varg(a,r);
        return r;
    }    
    template<typename T>
    Matrix<T> CIS(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vCIS(a,r);
        return r;
    }    
    template<typename T>
    Matrix<T> cospi(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vcospi(a,r);
        return r;
    }    
    template<typename T>
    Matrix<T> sinpi(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vsinpi(a,r);
        return r;
    }    
    template<typename T>
    Matrix<T> tanpi(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vtanpi(a,r);
        return r;
    }    
    template<typename T>
    Matrix<T> acospi(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vacospi(a,r);
        return r;
    }    
    template<typename T>
    Matrix<T> asinpi(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vasinpi(a,r);
        return r;
    }    
    template<typename T>
    Matrix<T> atanpi(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vatanpi(a,r);
        return r;
    }    
    template<typename T>
    Matrix<T> atan2pi(const Matrix<T> & a, Matrix<T> & b) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vatan2pi(a,b,r);
        return r;
    }    
    template<typename T>
    Matrix<T> cosd(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vcosd(a.size(),a.data(),r.data());
        return r;
    }    
    template<typename T>
    Matrix<T> sind(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vsind(a.size(),a.data(),r.data());
        return r;
    }    
    template<typename T>
    Matrix<T> tand(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vtand(a.size(),a.data(),r.data());
        return r;
    }       
    template<typename T>
    Matrix<T> lgamma(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vlgamma(a,r);
        return r;
    }       
    template<typename T>
    Matrix<T> tgamma(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vtgamma(a,r);
        return r;
    }       
    template<typename T>
    Matrix<T> expint1(const Matrix<T> & a) {
        Matrix<T> r(a.rows(),a.cols());
        cppmkl::vexpint1(a,r);
        return r;
    }       

    template<typename T>
    Vector<T> copy(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::cblas_copy(a.size(),a.data(),1,r.data(),1);
        return r;
    }       

    template<typename T> T sum(const Vector<T> & a) {        
        return cppmkl::cblas_asum(a.size(), a.data(),1);        
    }       

    template<typename T>
    Vector<T> add(const Vector<T> & a, const Vector<T> & b) {
        Vector<T> r(b);
        cppmkl::cblas_axpy(a.size(),1.0,a.data(),1,r.data(),1);
        return r;
    }       
    template<typename T>
    Vector<T> sub(const Vector<T> & a, const Vector<T> & b) {
        Vector<T> r(a);
        cppmkl::cblas_axpy(a.size(),-1.0,b.data(),1,r.data(),1);
        return r;
    }       
    template<typename T>
    T dot(const Vector<T> & a, Vector<T> & b) {
        return cppmkl::cblas_dot(a,b);
    }       
    template<typename T>
    T nrm2(const Vector<T> & a) {
        Vector<T> r(a);
        return cppmkl::cblas_nrm2(a);        
    }       
    
    template<typename T>
    void scale(Vector<T> & x, T alpha) {
        cppmkl::cblas_scal(x.size(),alpha,x.data(),1);
    }


    template<typename T>
    size_t min_index(const Matrix<T> & m) { return cppmkl::cblas_iamin(m.size(),m.data(),1); }
    template<typename T>
    size_t max_index(const Matrix<T> & m) { return cppmkl::cblas_iamax(m.size(),m.data(),1); }
    
    template<typename T>
    size_t min_index(const Vector<T> & v) { return cppmkl::cblas_iamin(v.size(),v.data(),1); }
    template<typename T>
    size_t max_index(const Vector<T> & v) { return cppmkl::cblas_iamax(v.size(),v.data(),1); }

    template<typename T>
    Matrix<T> deinterleave(size_t ch, size_t n, T * samples)
    {
        Matrix<T> r(ch,n);
        for(size_t i = 0; i < ch; i++)        
            for(size_t j = 0; j < n; j++)
                r(i,j) = samples[j*ch + i];
        return r;
    }
    template<typename T>
    std::ostream& operator << (std::ostream & o, const Matrix<T> & m )
    {
        for(size_t i = 0; i < m.rows(); i++)
        {
            for(size_t j = 0; j < m.cols(); j++)
                o << m(i,j) << ",";
            o << std::endl;
        }
        return o;
    }
    template<typename T>
    Vector<T> interleave(const Matrix<T> & m)
    {
        Vector<T> r(m.rows()*m.cols());
        int ch = m.rows();
        for(size_t i = 0; i < m.rows(); i++)
        for(size_t j = 0; j < m.cols(); j++)
            r[j*ch+i] = m(i,j);
        return r;
    }
    template<typename T>
    Vector<T> channel(size_t c, const Matrix<T> & m) {
        Vector<T> r(m.cols());
        for(size_t i = 0; i < m.cols(); i++) r[i] = m(c,i);
        return r;
    }

    template<typename T>
    Matrix<T> matmul(const Matrix<T> & a, const Matrix<T> & b)
    {            
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

}