// I have scultped out most of the API.
// I don't really use alot of it and don't even know how alot of it is even used.
// Your on your own with that stuff but it should at least work enough for you to improve it.
// Ideally, it would just wrap armadillo in the namespace.
// I am doing it this way as I learned with Eigen it was more conveneint to just do it this way.
// Eigen was severely difficult to wrap this is much easier.

#pragma once

#include<armadillo>
#include<vector>
#include<complex>

namespace Matlab {

    template<typename T>
    struct Vector {
        arma::Col<T> vector;

        Vector() = default;
        Vector(size_t n) { vector = arma::Col<T>(n); }
        Vector(const Vector<T> & v) { vector = v.vector; }
        Vector(const arma::Col<T> & v) { vector = v; }
        //Vector(std::initializer_list<T> & l) { vector = arma::Col<T>(l.begin(),l.end()}; }
        ~Vector() = default;

        Vector<T>& operator = (const Vector<T> & v) {
            vector = v.vector;
            return *this;
        }
        Vector<T>& operator = (const arma::Col<T> & v) {
            vector = v;
            return *this;
        }

        Vector<T> operator - () {
            Vector<T> r(-vector);
            return r;
        }
        Vector<T> operator + (const Vector<T> & b) {
            Vector<T> r(*this);
            r.vector = vector + b.vector;
            return r;
        }
        Vector<T> operator - (const Vector<T> & b) {
            Vector<T> r(*this);
            r.vector = vector - b.vector;
            return r;
        }
        Vector<T> operator * (const Vector<T> & b) {
            Vector<T> r(*this);
            r.vector = vector * b.vector;
            return r;
        }
        Vector<T> operator / (const Vector<T> & b) {
            Vector<T> r(*this);
            r.vector = vector / b.vector;
            return r;
        }
        Vector<T> operator % (const Vector<T> & b) {
            Vector<T> r(*this);
            r.vector = vector % b.vector;
            return r;
        }

        Vector<T> operator + (const T b) {
            Vector<T> r(*this);
            r.vector = vector + b;
            return r;
        }
        Vector<T> operator - (const T b) {
            Vector<T> r(*this);
            r.vector = vector - b;
            return r;
        }
        Vector<T> operator * (const T b) {
            Vector<T> r(*this);
            r.vector = vector * b;
            return r;
        }
        Vector<T> operator / (const T b) {
            Vector<T> r(*this);
            r.vector = vector / b;
            return r;        
        }
        Vector<T> operator % (const T b) {
            Vector<T> r(*this);
            r.vector = vector % b;
            return r;
        }
        
        T& at(size_t r) { return vector.at(r); }

        T& operator()(size_t i) { return vector(i); } 
        T& operator[](size_t i) { return vector(i); } 
        T __gettiem(size_t i) { return vector(i); }
        void __setitem(size_t i, T v) { vector(i) = v; }

        size_t size() const { return vector.size(); }
        void zeros() { vector.zeros(); }  
        void zeros(size_t r) { vector.zeros(r); }  
        void ones()  { vector.ones(); }
        void ones(size_t r)  {vector.ones(r); }        
        void randu() { vector.randu(); }        
        void random() { vector.randu(); }
        void random(const T min, const T max) { vector.randu(); vector = min + (max-min)*vector; }        
        void randu(size_t r) { vector.randu(r); }        
        void randn() { vector.randn(); }
        void randn(size_t r) { vector.randn(r); }
        void fill(const T val) { vector.fill(val); }        
        //void clean(const T v) { vector.clean(v); }
        void replace(const T oldv, const T newv) { vector.replace(oldv,newv); }
        void clamp(const T min, const T max) { vector.clamp(min,max); }  
        void set_size(size_t r) { vector.set_size(r); }        
        void resize(size_t n) { vector.resize(n); }
        void copy_size(const Vector<T> & m) { vector.copy_size(m.vector); }
        void reset() { vector.reset(); }
        T*   data() { return vector.memptr(); }    

        void print(const std::string & header = "") { vector.print(header); }
        void raw_print(const std::string & header = "") { vector.raw_print(header); }
        void brief_print(const std::string & header = "") { vector.brief_print(header); }

    };

    template<typename T>
    struct RowVector {
        arma::Row<T> vector;

        RowVector() = default;
        RowVector(size_t n) { vector = arma::Row<T>(n); }
        RowVector(const RowVector<T> & v) { vector = v.vector; }
        RowVector(const arma::Row<T> & v) { vector = v; }
        //RowVector(std::initializer_list<T> & l) { vector = arma::Row<T>(l}; }
        ~RowVector() = default;

        RowVector<T>& operator = (const RowVector<T> & v) {
            vector = v.vector;
            return *this;
        }
        RowVector<T>& operator = (const arma::Row<T> & v) {
            vector = v;
            return *this;
        }

        RowVector<T> operator - () {
            RowVector<T> r(-vector);
            return r;
        }

        RowVector<T> operator + (const RowVector<T> & b) {
            RowVector<T> r(*this);
            r.vector = vector + b.vector;
            return r;
        }
        RowVector<T> operator - (const RowVector<T> & b) {
            RowVector<T> r(*this);
            r.vector = vector - b.vector;
            return r;
        }
        RowVector<T> operator * (const RowVector<T> & b) {
            RowVector<T> r(*this);
            r.vector = vector * b.vector;
            return r;
        }
        RowVector<T> operator / (const RowVector<T> & b) {
            RowVector<T> r(*this);
            r.vector = vector / b.vector;
            return r;
        }
        RowVector<T> operator % (const RowVector<T> & b) {
            RowVector<T> r(*this);
            r.vector = vector % b.vector;
            return r;
        }

        RowVector<T> operator + (const T b) {
            RowVector<T> r(*this);
            r.vector = vector + b;
            return r;
        }
        RowVector<T> operator - (const T b) {
            RowVector<T> r(*this);
            r.vector = vector - b;
            return r;
        }
        RowVector<T> operator * (const T b) {
            RowVector<T> r(*this);
            r.vector = vector * b;
            return r;
        }
        RowVector<T> operator / (const T b) {
            RowVector<T> r(*this);
            r.vector = vector / b;
            return r;        
        }
        RowVector<T> operator % (const T b) {
            RowVector<T> r(*this);
            r.vector = vector % b;
            return r;
        }

        
        T& at(size_t r) { return vector.at(r); }
        size_t size() const { return vector.size(); }
        T& operator()(size_t i) { return vector(i); } 
        T& operator[](size_t i) { return vector(i); } 
        T __gettiem(size_t i) { return vector(i); }
        void __setitem(size_t i, T v) { vector(i) = v; }

        void zeros() { vector.zeros(); }  
        void zeros(size_t r) { vector.zeros(r); }  
        void ones()  { vector.ones(); }
        void ones(size_t r)  {vector.ones(r); }        
        void randu() { vector.randu(); }
        void random() { vector.randu(); }
        void random(const T min, const T max) { vector.ranu(); vector = min + (max-min)*vector; }        
        void randu(size_t r) { vector.randu(r); }
        void randn() { vector.randn(); }
        void randn(size_t r) { vector.randn(r); }
        void fill(const T val) { vector.fill(val); }        
        void clean(const T v) { vector.clean(v); }
        void replace(const T oldv, const T newv) { vector.replace(oldv,newv); }
        void clamp(const T min, const T max) { vector.clamp(min,max); }  
        void set_size(size_t r) { vector.set_size(r); }        
        void resize(size_t n) { vector.resize(n); }
        void copy_size(const RowVector<T> & m) { vector.copy_size(m.vector); }
        void reset() { vector.reset(); }
        T*   data() { return vector.mem_ptr(); }    

        void print(const std::string & header = "") { vector.print(header); }
        void raw_print(const std::string & header = "") { vector.raw_print(header); }
        void brief_print(const std::string & header = "") { vector.brief_print(header); }
    };

    template<typename T>
    struct MatrixView {
        arma::Mat<T> * m;
        size_t row;

        MatrixView(arma::Mat<T> * x, size_t r) : row(r),m(x) {}

        T& operator [](size_t i) { return (*m)(row,i); }
        T  __getitem(size_t i) { return (*m)(row,i); }
        void __setitem(size_t i, const T val) { (*m)(row,i) = val; }
    };

    template<typename T>
    struct Matrix {
        arma::Mat<T> matrix;

        Matrix() = default;
        Matrix(size_t rows, size_t cols) { matrix = arma::Mat<T>(rows,cols); }
        Matrix(const Matrix<T> & m) { matrix = m.matrix; }
        Matrix(const arma::Mat<T> & m) { matrix = m; }
        //Matrix(const std::vector<T> & v, size_t rows, size_t cols, bool copy_mem = true, bool strict = false) { matrix = arma::Mat<T>(v.data(),rows,cols,copy_mem,strict); }
        Matrix(const std::initializer_list<T> & l) { matrix = arma::Mat<T>(l); }
        ~Matrix() = default;

        Matrix<T>& operator = (const Matrix<T> & m) {
            matrix = m.matrix;
            return *this;
        }

        Matrix<T> operator - () {
            Matrix<T> r(-matrix);
            return r;
        }
        Matrix<T> operator + (const Matrix<T> & b) {
            Matrix<T> r(*this);
            r.matrix = matrix + b.matrix;
            return r;
        }
        Matrix<T> operator - (const Matrix<T> & b) {
            Matrix<T> r(*this);
            r.matrix = matrix - b.matrix;
            return r;
        }
        Matrix<T> operator * (const Matrix<T> & b) {
            Matrix<T> r(*this);
            r.matrix = matrix * b.matrix;
            return r;
        }
        Matrix<T> operator / (const Matrix<T> & b) {
            Matrix<T> r(*this);
            r.matrix = matrix / b.matrix;
            return r;
        }
        Matrix<T> operator % (const Matrix<T> & b) {
            Matrix<T> r(*this);
            r.matrix = matrix % b.matrix;
            return r;
        }

        Matrix<T> operator + (const T b) {
            Matrix<T> r(*this);
            r.matrix = matrix + b;
            return r;
        }
        Matrix<T> operator - (const T b) {
            Matrix<T> r(*this);
            r.matrix = matrix - b;
            return r;
        }
        Matrix<T> operator * (const T b) {
            Matrix<T> r(*this);
            r.matrix = matrix * b;
            return r;
        }
        Matrix<T> operator / (const T b) {
            Matrix<T> r(*this);
            r.matrix = matrix / b;
            return r;
        }
        Matrix<T> operator % (const T b) {
            Matrix<T> r(*this);
            r.matrix = matrix % b;
            return r;
        }

        size_t rows() const { return matrix.n_rows; }
        size_t cols() const { return matrix.n_cols; }

        
        void set_row(size_t r, const RowVector<T> & v) {
            matrix.row(r) = v.vector;
        }
        void set_col(size_t c, const Vector<T> & v) {
            matrix.col(c) = v.vector;
        }
        Vector<T>    get_col(size_t c) { return Vector<T>(matrix.col(c)); }
        RowVector<T> get_row(size_t r) { return RowVector<T>(matrix.row(r)); }

        T& operator()(size_t i, size_t j) { return matrix(i,j); } 
        MatrixView<T> operator[](size_t row) { return MatrixView<T>(&matrix,row); }    
        MatrixView<T> __getitem(size_t row) { return MatrixView<T>(&matrix,row); }
        void __setitem(size_t col, const Vector<T> & val) { matrix.col(col) = val.vector; }
        void __setitem(size_t row, const RowVector<T> & val) { matrix.row(row) = val.vector; }
        T& at(size_t r, size_t c) { return matrix.at(r,c);}        

        void zeros() { matrix.zeros(); }  
        void zeros(size_t r, size_t c) { matrix.zeros(r,c); }  
        void ones()  { matrix.ones(); }
        void ones(size_t r, size_t c)  { matrix.ones(r,c); }
        void eye()   { matrix.eye(); }
        void eye(size_t r, size_t c)   { matrix.eye(r,c); }
        void randu() { matrix.randu(); }
        void randu(size_t r, size_t c) { matrix.randu(r,c); }
        void randn() { matrix.randn(); }        
        void randn(size_t r, size_t c) { matrix.randn(r,c); }
        void fill(const T val) { matrix.fill(val); }
        //void clean( T v) { matrix.clean(v); }
        void replace(const T oldv, const T newv) { matrix.replace(oldv,newv); }
        void clamp(const T min, const T max) { matrix.clamp(min,max); }  
        void set_size(size_t r, size_t c) { matrix.set_size(r,c); }
        void reshape(size_t r, size_t c) { matrix.reshape(r,c); }
        void resize(size_t r, size_t c) { matrix.resize(r,c); }
        void copy_size(const Matrix<T> & m) { matrix.copy_size(m.matrix); }
        void reset() { matrix.reset(); }
        T*   data() { return matrix.memptr(); }    
        T*   colptr(size_t n) { return matrix.colptr(n); }

        Vector<T> as_col() { return Vector<T>(matrix.as_col()); }
        RowVector<T> as_row() { return RowVector<T>(matrix.as_row()); }

        Matrix<T> t() { return Matrix<T>(matrix.t()); }
        //Matrix<T> i() { return Matrix<T>(matrix.i()); }

        T min() { return matrix.min(); }
        T max() { return matrix.max(); }

        size_t index_min() { return matrix.index_min(); }
        size_t index_max() { return matrix.index_max(); }
        
        Matrix<T> eval() { return Matrix<T>(matrix.eval()); }

        bool is_empty() { return matrix.is_empty(); }
        bool is_trimatu() { return matrix.is_trimatu(); }
        bool is_trimatl() { return matrix.is_trimatl(); }
        bool is_diagmat() { return matrix.is_diagmat(); }
        bool is_square()  { return matrix.is_square(); }
        bool is_symmetric() { return matrix.is_symmetric(); }
        bool is_zero() { return matrix.is_zero(); }
        bool is_finite() { return matrix.is_finite(); }
        bool has_inf() { return matrix.has_inf(); }
        bool has_nan() { return matrix.has_nan(); }

        void print(const std::string & header = "") { matrix.print(header); }
        void raw_print(const std::string & header = "") { matrix.raw_print(header); }
        void brief_print(const std::string & header = "") { matrix.brief_print(header); }

        Vector<T> diag(size_t k=0) { return Vector<T>(matrix.diag(k)); }
        //void diag(const Vector<T> & v) { matrix.diag(v.vector); }
        
        void addToEachRow(const Vector<T> & v) { matrix.each_row() += v.vector; }

        void insert_rows(size_t r, size_t n) { matrix.insert_rows(r,n); }
        void insert_cols(size_t c, size_t n) { matrix.insert_cols(c,n); }
        
        void remove_row(size_t n) { matrix.shed_row(n); }
        void remove_rows(size_t f, size_t l) { matrix.shed_rows(f,l); }

        void remove_col(size_t n) { matrix.shed_col(n); }
        void remove_cols(size_t f, size_t l) { matrix.shed_cols(f,l); }
        
        void swap_rows(size_t r1, size_t r2) { matrix.swap_rows(r1,r2); }
        void swap_cols(size_t c1, size_t c2) { matrix.swap_cols(c1,c2); }

        void swap(Matrix<T> & x) { matrix.swap(x.matrix); }

        
        Matrix<T> submat(size_t r1, size_t c1, size_t r2, size_t c2) { return Matrix<T>(matrix.submat(r1,c1,r2,c2)); }
        

        void save(const std::string& filename) { matrix.save(filename); }
        void load(const std::string& filename) { matrix.load(filename); }
        
        //void save_csv(const std::string& filename, const std::string& header = "") { matrix.save(arma::csv_name(filename,header)); }
        //void load_csv(const std::string& filename, const std::string& header = "") { matrix.load(arma::csv_name(filename,header)); }
    };

    template<typename T>    
    Vector<T> operator + (const T a, const Vector<T> & b) {
        Vector<T> r(b.size());
        r.vector = a + b.vector;
        return r;
    }
    template<typename T>    
    Vector<T> operator - (const T a, const Vector<T> & b) {
        Vector<T> r(b.size());
        r.vector = a - b.vector;
        return r;
    }
    template<typename T>    
    Vector<T> operator * (const T a, const Vector<T> & b) {
        Vector<T> r(b.size());
        r.vector = a * b.vector;
        return r;
    }
    template<typename T>    
    Vector<T> operator / (const T a, const Vector<T> & b) {
        Vector<T> r(b.size());
        r.vector = a * b.vector;
        return r;
    }

    template<typename T>    
    Matrix<T> operator + (const T a, const Matrix<T> & b) {
        Matrix<T> r(b.size());
        r.vector = a + b.matrix;
        return r;
    }
    template<typename T>    
    Matrix<T> operator - (const T a, const Matrix<T> & b) {
        Matrix<T> r(b.size());
        r.vector = a - b.matrix;
        return r;
    }
    template<typename T>    
    Matrix<T> operator * (const T a, const Matrix<T> & b) {
        Matrix<T> r(b.size());
        r.vector = a * b.matrix;
        return r;
    }
    
    template<typename T>    
    Vector<T> linspace(T start, T end, size_t N=100) { return arma::linspace(start,end,N); }

    template<typename T>    
    Vector<T> logspace(T start, T end, size_t N=100) { return arma::logspace(start,end,N); }

    template<typename T>    
    Vector<T> regspace(T start, T end, size_t N=100) { return arma::regspace(start,end,N); }

    
    template<typename T>
    Vector<T> randperm(size_t n) { return Vector<T>(arma::randperm(n)); }
    
    template<typename T>
    Matrix<T> randperm(size_t r, size_t c) { return Matrix<T>(arma::randperm(r,c)); }

    template<typename T>
    Matrix<T> eye(size_t r, size_t c) { return Matrix<T>(arma::eye(r,c)); }

    template<typename T>
    T zero(const T & x) {
        T r(x);
        r.zeros();
        return r;
    }
    template<typename T>
    T ones(const T & x) {
        T r(x);
        r.ones();
        return r;
    }
    template<typename T>
    T random(const T & x) {
        T r(x);
        r.random();
        return r;
    }
    template<typename T>
    T randu(const T & x) {
        T r(x);
        r.randu();
        return r;
    }
    

    template<typename T>
    Vector<T> random_vector(size_t n, T min=(T)0, T max = (T)0) {
        Vector<T> r(n);
        r.random();
        r = min + (max-min)*r;
        return r;
    }

    template<typename T>
    Matrix<T> toeplitz(const Vector<T> & a) { return Matrix<T>(arma::toeplitz(a.vector)); }

    template<typename T>
    Matrix<T> circ_toeplitz(const Vector<T> & a) { return Matrix<T>(arma::circ_toeplitz(a.vector)); }

    template<typename T>
    Matrix<T> toeplitz(const Vector<T> & a, const Vector<T> & b) { return Matrix<T>(arma::toeplitz(a.vector,b.vector)); }

    
    template<typename T>
    T accu(const Vector<T> & a) { return arma::accu(a.vector); }
    template<typename T>
    T accu(const RowVector<T> & a) { return arma::accu(a.vector); }
    template<typename T>
    T accu(const Matrix<T> & a) { return arma::accu(a.matrix); }

    template<typename T>
    T arg(const Vector<T> & a) { return arma::arg(a.vector); }
    template<typename T>
    T arg(const RowVector<T> & a) { return arma::arg(a.vector); }
    template<typename T>
    T arg(const Matrix<T> & a) { return arma::arg(a.matrix); }

    template<typename T>
    T sum(const Vector<T> & a) { return arma::sum(a.vector); }
    template<typename T>
    T sum(const RowVector<T> & a) { return arma::sum(a.vector); }
    template<typename T>
    T sum(const Matrix<T> & a, size_t dim=0) { return arma::sum(a.matrix,dim); }

    template<typename T>
    Matrix<T> affmul(const Matrix<T> & a, const Matrix<T> & b) { return Matrix<T>(arma::affmul(a.matrix,b.matrix)); }
    template<typename T>
    Matrix<T> affmul(const Matrix<T> & a, const Vector<T> & b) { return Matrix<T>(arma::affmul(a.matrix,b.vector)); }

    template<typename T>
    T cdot(const Vector<T> & a) { return arma::cdot(a.vector); }
    template<typename T>
    T cdot(const RowVector<T> & a) { return arma::cdot(a.vector); }

    enum ApproxType
    {
        APPROX_ABSDIFF=0,
        APPROX_RELDIFF,
        APPROX_BOTH
    };
    template<typename T>
    bool approx_equal(const Matrix<T> & a, const Matrix<T> & b, ApproxType type, T tol) {
        std::string foo;
        if(type == APPROX_ABSDIFF) foo = "absdiff";
        else if(type == APPROX_RELDIFF) foo = "reldiff";
        else foo = "both";
        return arma::approx_equal(a.matrix,b.matrix, foo, tol);
    }

    template<typename T>
    Matrix<T> as_scalar(T expr) { return Matrix<T>(arma::as_scalar(expr)); }

    template<typename T>
    Matrix<T> clamp(const Matrix<T> & a, T min, T max) { return Matrix<T>(arma::clamp(a.matrix,min,max)); }
    
    template<typename T>
    T cond(const Matrix<T> & a) { return arma::cond(a.matrix); }

    template<typename T>
    Vector<T> cumsum(const Vector<T> & v) { return Vector<T>(arma::cumsum(v.vector)); }
    template<typename T>
    RowVector<T> cumsum(const RowVector<T> & v) { return Vector<T>(arma::cumsum(v.vector)); }
    template<typename T>
    Matrix<T> cumsum(const Matrix<T> & m) { return Matrix<T>(arma::cumsum(m.matrix)); }

    template<typename T>
    Vector<T> cumprod(const Vector<T> & v) { return Vector<T>(arma::cumprod(v.vector)); }
    template<typename T>
    RowVector<T> cumprod(const RowVector<T> & v) { return Vector<T>(arma::cumprod(v.vector)); }
    template<typename T>
    Matrix<T> cumprod(const Matrix<T> & m) { return Matrix<T>(arma::cumprod(m.matrix)); }

    template<typename T>
    Matrix<T> det(const Matrix<T> & m) { return Matrix<T>(arma::det(m.matrix)); }

    template<typename T>
    Matrix<T> diagmat(const Vector<T> & v) { return Matrix<T>(arma::diagmat(v.vector)); }
    template<typename T>
    Matrix<T> diagmat(const Matrix<T> & m) { return Matrix<T>(arma::diagmat(m.matrix)); }
    
    template<typename T>
    Vector<T> diagvec(const Matrix<T> & a, size_t k = 0) { return Vector<T>(arma::diagvec(a.matrix)); }

    template<typename T>
    Vector<T> diff(const Vector<T> & v, size_t k=1) { return Vector<T>(arma::diff(v.vector,k)); }

    template<typename T>
    Matrix<T> diff(const Matrix<T> & m, size_t k=1, size_t dim = 0) { return Matrix<T>(arma::diff(m.matrix,k,dim)); }

    template<typename T>
    T dot(const Vector<T> & a, const Vector<T> & b) { return arma::dot(a.vector,b.vector); }
    template<typename T>
    T dot(const RowVector<T> & a, const RowVector<T> & b) { return arma::dot(a.vector,b.vector); }

    template<typename T>
    T norm_dot(const Vector<T> & a, const Vector<T> & b) { return arma::norm_dot(a.vector,b.vector); }
    template<typename T>
    T norm_dot(const RowVector<T> & a, const RowVector<T> & b) { return arma::norm_dot(a.vector,b.vector); }

    template<typename T>
    void inplace_trans(Matrix<T> & m) { arma::inplace_trans(m.matrix); }

    template<typename T>
    T min(const Vector<T> & x) { return arma::min(x.vector); }
    template<typename T>
    T min(const RowVector<T> & x) { return arma::min(x.vector); }
    template<typename T>
    T min(const Matrix<T> & x,size_t dim=0) { return arma::min(x.matrix,dim); }

    template<typename T>
    T max(const Vector<T> & x) { return arma::max(x.vector); }
    template<typename T>
    T max(const RowVector<T> & x) { return arma::max(x.vector); }
    template<typename T>
    T max(const Matrix<T> & x, size_t dim=0) { return arma::max(x.matrix,dim); }

    template<typename T>
    Vector<T> normalize(const Vector<T> & x) { return Vector<T>(arma::normalise(x.vector)); }
    template<typename T>
    RowVector<T> normalize(const RowVector<T> & x) { return Vector<T>(arma::normalise(x.vector)); }
    template<typename T>
    Matrix<T> normalize(const Matrix<T> & x) { return Matrix<T>(arma::normalise(x.matrix)); }

    template<typename T>
    T prod(const Vector<T> & x) { return arma::prod(x.vector); }
    template<typename T>
    T prod(const RowVector<T> & x) { return arma::prod(x.vector); }
    template<typename T>
    T prod(const Matrix<T> & x, size_t dim=0) { return arma::prod(x.matrix,dim); }

    template<typename T>
    Matrix<T> powmat(const Matrix<T> & x, T n) { return Matrix<T>(arma::powmat(x.matrix,n)); }

    template<typename T>
    size_t rank(const Matrix<T> & x) { return arma::rank(x); }
    template<typename T>
    size_t rank(const Matrix<T> & x, T tol) { return arma::rank(x,tol); }

    template<typename T>
    Matrix<T> repelem(const Matrix<T> & a, size_t r, size_t c) { return Matrix<T>(arma::repelem(a.matrix,r,c)); }

    template<typename T>
    Matrix<T> repmat(const Matrix<T> & a, size_t r, size_t c) { return Matrix<T>(arma::repmat(a.matrix,r,c)); }

    template<typename T>
    Matrix<T> reshape(const Matrix<T> & a, size_t r, size_t c) { return Matrix<T>(arma::reshape(a.matrix,r,c)); }

    template<typename T>
    Vector<T> reverse(const Vector<T> & a) { return Vector<T>(arma::reverse(a.vector)); }
    template<typename T>
    Matrix<T> reverse(const Matrix<T> & a, size_t dim=0) { return Matrix<T>(arma::reverse(a.matrix,dim)); }

    template<typename T>
    Vector<T> shift(const Vector<T> & a, size_t n) { return Vector<T>(arma::shift(a.vector,n)); }
    template<typename T>
    Matrix<T> shift(const Matrix<T> & a, size_t n, size_t dim = 0) { return Matrix<T>(arma::shift(a.matrix,n,dim)); }

    template<typename T>
    Vector<T> shuffle(const Vector<T> & a, size_t n) { return Vector<T>(arma::shuffle(a.vector,n)); }
    template<typename T>
    Matrix<T> shuffle(const Matrix<T> & a, size_t n, size_t dim = 0) { return Matrix<T>(arma::shuffle(a.matrix,n,dim)); }

    template<typename T>
    size_t size(const Vector<T> & v) { return arma::size(v.vector); }
    template<typename T>
    size_t size(const RowVector<T> & v) { return arma::size(v.vector); }
    template<typename T>
    size_t size(const Matrix<T> & m) { return arma::size(m.matrix); }

    template<typename T>
    void sort(const Vector<T> & v, int direction=1) { return arma::sort(v.vector, direction==1? "ascend":"descend"); }
    template<typename T>
    void sort(const RowVector<T> & v, int direction=1) { return arma::sort(v.vector, direction==1? "ascend":"descend"); }
    template<typename T>
    void sort(const Matrix<T> & m, int direction=1, int dim=0) { return arma::sort(m.matrix, direction==1? "ascend":"descend",dim); }

    template<typename T>
    Matrix<T> symmatu(const Matrix<T> & m) { return Matrix<T>(arma::symmatu(m.matrix)); }

    template<typename T>
    T trace(const Matrix<T> & m) { return arma::trace(m.matrix); }

    template<typename T>
    Matrix<T> trans(const Matrix<T> & m) { return Matrix<T>(arma::trans(m.mnatrix)); }

    template<typename T>
    Matrix<T> trimatu(const Matrix<T> & m, size_t k = 0) { return Matrix<T>(arma::symmatu(m.matrix,k)); }

    template<typename T>
    Vector<T> vectorize(const Matrix<T> & m) { return Matrix<T>(arma::vectorise(m.matrix)); }


    template<typename T>
    Vector<T> abs(const Vector<T> & a) { return Vector<T>(arma::abs(a.vector)); }
    template<typename T>
    RowVector<T> abs(const RowVector<T> & a) { return RowVector<T>(arma::abs(a.vector)); }
    template<typename T>
    Matrix<T> abs(const Matrix<T> & a) { return Matrix<T>(arma::abs(a.matrix)); }

    template<typename T>
    Vector<T> exp(const Vector<T> & a) { return Vector<T>(arma::exp(a.vector)); }
    template<typename T>
    RowVector<T> exp(const RowVector<T> & a) { return RowVector<T>(arma::exp(a.vector)); }
    template<typename T>
    Matrix<T> exp(const Matrix<T> & a) { return Matrix<T>(arma::expmat(a.matrix)); }

    template<typename T>
    Vector<T> exp2(const Vector<T> & a) { return Vector<T>(arma::exp2(a.vector)); }
    template<typename T>
    RowVector<T> exp2(const RowVector<T> & a) { return RowVector<T>(arma::exp2(a.vector)); }
    template<typename T>
    Matrix<T> exp2(const Matrix<T> & a) { return Matrix<T>(arma::exp2(a.matrix)); }

    template<typename T>
    Vector<T> exp10(const Vector<T> & a) { return Vector<T>(arma::exp10(a.vector)); }
    template<typename T>
    RowVector<T> exp10(const RowVector<T> & a) { return RowVector<T>(arma::exp10(a.vector)); }
    template<typename T>
    Matrix<T> exp10(const Matrix<T> & a) { return Matrix<T>(arma::exp10(a.matrix)); }

    template<typename T>
    Vector<T> expm1(const Vector<T> & a) { return Vector<T>(arma::expm1(a.vector)); }
    template<typename T>
    RowVector<T> expm1(const RowVector<T> & a) { return RowVector<T>(arma::expm1(a.vector)); }
    template<typename T>
    Matrix<T> expm1(const Matrix<T> & a) { return Matrix<T>(arma::expm1(a.matrix)); }

    template<typename T>
    Vector<T> trunc_exp(const Vector<T> & a) { return Vector<T>(arma::trunc_exp(a.vector)); }
    template<typename T>
    RowVector<T> trunc_exp(const RowVector<T> & a) { return RowVector<T>(arma::trunc_exp(a.vector)); }
    template<typename T>
    Matrix<T> trunc_exp(const Matrix<T> & a) { return Matrix<T>(arma::trunc_exp(a.matrix)); }

    template<typename T>
    Vector<T> log(const Vector<T> & a) { return Vector<T>(arma::log(a.vector)); }
    template<typename T>
    RowVector<T> log(const RowVector<T> & a) { return RowVector<T>(arma::log(a.vector)); }
    template<typename T>
    Matrix<T> log(const Matrix<T> & a) { return Matrix<T>(arma::logmat(a.matrix)); }

    template<typename T>
    Vector<T> log2(const Vector<T> & a) { return Vector<T>(arma::log2(a.vector)); }
    template<typename T>
    RowVector<T> log2(const RowVector<T> & a) { return RowVector<T>(arma::log2(a.vector)); }
    template<typename T>
    Matrix<T> log2(const Matrix<T> & a) { return Matrix<T>(arma::log2(a.matrix)); }

    template<typename T>
    Vector<T> log1p(const Vector<T> & a) { return Vector<T>(arma::log1p(a.vector)); }
    template<typename T>
    RowVector<T> log1p(const RowVector<T> & a) { return RowVector<T>(arma::log1p(a.vector)); }
    template<typename T>
    Matrix<T> log1p(const Matrix<T> & a) { return Matrix<T>(arma::log1p(a.matrix)); }

    template<typename T>
    Vector<T> log10(const Vector<T> & a) { return Vector<T>(arma::log10(a.vector)); }
    template<typename T>
    RowVector<T> log10(const RowVector<T> & a) { return RowVector<T>(arma::log10(a.vector)); }
    template<typename T>
    Matrix<T> log10(const Matrix<T> & a) { return Matrix<T>(arma::log10(a.matrix)); }

    template<typename T>
    Vector<T> trunc_log(const Vector<T> & a) { return Vector<T>(arma::trunc_log(a.vector)); }
    template<typename T>
    RowVector<T> trunc_log(const RowVector<T> & a) { return RowVector<T>(arma::trunc_log(a.vector)); }
    template<typename T>
    Matrix<T> trunc_log(const Matrix<T> & a) { return Matrix<T>(arma::trunc_log(a.matrix)); }

    template<typename T>
    Vector<T> pow(const Vector<T> & a, T x) { return Vector<T>(arma::pow(a.vector,x)); }
    template<typename T>
    RowVector<T> pow(const RowVector<T> & a, T x) { return RowVector<T>(arma::pow(a.vector,x)); }
    template<typename T>
    Matrix<T> pow(const Matrix<T> & a, T x) { return Matrix<T>(arma::powmat(a.matrix,x)); }

    template<typename T>
    Vector<T> square(const Vector<T> & a) { return Vector<T>(arma::square(a.vector)); }
    template<typename T>
    RowVector<T> square(const RowVector<T> & a) { return RowVector<T>(arma::square(a.vector)); }
    template<typename T>
    Matrix<T> square(const Matrix<T> & a) { return Matrix<T>(arma::square(a.matrix)); }

    template<typename T>
    Vector<T> sqrt(const Vector<T> & a) { return Vector<T>(arma::sqrt(a.vector)); }
    template<typename T>
    RowVector<T> sqrt(const RowVector<T> & a) { return RowVector<T>(arma::sqrt(a.vector)); }
    template<typename T>
    Matrix<T> sqrt(const Matrix<T> & a) { return Matrix<T>(arma::sqrtmat(a.matrix)); }

    template<typename T>
    Vector<T> floor(const Vector<T> & a) { return Vector<T>(arma::floor(a.vector)); }
    template<typename T>
    RowVector<T> floor(const RowVector<T> & a) { return RowVector<T>(arma::floor(a.vector)); }
    template<typename T>
    Matrix<T> floor(const Matrix<T> & a) { return Matrix<T>(arma::floor(a.matrix)); }

    template<typename T>
    Vector<T> ceil(const Vector<T> & a) { return Vector<T>(arma::ceil(a.vector)); }
    template<typename T>
    RowVector<T> ceil(const RowVector<T> & a) { return RowVector<T>(arma::ceil(a.vector)); }
    template<typename T>
    Matrix<T> ceil(const Matrix<T> & a) { return Matrix<T>(arma::ceil(a.matrix)); }

    template<typename T>
    Vector<T> round(const Vector<T> & a) { return Vector<T>(arma::round(a.vector)); }
    template<typename T>
    RowVector<T> round(const RowVector<T> & a) { return RowVector<T>(arma::round(a.vector)); }
    template<typename T>
    Matrix<T> round(const Matrix<T> & a) { return Matrix<T>(arma::round(a.matrix)); }

    template<typename T>
    Vector<T> trunc(const Vector<T> & a) { return Vector<T>(arma::trunc(a.vector)); }
    template<typename T>
    RowVector<T> trunc(const RowVector<T> & a) { return RowVector<T>(arma::trunc(a.vector)); }
    template<typename T>
    Matrix<T> trunc(const Matrix<T> & a) { return Matrix<T>(arma::trunc(a.matrix)); }

    template<typename T>
    Vector<T> erf(const Vector<T> & a) { return Vector<T>(arma::erf(a.vector)); }
    template<typename T>
    RowVector<T> erf(const RowVector<T> & a) { return RowVector<T>(arma::erf(a.vector)); }
    template<typename T>
    Matrix<T> erf(const Matrix<T> & a) { return Matrix<T>(arma::erf(a.matrix)); }

    template<typename T>
    Vector<T> erfc(const Vector<T> & a) { return Vector<T>(arma::erfc(a.vector)); }
    template<typename T>
    RowVector<T> erfc(const RowVector<T> & a) { return RowVector<T>(arma::erfc(a.vector)); }
    template<typename T>
    Matrix<T> erfc(const Matrix<T> & a) { return Matrix<T>(arma::erfc(a.matrix)); }

    template<typename T>
    Vector<T> tgamma(const Vector<T> & a) { return Vector<T>(arma::tgamma(a.vector)); }
    template<typename T>
    RowVector<T> tgamma(const RowVector<T> & a) { return RowVector<T>(arma::tgamma(a.vector)); }
    template<typename T>
    Matrix<T> tgamma(const Matrix<T> & a) { return Matrix<T>(arma::tgamma(a.matrix)); }

    template<typename T>
    Vector<T> lgamma(const Vector<T> & a) { return Vector<T>(arma::lgamma(a.vector)); }
    template<typename T>
    RowVector<T> lgamma(const RowVector<T> & a) { return RowVector<T>(arma::lgamma(a.vector)); }
    template<typename T>
    Matrix<T> lgamma(const Matrix<T> & a) { return Matrix<T>(arma::lgamma(a.matrix)); }

    template<typename T>
    Vector<T> sign(const Vector<T> & a) { return Vector<T>(arma::sign(a.vector)); }
    template<typename T>
    RowVector<T> sign(const RowVector<T> & a) { return RowVector<T>(arma::sign(a.vector)); }
    template<typename T>
    Matrix<T> sign(const Matrix<T> & a) { return Matrix<T>(arma::sign(a.matrix)); }

    template<typename T>
    Vector<T> cos(const Vector<T> & a) { return Vector<T>(arma::cos(a.vector)); }
    template<typename T>
    RowVector<T> cos(const RowVector<T> & a) { return RowVector<T>(arma::cos(a.vector)); }
    template<typename T>
    Matrix<T> cos(const Matrix<T> & a) { return Matrix<T>(arma::cos(a.matrix)); }

    template<typename T>
    Vector<T> sin(const Vector<T> & a) { return Vector<T>(arma::sin(a.vector)); }
    template<typename T>
    RowVector<T> sin(const RowVector<T> & a) { return RowVector<T>(arma::sin(a.vector)); }
    template<typename T>
    Matrix<T> sin(const Matrix<T> & a) { return Matrix<T>(arma::sin(a.matrix)); }

    template<typename T>
    Vector<T> tan(const Vector<T> & a) { return Vector<T>(arma::tan(a.vector)); }
    template<typename T>
    RowVector<T> tan(const RowVector<T> & a) { return RowVector<T>(arma::tan(a.vector)); }
    template<typename T>
    Matrix<T> tan(const Matrix<T> & a) { return Matrix<T>(arma::tan(a.matrix)); }

    template<typename T>
    Vector<T> acos(const Vector<T> & a) { return Vector<T>(arma::acos(a.vector)); }
    template<typename T>
    RowVector<T> acos(const RowVector<T> & a) { return RowVector<T>(arma::acos(a.vector)); }
    template<typename T>
    Matrix<T> acos(const Matrix<T> & a) { return Matrix<T>(arma::acos(a.matrix)); }

    template<typename T>
    Vector<T> asin(const Vector<T> & a) { return Vector<T>(arma::asin(a.vector)); }
    template<typename T>
    RowVector<T> asin(const RowVector<T> & a) { return RowVector<T>(arma::asin(a.vector)); }
    template<typename T>
    Matrix<T> asin(const Matrix<T> & a) { return Matrix<T>(arma::asin(a.matrix)); }

    template<typename T>
    Vector<T> atan(const Vector<T> & a) { return Vector<T>(arma::atan(a.vector)); }
    template<typename T>
    RowVector<T> atan(const RowVector<T> & a) { return RowVector<T>(arma::atan(a.vector)); }
    template<typename T>
    Matrix<T> atan(const Matrix<T> & a) { return Matrix<T>(arma::atan(a.matrix)); }


    template<typename T>
    Vector<T> cosh(const Vector<T> & a) { return Vector<T>(arma::cosh(a.vector)); }
    template<typename T>
    RowVector<T> cosh(const RowVector<T> & a) { return RowVector<T>(arma::cosh(a.vector)); }
    template<typename T>
    Matrix<T> cosh(const Matrix<T> & a) { return Matrix<T>(arma::cosh(a.matrix)); }

    template<typename T>
    Vector<T> sinh(const Vector<T> & a) { return Vector<T>(arma::sinh(a.vector)); }
    template<typename T>
    RowVector<T> sinh(const RowVector<T> & a) { return RowVector<T>(arma::sinh(a.vector)); }
    template<typename T>
    Matrix<T> sinh(const Matrix<T> & a) { return Matrix<T>(arma::sinh(a.matrix)); }

    template<typename T>
    Vector<T> tanh(const Vector<T> & a) { return Vector<T>(arma::tanh(a.vector)); }
    template<typename T>
    RowVector<T> tanh(const RowVector<T> & a) { return RowVector<T>(arma::tanh(a.vector)); }
    template<typename T>
    Matrix<T> tanh(const Matrix<T> & a) { return Matrix<T>(arma::tanh(a.matrix)); }

    template<typename T>
    Vector<T> acosh(const Vector<T> & a) { return Vector<T>(arma::acosh(a.vector)); }
    template<typename T>
    RowVector<T> acosh(const RowVector<T> & a) { return RowVector<T>(arma::acosh(a.vector)); }
    template<typename T>
    Matrix<T> acosh(const Matrix<T> & a) { return Matrix<T>(arma::acosh(a.matrix)); }

    template<typename T>
    Vector<T> asinh(const Vector<T> & a) { return Vector<T>(arma::asinh(a.vector)); }
    template<typename T>
    RowVector<T> asinh(const RowVector<T> & a) { return RowVector<T>(arma::asinh(a.vector)); }
    template<typename T>
    Matrix<T> asinh(const Matrix<T> & a) { return Matrix<T>(arma::asinh(a.matrix)); }

    template<typename T>
    Vector<T> atanh(const Vector<T> & a) { return Vector<T>(arma::atanh(a.vector)); }
    template<typename T>
    RowVector<T> atanh(const RowVector<T> & a) { return RowVector<T>(arma::atanh(a.vector)); }
    template<typename T>
    Matrix<T> atanh(const Matrix<T> & a) { return Matrix<T>(arma::atanh(a.matrix)); }

    template<typename T>
    Vector<T> atan2(const Vector<T> & a, const Vector<T> & b) { return Vector<T>(arma::atan2(a.vector,b.vector)); }
    template<typename T>
    RowVector<T> atan2(const RowVector<T> & a, const RowVector<T> & b) { return RowVector<T>(arma::atan2(a.vector,b.vector)); }
    template<typename T>
    Matrix<T> atan2(const Matrix<T> & a, const Matrix<T> & b) { return Matrix<T>(arma::atan2(a.matrix,b.matrix)); }

    template<typename T>
    Vector<T> hypot(const Vector<T> & a, const Vector<T> & b) { return Vector<T>(arma::hypot(a.vector,b.vector)); }
    template<typename T>
    RowVector<T> hypot(const RowVector<T> & a, const RowVector<T> & b) { return RowVector<T>(arma::hypot(a.vector,b.vector)); }
    template<typename T>
    Matrix<T> hypot(const Matrix<T> & a, const Matrix<T> & b) { return Matrix<T>(arma::hypot(a.matrix,b.matrix)); }

    template<typename T>
    Vector<T> sinc(const Vector<T> & a) { return Vector<T>(arma::sinc(a.vector)); }
    template<typename T>
    RowVector<T> sinc(const RowVector<T> & a) { return RowVector<T>(arma::sinc(a.vector)); }
    template<typename T>
    Matrix<T> sinc(const Matrix<T> & a) { return Matrix<T>(arma::sinc(a.matrix)); }

    /*
    template<typename T>
    Vector<T> col(const Matrix<T> & a) { return Vector<T>(a.matrix)); }
    
    template<typename T>
    RowVector<T> row(const Matrix<T> & a) { return Vector<T>(arma::row(a.matrix)); }
    */
   
    template<typename T>
    void addToEachRow(Matrix<T> & m, const RowVector<T> & v) {
        for(size_t i = 0; i < m.rows(); i++)
            m.matrix.row(i) += v.vector;
    }


    template<typename T>
    void addToEachCol(Matrix<T> & m, const Vector<T> & v) {
        for(size_t i = 0; i < m.cols(); i++)
            m.matrix.col(i) += v.vector;
    }

// complex

    template<typename T>
    void set_real(Vector<std::complex<T>> & v, const T & value) { v.vector.set_real(value); }
    template<typename T>
    void set_real(RowVector<std::complex<T>> & v, const T & value) { v.vector.set_real(value); }
    template<typename T>
    void set_real(Matrix<std::complex<T>> & m, const T & value) { m.matrix.set_real(value); }


    template<typename T>
    void set_imag(Vector<std::complex<T>> & v, const T & value) { v.vector.set_imag(value); }
    template<typename T>
    void set_imag(RowVector<std::complex<T>> & v, const T & value) { v.vector.set_imag(value); }
    template<typename T>
    void set_imag(Matrix<std::complex<T>> & m, const T & value) { m.matrix.set_imag(value); }

    template<typename T>
    Matrix<std::complex<T>> st(Matrix<std::complex<T>> & m) { return Matrix<T>(m.matrix.st()); }

    template<typename T>
    bool is_hermitian(const Matrix<std::complex<T>>& m) { return m.matrix.is_hermitian(); }

    template<typename T>
    Matrix<T> real(Matrix<std::complex<T>> & m) { return Matrix<T>(arma::real(m.matrix)); }

    template<typename T>
    Matrix<T> imag(Matrix<std::complex<T>> & m) { return Matrix<T>(arma::imag(m.matrix)); }
  
    template<typename T>
    std::complex<double> log_det(Matrix<T> & m) { return std::complex<double>(arma::log_det(m.matrix)); }

    template<typename T>
    bool log_det(Matrix<T> & m, double val, double sign) { return arma::log_det(val,sign,m.matrix); }

    template<typename T>
    double log_det_sympd(Matrix<T> & m) { return arma::log_det_sympd(m.matrix); }
    
    template<typename T>
    bool   log_det_sympd(Matrix<T> & m, double result) { return arma::log_det(result,m.matrix); }
  
    template<typename T>
    Matrix<std::complex<T>> logmat(Matrix<T> & m) { return Matrix<T>(arma::logmat(m.matrix)); }

    template<typename T>
    Matrix<T> logmat_sympd(Matrix<T> & m) { return Matrix<T>(arma::log_det_sympd(m.matrix)); }

    template<typename T>
    Matrix<std::complex<T>> sqrtmat(Matrix<T> & m) { return Matrix<std::complex<T>>(arma::sqrtmat(m.matrix)); }
    
    template<typename T>
    Matrix<T> sqrtmat_sympd(Matrix<T> & m) { return Matrix<T>(arma::sqrtmat_sympd(m.matrix)); }


// solvers
    template<typename T>
    Matrix<T> chol(Matrix<T> & m, std::string & layout = "upper") { return Matrix<T>(arma::chol(m.matrix,layout)); }

    template<typename T>
    Vector<T> eig_sym(Matrix<T> & m) { return Vector<T>(arma::chol(m.matrix)); }
    
    template<typename T>
    Vector<std::complex<T>> eig_gen(Matrix<T> & m) { return Vector<std::complex<T>>(arma::eig_gen(m.matrix)); }

    template<typename T>
    Vector<std::complex<T>> eig_pair(Matrix<T> & a, Matrix<T> & b) { return Vector<std::complex<T>>(arma::eig_pair(a.matrix,b.matrix)); }

    template<typename T>
    Matrix<T> hess(Matrix<T> & x) { return Matrix<T>(arma::hess(x.matrix)); }

    template<typename T>
    Matrix<T> inv(Matrix<T> & x) { return Matrix<T>(arma::inv(x.matrix)); }

    template<typename T>
    Matrix<T> inv_sympd(Matrix<T> & x) { return Matrix<T>(arma::inv_sympd(x.matrix)); }

    template<typename T>
    void lu(Matrix<T> & L, Matrix<T> & U, Matrix<T>& P, Matrix<T> & X) { arma::lu(L.matrix,U.matrix,P.matrix,X.matrix); }

    template<typename T>
    void lu(Matrix<T> & L, Matrix<T> & U, Matrix<T> & X) { arma::lu(L.matrix,U.matrix,X.matrix); }

    // I don't really know what this does, I don't really use the solvers much.
    template<typename T>
    Matrix<T> null(Matrix<T> & A) { return Matrix<T>(arma::null(A.matrix)); }

    template<typename T>
    Matrix<T> orth(Matrix<T> & A) { return Matrix<T>(arma::orth(A.matrix)); }

    template<typename T>
    Matrix<T> pinv(Matrix<T> & A) { return Matrix<T>(arma::pinv(A.matrix)); }

    template<typename T>
    Matrix<T> qr(Matrix<T> & Q, Matrix<T> & R, Matrix<T> & X) { return Matrix<T>(arma::qr(Q.matrix, R.matrix, X.matrix)); }

    template<typename T>
    Matrix<T> qr_econ(Matrix<T> & Q, Matrix<T> & R, Matrix<T> & X) { return Matrix<T>(arma::qr_econ(Q.matrix, R.matrix, X.matrix)); }

    template<typename T>
    Matrix<T> schur(Matrix<T> & A) { return Matrix<T>(arma::schur(A.matrix)); }

    template<typename T>
    Vector<T> solve(Matrix<T> & A, Vector<T> & b) { return Matrix<T>(arma::pinv(A.matrix,b.vector)); }

    template<typename T>
    Vector<T> svd(Matrix<T> & X) { return Matrix<T>(arma::svd(X.matrix)); }

// DSP
    template<typename T>
    Vector<T> conv(Vector<T> & A, Vector<T> & b) { return Vector<T>(arma::conv(A.vector, b.vector)); }

    template<typename T>
    Matrix<T> conv2(Matrix<T> & A, Matrix<T> & b) { return Matrix<T>(arma::conv2(A.matrix, b.matrix)); }

    template<typename T>
    Vector<std::complex<T>> fft(Vector<T> & A) { return Vector<std::complex<T>>(arma::fft(A.vector)); }

    template<typename T>
    Vector<std::complex<T>> fft(Vector<std::complex<T>> & A) { return Vector<std::complex<T>>(arma::fft(A.vector)); }
    
    template<typename T>
    Vector<std::complex<T>> ifft(Vector<std::complex<T>> & A) { return Vector<std::complex<T>>(arma::ifft(A.vector)); }

    template<typename T>
    Matrix<std::complex<T>> fft2(Matrix<T> & A) { return Matrix<T>(arma::fft(A.matrix)); }

    template<typename T>
    Matrix<std::complex<T>> fft2(Matrix<std::complex<T>> & A) { return Matrix<std::complex<T>>(arma::fft(A.matrix)); }

    template<typename T>
    Matrix<std::complex<T>> ifft2(Matrix<std::complex<T>> & A) { return Matrix<std::complex<T>>(arma::ifft2(A.matrix)); }

    template<typename T>
    void interp1(Vector<T> & X, Vector<T> & Y, Vector<T> & XI, Vector<T> & YI, std::string& method = "linear", double extrap = arma::datum::nan ) { (arma::interp1(X.vector, Y.vector, XI.vector, YI.vector)); }

    template<typename T>
    void interp2(Matrix<T> & X, Matrix<T> & Y, Matrix<T> & XI, Matrix<T> & YI, std::string& method = "linear", double extrap = arma::datum::nan ) { (arma::interp2(X.matrix, Y.matrix, XI.matrix, YI.matrix)); }


    template<typename T>
    Vector<T> polyfit(Vector<T> & X, Vector<T> & Y, Vector<T> & N) { return Vector<T>(arma::polyfit(X.vector, Y.vector, N.vector)); }

    template<typename T>
    Vector<T> polyval(Vector<T> & P, Vector<T> & X) { return Vector<T>(arma::polyval(P.vector, X.vector)); }

// Statistics

    template<typename T>
    double mean(Vector<T> & V) { return Vector<T>(arma::mean(V.vector)); }
    template<typename T>
    double median(Vector<T> & V) { return Vector<T>(arma::median(V.vector)); }
    template<typename T>
    double stddev(Vector<T> & V) { return Vector<T>(arma::stddev(V.vector)); }
    template<typename T>
    double var(Vector<T> & V) { return Vector<T>(arma::var(V.vector)); }
    template<typename T>
    double range(Vector<T> & V) { return Vector<T>(arma::range(V.vector)); }

    template<typename T>
    double mean(Matrix<T> & V) { return Matrix<T>(arma::mean(V.matrix)); }
    template<typename T>
    double median(Matrix<T> & V) { return Matrix<T>(arma::median(V.matrix)); }
    template<typename T>
    double stddev(Matrix<T> & V) { return Matrix<T>(arma::stddev(V.matrix)); }
    template<typename T>
    double var(Matrix<T> & V) { return Matrix<T>(arma::var(V.matrix)); }
    template<typename T>
    double range(Matrix<T> & V) { return Matrix<T>(arma::range(V.matrix)); }

    template<typename T>
    Vector<T> cov(Vector<T> & X, Vector<T>& Y) { return Vector<T>(arma::cov(X.vector,Y.vector)); }
    template<typename T>
    Vector<T> cov(Vector<T> & X) { return Vector<T>(arma::cov(X.vector)); }
    template<typename T>
    Matrix<T> cov(Matrix<T> & X, Matrix<T>& Y) { return Matrix<T>(arma::cov(X.matrix,Y.matrix)); }
    template<typename T>
    Matrix<T> cov(Matrix<T> & X) { return Matrix<T>(arma::cov(X.Matrix)); }

    template<typename T>
    Vector<T> cor(Vector<T> & X, Vector<T>& Y) { return Vector<T>(arma::cor(X.vector,Y.vector)); }
    template<typename T>
    Vector<T> cor(Vector<T> & X) { return Vector<T>(arma::cor(X.vector)); }
    template<typename T>
    Matrix<T> cor(Matrix<T> & X, Matrix<T>& Y) { return Matrix<T>(arma::cor(X.matrix,Y.matrix)); }
    template<typename T>
    Matrix<T> cor(Matrix<T> & X) { return Matrix<T>(arma::cor(X.matrix)); }

    
    template<typename T>
    Vector<unsigned> hist(Vector<T> & X) { return Vector<unsigned>(arma::hist(X.vector)); }    
    template<typename T>
    Vector<unsigned> hist(Vector<T> & X, size_t n_bins) { return Vector<unsigned>(arma::cov(X.vector,n_bins)); }
    template<typename T>
    Matrix<unsigned> hist(Matrix<T> & X) { return Matrix<unsigned>(arma::hist(X.Matrix)); }
    template<typename T>
    Matrix<unsigned> hist(Matrix<T> & X, size_t n_bins) { return Matrix<unsigned>(arma::cov(X.matrix,n_bins)); }

    template<typename T>
    Vector<unsigned> histc(Vector<T> & X) { return Vector<unsigned>(arma::hist(X.vector)); }
    template<typename T>
    Vector<unsigned> histc(Vector<T> & X, size_t n_bins) { return Vector<unsigned>(arma::cov(X.vector,n_bins)); }
    template<typename T>
    Matrix<unsigned> histc(Matrix<T> & X) { return Matrix<unsigned>(arma::hist(X.Matrix)); }
    template<typename T>
    Matrix<unsigned> histc(Matrix<T> & X, size_t n_bins) { return Matrix<unsigned>(arma::cov(X.matrix,n_bins)); }
    
    template<typename T>
    Vector<T> quantile(Vector<T> & V, Vector<T> & P) { return Vector<T>(arma::quantile(V.vector,P.vector));}
    template<typename T>
    Vector<T> quantile(Matrix<T> & V, Vector<T> & P) { return Matrix<T>(arma::quantile(V.matrix,P.matrix));}

    template<typename T>
    Matrix<T> princomp(Matrix<T> & x) { return Matrix<T>(arma::princomp(x.matrix)); }

    template<typename T>
    void princomp(Matrix<T> & MC, Matrix<T> & X) { (arma::princomp(MC.matrix, X.matrix));}

    template<typename T>
    void princomp(Matrix<T> & MC, Matrix<T> & mat_score, Matrix<T> & X) { (arma::princomp(MC.matrix,mat_score.matrix, X.matrix));}

    template<typename T>
    void princomp(Matrix<T> & MC, Matrix<T> & mat_score, Vector<T> & latent, Matrix<T> & X) { (arma::princomp(MC.matrix,mat_score.matrix, latent.vector, X.matrix)); }

    template<typename T>
    void princomp(Matrix<T> & MC, Matrix<T> & mat_score, Vector<T> & latent, Vector<T> & tsqr, Matrix<T> & X) { (arma::princomp(MC.matrix,mat_score.matrix, latent.vector, tsqr.vector, X.matrix)); }

    template<typename T>
    void princomp(Matrix<std::complex<T>> & MC, Matrix<std::complex<T>> & mat_score, Vector<std::complex<T>> & latent, Vector<std::complex<T>> & tsqr, Matrix<std::complex<T>> & X) { (arma::princomp(MC.matrix,mat_score.matrix, latent.vector, tsqr.vector, X.matrix)); }


};