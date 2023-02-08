#pragma once

namespace Octopus
{
    struct OctopusMatrixXf;

    struct MatrixViewXf
    {
        OctopusMatrixXf * matrix;
        size_t row;

        MatrixViewXf(OctopusMatrixXf * m, size_t r) : matrix(m),row(r) {}

        float& operator[](size_t i);
        float  operator[](size_t i) const;

        float __getitem__(size_t i);
        void __setitem__(size_t i, const float & v);

    };
        

    struct OctopusMatrixXf : public FloatMatrix
    {
        OctopusMatrixXf() = default;
        OctopusMatrixXf(const FloatMatrix &v) : FloatMatrix(v) {}
        OctopusMatrixXf(size_t i,size_t j) : FloatMatrix(i,j) {}

        using FloatMatrix::operator =;
        using FloatMatrix::operator ();
        using FloatMatrix::insert;
        using FloatMatrix::append;
        using FloatMatrix::fill;
        using FloatMatrix::extract;
        using FloatMatrix::extract_n;
        using FloatMatrix::transpose;
        using FloatMatrix::rows;
        using FloatMatrix::cols;
        using FloatMatrix::row;
        //using FloatMatrix::col;
        
       

        OctopusMatrixXf operator + (const OctopusMatrixXf & b) {
            OctopusMatrixXf r;
            ValueList l;
            l(0) = "plus";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).float_matrix_value();
            return r;
        }
        OctopusMatrixXf operator - (const OctopusMatrixXf & b) {
            OctopusMatrixXf r;
            ValueList l;
            l(0) = "minus";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).float_matrix_value();
            return r;
        }
        OctopusMatrixXf operator * (const OctopusMatrixXf & b) {
            OctopusMatrixXf r;
            ValueList l;
            l(0) = "times";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).float_matrix_value();
            return r;
        }
        OctopusMatrixXf operator / (const OctopusMatrixXf & b) {
            OctopusMatrixXf r;
            ValueList l;
            l(0) = "divide";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).float_matrix_value();
            return r;
        }
        OctopusMatrixXf operator + (const float b) {
            OctopusMatrixXf r;
            ValueList l;
            l(0) = "plus";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).float_matrix_value();
            return r;
        }
        OctopusMatrixXf operator - (const float b) {
            OctopusMatrixXf r;
            ValueList l;
            l(0) = "minus";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).float_matrix_value();
            return r;
        }
        OctopusMatrixXf operator * (const float b) {
            OctopusMatrixXf r;
            ValueList l;
            l(0) = "times";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).float_matrix_value();
            return r;
        }
        OctopusMatrixXf operator / (const float b) {
            OctopusMatrixXf r;
            ValueList l;
            l(0) = "divide";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).float_matrix_value();
            return r;
        }   
        void print()
        {
            ValueList l;
            l(0) = *this;
            octave::feval("display",l,0);
        }
        OctopusMatrixXf& operator += (const OctopusMatrixXf & b) {
            *this = *this + b;
            return *this;
        }
        OctopusMatrixXf& operator -= (const OctopusMatrixXf & b) {
            *this = *this - b;
            return *this;
        }
        OctopusMatrixXf& operator *= (const OctopusMatrixXf & b) {
            *this = *this * b;
            return *this;
        }
        OctopusMatrixXf& operator /= (const OctopusMatrixXf & b) {
            *this = *this / b;
            return *this;
        }
        OctopusMatrixXf& operator += (const float b) {
            *this = *this + b;
            return *this;
        }
        OctopusMatrixXf& operator -= (const float b) {
            *this = *this - b;
            return *this;
        }
        OctopusMatrixXf& operator *= (const float b) {
            *this = *this * b;
            return *this;
        }
        OctopusMatrixXf& operator /= (const float b) {
            *this = *this / b;
            return *this;
        }
        

        OctopusMatrixXf addToEachRow(OctopusRowVectorXf & v) {
            OctopusMatrixXf r(*this);
            
            for(size_t i = 0; i < rows(); i++)
            {
                for(size_t j = 0; j < cols(); j++)
                {
                    r(i,j) += v(j);
                }
            }
            return r;
        }
        OctopusMatrixXf addToEachRow(OctopusMatrixXf & v) {
            OctopusMatrixXf r(*this);
            for(size_t i = 0; i < rows(); i++)
            {
                for(size_t j = 0; j < cols(); j++)
                {
                    r(i,j) += v(0,j);
                }
            }
            return r;
        }
        OctopusMatrixXf eval() {
            OctopusMatrixXf r(*this);
            return r;
        }
        void printrowscols() const {
            std::cout << "rows=" << rows() << " cols=" << cols() << std::endl;
        }
        
        OctopusMatrixXf matmul(const OctopusMatrixXf & b)
        {            
            return *this * b;
        }
        OctopusMatrixXf hadamard(const OctopusMatrixXf & b)
        {            
            OctopusMatrixXf r;
            ValueList l;
            l(0) = "times";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).float_matrix_value();
            return r;
        }            

        #ifdef SWIG
        %extend
        {
            MatrixViewXf __getitem__(size_t i) { return MatrixViewXf($self,i-1); }            

            size_t rows() { return $self->rows(); }
            size_t cols() { return $self->cols(); }

            void fill(float v) { $self->fill(v); }
            
            //float min() { return $self->min(); }
            //float max() { return $self->max(); }
            
            const float* data() { return $self->data(); }

            void copy(size_t m, size_t n, double * p) {
                $self->resize(m,n);
                for(size_t i = 0; i < $self->rows(); i++)
                for(size_t j = 0; j < $self->cols(); j++)
                    (*$self)(i,j) = p[i*n + j];
            }

            void copy(size_t m, size_t n, float * p) {
                $self->resize(m,n);
                for(size_t i = 0; i < $self->rows(); i++)
                for(size_t j = 0; j < $self->cols(); j++)
                    (*$self)(i,j) = p[i*n + j];
            }
            OctopusMatrixXf __add__ (const OctopusMatrixXf & b) {
                assert($self->rows() == b.rows() && $self->cols() == b.cols());
                OctopusMatrixXf r($self->size());
                for(size_t i = 0; i < $self->rows(); i++)
                for(size_t j = 0; j < $self->cols(); j++)
                    r(i,j) = (*$self)(i,j) + b(i,j);
                return r;
            }
            OctopusMatrixXf __sub__ (const OctopusMatrixXf & b) {
                assert($self->rows() == b.rows() && $self->cols() == b.cols());
                OctopusMatrixXf r($self->size());
                for(size_t i = 0; i < $self->rows(); i++)
                for(size_t j = 0; j < $self->cols(); j++)
                    r(i,j) = (*$self)(i,j) - b(i,j);                    
                return r;
            }
            OctopusMatrixXf __mul__ (const OctopusMatrixXf & b) {
                return $self->matmul((*$self),b);
            }
            OctopusMatrixXf __div__ (const OctopusMatrixXf & b) {
                assert($self->rows() == b.rows() && $self->cols() == b.cols());
                OctopusMatrixXf r($self->size());
                for(size_t i = 0; i < $self->rows(); i++)
                for(size_t j = 0; j < $self->cols(); j++)
                    r(i,j) = (*$self)(i,j) / b(i,j);
                return r;
            } 

            OctopusMatrixXf __add__ (const float b) {
                OctopusMatrixXf r($self->size());
                for(size_t i = 0; i < $self->rows(); i++)
                for(size_t j = 0; j < $self->cols(); j++)
                    r(i,j) = (*$self)(i,j) + b;
                return r;
            }
            OctopusMatrixXf __sub__ (const float b) {
                OctopusMatrixXf r($self->size());
                for(size_t i = 0; i < $self->rows(); i++)
                for(size_t j = 0; j < $self->cols(); j++)
                    r(i,j) = (*$self)(i,j) - b;                    
                return r;
            }
            OctopusMatrixXf __mul__ (const float b) {
                OctopusMatrixXf r($self->size());
                for(size_t i = 0; i < $self->rows(); i++)
                for(size_t j = 0; j < $self->cols(); j++)
                    r(i,j) = (*$self)(i,j) * b;
                return r;
            }
            OctopusMatrixXf __div__ (const float b) {
                OctopusMatrixXf r($self->size());
                for(size_t i = 0; i < $self->rows(); i++)
                for(size_t j = 0; j < $self->cols(); j++)
                    r(i,j) = (*$self)(i,j) / b;
                return r;
            } 
        }
        #endif
    };   
    inline float& MatrixViewXf::operator[](size_t i) {
            return (*matrix)(row,i);
    }
    inline float MatrixViewXf::operator[](size_t i) const {
        return (*matrix)(row,i);
    }
    inline float MatrixViewXf::__getitem__(size_t i) { return (*this)[i-1]; }
    inline void MatrixViewXf::__setitem__(size_t i, const float & v) { (*this)[i-1] = v; }
}