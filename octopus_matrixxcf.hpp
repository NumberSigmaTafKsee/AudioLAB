#pragma once

namespace Octopus
{
    struct OctopusMatrixXcf;
    struct MatrixViewXcf
    {
        OctopusMatrixXcf * matrix;
        size_t row;

        MatrixViewXcf(OctopusMatrixXcf * m, size_t r) : matrix(m),row(r) {}

        std::complex<float>& operator[](size_t i);
        std::complex<float>  operator[](size_t i) const;
    };
    struct OctopusMatrixXcf : public FloatComplexMatrix
    {
        OctopusMatrixXcf() = default;
        OctopusMatrixXcf(const FloatComplexMatrix &v) : FloatComplexMatrix(v) {}
        OctopusMatrixXcf(size_t i,size_t j) : FloatComplexMatrix(i,j) {}

        using FloatComplexMatrix::operator =;
        using FloatComplexMatrix::operator ();
        using FloatComplexMatrix::insert;
        using FloatComplexMatrix::append;
        using FloatComplexMatrix::fill;
        using FloatComplexMatrix::extract;
        using FloatComplexMatrix::extract_n;
        using FloatComplexMatrix::transpose;
        using FloatComplexMatrix::size;
        using FloatComplexMatrix::min;
        using FloatComplexMatrix::max;
        using FloatComplexMatrix::resize;
        using FloatComplexMatrix::clear;
        using FloatComplexMatrix::rows;
        using FloatComplexMatrix::cols;
        using FloatComplexMatrix::row;
        
        void print() {
            std::cout << "Matrix(" << rows() << "," << cols() << ")\n";
            for(size_t i = 0; i < rows(); i++)
            for(size_t j = 0; j < cols(); j++)
                std::cout << (*this)(i,j) << ",";
            std::cout << "\n";
        }
        OctopusMatrixXcf operator + (const OctopusMatrixXcf & b) {
            OctopusMatrixXcf r;
            ValueList l;
            l(0) = "plus";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).float_complex_matrix_value();
            return r;
        }
        OctopusMatrixXcf operator - (const OctopusMatrixXcf & b) {
            OctopusMatrixXcf r;
            ValueList l;
            l(0) = "minus";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).float_complex_matrix_value();
            return r;
        }
        OctopusMatrixXcf operator * (const OctopusMatrixXcf & b) {
            OctopusMatrixXcf r;
            ValueList l;
            l(0) = "mul";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).float_complex_matrix_value();
            return r;
        }
        OctopusMatrixXcf operator / (const OctopusMatrixXcf & b) {
            OctopusMatrixXcf r;
            ValueList l;
            l(0) = "div";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).float_complex_matrix_value();
            return r;
        }
        OctopusMatrixXcf operator + (const std::complex<float> b) {
            OctopusMatrixXcf r;
            ValueList l;
            l(0) = "plus";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).float_complex_matrix_value();
            return r;
        }
        OctopusMatrixXcf operator - (const std::complex<float> b) {
            OctopusMatrixXcf r;
            ValueList l;
            l(0) = "minus";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).float_complex_matrix_value();
            return r;
        }
        OctopusMatrixXcf operator * (const std::complex<float> b) {
            OctopusMatrixXcf r;
            ValueList l;
            l(0) = "times";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).float_complex_matrix_value();
            return r;
        }
        OctopusMatrixXcf operator / (const std::complex<float> b) {
            OctopusMatrixXcf r;
            ValueList l;
            l(0) = "divide";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).float_complex_matrix_value();
            return r;
        }   

        #ifdef SWIG
        %extend
        {
            MatrixViewXcf __getitem__(size_t i) { return MatrixViewXcf($self,i); }            
            
            size_t rows() { return $self->rows(); }
            size_t cols() { return $self->cols(); }

            void fill(std::complex<float> v) { $self->fill(v); }
            
            //std::complex<float> min() { return $self->min(); }
            //std::complex<float> max() { return $self->max(); }
            
            const std::complex<float>* data() { return $self->data(); }

            void copy(size_t m, size_t n, std::complex<float> * p) {
                $self->resize(m,n);
                for(size_t i = 0; i < $self->rows(); i++)
                for(size_t j = 0; j < $self->cols(); j++)
                    (*$self)(i,j) = p[i*n + j];
            }

            void copy(size_t m, size_t n, std::complex<double> * p) {
                $self->resize(m,n);
                for(size_t i = 0; i < $self->rows(); i++)
                for(size_t j = 0; j < $self->cols(); j++)
                    (*$self)(i,j) = p[i*n + j];
            }

            OctopusMatrixXcf __add__ (const OctopusMatrixXcf & b) {
                assert($self->rows() == b.rows() && $self->cols() == b.cols());
                OctopusMatrixXcf r($self->size());
                for(size_t i = 0; i < $self->rows(); i++)
                for(size_t j = 0; j < $self->cols(); j++)
                    r(i,j) = (*$self)(i,j) + b(i,j);
                return r;
            }
            OctopusMatrixXcf __sub__ (const OctopusMatrixXcf & b) {
                assert($self->rows() == b.rows() && $self->cols() == b.cols());
                OctopusMatrixXcf r($self->size());
                for(size_t i = 0; i < $self->rows(); i++)
                for(size_t j = 0; j < $self->cols(); j++)
                    r(i,j) = (*$self)(i,j) - b(i,j);                    
                return r;
            }
            OctopusMatrixXcf __mul__ (const OctopusMatrixXcf & b) {
                return $self->matmul((*$self),b);
            }
            OctopusMatrixXcf __div__ (const OctopusMatrixXcf & b) {
                assert($self->rows() == b.rows() && $self->cols() == b.cols());
                OctopusMatrixXcf r($self->size());
                for(size_t i = 0; i < $self->rows(); i++)
                for(size_t j = 0; j < $self->cols(); j++)
                    r(i,j) = (*$self)(i,j) / b(i,j);
                return r;
            } 

            OctopusMatrixXcf __add__ (const std::complex<float> b) {
                OctopusMatrixXcf r($self->size());
                for(size_t i = 0; i < $self->rows(); i++)
                for(size_t j = 0; j < $self->cols(); j++)
                    r(i,j) = (*$self)(i,j) + b;
                return r;
            }
            OctopusMatrixXcf __sub__ (const std::complex<float> b) {
                OctopusMatrixXcf r($self->size());
                for(size_t i = 0; i < $self->rows(); i++)
                for(size_t j = 0; j < $self->cols(); j++)
                    r(i,j) = (*$self)(i,j) - b;                    
                return r;
            }
            OctopusMatrixXcf __mul__ (const std::complex<float> b) {
                OctopusMatrixXcf r($self->size());
                for(size_t i = 0; i < $self->rows(); i++)
                for(size_t j = 0; j < $self->cols(); j++)
                    r(i,j) = (*$self)(i,j) * b;
                return r;
            }
            OctopusMatrixXcf __div__ (const std::complex<float> b) {
                OctopusMatrixXcf r($self->size());
                for(size_t i = 0; i < $self->rows(); i++)
                for(size_t j = 0; j < $self->cols(); j++)
                    r(i,j) = (*$self)(i,j) / b;
                return r;
            } 
        }
        #endif
    };    

    inline std::complex<float>& MatrixViewXcf::operator[](size_t i) {
            return (*matrix)(row,i);
    }
    inline std::complex<float> MatrixViewXcf::operator[](size_t i) const {
        return (*matrix)(row,i);
    }
}