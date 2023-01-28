#pragma once

namespace Octopus
{
    struct OctopusMatrixXcd;
    struct MatrixViewXcd
    {
        OctopusMatrixXcd * matrix;
        size_t row;

        MatrixViewXcd(OctopusMatrixXcd * m, size_t r) : matrix(m),row(r) {}

        std::complex<double>& operator[](size_t i);
        std::complex<double>  operator[](size_t i) const;
    };
    struct OctopusMatrixXcd : public ComplexMatrix
    {
        OctopusMatrixXcd() = default;
        OctopusMatrixXcd(const ComplexMatrix &v) : ComplexMatrix(v) {}
        OctopusMatrixXcd(size_t i,size_t j) : ComplexMatrix(i,j) {}

        using ComplexMatrix::operator =;
        using ComplexMatrix::operator ();
        using ComplexMatrix::insert;
        using ComplexMatrix::append;
        using ComplexMatrix::fill;
        using ComplexMatrix::extract;
        using ComplexMatrix::extract_n;
        using ComplexMatrix::transpose;
        using ComplexMatrix::size;
        using ComplexMatrix::min;
        using ComplexMatrix::max;
        using ComplexMatrix::resize;        
        using ComplexMatrix::clear;
        using ComplexMatrix::rows;
        using ComplexMatrix::cols;
        using ComplexMatrix::row;

        void print() {
            std::cout << "Matrix(" << rows() << "," << cols() << ")\n";
            for(size_t i = 0; i < rows(); i++)
            for(size_t j = 0; j < cols(); j++)
                std::cout << (*this)(i,j) << ",";
            std::cout << "\n";
        }
        OctopusMatrixXcd operator + (const OctopusMatrixXcd & b) {
            OctopusMatrixXcd r;
            ValueList l;
            l(0) = "plus";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).complex_matrix_value();
            return r;
        }
        OctopusMatrixXcd operator - (const OctopusMatrixXcd & b) {
            OctopusMatrixXcd r;
            ValueList l;
            l(0) = "minus";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).complex_matrix_value();
            return r;
        }
        OctopusMatrixXcd operator * (const OctopusMatrixXcd & b) {
            OctopusMatrixXcd r;
            ValueList l;
            l(0) = "mul";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).complex_matrix_value();
            return r;
        }
        OctopusMatrixXcd operator / (const OctopusMatrixXcd & b) {
            OctopusMatrixXcd r;
            ValueList l;
            l(0) = "div";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).complex_matrix_value();
            return r;
        }
        OctopusMatrixXcd operator + (const std::complex<double> b) {
            OctopusMatrixXcd r;
            ValueList l;
            l(0) = "plus";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).complex_matrix_value();
            return r;
        }
        OctopusMatrixXcd operator - (const std::complex<double> b) {
            OctopusMatrixXcd r;
            ValueList l;
            l(0) = "minus";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).complex_matrix_value();
            return r;
        }
        OctopusMatrixXcd operator * (const std::complex<double> b) {
            OctopusMatrixXcd r;
            ValueList l;
            l(0) = "times";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).complex_matrix_value();
            return r;
        }
        OctopusMatrixXcd operator / (const std::complex<double> b) {
            OctopusMatrixXcd r;
            ValueList l;
            l(0) = "divide";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).complex_matrix_value();
            return r;
        }   

        #ifdef SWIG
        %extend
        {
            MatrixViewXcd __getitem__(size_t i) { return MatrixViewXcd($self,i); }            

            size_t rows() { return $self->rows(); }
            size_t cols() { return $self->cols(); }

            void fill(std::complex<double> v) { $self->fill(v); }
            
            //std::complex<double> min() { return $self->min(); }
            //std::complex<double> max() { return $self->max(); }
            
            const std::complex<double>* data() { return $self->data(); }

            void copy(size_t m, size_t n, std::complex<double> * p) {
                $self->resize(m,n);
                for(size_t i = 0; i < $self->rows(); i++)
                for(size_t j = 0; j < $self->cols(); j++)
                    (*$self)(i,j) = p[i*n + j];
            }

            void copy(size_t m, size_t n, std::complex<float> * p) {
                $self->resize(m,n);
                for(size_t i = 0; i < $self->rows(); i++)
                for(size_t j = 0; j < $self->cols(); j++)
                    (*$self)(i,j) = p[i*n + j];
            }

            OctopusMatrixXcd __add__ (const OctopusMatrixXcd & b) {                
                assert($self->rows() == b.rows() && $self->cols() == b.cols());
                OctopusMatrixXcd r($self->size());
                for(size_t i = 0; i < $self->rows(); i++)
                for(size_t j = 0; j < $self->cols(); j++)
                    r(i,j) = (*$self)(i,j) + b(i,j);
                return r;
            }
            OctopusMatrixXcd __sub__ (const OctopusMatrixXcd & b) {            
                assert($self->rows() == b.rows() && $self->cols() == b.cols());    
                OctopusMatrixXcd r($self->size());
                for(size_t i = 0; i < $self->rows(); i++)
                for(size_t j = 0; j < $self->cols(); j++)
                    r(i,j) = (*$self)(i,j) - b(i,j);                    
                return r;
            }
            OctopusMatrixXcd __mul__ (const OctopusMatrixXcd & b) {
                return $self->matmul((*$self),b);
            }
            OctopusMatrixXcd __div__ (const OctopusMatrixXcd & b) {            
                assert($self->rows() == b.rows() && $self->cols() == b.cols());    
                OctopusMatrixXcd r($self->size());
                for(size_t i = 0; i < $self->rows(); i++)
                for(size_t j = 0; j < $self->cols(); j++)
                    r(i,j) = (*$self)(i,j) / b(i,j);
                return r;
            } 

            OctopusMatrixXcd __add__ (const std::complex<double> b) {                
                OctopusMatrixXcd r($self->size());
                for(size_t i = 0; i < $self->rows(); i++)
                for(size_t j = 0; j < $self->cols(); j++)
                    r(i,j) = (*$self)(i,j) + b;
                return r;
            }
            OctopusMatrixXcd __sub__ (const std::complex<double> b) {            
                OctopusMatrixXcd r($self->size());
                for(size_t i = 0; i < $self->rows(); i++)
                for(size_t j = 0; j < $self->cols(); j++)
                    r(i,j) = (*$self)(i,j) - b;                    
                return r;
            }
            OctopusMatrixXcd __mul__ (const std::complex<double> b) {                
                OctopusMatrixXcd r($self->size());
                for(size_t i = 0; i < $self->rows(); i++)
                for(size_t j = 0; j < $self->cols(); j++)
                    r(i,j) = (*$self)(i,j) * b;
                return r;
            }
            OctopusMatrixXcd __div__ (const std::complex<double> b) {                
                OctopusMatrixXcd r($self->size());
                for(size_t i = 0; i < $self->rows(); i++)
                for(size_t j = 0; j < $self->cols(); j++)
                    r(i,j) = (*$self)(i,j) / b;
                return r;
            } 
        }
        #endif
 
    };        

    inline std::complex<double>& MatrixViewXcd::operator[](size_t i) {
            return (*matrix)(row,i);
    }
    inline std::complex<double> MatrixViewXcd::operator[](size_t i) const {
        return (*matrix)(row,i);
    }
}