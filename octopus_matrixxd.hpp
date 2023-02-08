#pragma once

namespace Octopus
{
    struct OctopusMatrixXd;
    struct MatrixViewXd
    {
        OctopusMatrixXd * matrix;
        size_t row;

        MatrixViewXd(OctopusMatrixXd * m, size_t r) : matrix(m),row(r) {}

        double& operator[](size_t i);
        double  operator[](size_t i) const;

        double __getitem__(size_t i);
        void __setitem__(size_t i, const double & v);

    };
    struct OctopusMatrixXd : public Matrix
    {
        OctopusMatrixXd() = default;
        OctopusMatrixXd(const Matrix &v) : Matrix(v) {}
        OctopusMatrixXd(size_t i,size_t j) : Matrix(i,j) {}

        using Matrix::operator =;
        using Matrix::operator ();
        using Matrix::insert;
        using Matrix::append;
        using Matrix::fill;
        using Matrix::extract;
        using Matrix::extract_n;
        using Matrix::transpose;
        using Matrix::size;
        using Matrix::min;
        using Matrix::max;
        using Matrix::resize;
        using Matrix::clear;
        using Matrix::rows;
        using Matrix::cols;
        using Matrix::row;

        void print() {
            std::cout << "Matrix(" << rows() << "," << cols() << ")\n";
            for(size_t i = 0; i < rows(); i++)
            for(size_t j = 0; j < cols(); j++)
                std::cout << (*this)(i,j) << ",";
            std::cout << "\n";
        }

        #ifdef SWIG
        %extend
        {
            double get(size_t i, size_t j) { return (*$self)(i,j); }
            void set(size_t i, size_t j, double v) { (*$self)(i,j) = v; }
            size_t rows() { return $self->rows(); }
            size_t cols() { return $self->cols(); }
        }
        #endif
        OctopusMatrixXd operator + (const OctopusMatrixXd & b) {
            OctopusMatrixXd r;
            ValueList l;
            l(0) = "plus";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).matrix_value();
            return r;
        }
        OctopusMatrixXd operator - (const OctopusMatrixXd & b) {
            OctopusMatrixXd r;
            ValueList l;
            l(0) = "minus";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).matrix_value();
            return r;
        }
        OctopusMatrixXd operator * (const OctopusMatrixXd & b) {
            OctopusMatrixXd r;
            ValueList l;
            l(0) = "mul";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).matrix_value();
            return r;
        }
        OctopusMatrixXd operator / (const OctopusMatrixXd & b) {
            OctopusMatrixXd r;
            ValueList l;
            l(0) = "div";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).matrix_value();
            return r;
        }
        OctopusMatrixXd operator + (const double b) {
            OctopusMatrixXd r;
            ValueList l;
            l(0) = "plus";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).matrix_value();
            return r;
        }
        OctopusMatrixXd operator - (const double b) {
            OctopusMatrixXd r;
            ValueList l;
            l(0) = "minus";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).matrix_value();
            return r;
        }
        OctopusMatrixXd operator * (const double b) {
            OctopusMatrixXd r;
            ValueList l;
            l(0) = "times";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).matrix_value();
            return r;
        }
        OctopusMatrixXd operator / (const double b) {
            OctopusMatrixXd r;
            ValueList l;
            l(0) = "divide";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).matrix_value();
            return r;
        }   

        #ifdef SWIG
        %extend
        {
            MatrixViewXd __getitem__(size_t i) { return MatrixViewXd($self,i-1); }            

            size_t rows() { return $self->rows(); }
            size_t cols() { return $self->cols(); }

            void fill(double v) { $self->fill(v); }

            //double min() { return $self->min(); }
            //double max() { return $self->max(); }

            const double* data() { return $self->data(); }

            void copy(size_t m, size_t n, float * p) {
                $self->resize(m,n);
                for(size_t i = 0; i < $self->rows(); i++)
                for(size_t j = 0; j < $self->cols(); j++)
                    (*$self)(i,j) = p[i*n + j];
            }

            void copy(size_t m, size_t n, double * p) {
                $self->resize(m,n);
                for(size_t i = 0; i < $self->rows(); i++)
                for(size_t j = 0; j < $self->cols(); j++)
                    (*$self)(i,j) = p[i*n + j];
            }

            OctopusMatrixXd __add__ (const OctopusMatrixXd & b) {                
                OctopusMatrixXd r($self->size());
                for(size_t i = 0; i < $self->rows(); i++)
                for(size_t j = 0; j < $self->cols(); j++)            
                    r(i,j) = (*$self)[i] + b(i,j);
                return r;
            }
            OctopusMatrixXd __sub__ (const OctopusMatrixXd & b) {
                assert($self->size() == b.size());
                OctopusMatrixXd r($self->size());
                for(size_t i = 0; i < $self->rows(); i++)
                for(size_t j = 0; j < $self->cols(); j++)            
                    r(i,j) = (*$self)[i] - b(i,j);
                return r;
            }
            OctopusMatrixXd __mul__ (const OctopusMatrixXd & b) {
                assert($self->size() == b.size());
                OctopusMatrixXd r($self->size());
                for(size_t i = 0; i < $self->rows(); i++)
                for(size_t j = 0; j < $self->cols(); j++)            
                    r(i,j) = (*$self)[i] * b(i,j);
                return r;
            }
            OctopusMatrixXd __div__ (const OctopusMatrixXd & b) {
                assert($self->size() == b.size());
                OctopusMatrixXd r($self->size());
                for(size_t i = 0; i < $self->rows(); i++)
                for(size_t j = 0; j < $self->cols(); j++)            
                    r(i,j) = (*$self)[i] / b(i,j);
                return r;
            } 

            
            OctopusMatrixXd __add__ (const double b) {
                assert($self->rows() == b.rows() && $self->cols() == b.cols());
                OctopusMatrixXd r($self->size());
                for(size_t i = 0; i < $self->rows(); i++)
                for(size_t j = 0; j < $self->cols(); j++)
                    r(i,j) = (*$self)(i,j) + b;
                return r;
            }
            OctopusMatrixXd __sub__ (const double b) {
                assert($self->rows() == b.rows() && $self->cols() == b.cols());
                OctopusMatrixXd r($self->size());
                for(size_t i = 0; i < $self->rows(); i++)
                for(size_t j = 0; j < $self->cols(); j++)
                    r(i,j) = (*$self)(i,j) - b;                    
                return r;
            }
            OctopusMatrixXd __mul__ (const double b) {
                assert($self->rows() == b.rows() && $self->cols() == b.cols());
                OctopusMatrixXd r($self->size());
                for(size_t i = 0; i < $self->rows(); i++)
                for(size_t j = 0; j < $self->cols(); j++)
                    r(i,j) = (*$self)(i,j) * b;
                return r;
            }
            OctopusMatrixXd __div__ (const double b) {
                assert($self->rows() == b.rows() && $self->cols() == b.cols());
                OctopusMatrixXd r($self->size());
                for(size_t i = 0; i < $self->rows(); i++)
                for(size_t j = 0; j < $self->cols(); j++)
                    r(i,j) = (*$self)(i,j) / b;
                return r;
            } 
        }
        #endif
    };    

    inline double& MatrixViewXd::operator[](size_t i) {
            return (*matrix)(row,i);
    }
    inline double MatrixViewXd::operator[](size_t i) const {
        return (*matrix)(row,i);
    }
    inline double MatrixViewXd::__getitem__(size_t i) { return (*this)[i-1]; }
    inline void MatrixViewXd::__setitem__(size_t i, const double & v) { (*this)[i-1] = v; }
}