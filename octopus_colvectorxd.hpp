#pragma once

namespace Octopus
{
        struct OctopusColVectorXd: public ColumnVector
    {
        OctopusColVectorXd() = default;
        OctopusColVectorXd(const ColumnVector &v) : ColumnVector(v) {}
        OctopusColVectorXd(size_t i) : ColumnVector(i) {}

        using ColumnVector::operator =;
        using ColumnVector::operator ();
        using ColumnVector::insert;
        //using ColumnVector::append;
        using ColumnVector::fill;
        using ColumnVector::extract;
        using ColumnVector::extract_n;
        using ColumnVector::transpose;
        using ColumnVector::size;
        using ColumnVector::min;
        using ColumnVector::max;
        using ColumnVector::resize;
        using ColumnVector::clear;

        OctopusColVectorXd operator + (const OctopusColVectorXd & b) {
            OctopusColVectorXd r;
            ValueList l;
            l(0) = "plus";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).column_vector_value();
            return r;
        }
        OctopusColVectorXd operator - (const OctopusColVectorXd & b) {
            OctopusColVectorXd r;
            ValueList l;
            l(0) = "minus";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).column_vector_value();
            return r;
        }
        OctopusColVectorXd operator * (const OctopusColVectorXd & b) {
            OctopusColVectorXd r;
            ValueList l;
            l(0) = "times";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).column_vector_value();
            return r;
        }
        OctopusColVectorXd operator / (const OctopusColVectorXd & b) {
            OctopusColVectorXd r;
            ValueList l;
            l(0) = "divide";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).column_vector_value();
            return r;
        }
        OctopusColVectorXd operator + (const double b) {
            OctopusColVectorXd r;
            ValueList l;
            l(0) = "plus";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).column_vector_value();
            return r;
        }
        OctopusColVectorXd operator - (const double b) {
            OctopusColVectorXd r;
            ValueList l;
            l(0) = "minus";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).column_vector_value();
            return r;
        }
        OctopusColVectorXd operator * (const double b) {
            OctopusColVectorXd r;
            ValueList l;
            l(0) = "times";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).column_vector_value();
            return r;
        }
        OctopusColVectorXd operator / (const double b) {
            OctopusColVectorXd r;
            ValueList l;
            l(0) = "divide";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).column_vector_value();
            return r;
        }   

        #ifdef SWIG
        %extend
        {
            double __getitem__(size_t i) { return (*$self)(i); }
            void __setitem__(size_t i, double v) { (*$self)(i) = v; }

            size_t size() const { return $self->size(1); }

            void fill(double v) { $self->fill(v); }
            
            // double min() { return $self->min(); }
            // double max() { return $self->max(); }
            
            const double* data() { return $self->data(); }
            
            void copy(size_t n, float * p) {
                $self->resize(n);
                for(size_t i = 0; i < n; i++) (*$self)(i) = p[i];
            }
            void copy(size_t n, double * p) {
                $self->resize(n);
                for(size_t i = 0; i < n; i++) (*$self)(i) = p[i];
            }

            OctopusColVectorXd __add__ (const OctopusColVectorXd & b) {
                assert($self->size(1) == b.size(1));
                OctopusColVectorXd r($self->size(1));
                for(size_t i = 0; i < $self->size(1); i++)
                    r(i) = (*$self)(i) + b(i);
                return r;
            }
            OctopusColVectorXd __sub__ (const OctopusColVectorXd & b) {
                assert($self->size(1) == b.size(1));
                OctopusColVectorXd r($self->size(1));
                for(size_t i = 0; i < $self->size(1); i++)
                    r(i) = (*$self)(i) - b(i);
                return r;
            }
            OctopusColVectorXd __mul__ (const OctopusColVectorXd & b) {
                assert($self->size(1) == b.size(1));
                OctopusColVectorXd r($self->size(1));
                for(size_t i = 0; i < $self->size(1); i++)
                    r(i) = (*$self)(i) * b(i);
                return r;
            }
            OctopusColVectorXd __div__ (const OctopusColVectorXd & b) {
                assert($self->size(1) == b.size(1));
                OctopusColVectorXd r($self->size(1));
                for(size_t i = 0; i < $self->size(1); i++)
                    r(i) = (*$self)(i) / b(i);
                return r;
            } 


            OctopusColVectorXd __add__ (double b) {
                OctopusColVectorXd r($self->size(1));
                for(size_t i = 0; i < $self->size(1); i++)
                    r(i) = (*$self)(i) + b;
                return r;
            }
            OctopusColVectorXd __sub__ (double b) {
                OctopusColVectorXd r($self->size(1));
                for(size_t i = 0; i < $self->size(1); i++)
                    r(i) = (*$self)(i) - b;
                return r;
            }
            OctopusColVectorXd __mul__ (double b) {
                OctopusColVectorXd r($self->size(1));
                for(size_t i = 0; i < $self->size(1); i++)
                    r(i) = (*$self)(i) * b;
                return r;
            }
            OctopusColVectorXd __div__ (double b) {
                OctopusColVectorXd r($self->size(1));
                for(size_t i = 0; i < $self->size(1); i++)
                    r(i) = (*$self)(i) / b;
                return r;
            } 

        }
        #endif

        void print()
        {
            ValueList l;
            l(0) = *this;
            octave::feval("display",l,0);
        }
    };
}