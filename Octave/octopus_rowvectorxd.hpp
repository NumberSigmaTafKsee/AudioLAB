#pragma once

namespace Octopus
{
    struct OctopusRowVectorXd: public RowVector
    {
        OctopusRowVectorXd() : RowVector() {}
        OctopusRowVectorXd(const RowVector & v) : RowVector(v) {} 
        OctopusRowVectorXd(size_t i) : RowVector(i) {}

        using RowVector::operator =;
        using RowVector::operator ();
        using RowVector::insert;
        using RowVector::append;
        using RowVector::fill;
        using RowVector::extract;
        using RowVector::extract_n;
        using RowVector::transpose;
        using RowVector::size;
        using RowVector::min;
        using RowVector::max;
        using RowVector::resize;
        using RowVector::clear;

        void print()
        {
            ValueList l;
            l(0) = *this;
            octave::feval("display",l,0);
        }

        OctopusRowVectorXd operator + (const OctopusRowVectorXd & b) {
            OctopusRowVectorXd r;
            ValueList l;
            l(0) = "add";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).row_vector_value();
            return r;
        }
        OctopusRowVectorXd operator - (const OctopusRowVectorXd & b) {
            OctopusRowVectorXd r;
            ValueList l;
            l(0) = "sub";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).row_vector_value();
            return r;
        }
        OctopusRowVectorXd operator * (const OctopusRowVectorXd & b) {
            OctopusRowVectorXd r;
            ValueList l;
            l(0) = "times";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).row_vector_value();
            return r;
        }
        OctopusRowVectorXd operator / (const OctopusRowVectorXd & b) {
            OctopusRowVectorXd r;
            ValueList l;
            l(0) = "divide";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).row_vector_value();
            return r;
        }
        OctopusRowVectorXd operator + (const double b) {
            OctopusRowVectorXd r;
            ValueList l;
            l(0) = "plus";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).row_vector_value();
            return r;
        }
        OctopusRowVectorXd operator - (const double b) {
            OctopusRowVectorXd r;
            ValueList l;
            l(0) = "minus";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).row_vector_value();
            return r;
        }
        OctopusRowVectorXd operator * (const double b) {
            OctopusRowVectorXd r;
            ValueList l;
            l(0) = "times";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).row_vector_value();
            return r;
        }
        OctopusRowVectorXd operator / (const double b) {
            OctopusRowVectorXd r;
            ValueList l;
            l(0) = "divide";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).row_vector_value();
            return r;
        }        

        #ifdef SWIG
        %extend
        {
            double __getitem__(size_t i) { return (*$self)(i); }
            void __setitem__(size_t i, double v) { (*$self)(i) = v; }
            
            size_t size() const { return $self->size(1); }
            
            void fill(double v) { $self->fill(v); }
            
            //double min() { return $self->min(); }
            //double max() { return $self->max(); }
            
            const double* data() { return $self->data(); }

            void copy(size_t n, float * p) {
                $self->resize(n);
                for(size_t i = 0; i < n; i++) (*$self)(i) = p[i];
            }
            void copy(size_t n, double * p) {
                $self->resize(n);
                for(size_t i = 0; i < n; i++) (*$self)(i) = p[i];
            }
            
            OctopusRowVectorXd __add__ (const OctopusRowVectorXd & b) {
                assert($self->size(1) == b.size(1));
                OctopusRowVectorXd r($self->size(1));
                for(size_t i = 0; i < $self->size(1); i++)
                    r(i) = (*$self)(i) + b(i);
                return r;
            }
            OctopusRowVectorXd __sub__ (const OctopusRowVectorXd & b) {
                assert($self->size(1) == b.size(1));
                OctopusRowVectorXd r($self->size(1));
                for(size_t i = 0; i < $self->size(1); i++)
                    r(i) = (*$self)(i) - b(i);
                return r;
            }
            OctopusRowVectorXd __mul__ (const OctopusRowVectorXd & b) {
                assert($self->size(1) == b.size(1));
                OctopusRowVectorXcd r($self->size(1));
                for(size_t i = 0; i < $self->size(1); i++)
                    r(i) = (*$self)(i) * b(i);
                return r;
            }
            OctopusRowVectorXd __div__ (const OctopusRowVectorXd & b) {                
                OctopusRowVectorXd r($self->size(1));
                for(size_t i = 0; i < $self->size(1); i++)
                    r(i) = (*$self)(i) / b;
                return r;
            } 

            OctopusRowVectorXd __add__ (const double b) {
                assert($self->size(1) == b.size(1));
                OctopusRowVectorXd r($self->size(1));
                for(size_t i = 0; i < $self->size(1); i++)
                    r(i) = (*$self)(i) + b(i);
                return r;
            }
            OctopusRowVectorXd __sub__ (const double b) {                
                OctopusRowVectorXd r($self->size(1));
                for(size_t i = 0; i < $self->size(1); i++)
                    r(i) = (*$self)(i) - b;
                return r;
            }
            OctopusRowVectorXd __mul__ (const double b) {                
                OctopusRowVectorXcd r($self->size(1));
                for(size_t i = 0; i < $self->size(1); i++)
                    r(i) = (*$self)(i) * b;
                return r;
            }
            OctopusRowVectorXd __div__ (const double & b) {                
                OctopusRowVectorXd r($self->size(1));
                for(size_t i = 0; i < $self->size(1); i++)
                    r(i) = (*$self)(i) / b;
                return r;
            } 
        }
        #endif
    };
}