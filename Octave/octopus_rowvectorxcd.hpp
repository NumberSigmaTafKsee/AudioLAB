#pragma once

namespace Octopus
{
    struct OctopusRowVectorXcd: public ComplexRowVector
    {
        OctopusRowVectorXcd() =default;
        OctopusRowVectorXcd(const ComplexRowVector & v) : ComplexRowVector(v) {}   
        OctopusRowVectorXcd(size_t i) : ComplexRowVector(i) {}   

        using ComplexRowVector::operator =;
        using ComplexRowVector::operator ();
        using ComplexRowVector::insert;
        using ComplexRowVector::append;
        using ComplexRowVector::fill;
        using ComplexRowVector::extract;
        using ComplexRowVector::extract_n;
        using ComplexRowVector::transpose;
        using ComplexRowVector::size;
        using ComplexRowVector::min;
        using ComplexRowVector::max;
        using ComplexRowVector::resize;
        using ComplexRowVector::clear;

        void print()
        {
            ValueList l;
            l(0) = *this;
            octave::feval("display",l,0);
        }
        OctopusRowVectorXcd operator + (const OctopusRowVectorXcd & b) {
            OctopusRowVectorXcd r;
            ValueList l;
            l(0) = "add";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).complex_row_vector_value();
            return r;
        }
        OctopusRowVectorXcd operator - (const OctopusRowVectorXcd & b) {
            OctopusRowVectorXcd r;
            ValueList l;
            l(0) = "sub";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).complex_row_vector_value();
            return r;
        }
        OctopusRowVectorXcd operator * (const OctopusRowVectorXcd & b) {
            OctopusRowVectorXcd r;
            ValueList l;
            l(0) = "times";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).complex_row_vector_value();
            return r;
        }
        OctopusRowVectorXcd operator / (const OctopusRowVectorXcd & b) {
            OctopusRowVectorXcd r;
            ValueList l;
            l(0) = "div";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).complex_row_vector_value();
            return r;
        }
        OctopusRowVectorXcd operator + (const std::complex<float> b) {
            OctopusRowVectorXcd r;
            ValueList l;
            l(0) = "plus";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).complex_row_vector_value();
            return r;
        }
        OctopusRowVectorXcd operator - (const std::complex<float> b) {
            OctopusRowVectorXcd r;
            ValueList l;
            l(0) = "divide";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).complex_row_vector_value();
            return r;
        }
        OctopusRowVectorXcd operator * (const std::complex<float> b) {
            OctopusRowVectorXcd r;
            ValueList l;
            l(0) = "times";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).complex_row_vector_value();
            return r;
        }
        OctopusRowVectorXcd operator / (const std::complex<float> b) {
            OctopusRowVectorXcd r;
            ValueList l;
            l(0) = "div";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).complex_row_vector_value();
            return r;
        }   

        #ifdef SWIG
        %extend
        {
            std::complex<double> __getitem__(size_t i) { return (*$self)(i); }
            void __setitem__(size_t i, std::complex<double> v) { (*$self)(i) = v; }
            size_t size() const { return $self->size(1); }
            void fill(std::complex<double> v) { $self->fill(v); }
            //std::complex<float> min() { return $self->min(); }
            //std::complex<float> max() { return $self->max(); }
            const std::complex<double>* data() { return $self->data(); }
            
            void copy(size_t n, std::complex<float> * p) {
                $self->resize(n);
                for(size_t i = 0; i < n; i++) (*$self)(i) = p[i];
            }
            void copy(size_t n, std::complex<double> * p) {
                $self->resize(n);
                for(size_t i = 0; i < n; i++) (*$self)(i) = p[i];
            }

            OctopusRowVectorXcd __add__ (const OctopusRowVectorXcd & b) {
                assert($self->size(1) == b.size(1));
                OctopusRowVectorXcd r($self->size(1));
                for(size_t i = 0; i < $self->size(1); i++)
                    r(i) = (*$self)(i) + b(i);
                return r;
            }
            OctopusRowVectorXcd __sub__ (const OctopusRowVectorXcd & b) {
                assert($self->size(1) == b.size(1));
                OctopusRowVectorXcd r($self->size(1));
                for(size_t i = 0; i < $self->size(1); i++)
                    r(i) = (*$self)(i) - b(i);
                return r;
            }
            OctopusRowVectorXcd __mul__ (const OctopusRowVectorXcd & b) {
                assert($self->size(1) == b.size(1));
                OctopusRowVectorXcd r($self->size(1));
                for(size_t i = 0; i < $self->size(1); i++)
                    r(i) = (*$self)(i) * b(i);
                return r;
            }
            OctopusRowVectorXcd __div__ (const OctopusRowVectorXcd & b) {
                assert($self->size(1) == b.size(1));
                OctopusRowVectorXcd r($self->size(1));
                for(size_t i = 0; i < $self->size(1); i++)
                    r(i) = (*$self)(i) / b(i);
                return r;
            }     

            OctopusRowVectorXcd __add__ (const std::complex<double> b) {                
                OctopusRowVectorXcd r($self->size(1));
                for(size_t i = 0; i < $self->size(1); i++)
                    r(i) = (*$self)(i) + b;
                return r;
            }
            OctopusRowVectorXcd __sub__ (const std::complex<double> b) {                
                OctopusRowVectorXcd r($self->size(1));
                for(size_t i = 0; i < $self->size(1); i++)
                    r(i) = (*$self)(i) - b;
                return r;
            }
            OctopusRowVectorXcd __mul__ (const std::complex<double> b) {                
                OctopusRowVectorXcd r($self->size(1));
                for(size_t i = 0; i < $self->size(1); i++)
                    r(i) = (*$self)(i) * b;
                return r;
            }
            OctopusRowVectorXcd __div__ (const std::complex<double> b) {                
                OctopusRowVectorXcd r($self->size(1));
                for(size_t i = 0; i < $self->size(1); i++)
                    r(i) = (*$self)(i) / b;
                return r;
            }            
        
        }
        #endif
    };
}