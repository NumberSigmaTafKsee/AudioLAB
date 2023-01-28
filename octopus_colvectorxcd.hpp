#pragma once

namespace Octopus
{
    struct OctopusColVectorXcd: public ComplexColumnVector
    {
        OctopusColVectorXcd() = default;
        OctopusColVectorXcd(const ComplexColumnVector &v) : ComplexColumnVector(v) {}
        OctopusColVectorXcd(size_t i) : ComplexColumnVector(i) {}

        using ComplexColumnVector::operator =;
        using ComplexColumnVector::operator ();
        using ComplexColumnVector::insert;
        //using ComplexColumnVector::append;
        using ComplexColumnVector::fill;
        using ComplexColumnVector::extract;
        using ComplexColumnVector::extract_n;
        using ComplexColumnVector::transpose;
        using ComplexColumnVector::size;
        using ComplexColumnVector::min;
        using ComplexColumnVector::max;
        using ComplexColumnVector::resize;
        using ComplexColumnVector::clear;

        void print()
        {
            ValueList l;
            l(0) = *this;
            octave::feval("display",l,0);
        }

        OctopusColVectorXcd operator + (const OctopusColVectorXcd & b) {
            OctopusColVectorXcd r;
            ValueList l;
            l(0) = "plus";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).complex_column_vector_value();
            return r;
        }
        OctopusColVectorXcd operator - (const OctopusColVectorXcd & b) {
            OctopusColVectorXcd r;
            ValueList l;
            l(0) = "minus";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).complex_column_vector_value();
            return r;
        }
        OctopusColVectorXcd operator * (const OctopusColVectorXcd & b) {
            OctopusColVectorXcd r;
            ValueList l;
            l(0) = "times";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).complex_column_vector_value();
            return r;
        }
        OctopusColVectorXcd operator / (const OctopusColVectorXcd & b) {
            OctopusColVectorXcd r;
            ValueList l;
            l(0) = "divide";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).complex_column_vector_value();
            return r;
        }
        OctopusColVectorXcd operator + (const std::complex<double> b) {
            OctopusColVectorXcd r;
            ValueList l;
            l(0) = "plus";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).complex_column_vector_value();
            return r;
        }
        OctopusColVectorXcd operator - (const std::complex<double> b) {
            OctopusColVectorXcd r;
            ValueList l;
            l(0) = "minus";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).complex_column_vector_value();
            return r;
        }
        OctopusColVectorXcd operator * (const std::complex<double> b) {
            OctopusColVectorXcd r;
            ValueList l;
            l(0) = "times";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).complex_column_vector_value();
            return r;
        }
        OctopusColVectorXcd operator / (const std::complex<double> b) {
            OctopusColVectorXcd r;
            ValueList l;
            l(0) = "divide";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).complex_column_vector_value();
            return r;
        }   

        #ifdef SWIG
        %extend {
            void __setitem__(size_t i, std::complex<double> v) { (*$self)(i) = v; }
            size_t size() const { return $self->size(1); }

            void fill(std::complex<double> v) { $self->fill(v); }
            
            //std::complex<double> min() { return $self->min(); }
            //std::complex<double> max() { return $self->max(); }

            const std::complex<double>* data() { return $self->data(); }
            void copy(size_t n, std::complex<double> * p) {
                $self->resize(n);
                for(size_t i = 0; i < n; i++) (*$self)(i) = p[i];
            }
            OctopusColVectorXcd __add__ (const OctopusColVectorXcd & b) {
                assert($self->size(1) == b.size(1));
                OctopusColVectorXcd r($self->size(1));
                for(size_t i = 0; i < $self->size(1); i++)
                    r(i) = (*$self)(i) + b(i);
                return r;
            }
            OctopusColVectorXcd __sub__ (const OctopusColVectorXcd & b) {
                assert($self->size(1) == b.size(1));
                OctopusColVectorXcd r($self->size(1));
                for(size_t i = 0; i < $self->size(1); i++)
                    r(i) = (*$self)(i) - b(i);
                return r;
            }
            OctopusColVectorXcd __mul__ (const OctopusColVectorXcd & b) {
                assert($self->size(1) == b.size(1));
                OctopusColVectorXcd r($self->size(1));
                for(size_t i = 0; i < $self->size(1); i++)
                    r(i) = (*$self)(i) * b(i);
                return r;
            }
            OctopusColVectorXcd __div__ (const OctopusColVectorXcd & b) {
                assert($self->size(1) == b.size(1));
                OctopusColVectorXcd r($self->size(1));
                for(size_t i = 0; i < $self->size(1); i++)
                    r(i) = (*$self)(i) / b(i);
                return r;
            } 

            OctopusColVectorXcd __add__ (std::complex<double> b) {
                OctopusColVectorXcd r($self->size(1));
                for(size_t i = 0; i < $self->size(1); i++)
                    r(i) = (*$self)(i) + b;
                return r;
            }
            OctopusColVectorXcd __sub__ (std::complex<double> b) {
                OctopusColVectorXcd r($self->size(1));
                for(size_t i = 0; i < $self->size(1); i++)
                    r(i) = (*$self)(i) - b;
                return r;
            }
            OctopusColVectorXcd __mul__ (std::complex<double> b) {
                OctopusColVectorXcd r($self->size(1));
                for(size_t i = 0; i < $self->size(1); i++)
                    r(i) = (*$self)(i) * b;
                return r;
            }
            OctopusColVectorXcd __div__ (std::complex<double> b) {
                OctopusColVectorXcd r($self->size(1));
                for(size_t i = 0; i < $self->size(1); i++)
                    r(i) = (*$self)(i) / b;
                return r;
            } 
        }
        #endif
    };
}