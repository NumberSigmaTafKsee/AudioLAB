#pragma once

namespace Octopus
{
    struct OctopusRowVectorXcf : public FloatComplexRowVector
    {
        OctopusRowVectorXcf() = default;
        OctopusRowVectorXcf(const FloatComplexRowVector & v) : FloatComplexRowVector(v) {}
        OctopusRowVectorXcf(size_t i) : FloatComplexRowVector(i) {}

        using FloatComplexRowVector::operator =;
        using FloatComplexRowVector::operator ();
        using FloatComplexRowVector::insert;
        using FloatComplexRowVector::append;
        using FloatComplexRowVector::fill;
        using FloatComplexRowVector::extract;
        using FloatComplexRowVector::extract_n;
        using FloatComplexRowVector::transpose;
        using FloatComplexRowVector::size;
        using FloatComplexRowVector::min;
        using FloatComplexRowVector::max;
        using FloatComplexRowVector::resize;
        using FloatComplexRowVector::clear;

        OctopusRowVectorXcf operator + (const OctopusRowVectorXcf & b) {
            OctopusRowVectorXcf r;
            ValueList l;
            l(0) = "add";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).float_complex_row_vector_value();
            return r;
        }
        OctopusRowVectorXcf operator - (const OctopusRowVectorXcf & b) {
            OctopusRowVectorXcf r;
            ValueList l;
            l(0) = "sub";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).float_complex_row_vector_value();
            return r;
        }
        OctopusRowVectorXcf operator * (const OctopusRowVectorXcf & b) {
            OctopusRowVectorXcf r;
            ValueList l;
            l(0) = "times";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).float_complex_row_vector_value();
            return r;
        }
        OctopusRowVectorXcf operator / (const OctopusRowVectorXcf & b) {
            OctopusRowVectorXcf r;
            ValueList l;
            l(0) = "divide";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).float_complex_row_vector_value();
            return r;
        }
        OctopusRowVectorXcf operator + (const std::complex<float> b) {
            OctopusRowVectorXcf r;
            ValueList l;
            l(0) = "plus";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).float_complex_row_vector_value();
            return r;
        }
        OctopusRowVectorXcf operator - (const std::complex<float> b) {
            OctopusRowVectorXcf r;
            ValueList l;
            l(0) = "minus";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).float_complex_row_vector_value();
            return r;
        }
        OctopusRowVectorXcf operator * (const std::complex<float> b) {
            OctopusRowVectorXcf r;
            ValueList l;
            l(0) = "times";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).float_complex_row_vector_value();
            return r;
        }
        OctopusRowVectorXcf operator / (const std::complex<float> b) {
            OctopusRowVectorXcf r;
            ValueList l;
            l(0) = "divide";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).float_complex_row_vector_value();
            return r;
        }        


        #ifdef SWIG
        %extend
        {
            std::complex<float> __getitem__(size_t i) { return (*$self)(i); }
            void __setitem__(size_t i, std::complex<float> v) { (*$self)(i) = v; }
            
            size_t size() const { return $self->size(1); }
            
            void fill(std::complex<float> v) { $self->fill(v); }
            
            // std::complex<float> min() { return $self->min(); }
            // std::complex<float> max() { return $self->max(); }
            
            const std::complex<float>* data() { return $self->data(); }

            void copy(size_t n, std::complex<float> * p) {
                $self->resize(n);
                for(size_t i = 0; i < n; i++) (*$self)(i) = p[i];
            }
            void copy(size_t n, std::complex<double> * p) {
                $self->resize(n);
                for(size_t i = 0; i < n; i++) (*$self)(i) = p[i];
            }
            
            OctopusRowVectorXcf __add__ (const OctopusRowVectorXcf & b) {
                assert($self->size(1) == b.size(1));
                OctopusRowVectorXcf r($self->size(1));
                for(size_t i = 0; i < $self->size(1); i++)
                    r(i) = (*$self)(i) + b(i);
                return r;
            }
            OctopusRowVectorXcf __sub__ (const OctopusRowVectorXcf & b) {
                assert($self->size(1) == b.size(1));
                OctopusRowVectorXcf r($self->size(1));
                for(size_t i = 0; i < $self->size(1); i++)
                    r(i) = (*$self)(i) - b(i);
                return r;
            }
            OctopusRowVectorXcf __mul__ (const OctopusRowVectorXcf & b) {
                assert($self->size(1) == b.size(1));
                OctopusRowVectorXcf r($self->size(1));
                for(size_t i = 0; i < $self->size(1); i++)
                    r(i) = (*$self)(i) * b(i);
                return r;
            }
            OctopusRowVectorXcf __div__ (const OctopusRowVectorXcf & b) {
                assert($self->size(1) == b.size(1));
                OctopusRowVectorXcf r($self->size(1));
                for(size_t i = 0; i < $self->size(1); i++)
                    r(i) = (*$self)(i) / b(i);
                return r;
            } 

            OctopusRowVectorXcf __add__ (const std::complex<float> b) {
                assert($self->size(1) == b.size(1));
                OctopusRowVectorXcf r($self->size(1));
                for(size_t i = 0; i < $self->size(1); i++)
                    r(i) = (*$self)(i) + b;
                return r;
            }
            OctopusRowVectorXcf __sub__ (const std::complex<float> b) {
                assert($self->size(1) == b.size(1));
                OctopusRowVectorXcf r($self->size(1));
                for(size_t i = 0; i < $self->size(1); i++)
                    r(i) = (*$self)(i) - b;
                return r;
            }
            OctopusRowVectorXcf __mul__ (const std::complex<float> b) {
                assert($self->size(1) == b.size(1));
                OctopusRowVectorXcf r($self->size(1));
                for(size_t i = 0; i < $self->size(1); i++)
                    r(i) = (*$self)(i) * b;
                return r;
            }
            OctopusRowVectorXcf __div__ (const std::complex<float> b) {
                assert($self->size(1) == b.size(1));
                OctopusRowVectorXcf r($self->size(1));
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