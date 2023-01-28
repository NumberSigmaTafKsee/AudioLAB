#pragma once

namespace Octopus
{
        struct OctopusColVectorXcf : public FloatComplexColumnVector
    {
        OctopusColVectorXcf() = default;
        OctopusColVectorXcf(const FloatComplexColumnVector &v) : FloatComplexColumnVector(v) {}
        OctopusColVectorXcf(size_t i) : FloatComplexColumnVector(i) {}

        using FloatComplexColumnVector::operator =;
        using FloatComplexColumnVector::operator ();
        using FloatComplexColumnVector::insert;
        //using FloatComplexColumnVector::append;
        using FloatComplexColumnVector::fill;
        using FloatComplexColumnVector::extract;
        using FloatComplexColumnVector::extract_n;
        using FloatComplexColumnVector::transpose;
        using FloatComplexColumnVector::size;
        using FloatComplexColumnVector::min;
        using FloatComplexColumnVector::max;
        using FloatComplexColumnVector::resize;
        using FloatComplexColumnVector::clear;

        OctopusColVectorXcf operator + (const OctopusColVectorXcf & b) {
            OctopusColVectorXcf r;
            ValueList l;
            l(0) = "plus";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).float_complex_column_vector_value();
            return r;
        }
        OctopusColVectorXcf operator - (const OctopusColVectorXcf & b) {
            OctopusColVectorXcf r;
            ValueList l;
            l(0) = "minus";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).float_complex_column_vector_value();
            return r;
        }
        OctopusColVectorXcf operator * (const OctopusColVectorXcf & b) {
            OctopusColVectorXcf r;
            ValueList l;
            l(0) = "times";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).float_complex_column_vector_value();
            return r;
        }
        OctopusColVectorXcf operator / (const OctopusColVectorXcf & b) {
            OctopusColVectorXcf r;
            ValueList l;
            l(0) = "divide";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).float_complex_column_vector_value();
            return r;
        }
        OctopusColVectorXcf operator + (const std::complex<float> b) {
            OctopusColVectorXcf r;
            ValueList l;
            l(0) = "plus";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).float_complex_column_vector_value();
            return r;
        }
        OctopusColVectorXcf operator - (const std::complex<float> b) {
            OctopusColVectorXcf r;
            ValueList l;
            l(0) = "minus";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).float_complex_column_vector_value();
            return r;
        }
        OctopusColVectorXcf operator * (const std::complex<float> b) {
            OctopusColVectorXcf r;
            ValueList l;
            l(0) = "times";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).float_complex_column_vector_value();
            return r;
        }
        OctopusColVectorXcf operator / (const std::complex<float> b) {
            OctopusColVectorXcf r;
            ValueList l;
            l(0) = "divide";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).float_complex_column_vector_value();
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

            OctopusColVectorXcf __add__ (const OctopusColVectorXcf & b) {
                assert($self->size(1) == b.size(1));
                OctopusColVectorXcf r($self->size(1));
                for(size_t i = 0; i < $self->size(1); i++)
                    r(i) = (*$self)(i) + b(i);
                return r;
            }
            OctopusColVectorXcf __sub__ (const OctopusColVectorXcf & b) {
                assert($self->size(1) == b.size(1));
                OctopusColVectorXcf r($self->size(1));
                for(size_t i = 0; i < $self->size(1); i++)
                    r(i) = (*$self)(i) - b(i);
                return r;
            }
            OctopusColVectorXcf __mul__ (const OctopusColVectorXcf & b) {
                assert($self->size(1) == b.size(1));
                OctopusColVectorXcf r($self->size(1));
                for(size_t i = 0; i < $self->size(1); i++)
                    r(i) = (*$self)(i) * b(i);
                return r;
            }
            OctopusColVectorXcf __div__ (const OctopusColVectorXcf & b) {
                assert($self->size(1) == b.size(1));
                OctopusColVectorXcf r($self->size(1));
                for(size_t i = 0; i < $self->size(1); i++)
                    r(i) = (*$self)(i) / b(i);
                return r;
            } 

            OctopusColVectorXcf __add__ (std::complex<float> b) {
                OctopusColVectorXcf r($self->size(1));
                for(size_t i = 0; i < $self->size(1); i++)
                    r(i) = (*$self)(i) + b;
                return r;
            }
            OctopusColVectorXcf __sub__ (std::complex<float> b) {
                OctopusColVectorXcf r($self->size(1));
                for(size_t i = 0; i < $self->size(1); i++)
                    r(i) = (*$self)(i) - b;
                return r;
            }
            OctopusColVectorXcf __mul__ (std::complex<float> b) {
                OctopusColVectorXcf r($self->size(1));
                for(size_t i = 0; i < $self->size(1); i++)
                    r(i) = (*$self)(i) * b;
                return r;
            }
            OctopusColVectorXcf __div__ (std::complex<float> b) {
                OctopusColVectorXcf r($self->size(1));
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