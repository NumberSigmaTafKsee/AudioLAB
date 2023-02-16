#pragma once

namespace Octopus
{
        struct OctopusColVectorXf : public FloatColumnVector
    {
        OctopusColVectorXf() = default;
        OctopusColVectorXf(const FloatColumnVector &v) : FloatColumnVector(v) {}
        OctopusColVectorXf(size_t i) : FloatColumnVector(i) {}

        using FloatColumnVector::operator =;
        using FloatColumnVector::operator ();
        using FloatColumnVector::insert;
        //using FloatColumnVector::append;
        using FloatColumnVector::fill;
        using FloatColumnVector::extract;
        using FloatColumnVector::extract_n;
        using FloatColumnVector::transpose;
        using FloatColumnVector::size;
        using FloatColumnVector::min;
        using FloatColumnVector::max;
        using FloatColumnVector::resize;
        using FloatColumnVector::clear;

        OctopusColVectorXf operator + (const OctopusColVectorXf & b) {
            OctopusColVectorXf r;
            ValueList l;
            l(0) = "plus";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).float_column_vector_value();
            return r;
        }
        OctopusColVectorXf operator - (const OctopusColVectorXf & b) {
            OctopusColVectorXf r;
            ValueList l;
            l(0) = "minus";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).float_column_vector_value();
            return r;
        }
        OctopusColVectorXf operator * (const OctopusColVectorXf & b) {
            OctopusColVectorXf r;
            ValueList l;
            l(0) = "times";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).float_column_vector_value();
            return r;
        }
        OctopusColVectorXf operator / (const OctopusColVectorXf & b) {
            OctopusColVectorXf r;
            ValueList l;
            l(0) = "divide";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).float_column_vector_value();
            return r;
        }
        OctopusColVectorXf operator + (const float b) {
            OctopusColVectorXf r;
            ValueList l;
            l(0) = "plus";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).float_column_vector_value();
            return r;
        }
        OctopusColVectorXf operator - (const float b) {
            OctopusColVectorXf r;
            ValueList l;
            l(0) = "minus";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).float_column_vector_value();
            return r;
        }
        OctopusColVectorXf operator * (const float b) {
            OctopusColVectorXf r;
            ValueList l;
            l(0) = "times";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).float_column_vector_value();
            return r;
        }
        OctopusColVectorXf operator / (const float b) {
            OctopusColVectorXf r;
            ValueList l;
            l(0) = "divide";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).float_column_vector_value();
            return r;
        }   

        #ifdef SWIG
        %extend
        {
            float __getitem__(size_t i) { return (*$self)(i); }
            void __setitem__(size_t i, float v) { (*$self)(i) = v; }
            size_t size() const { return $self->size(1); }

            void fill(float v) { $self->fill(v); }

            // float min() { return $self->min(); }
            // float max() { return $self->max(); }

            const float* data() { return $self->data(); }

            void copy(size_t n, float * p) {
                $self->resize(n);
                for(size_t i = 0; i < n; i++) (*$self)(i) = p[i];
            }
            void copy(size_t n, double * p) {
                $self->resize(n);
                for(size_t i = 0; i < n; i++) (*$self)(i) = p[i];
            }

            OctopusColVectorXf __add__ (const OctopusColVectorXf & b) {
                assert($self->size(1) == b.size(1));
                OctopusColVectorXf r($self->size(1));
                for(size_t i = 0; i < $self->size(1); i++)
                    r(i) = (*$self)(i) + b(i);
                return r;
            }
            OctopusColVectorXf __sub__ (const OctopusColVectorXf & b) {
                assert($self->size(1) == b.size(1));
                OctopusColVectorXf r($self->size(1));
                for(size_t i = 0; i < $self->size(1); i++)
                    r(i) = (*$self)(i) - b(i);
                return r;
            }
            OctopusColVectorXf __mul__ (const OctopusColVectorXf & b) {
                assert($self->size(1) == b.size(1));
                OctopusColVectorXf r($self->size(1));
                for(size_t i = 0; i < $self->size(1); i++)
                    r(i) = (*$self)(i) * b(i);
                return r;
            }
            OctopusColVectorXf __div__ (const OctopusColVectorXf & b) {
                assert($self->size(1) == b.size(1));
                OctopusColVectorXf r($self->size(1));
                for(size_t i = 0; i < $self->size(1); i++)
                    r(i) = (*$self)(i) / b(i);
                return r;
            } 

            OctopusColVectorXf __add__ (float b) {
                OctopusColVectorXf r($self->size(1));
                for(size_t i = 0; i < $self->size(1); i++)
                    r(i) = (*$self)(i) + b;
                return r;
            }
            OctopusColVectorXf __sub__ (float b) {
                OctopusColVectorXf r($self->size(1));
                for(size_t i = 0; i < $self->size(1); i++)
                    r(i) = (*$self)(i) - b;
                return r;
            }
            OctopusColVectorXf __mul__ (float b) {
                OctopusColVectorXf r($self->size(1));
                for(size_t i = 0; i < $self->size(1); i++)
                    r(i) = (*$self)(i) * b;
                return r;
            }
            OctopusColVectorXf __div__ (float b) {
                OctopusColVectorXf r($self->size(1));
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