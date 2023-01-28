#pragma once

namespace Octopus
{
    struct OctopusRowVectorXf : public FloatRowVector
    {
        OctopusRowVectorXf() : FloatRowVector() {}
        OctopusRowVectorXf(const FloatRowVector & v) : FloatRowVector(v) {}
        OctopusRowVectorXf(size_t i) : FloatRowVector(i) {}

        using FloatRowVector::operator =;
        using FloatRowVector::operator ();
        using FloatRowVector::insert;
        using FloatRowVector::append;
        using FloatRowVector::fill;
        using FloatRowVector::extract;
        using FloatRowVector::extract_n;
        using FloatRowVector::transpose;
        using FloatRowVector::size;
        using FloatRowVector::min;
        using FloatRowVector::max;
        using FloatRowVector::resize;
        using FloatRowVector::clear;

        void print()
        {
            ValueList l;
            l(0) = *this;
            octave::feval("display",l,0);
        }
        OctopusRowVectorXf operator + (const OctopusRowVectorXf & b) {
            OctopusRowVectorXf r;
            ValueList l;
            l(0) = "plus";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).float_row_vector_value();
            return r;
        }
        OctopusRowVectorXf operator - (const OctopusRowVectorXf & b) {
            OctopusRowVectorXf r;
            ValueList l;
            l(0) = "minus";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).float_row_vector_value();
            return r;
        }
        OctopusRowVectorXf operator * (const OctopusRowVectorXf & b) {
            OctopusRowVectorXf r;
            ValueList l;
            l(0) = "times";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).float_row_vector_value();
            return r;
        }
        OctopusRowVectorXf operator / (const OctopusRowVectorXf & b) {
            OctopusRowVectorXf r;
            ValueList l;
            l(0) = "divide";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).float_row_vector_value();
            return r;
        }
        OctopusRowVectorXf operator + (const float b) {
            OctopusRowVectorXf r;
            ValueList l;
            l(0) = "plus";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).float_row_vector_value();
            return r;
        }
        OctopusRowVectorXf operator - (const float b) {
            OctopusRowVectorXf r;
            ValueList l;
            l(0) = "minus";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).float_row_vector_value();
            return r;
        }
        OctopusRowVectorXf operator * (const float b) {
            OctopusRowVectorXf r;
            ValueList l;
            l(0) = "times";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).float_row_vector_value();
            return r;
        }
        OctopusRowVectorXf operator / (const float b) {
            OctopusRowVectorXf r;
            ValueList l;
            l(0) = "divide";
            l(1) = *this;
            l(2) = b;
            l = octave::feval("cwiseops",l,1);
            r = l(0).float_row_vector_value();
            return r;
        }        

        #ifdef SWIG
        %extend
        {
            float __getitem__(size_t i) { return (*$self)(i-1); }
            void __setitem__(size_t i, float v) { (*$self)(i-1) = v; }
            
            size_t size() const { return $self->size(1); }
            
            void fill(float v) { $self->fill(v); }
            
            //float min() { return $self->min(); }
            //float max() { return $self->max(); }
            
            const float* data() { return $self->data(); }
            
            void copy(size_t n, float * p) {
                $self->resize(n);
                for(size_t i = 0; i < n; i++) (*$self)(i) = p[i];
            }
            void copy(size_t n, double * p) {
                $self->resize(n);
                for(size_t i = 0; i < n; i++) (*$self)(i) = p[i];
            }

            OctopusRowVectorXf __add__ (const OctopusRowVectorXf & b) {
                assert($self->size(1) == b.size(1));
                OctopusRowVectorXf r($self->size(1));
                for(size_t i = 0; i < $self->size(1); i++)
                    r(i) = (*$self)(i) + b(i);
                return r;
            }
            OctopusRowVectorXf __sub__ (const OctopusRowVectorXf & b) {
                assert($self->size(1) == b.size(1));
                OctopusRowVectorXf r($self->size(1));
                for(size_t i = 0; i < $self->size(1); i++)
                    r(i) = (*$self)(i) - b(i);
                return r;
            }
            OctopusRowVectorXf __mul__ (const OctopusRowVectorXf & b) {
                assert($self->size(1) == b.size(1));
                OctopusRowVectorXf r($self->size(1));
                for(size_t i = 0; i < $self->size(1); i++)
                    r(i) = (*$self)(i) * b(i);
                return r;
            }
            OctopusRowVectorXcf __div__ (const OctopusRowVectorXf & b) {
                assert($self->size(1) == b.size(1));
                OctopusRowVectorXf r($self->size(1));
                for(size_t i = 0; i < $self->size(1); i++)
                    r(i) = (*$self)(i) / b(i);
                return r;
            } 

            OctopusRowVectorXf __add__ (const float b) {
                assert($self->size(1) == b.size(1));
                OctopusRowVectorXf r($self->size(1));
                for(size_t i = 0; i < $self->size(1); i++)
                    r(i) = (*$self)(i) + b;
                return r;
            }
            OctopusRowVectorXf __sub__ (const float b) {
                assert($self->size(1) == b.size(1));
                OctopusRowVectorXf r($self->size(1));
                for(size_t i = 0; i < $self->size(1); i++)
                    r(i) = (*$self)(i) - b;
                return r;
            }
            OctopusRowVectorXf __mul__ (const float b) {
                assert($self->size(1) == b.size(1));
                OctopusRowVectorXf r($self->size(1));
                for(size_t i = 0; i < $self->size(1); i++)
                    r(i) = (*$self)(i) * b;
                return r;
            }
            OctopusRowVectorXcf __div__ (const float b) {
                assert($self->size(1) == b.size(1));
                OctopusRowVectorXf r($self->size(1));
                for(size_t i = 0; i < $self->size(1); i++)
                    r(i) = (*$self)(i) / b;
                return r;
            }             
        }
        #endif
    };
}