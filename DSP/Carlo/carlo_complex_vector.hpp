#pragma once

namespace Casino::MKL
{
    template<typename T>
    struct ComplexVector : public vector_base<std::complex<T>>
    {                
        using vecbase = vector_base<std::complex<T>>;
        

        using vecbase::size;
        using vecbase::resize;
        using vecbase::data;
        using vecbase::push_back;
        using vecbase::pop_back;
        using vecbase::front;
        using vecbase::back;
        using vecbase::at;
        using vecbase::operator [];
        using vecbase::operator =;

        ComplexVector() = default;
        ComplexVector(size_t i) : vecbase(i) {}
        ComplexVector(const vecbase& v) : vecbase(v) {}
        ComplexVector(const ComplexVector<T> & v) : vecbase(v) {}

        void fill(const std::complex<T>& c) {
            for(size_t i = 0; i < size(); i++) (*this)[i] = c;
        }

        ComplexVector<T>& operator +=  (const ComplexVector<T> & v) { 
            cppmkl::vadd(*this,v,*this);
            return *this;
        }
        ComplexVector<T>& operator -=  (const ComplexVector<T> & v) { 
            cppmkl::vsub(*this,v,*this);
            return *this;
        }
        ComplexVector<T>& operator *=  (const ComplexVector<T> & v) { 
            cppmkl::vmul(*this,v,*this);
            return *this;
        }
        ComplexVector<T>& operator /=  (const ComplexVector<T> & v) { 
            cppmkl::vdiv(*this,v,*this);
            return *this;        
	    }

        ComplexVector<T>& operator +=  (const std::complex<T> & x) { 
            ComplexVector<T> v(size());
            v.fill(x);
            cppmkl::vadd(*this,v,*this);
            return *this;
        }
        ComplexVector<T>& operator -=  (const std::complex<T> & x) { 
            ComplexVector<T> v(size());
            v.fill(x);
            cppmkl::vsub(*this,v,*this);
            return *this;
        }
        ComplexVector<T>& operator *=  (const std::complex<T> & x) { 
            ComplexVector<T> v(size());
            v.fill(x);
            cppmkl::vmul(*this,v,*this);
            return *this;
        }
        ComplexVector<T>& operator /=  (const std::complex<T> & x) { 
            ComplexVector<T> v(size());
            v.fill(x);
            cppmkl::vdiv(*this,v,*this);
            return *this;        
	    }

        ComplexVector<T> operator - () {
            ComplexVector<T> r(*this);
            for(size_t i = 0; i < size(); i++) r[i] = std::conj(r[i]);
            return r;
        }
        ComplexVector<T>& operator = (const std::complex<T>& v)
        {
            fill(v);
            return *this;
        }

        Vector<T> real() {
            Vector<T> r(size());
            for(size_t i = 0; i < size(); i++) r[i] = (*this)[i].real();
            return r;
        }
        Vector<T> imag() {
            Vector<T> r(size());
            for(size_t i = 0; i < size(); i++) r[i] = (*this)[i].imag();
            return r;
        }
        void real(const Vector<T> & r) {
            for(size_t i = 0; i < size(); i++) (*this)[i].real(r[i]);
        }
        void imag(const Vector<T> & r) {
            for(size_t i = 0; i < size(); i++) (*this)[i].imag(r[i]);
        }

        Vector<T> abs() {
            Vector<T> r(size());
            for(size_t i = 0; i < size(); i++) r[i] = std::abs((*this)[i]);
            return r;
        }
        Vector<T> arg() {
            Vector<T> r(size());
            for(size_t i = 0; i < size(); i++) r[i] = std::abs((*this)[i]);
            return r;
        }
        ComplexVector<T> conj() {
            ComplexVector<T> r(size());
            for(size_t i = 0; i < size(); i++) r[i] = std::conj((*this)[i]);
            return r;
        }
        ComplexVector<T> proj() {
            ComplexVector<T> r(size());
            for(size_t i = 0; i < size(); i++) r[i] = std::proj((*this)[i]);
            return r;
        }
        Vector<T> norm() {
            Vector<T> r(size());
            for(size_t i = 0; i < size(); i++) r[i] = std::norm((*this)[i]);
            return r;
        }

        void print() {
            std::cout << "Vector[" << size() << "]=";
            for(size_t i = 0; i < size()-1; i++) std::cout << (*this)[i] << ",";
            std::cout << (*this)[size()-1] << std::endl;
        }        

        #ifdef SWIG
        %extend {
            size_t size() const { return $self->size(); }
            void resize(size_t i) { $self->resize(i); }
            std::complex<T>* data() { return $self->data(); }
            void push_back(const std::complex<T> v) { $self->push_back(v); }
            void pop_back() { $self->pop_back(); }
            std::complex<T> front() { return $self->front(); }
            std::complex<T> back() { return $self->back(); }
            std::complex<T> at(size_t i) { return $self->at(i); }

            std::complex<T> __getitem__(size_t i) { return (*$self)[i-1]; }
            void __setitem__(size_t i, std::complex<T> v) { (*$self)[i-1] = v; }
            std::complex<T>* data() { return $self->data(); }
            ComplexVector<T> __add__(const ComplexVector<T> & m) {
                ComplexVector<T> r(*$self);
                r += m;
                return r;
            }
            ComplexVector<T> __add__(const std::complex<T> & m) {
                ComplexVector<T> r(*$self);
                r += m;
                return r;
            }
            ComplexVector<T> __sub__(const ComplexVector<T> & m) {
                ComplexVector<T> r(*$self);
                r -= m;
                return r;
            }
            ComplexVector<T> __sub__(const std::complex<T> & m) {
                ComplexVector<T> r(*$self);
                r -= m;
                return r;
            }
            ComplexVector<T> __mul__(const ComplexVector<T> & m) {
                ComplexVector<T> r(*$self);
                r *= m;
                return r;
            }
            ComplexVector<T> __mul__(const std::complex<T> & m) {
                ComplexVector<T> r(*$self);
                r *= m;
                return r;
            }
            ComplexVector<T> __div__(const ComplexVector<T> & m) {
                ComplexVector<T> r(*$self);
                r /= m;
                return r;
            }
            ComplexVector<T> __div__(const std::complex<T> & m) {
                ComplexVector<T> r(*$self);
                r /= m;
                return r;
            }
        }
        #endif
    };

    /////////////////////////////////////
// ComplexVector
/////////////////////////////////////
    template<typename T>
    Vector<T> abs(const ComplexVector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vabs(a,r);
        return r;
    }
    template<typename T>
    ComplexVector<T> sqrt(const ComplexVector<T> & a) {
        ComplexVector<T> r(a.size());
        cppmkl::vsqrt(a,r);
        return r;
    }
    template<typename T>
    ComplexVector<T> pow(const ComplexVector<T> & a,const ComplexVector<T> & b) {
        ComplexVector<T> r(a.size());
        cppmkl::vpow(a,b,r);
        return r;
    }
    template<typename T>
    ComplexVector<T> exp(const ComplexVector<T> & a) {
        ComplexVector<T> r(a.size());
        cppmkl::vexp(a,r);
        return r;
    }    
    template<typename T>
    ComplexVector<T> log(const ComplexVector<T> & a) {
        ComplexVector<T> r(a.size());
        cppmkl::vln(a,r);
        return r;        
    }
    template<typename T>
    ComplexVector<T> cos(const ComplexVector<T> & a) {
        ComplexVector<T> r(a.size());
        cppmkl::vcos(a,r);
        return r;
    }
    template<typename T>
    ComplexVector<T> sin(const ComplexVector<T> & a) {
        ComplexVector<T> r(a.size());
        cppmkl::vsin(a,r);
        return r;            
    }
    template<typename T>
    ComplexVector<T> tan(const ComplexVector<T> & a) {
        ComplexVector<T> r(a.size());
        cppmkl::vtan(a,r);
        return r;
    }
    template<typename T>
    ComplexVector<T> cosh(const ComplexVector<T> & a) {
        ComplexVector<T> r(a.size());
        cppmkl::vcosh(a,r);
        return r;
    }
    template<typename T>
    ComplexVector<T> sinh(const ComplexVector<T> & a) {
        ComplexVector<T> r(a.size());
        cppmkl::vsinh(a,r);
        return r;            
    }
    template<typename T>
    ComplexVector<T> tanh(const ComplexVector<T> & a) {
        ComplexVector<T> r(a.size());
        cppmkl::vtanh(a,r);
        return r;            
    }
    template<typename T>
    ComplexVector<T> acos(const ComplexVector<T> & a) {
        ComplexVector<T> r(a.size());
        cppmkl::vacos(a,r);
        return r;
    }
    template<typename T>
    ComplexVector<T> asin(const ComplexVector<T> & a) {
        ComplexVector<T> r(a.size());
        cppmkl::vasin(a,r);
        return r;            
    }
    template<typename T>
    ComplexVector<T> atan(const ComplexVector<T> & a) {
        ComplexVector<T> r(a.size());
        cppmkl::vatan(a,r);
        return r;
    }
    template<typename T>
    ComplexVector<T> acosh(const ComplexVector<T> & a) {
        ComplexVector<T> r(a.size());
        cppmkl::vacosh(a,r);
        return r;
    }
    template<typename T>
    ComplexVector<T> asinh(const ComplexVector<T> & a) {
        ComplexVector<T> r(a.size());
        cppmkl::vasinh(a,r);
        return r;            
    }
    template<typename T>
    ComplexVector<T> atanh(const ComplexVector<T> & a) {
        ComplexVector<T> r(a.size());
        cppmkl::vatanh(a,r);
        return r;
    }        

/////////////////////////////////////
// ComplexVector
/////////////////////////////////////
    template<typename T>
    ComplexVector<T> operator * (std::complex<T> a, const ComplexVector<T> & b) {
        ComplexVector<T> x(b.size());
        x.fill(a);
        ComplexVector<T> r(b.size());
        cppmkl::vmul(x,b,r);            
        return r;
    }
    template<typename T>
    ComplexVector<T> operator / (std::complex<T> a, const ComplexVector<T> & b) {
        ComplexVector<T> x(b.size());
        x.fill(a);
        ComplexVector<T> r(b.size());
        cppmkl::vdiv(x,b,r);            
        return r;
    }
    template<typename T>
    ComplexVector<T> operator + (std::complex<T> a, const ComplexVector<T> & b) {
        ComplexVector<T> x(b.size());
        x.fill(a);
        ComplexVector<T> r(b.size());
        cppmkl::vadd(x,b,r);            
        return r;
    }
    template<typename T>
    ComplexVector<T> operator - (std::complex<T> a, const ComplexVector<T> & b) {
        ComplexVector<T> x(b.size());
        x.fill(a);
        ComplexVector<T> r(b.size());
        cppmkl::vsub(x,b,r);            
        return r;
    }

    template<typename T>
    ComplexVector<T> operator * (const ComplexVector<T> & a, std::complex<T> b) {
        ComplexVector<T> x(a.size());
        x.fill(b);
        ComplexVector<T> r(a.size());
        cppmkl::vmul(a,x,r);            
        return r;
    }
    template<typename T>
    ComplexVector<T> operator / (const ComplexVector<T> & a , std::complex<T> b) {
        ComplexVector<T> x(a.size());
        x.fill(b);
        ComplexVector<T> r(a.size());
        cppmkl::vdiv(a,x,r);            
        return r;
    }
    template<typename T>
    ComplexVector<T> operator + (const ComplexVector<T> & a, std::complex<T> b) {
        ComplexVector<T> x(a.size());
        x.fill(b);
        ComplexVector<T> r(a.size());
        cppmkl::vadd(a,x,r);            
        return r;
    }
    template<typename T>
    ComplexVector<T> operator - (const ComplexVector<T> & a, std::complex<T> b) {
        ComplexVector<T> x(a.size());
        x.fill(b);
        ComplexVector<T> r(a.size());
        cppmkl::vsub(a,x,r);            
        return r;
    }
    
    template<typename T>
    ComplexVector<T> operator - (const ComplexVector<T> & a, const ComplexVector<T> & b) {        
        ComplexVector<T> r(a.size());
        cppmkl::vsub(a,b,r);            
        return r;
    }
    template<typename T>
    ComplexVector<T> operator + (const ComplexVector<T> & a, const ComplexVector<T> & b) {        
        ComplexVector<T> r(a.size());
        cppmkl::vadd(a,b,r);                    
        return r;
    }
    template<typename T>
    ComplexVector<T> operator * (const ComplexVector<T> & a, const ComplexVector<T> & b) {        
        ComplexVector<T> r(a.size());
        cppmkl::vmul(a,b,r);            
        return r;
    }
    template<typename T>
    ComplexVector<T> operator / (const ComplexVector<T> & a, const ComplexVector<T> & b) {        
        ComplexVector<T> r(a.size());
        cppmkl::vdiv(a,b,r);            
        return r;
    }

}