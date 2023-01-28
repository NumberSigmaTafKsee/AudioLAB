#pragma once

namespace Casino::MKL
{
    template<typename T>
    struct Vector : public vector_base<T>
    {                
        using vector_base<T>::size;
        using vector_base<T>::resize;
        using vector_base<T>::data;
        using vector_base<T>::push_back;
        using vector_base<T>::pop_back;
        using vector_base<T>::front;
        using vector_base<T>::back;
        using vector_base<T>::at;
        using vector_base<T>::operator [];
        using vector_base<T>::operator =;

        Vector() = default;
        Vector(size_t n) { assert(n > 0); resize(n); zero(); }
        Vector(const Vector<T> & v) { (*this) = v; }
        Vector(const std::vector<T> & v) {
            resize(v.size());
            memcpy(this->data(),v.data(),v.size()*sizeof(T));
        }
        /*
        Vector(const VectorXf& v) {
            resize(v.size());
            for(size_t i = 0; i < v.size(); i++) (*this)[i] = v[i];
        }
        Vector(const VectorXd& v) {
            resize(v.size());
            for(size_t i = 0; i < v.size(); i++) (*this)[i] = v[i];
        }
        */
        T min() { return *std::min_element(this->begin(),this->end()); }
        T max() { return *std::max_element(this->begin(),this->end()); }

        T& operator()(size_t i) { return (*this)[i]; }
        T operator()(size_t i) const { return (*this)[i]; }

        size_t min_index() { return std::distance(this->begin(), std::min_element(this->begin(),this->end())); }
        size_t max_index() { return std::distance(this->begin(), std::max_element(this->begin(),this->end())); }
                
        Vector<T>& operator +=  (const Vector<T> & v) { 
            cppmkl::vadd(*this,v,*this);
            return *this;
        }
        Vector<T>& operator -=  (const Vector<T> & v) { 
            cppmkl::vsub(*this,v,*this);
            return *this;
        }
        Vector<T>& operator *=  (const Vector<T> & v) { 
            cppmkl::vmul(*this,v,*this);
            return *this;
        }
        Vector<T>& operator /=  (const Vector<T> & v) { 
            cppmkl::vdiv(*this,v,*this);
            return *this;        
	    }
        Vector<T>& operator +=  (const T & x) { 
            Vector<T> v(size());
            v.fill(x);
            cppmkl::vadd(*this,v,*this);
            return *this;
        }
        Vector<T>& operator -=  (const T & x) { 
            Vector<T> v(size());
            v.fill(x);
            cppmkl::vsub(*this,v,*this);
            return *this;
        }
        Vector<T>& operator *=  (const T & x) { 
            Vector<T> v(size());
            v.fill(x);
            cppmkl::vmul(*this,v,*this);
            return *this;
        }
        Vector<T>& operator /=  (const T & x) { 
            Vector<T> v(size());
            v.fill(x);
            cppmkl::vdiv(*this,v,*this);
            return *this;        
	    }

        Vector<T> operator - () {
            Vector<T> r(*this);
            return (T)-1.0 * r;
        }

        Vector<T>& operator = (const T v)
        {
            fill(v);
            return *this;
        }
                
        void zero() {
            memset(data(),0x00,size()*sizeof(T));
        }
        void zeros() { zero(); }
        void fill(T x) {
            for(size_t i = 0; i < size(); i++) (*this)[i] = x;
        }
        void ones() {
            fill((T)1);
        }
        void random(T min = T(0), T max = T(1)) {
            Default noise;
            for(size_t i = 0; i < size(); i++) (*this)[i] = noise.random(min,max);
        }
        void randu(T min = T(0), T max = T(1)) { random(min,max); }
        void randn(T min = T(0), T max = T(1)) { random(min,max); }

        
        void clamp(T min = T(-1), T max = T(1)) {
            for(size_t i = 0; i < size(); i++)
            {
                if((*this)[i] < min) (*this)[i] = min;
                if((*this)[i] < max) (*this)[i] = max;
            }
        }

        void set_size(size_t n) { resize(n); }


        Vector<T> eval() {
            return Vector<T>(*this); 
        }
        Vector<T> slice(size_t start, size_t len) {
            Vector<T> x(len);
            memcpy(x.data(),data()+start,len*sizeof(T));
            return x;
        }

        void print() {
            std::cout << "Vector[" << size() << "]=";
            for(size_t i = 0; i < size()-1; i++) std::cout << (*this)[i] << ",";
            std::cout << (*this)[size()-1] << std::endl;
        }
        
        // eigen compatibility
        void setZero() { zero(); }
        void setOnes() { ones(); }
        void setRandom() { random(); }                



        Vector<T> cwiseMin(T min) {
            Vector<T> r(size());
            for(size_t i = 0; i < size(); i++)
                if((*this)[i] > min) r[i] = min;
                else r[i] = (*this)[i];
            return r;
        }
        Vector<T> cwiseMax(T max) {
            Vector<T> r(size());
            for(size_t i = 0; i < size(); i++)
                if((*this)[i] < max) r[i] = max;
                else r[i] = (*this)[i];
            return r;
        }
        Vector<T> cwiseProd(Vector<T> & x) {
            Vector<T> r(*this);
            r *= x;
            return r;
        }
        Vector<T> cwiseAdd(Vector<T> & x) {
            Vector<T> r(*this);
            r += x;
            return r;
        }
        Vector<T> cwiseSub(Vector<T> & x) {
            Vector<T> r(*this);
            r -= x;
            return r;
        }
        Vector<T> cwiseDiv(Vector<T> & x) {
            Vector<T> r(*this);
            r /= x;
            return r;
        }

        T sum() {        
            return cppmkl::cblas_asum(size(), data(),1);                
        }
        T prod() {
            T p = (T)1.0;
            for(size_t i = 0; i < size(); i++) p *= (*this)[i];
            return p;
        }
        T mean() {
            return sum()/(T)size();
        }
        T geometric_mean() {
            T r = prod();
            return std::pow(r,(T)1.0/(T)size());
        }
        T harmonic_mean() {
            T r = 1.0/sum();
            return size()*std::pow(r,(T)-1.0);
        }
        T stddev() {
            T u = sum();
            T r = 0;
            for(size_t i = 0; i < size(); i++)
                r += std::pow((*this)[i]-u,2.0);
            return r/(T)size();
        }

        #ifdef SWIG
        %extend {
            T __getitem__(size_t i) { return (*$self)[i-1]; }
            void __setitem__(size_t i, T v) { (*$self)[i-1] = v; }
            Vector<T> __add__(const Vector<T> & m) {
                Vector<T> r(*$self);
                r += m;
                return r;
            }
            T* data() { return $self->data(); }
            Vector<T> __add__(const T m) {
                Vector<T> r(*$self);
                r += m;
                return r;
            }
            Vector<T> __sub__(const Vector<T> & m) {
                Vector<T> r(*$self);
                r -= m;
                return r;
            }
            Vector<T> __sub__(const T m) {
                Vector<T> r(*$self);
                r -= m;
                return r;
            }
            Vector<T> __mul__(const Vector<T> & m) {
                Vector<T> r(*$self);
                r *= m;
                return r;
            }
            Vector<T> __mul__(const T m) {
                Vector<T> r(*$self);
                r *= m;
                return r;
            }
            Vector<T> __div__(const Vector<T> & m) {
                Vector<T> r(*$self);
                r /= m;
                return r;
            }
            Vector<T> __div__(const T m) {
                Vector<T> r(*$self);
                r /= m;
                return r;
            }
        }
        #endif
    };

    template<typename T>
    std::ostream& operator << (std::ostream & o, const Vector<T> & v )
    {
        for(size_t i = 0; i < v.size(); i++)
        {               
            o << v[i] << ",";            
        }
        o << std::endl;
        return o;
    }

    /////////////////////////////////////
    // Vector
    /////////////////////////////////////



    template<typename T>
    Vector<T> sqr(Vector<T> & a) {
        Vector<T> r(a.size());                
        cppmkl::vsqr(a,r);
        return r;
    }
    template<typename T>
    Vector<T> abs(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vabs(a,r);
        return r;
    }    
    template<typename T>
    Vector<T> inv(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vinv(a,r);
        return r;            
    }    
    template<typename T>
    Vector<T> sqrt(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vsqrt(a,r);
        return r;
    }
    template<typename T>
    Vector<T> rsqrt(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vinvsqrt(a,r);
        return r;
    }
    template<typename T>
    Vector<T> cbrt(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vcbrt(a,r);
        return r;
    }
    template<typename T>
    Vector<T> rcbrt(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vinvcbrt(a,r);
        return r;
    }
    
    template<typename T>
    Vector<T> pow(const Vector<T> & a,const Vector<T> & b) {
        Vector<T> r(a.size());
        cppmkl::vpow(a,b,r);
        return r;
    }
    template<typename T>
    Vector<T> pow2o3(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vpow2o3(a,r);
        return r;
    }
    template<typename T>
    Vector<T> pow3o2(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vpow3o2(a,r);
        return r;
    }
    template<typename T>
    Vector<T> pow(const Vector<T> & a,const T b) {
        Vector<T> r(a.size());
        cppmkl::vpowx(a,b,r);
        return r;
    }
    template<typename T>
    Vector<T> hypot(const Vector<T> & a,const Vector<T> & b) {
        Vector<T> r(a.size());
        cppmkl::vhypot(a,b,r);
        return r;
    }    
    template<typename T>
    Vector<T> exp(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vexp(a,r);
        return r;
    }
    template<typename T>
    Vector<T> exp2(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vexp2(a,r);
        return r;
    }
    template<typename T>
    Vector<T> exp10(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vexp10(a,r);
        return r;
    }
    template<typename T>
    Vector<T> expm1(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vexpm1(a,r);
        return r;
    }
    template<typename T>
    Vector<T> log(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vln(a,r);
        return r;        
    }    
    template<typename T>
    Vector<T> log10(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vlog10(a,r);
        return r;
    }
    template<typename T>
    Vector<T> log2(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vlog2(a,r);
        return r;
    }
    template<typename T>
    Vector<T> logb(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vlogb(a,r);
        return r;
    }
    template<typename T>
    Vector<T> log1p(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vlog1p(a,r);
        return r;
    }
    
    template<typename T>
    Vector<T> cos(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vcos(a,r);
        return r;
    }
    template<typename T>
    Vector<T> sin(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vsin(a,r);
        return r;            
    }
    template<typename T>
    Vector<T> tan(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vtan(a,r);
        return r;
    }
    template<typename T>
    Vector<T> cosh(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vcosh(a,r);
        return r;
    }
    template<typename T>
    Vector<T> sinh(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vsinh(a,r);
        return r;            
    }
    template<typename T>
    Vector<T> tanh(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vtanh(a,r);
        return r;            
    }
    template<typename T>
    Vector<T> acos(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vacos(a,r);
        return r;
    }
    template<typename T>
    Vector<T> asin(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vasin(a,r);
        return r;            
    }
    template<typename T>
    Vector<T> atan(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vatan(a,r);
        return r;
    }
    template<typename T>
    Vector<T> atan2(const Vector<T> & a,const Vector<T> &n) {
        Vector<T> r(a.size());
        cppmkl::vatan2(a,n,r);
        return r;
    }
    template<typename T>
    Vector<T> acosh(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vacosh(a,r);
        return r;
    }
    template<typename T>
    Vector<T> asinh(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vasinh(a,r);
        return r;            
    }
    template<typename T>
    Vector<T> atanh(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vatanh(a,r);
        return r;
    }        
    template<typename T>
    void sincos(const Vector<T> & a, Vector<T> & b, Vector<T> & r) {        
        cppmkl::vsincos(a,b,r);        
    }
    template<typename T>
    Vector<T> erf(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::verf(a,r);
        return r;
    }
    template<typename T>
    Vector<T> erfinv(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::verfinv(a,r);
        return r;
    }
    template<typename T>
    Vector<T> erfc(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::verfc(a,r);
        return r;
    }
    template<typename T>
    Vector<T> cdfnorm(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vcdfnorm(a,r);
        return r;
    }    
    template<typename T>
    Vector<T> cdfnorminv(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vcdfnorminv(a,r);
        return r;
    }    
    template<typename T>
    Vector<T> floor(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vfloor(a,r);
        return r;
    }    
    template<typename T>
    Vector<T> ceil(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vceil(a,r);
        return r;
    }    
    template<typename T>
    Vector<T> trunc(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vtrunc(a,r);
        return r;
    }    
    template<typename T>
    Vector<T> round(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vround(a,r);
        return r;        
    }    
    template<typename T>
    Vector<T> nearbyint(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vnearbyint(a,r);
        return r;
    }    
    template<typename T>
    Vector<T> rint(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vrint(a,r);
        return r;
    }    
    template<typename T>
    Vector<T> fmod(const Vector<T> & a, Vector<T> & b) {
        Vector<T> r(a.size());
        cppmkl::vmodf(a,b,r);
        return r;
    }    
    
    template<typename T>
    Vector<T> CIS(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vCIS(a,r);
        return r;
    }    
    template<typename T>
    Vector<T> cospi(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vcospi(a,r);
        return r;
    }    
    template<typename T>
    Vector<T> sinpi(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vsinpi(a,r);
        return r;
    }    
    template<typename T>
    Vector<T> tanpi(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vtanpi(a,r);
        return r;
    }    
    template<typename T>
    Vector<T> acospi(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vacospi(a,r);
        return r;
    }    
    template<typename T>
    Vector<T> asinpi(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vasinpi(a,r);
        return r;
    }    
    template<typename T>
    Vector<T> atanpi(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vatanpi(a,r);
        return r;
    }    
    template<typename T>
    Vector<T> atan2pi(const Vector<T> & a, Vector<T> & b) {
        Vector<T> r(a.size());
        cppmkl::vatan2pi(a,b,r);
        return r;
    }    
    template<typename T>
    Vector<T> cosd(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vcosd(a.size(),a.data(),r.data());
        return r;
    }    
    template<typename T>
    Vector<T> sind(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vsind(a.size(),a.data(),r.data());
        return r;
    }    
    template<typename T>
    Vector<T> tand(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vtand(a.size(),a.data(),r.data());
        return r;
    }       
    template<typename T>
    Vector<T> lgamma(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vlgamma(a,r);
        return r;
    }       
    template<typename T>
    Vector<T> tgamma(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vtgamma(a,r);
        return r;
    }       
    template<typename T>
    Vector<T> expint1(const Vector<T> & a) {
        Vector<T> r(a.size());
        cppmkl::vexpint1(a,r);
        return r;
    }       

/////////////////////////////////////
// Vector
/////////////////////////////////////
    template<typename T>
    Vector<T> linspace(size_t n, T start, T inc=1) {
        Vector<T> r(n);
        for(size_t i = 0; i < n; i++) {
            r[i] = start + i*inc;
        }
        return r;
    }
    template<typename T>
    Vector<T> linspace(T start, T end, T inc=1) {
        size_t n = (end - start)/inc;
        Vector<T> r(n);
        for(size_t i = 0; i < n; i++) {
            r[i] = start + i*inc;
        }
        return r;
    }

    template<typename T>
    Vector<T> operator * (T a, const Vector<T> & b) {
        Vector<T> x(b.size());
        x.fill(a);
        Vector<T> r(b.size());
        cppmkl::vmul(x,b,r);            
        return r;
    }
    template<typename T>
    Vector<T> operator / (T a, const Vector<T> & b) {
        Vector<T> x(b.size());
        x.fill(a);
        Vector<T> r(b.size());
        cppmkl::vdiv(x,b,r);            
        return r;
    }
    template<typename T>
    Vector<T> operator + (T a, const Vector<T> & b) {
        Vector<T> x(b.size());
        x.fill(a);
        Vector<T> r(b.size());
        cppmkl::vadd(x,b,r);            
        return r;
    }
    template<typename T>
    Vector<T> operator - (T a, const Vector<T> & b) {
        Vector<T> x(b.size());
        x.fill(a);
        Vector<T> r(b.size());
        cppmkl::vsub(x,b,r);            
        return r;
    }

    template<typename T>
    Vector<T> operator * (const Vector<T> & a, T b) {
        Vector<T> x(a.size());
        x.fill(b);
        Vector<T> r(a.size());
        cppmkl::vmul(a,x,r);            
        return r;
    }
    template<typename T>
    Vector<T> operator / (const Vector<T> & a , T b) {
        Vector<T> x(a.size());
        x.fill(b);
        Vector<T> r(a.size());
        cppmkl::vdiv(a,x,r);            
        return r;
    }
    template<typename T>
    Vector<T> operator + (const Vector<T> & a, T b) {
        Vector<T> x(a.size());
        x.fill(b);
        Vector<T> r(a.size());
        cppmkl::vadd(a,x,r);            
        return r;
    }
    template<typename T>
    Vector<T> operator - (const Vector<T> & a, T b) {
        Vector<T> x(a.size());
        x.fill(b);
        Vector<T> r(a.size());
        cppmkl::vsub(a,x,r);            
        return r;
    }
    
    template<typename T>
    Vector<T> operator - (const Vector<T> & a, const Vector<T> & b) {        
        Vector<T> r(a.size());
        cppmkl::vsub(a,b,r);            
        return r;
    }
    template<typename T>
    Vector<T> operator + (const Vector<T> & a, const Vector<T> & b) {        
        Vector<T> r(a.size());
        cppmkl::vadd(a,b,r);                    
        return r;
    }
    template<typename T>
    Vector<T> operator * (const Vector<T> & a, const Vector<T> & b) {        
        Vector<T> r(a.size());
        cppmkl::vmul(a,b,r);            
        return r;
    }
    template<typename T>
    Vector<T> operator / (const Vector<T> & a, const Vector<T> & b) {        
        Vector<T> r(a.size());
        cppmkl::vdiv(a,b,r);            
        return r;
    }

}