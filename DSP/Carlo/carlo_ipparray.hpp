#pragma once

namespace Casino::IPP
{
    template<typename T>
    struct IPPArray
    {
    private:
        std::shared_ptr<T> ptr;
        T * array;
        size_t   len;
        int      r,w;
    public:
    
        IPPArray(const IPPArray<T> & a) {
            *this = a;
        }
        virtual ~IPPArray() = default;

        IPPArray(size_t n) {
            array = Malloc<T>(len = n);
            ptr = std::shared_ptr<T>(array,[](T* p) { Free<T>(p); });
            assert(array != NULL);
            Zero<T>(array,n);
            r = 0;
            w = 0;
        }
                
        
        void resize(size_t n) {
            T * p = Malloc<T>(n);
            Move<T>(array,p,n);
            Free<T>(array);
            array  = p;
            len    = n;
        }
        void fill(T value) {
            if(array == NULL) return;
            Set<T>(value,array,len);
        }
        T sum() {
            T r = 0;
            Sum<T>(array,len,&r);
            return r;
        }
        
        T& operator[] (size_t i) { return array[i]; }

        T      __getitem__(size_t i) { return array[i]; }
        void   __setitem__(size_t i, T v) { array[i] = v; }

        IPPArray<T>& operator = (const IPPArray & x) {
            ptr.reset();
            ptr = x.ptr;
            array = x.array;
            len = x.len;
            return *this;
        }
        
        void ring_push(const T& value) {
            array[w++] = value;
            w = w % len;
        }
        T ring_pop() {
            T v = array[r++];
            r = r % len;
            return v;
        }
        T ring_linear_pop() {
            T v1 = array[r];
            r = (r+1) % len;
            T v2 = array[r];            
            T frac = v1 - std::floor(v1);
            return v1 + frac*(v2-v1);
        }

        IPPArray<T> operator + (const T& value) {
            IPPArray<T> r(*this);
            AddC<T>(array,value,r.array,len);
            return r;
        }
        IPPArray<T> operator + (const IPPArray<T> & b) {
            IPPArray<T> r(*this);
            assert(len == b.len);
            Add<T>(array,b.array,r.array,len);
            return r;
        }
        IPPArray<T> operator - (const T& value) {
            IPPArray<T> r(*this);
            SubC<T>(array,value,r.array,len);
            return r;
        }
        IPPArray<T> operator - (const IPPArray<T> & b) {
            IPPArray<T> r(*this);
            assert(len == b.len);
            Sub<T>(array,b.array,r.array,len);
            return r;
        }        
        IPPArray<T> operator * (const T& value) {
            IPPArray<T> r(*this);
            MulC<T>(array,value,r.array,len);
            return r;
        }        
        IPPArray<T> operator * (const IPPArray<T> & b) {
            IPPArray<T> r(*this);
            assert(len == b.len);
            Mul<T>(array,b.array,r.array,len);
            return r;
        }
        IPPArray<T> operator / (const T& value) {
            IPPArray<T> r(*this);
            DivC<T>(array,value,r.array,len);
            return r;
        }        
        IPPArray<T> operator / (const IPPArray<T>& b) {
            IPPArray<T> r(*this);
            assert(len == b.len);
            Div<T>(array,b.array,r.array,len);
            return r;
        }        

        void print() {
            std::cout << "Array[" << len << "]=";
            for(size_t i = 0; i < len-1; i++) std::cout << array[i] << ",";
            std::cout << array[len-1] << std::endl;
        }
        IPPArray<T>& copy(const IPPArray<T> & a) {
            ptr.reset();
            array = Malloc<T>(a.len);
            memcpy(array,a.array,a.len*sizeof(T));
            return *this; 
       }        
    };
    
    template<typename T>
    void move(IPPArray<T> & dst, const IPPArray<T> & src)
    {
        if(src.len <= dst.len)
            Move(src.array,dst.array,src.len);
        else
            Move(src.array,dst.array,dst.len);
    }

    template<typename T>
    void copy(IPPArray<T> & dst, const IPPArray<T> & src)
    {
        if(src.len <= dst.len)
            Move(src.array,dst.array,src.len);
        else
            Move(src.array,dst.array,dst.len);
    }

    template<typename T>
    void fill(IPPArray<T> & a, const T val)
    {
        Set(val,a.array,a.len);
    }

    template<typename T>
    void zero(IPPArray<T> & a) {
        Zero(a.array, a.len);
    }

    template<typename T>
    void sinc(const IPPArray<T> & src, IPPArray<T> & dst) {
        assert(src.len == dst.len);
        if(src.len <= dst.len)
            Sinc(src.array,dst.array,src.len);
        else
            Sinc(src.array,dst.array,dst.len);
    }

    template<typename T1, typename T2>
    void complex(const IPPArray<T1> & real, const IPPArray<T1> & imag, IPPArray<T2> & dst)
    {
        assert(real.len == imag.len && real.len == dst.len);
        RealToComplex(real.array,imag.array,dst.array,dst.len);
    }
    template<typename T1, typename T2>
    void real(const IPPArray<T2> & src, IPPArray<T1> & real, IPPArray<T1> & imag)
    {
        assert(real.len == imag.len && real.len == src.len);
        ComplexToReal(src.array,real.array,imag.array,src.len);
    }
    template<typename T>
    void magnitude(const IPPArray<T> & real, const IPPArray<T> & imag, IPPArray<T> & dst)
    {
        assert(real.len == imag.len && real.len == dst.len);
        Magnitude(real.array,imag.array,dst.array,dst.len);        
    }
    template<typename T>
    void phase(const IPPArray<T> & real, const IPPArray<T> & imag, IPPArray<T> & dst)
    {
        assert(real.len == imag.len && real.len == dst.len);
        Phase(real.array,imag.array,dst.array,dst.len);
    }
    template<typename T>
    void split_complex(const std::complex<T> & in, IPPArray<T> & real, IPPArray<T> & imag)
    {
        assert(real.len == imag.len && real.len == in.len);
        for(size_t i = 0; i < in.size(); i++)
        {
            real[i] = in[i].real();
            imag[i] = in[i].imag();
        }
    }
    template<typename T>
    void merge_complex(const IPPArray<T> & real, const IPPArray<T> & imag, IPPArray<T> & dst)
    {
        assert(real.len == imag.len && real.len == dst.len);
        for(size_t i = 0; i < dst.size(); i++)
        {
            dst[i].real(real[i]);
            dst[i].imag(imag[i]);
        }
    }
    template<typename T>
    void cart2polar(const IPPArray<T> & real, const IPPArray<T> & imag, IPPArray<T> & mag, IPPArray<T> & phase)
    {
        assert(real.len == imag.len && real.len == mag.len);
        assert(real.len == mag.len && real.len == phase.len);
        Cart2Polar(real.array,imag.array,mag.array,phase.array,real.len);
    }   
    template<typename T>
    void polar2cart(const IPPArray<T> & mag, const IPPArray<T> & phase, IPPArray<T> & real, IPPArray<T> & imag)
    {
        assert(real.len == imag.len && real.len == mag.len);
        assert(real.len == mag.len && real.len == phase.len);
        Polar2Cart(mag.array,phase.array,real.array,imag.array,mag.len);
    }


}