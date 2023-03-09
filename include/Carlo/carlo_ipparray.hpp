#pragma once

#include <iterator>

namespace Casino::IPP
{
    template<typename T>
    struct IPPArray
    {
    private:
        std::shared_ptr<T> ptr;
        T * array;
        size_t   len;
        
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
        }
         
        T* data() { return array; }
        size_t size() const { return len; }
                       
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
        
        IPPArray<T>& copy(const IPPArray<T> & a) {
            ptr.reset();
            array = Malloc<T>(a.len);
            memcpy(array,a.array,a.len*sizeof(T));
            return *this; 
       }        
       
       void print() {
            std::cout << "Array[" << len << "]=";
            for(size_t i = 0; i < len-1; i++) std::cout << array[i] << ",";
            std::cout << array[len-1] << std::endl;
        }

		T* begin() { return array; }
		T* end() { return array+len; }
		
    };

	template<typename T>
	struct RingBuffer
	{
		IPPArray<T> array;
		int      r,w;
		
		RingBuffer(size_t n) : array(n) {
			r = w = 0;
		}
		void setDelay(size_t n) {
			w += n % array.size();
		}
		
		void push(const T& value) {
            array[w++] = value;
            w = w % array.size();
        }
        T pop() {
            T v = array[r++];
            r = r % array.size();
            return v;
        }
        T get() {
			return array[r];
		}
		void set(const T& v) {
			array[w] = v;
		} 
        T lipop() {
            T v1 = array[r];
            r = (r+1) % array.size();
            T v2 = array[r];            
            T frac = v1 - std::floor(v1);
            return v1 + frac*(v2-v1);
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
    void merge_complex(const IPPArray<T> & real, const IPPArray<std::complex<T>> & imag, IPPArray<T> & dst)
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

	template<typename T>
	std::ostream& operator << (std::ostream & o, IPPArray<T> & a)
	{		
		for(size_t i = 0; i < a.size(); i++) {
			o << a[i] << ",";
		}		
		return o;
	}
}
