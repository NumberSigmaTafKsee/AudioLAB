#pragma once

namespace Casino
{
    template<typename T>
    // r = frac
    // x = [i]
    // y = [i+1]
    T linear_interpolate(T x, T y, T r)
    {        
        return x + r*(y-x);
    }
    template<typename T>
    T cubic_interpolate(T finpos, T xm1, T x0, T x1, T x2)
    {
        //T xm1 = x [inpos - 1];
        //T x0  = x [inpos + 0];
        //T x1  = x [inpos + 1];
        //T x2  = x [inpos + 2];
        T a = (3 * (x0-x1) - xm1 + x2) / 2;
        T b = 2*x1 + xm1 - (5*x0 + x2) / 2;
        T c = (x1 - xm1) / 2;
        return (((a * finpos) + b) * finpos + c) * finpos + x0;
    }
    // original
    template<typename T>
    // x = frac
    // y0 = [i-1]
    // y1 = [i]
    // y2 = [i+1]
    // y3 = [i+2]
    T hermite1(T x, T y0, T y1, T y2, T y3)
    {
        // 4-point, 3rd-order Hermite (x-form)
        T c0 = y1;
        T c1 = 0.5f * (y2 - y0);
        T c2 = y0 - 2.5f * y1 + 2.f * y2 - 0.5f * y3;
        T c3 = 1.5f * (y1 - y2) + 0.5f * (y3 - y0);
        return ((c3 * x + c2) * x + c1) * x + c0;
    }

    // james mccartney
    template<typename T>
    // x = frac
    // y0 = [i-1]
    // y1 = [i]
    // y2 = [i+1]
    // y3 = [i+2]
    T hermite2(T x, T y0, T y1, T y2, T y3)
    {
        // 4-point, 3rd-order Hermite (x-form)
        T c0 = y1;
        T c1 = 0.5f * (y2 - y0);
        T c3 = 1.5f * (y1 - y2) + 0.5f * (y3 - y0);
        T c2 = y0 - y1 + c1 - c3;
        return ((c3 * x + c2) * x + c1) * x + c0;
    }

    // james mccartney
    template<typename T>
    // x = frac
    // y0 = [i-1]
    // y1 = [i]
    // y2 = [i+1]
    // y3 = [i+2]
    T hermite3(T x, T y0, T y1, T y2, T y3)
    {
            // 4-point, 3rd-order Hermite (x-form)
            T c0 = y1;
            T c1 = 0.5f * (y2 - y0);
            T y0my1 = y0 - y1;
            T c3 = (y1 - y2) + 0.5f * (y3 - y0my1 - y2);
            T c2 = y0my1 + c1 - c3;

            return ((c3 * x + c2) * x + c1) * x + c0;
    }

    // laurent de soras
    template<typename T>
    // x[i-1]
    // x[i]
    // x[i+1]
    // x[i+2]    
    inline T hermite4(T frac_pos, T xm1, T x0, T x1, T x2)
    {
        const T    c     = (x1 - xm1) * 0.5f;
        const T    v     = x0 - x1;
        const T    w     = c + v;
        const T    a     = w + v + (x2 - x0) * 0.5f;
        const T    b_neg = w + a;

        return ((((a * frac_pos) - b_neg) * frac_pos + c) * frac_pos + x0);
    }

    template<typename T>
    sample_vector<T> mix(const sample_vector<T> & a, const sample_vector<T> & b, T f)
    {
        assert(a.size() == b.size());
        sample_vector<T> r(a.size());
        #pragma omp simd
        for(size_t i = 0; i < r.size(); i++) r[i] = a[i] + f*(b[i]-a[i]);
        return r;
    }

    template<typename T>
    sample_vector<T> interp2x(const sample_vector<T> & a)
    {    
        sample_vector<T> r(a.size()*2);
        size_t n=0;
        #pragma omp simd
        for(size_t i = 0; i < a.size(); i++) 
        {
            r[n++] = a[i];
            r[n++] = cubic_interpolate(T(0.5),a[i==0? 0:i-1],a[i],a[(i+1) % a.size()],a[(i+2) % a.size()]);
        }
        return r;
    }

    template<typename T>
    sample_vector<T> interp4x(const sample_vector<T> & a)
    {    
        sample_vector<T> r(a.size()*4);
        size_t n=0;
        #pragma omp simd
        for(size_t i = 0; i < a.size(); i++) 
        {
            r[n++] = a[i];
            r[n++] = cubic_interpolate(T(0.25),a[i==0? 0:i-1],a[i],a[(i+1) % a.size()],a[(i+2) % a.size()]);
            r[n++] = cubic_interpolate(T(0.5),a[i==0? 0:i-1],a[i],a[(i+1) % a.size()],a[(i+2) % a.size()]);
            r[n++] = cubic_interpolate(T(0.75),a[i==0? 0:i-1],a[i],a[(i+1) % a.size()],a[(i+2) % a.size()]);
        }
        return r;
    }

    template<typename T>
    struct RingBuffer : public sample_vector<T>
    {
        size_t r=0;
        size_t w=0;

        RingBuffer(size_t n) {
            sample_vector<T>::resize(n);
        }

        void set_write_position(size_t n) {
            w = (n % sample_vector<T>::size());
        }  
        T    get() {
            return (*this)[r++];
        }
        void push(T x) {
            (*this)[w++] = x;
            w = (w % sample_vector<T>::size());
        }
        T linear() {
            T x = (*this)[r];
            T x1= (*this)[r++];
            T f = x - floor(x);
            r = r % sample_vector<T>::size();
            return linear_interpolate(x,x1,f);        
        }
        T cubic() {
            T xm1= (*this)[(r-1) % sample_vector<T>::size()];
            T x = (*this)[r];
            T x1= (*this)[(r+1) % sample_vector<T>::size()];
            T x2= (*this)[(r+2) % sample_vector<T>::size()];
            T f = x - floor(x);
            r++;
            r = r % sample_vector<T>::size();
            return cubic_interpolate(f,xm1,x,x1,x2);        
        }
    };
}