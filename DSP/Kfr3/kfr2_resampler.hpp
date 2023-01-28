#pragma once

#include <cmath>

namespace DSP1
{

    ///////////////////////////////////////////////////////////////
    // Interpolation
    ///////////////////////////////////////////////////////////////

    template<typename T>
    // r = frac
    // x = [i]
    // y = [i+1]
    T linear_interpolate(T x, T y, T r)
    {        
        return r +*x (1.0-r)*y;
        
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
    class Decimator5
    {
    private:
    float R1,R2,R3,R4,R5;
    const float h0;
    const float h1;
    const float h3;
    const float h5;
    public:
    
    Decimator5():h0(346/692.0f),h1(208/692.0f),h3(-44/692.0f),h5(9/692.0f)
    {
        R1=R2=R3=R4=R5=0.0f;
    }
    float Calc(const float x0,const float x1)
    {
        float h5x0=h5*x0;
        float h3x0=h3*x0;
        float h1x0=h1*x0;
        float R6=R5+h5x0;
        R5=R4+h3x0;
        R4=R3+h1x0;
        R3=R2+h1x0+h0*x1;
        R2=R1+h3x0;
        R1=h5x0;
        return R6;
    }
    };


    template<typename T>
    class Decimator7
    {
    private:
    float R1,R2,R3,R4,R5,R6,R7;
    const float h0,h1,h3,h5,h7;
    public:
    Decimator7():h0(802/1604.0f),h1(490/1604.0f),h3(-116/1604.0f),h5(33/1604.0f),h7(-6/1604.0f)
    {
        R1=R2=R3=R4=R5=R6=R7=0.0f;
    }
    float Calc(const float x0,const float x1)
    {
        float h7x0=h7*x0;
        float h5x0=h5*x0;
        float h3x0=h3*x0;
        float h1x0=h1*x0;
        float R8=R7+h7x0;
        R7=R6+h5x0;
        R6=R5+h3x0;
        R5=R4+h1x0;
        R4=R3+h1x0+h0*x1;
        R3=R2+h3x0;
        R2=R1+h5x0;
        R1=h7x0;
        return R8;
    }
    };

    template<typename T>
    class Decimator9
    {
    private:
    float R1,R2,R3,R4,R5,R6,R7,R8,R9;
    const float h0,h1,h3,h5,h7,h9;
    public:
    Decimator9():h0(8192/16384.0f),h1(5042/16384.0f),h3(-1277/16384.0f),h5(429/16384.0f),h7(-116/16384.0f),h9(18/16384.0f)
    {
        R1=R2=R3=R4=R5=R6=R7=R8=R9=0.0f;
    }
    float Calc(const float x0,const float x1)
    {
        float h9x0=h9*x0;
        float h7x0=h7*x0;
        float h5x0=h5*x0;
        float h3x0=h3*x0;
        float h1x0=h1*x0;
        float R10=R9+h9x0;
        R9=R8+h7x0;
        R8=R7+h5x0;
        R7=R6+h3x0;
        R6=R5+h1x0;
        R5=R4+h1x0+h0*x1;
        R4=R3+h3x0;
        R3=R2+h5x0;
        R2=R1+h7x0;
        R1=h9x0;
        return R10;
    }
    };


    ///////////////////////////////////////////////////////////////
    // Resampler/Upsample/Downsample
    ///////////////////////////////////////////////////////////////

    // this is interpolator/decimator
    template<typename T>
    struct Resampler
    {
        kfr::samplerate_converter<T> *resampler;

        Resampler(double insr, double outsr)
        {
            resampler = new kfr::samplerate_converter<T>(kfr::resample_quality::normal,outsr,insr);
        }
        ~Resampler() {
            if(resampler) delete resampler;
        }
        void Process(kfr::univector<T> & out, kfr::univector<T> & in) {
            resampler->process(out,in);        
        }
    };

    template<typename T>
    kfr::univector<T> upsample2x(kfr::univector<T> in)
    {
        kfr::univector<T> out(in.size()*2);
        zeros(out);
        for(size_t i = 0; i < in.size(); i++)
            out[i*2] = in[i];
        return out;
    }
    template<typename T>
    kfr::univector<T> upsample4x(kfr::univector<T> in)
    {
        kfr::univector<T> out(in.size()*4);
        zeros(out);
        for(size_t i = 0; i < in.size(); i++)
            out[i*4] = in[i];
        return out;
    }
    template<typename T>
    kfr::univector<T> upsample2N(size_t n, kfr::univector<T> in)
    {
        kfr::univector<T> out(in.size()*2*n);
        zeros(out);
        for(size_t i = 0; i < in.size(); i++)
            out[i*2*n] = in[i];
        return out;
    }
    template<typename T>
    kfr::univector<T> downsample2x(kfr::univector<T> in)
    {
        kfr::univector<T> out(in.size()/2);
        zeros(out);
        for(size_t i = 0; i < in.size()/2; i++)
            out[i] = in[i*2];
        return out;
    }
    template<typename T>
    kfr::univector<T> downsample4x(kfr::univector<T> in)
    {
        kfr::univector<T> out(in.size()/4);
        zeros(out);
        for(size_t i = 0; i < in.size()/4; i++)
            out[i] = in[i*4];
        return out;
    }
    template<typename T>
    kfr::univector<T> downsample2N(size_t n, kfr::univector<T> in)
    {
        kfr::univector<T> out(in.size()/(2*n));
        zeros(out);
        for(size_t i = 0; i < in.size()/(2*n); i++)
            out[i] = in[i*2*n];
        return out;
    }

}