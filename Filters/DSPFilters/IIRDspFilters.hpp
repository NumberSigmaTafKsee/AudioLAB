#pragma once
#include <iostream>
#include <vector>
#include <Undenormal.hpp>

#define MAX_ORDER 64

namespace Filters
{
    typedef float DspFloatType;

    struct FilterCoefficients
    {
        DspFloatType a[2];
        DspFloatType b[3];
    };

    struct BiquadSection
    {
        DspFloatType z[3];
        DspFloatType p[3];

        BiquadSection()
        {
            memset(z, 0, sizeof(z));
            memset(p, 0, sizeof(p));
        }
        BiquadSection(const FilterCoefficients &c)
        {
            z[0] = c.b[0];
            z[1] = c.b[1];
            z[2] = c.b[2];
            p[0] = c.a[0];
            p[1] = c.a[1];
        }
        BiquadSection(DspFloatType z1, DspFloatType z2, DspFloatType z3, DspFloatType p1, DspFloatType p2)
        {
            z[0] = z1;
            z[1] = z2;
            z[2] = z3;
            p[0] = p1;
            p[1] = p2;
        }
        BiquadSection(const BiquadSection &b)
        {
            memcpy(z, b.z, sizeof(z));
            memcpy(p, b.p, sizeof(p));
        }
        void setCoefficients(DspFloatType z1, DspFloatType z2, DspFloatType z3, DspFloatType p1, DspFloatType p2)
        {
            z[0] = z1;
            z[1] = z2;
            z[2] = z3;
            p[0] = p1;
            p[1] = p2;
        }
        void setCoefficients(DspFloatType n[3], DspFloatType d[2])
        {
            memcpy(z, n, sizeof(z));
            memcpy(p, d, sizeof(p));
        }
        void setCoefficients(const FilterCoefficients &c)
        {
            z[0] = c.b[0];
            z[1] = c.b[1];
            z[2] = c.b[2];
            p[0] = c.a[0];
            p[1] = c.a[1];
        }
        BiquadSection &operator=(const BiquadSection &b)
        {
            memcpy(z, b.z, sizeof(z));
            memcpy(p, b.p, sizeof(p));
            return *this;
        }

        void print()
        {
            std::cout << z[0] << " + " << z[1] << " z^-1 + " << z[2] << " z^-1\n";
            std::cout << "-------------------------------------------------------------\n";
            std::cout << " 1 + " << p[0] << +" z^-1 + " << p[1] << " z^-2\n";
        }
    };


    using BiquadSOS = std::vector<BiquadSection>;
    
    struct BiquadTransposedTypeII 
    {
        BiquadSection biquad;
        DspFloatType x, y, d1, d2;

        BiquadTransposedTypeII() 
        {
            x = y = 0;
            d1 = d2 = 0;
        }
        BiquadTransposedTypeII(const BiquadSection &b) : biquad(b)
        {
            x = y = 0;
            d1 = d2 = 0;
        }
        BiquadTransposedTypeII &operator=(const BiquadTransposedTypeII &b)
        {
            biquad = b.biquad;
            x = b.x;
            y = b.y;
            d1 = b.d1;
            d2 = b.d2;
            return *this;
        }
        void setCoefficients(const BiquadSection &b)
        {
            biquad = b;
        }        
        void setBiquad(const BiquadSection &b)
        {
            biquad = b;
        }

        // transposed is just flip - to +
        DspFloatType Tick(DspFloatType I, DspFloatType A = 1, DspFloatType X = 0, DspFloatType Y = 0)
        {
            Undenormal denormal;
            x = I;
            y = biquad.z[0] * x + d1;
            d1 = biquad.z[1] * x - biquad.p[0] * y + d2;
            d2 = biquad.z[2] * x - biquad.p[1] * y;
            return A * y;
        }
    };

    struct FilterBase
    {
        
        std::vector<BiquadTransposedTypeII> biquads;
        size_t order;
        DspFloatType fc,sr,R,q,bw,g,ripple,rolloff,stop,pass;
        bool init = false;

        FilterBase() 
        {
        }
        FilterBase(const BiquadSOS &s)
        {
            setCoefficients(s);
        }
        void setCoefficients(const BiquadSOS &s)
        {        
            if(s.size() == 0) {
                init = false;
                return;
            }
            biquads.resize(s.size());
            for (size_t i = 0; i < s.size(); i++)
            {
                biquads[i].setCoefficients(s[i]);
            }
            init = true;
        }
        DspFloatType Tick(DspFloatType I, DspFloatType A = 1, DspFloatType X = 0, DspFloatType Y = 0)
        {
            if(!init) return 0;
            DspFloatType o = biquads[0].Tick(I, A, X, Y);
            for (size_t i = 1; i < biquads.size(); i++)
                o = biquads[i].Tick(o, A, X, Y);
            return A * o;
        }

        void ProcessBlock(size_t n, DspFloatType * in, DspFloatType * out) {
            for(size_t i = 0; i < n; i++) out[i] = Tick(in[i]);
        }
        void ProcessBlock(size_t n, DspFloatType * in, DspFloatType * out, DspFloatType * A) {
            for(size_t i = 0; i < n; i++) out[i] = Tick(in[i],A[i]);
        }
        void ProcessBlock(size_t n, DspFloatType * in, DspFloatType * out, DspFloatType * A, DspFloatType * X) {
            for(size_t i = 0; i < n; i++) out[i] = Tick(in[i],A[i],X[i]);
        }
        void ProcessBlock(size_t n, DspFloatType * in, DspFloatType * out, DspFloatType * A, DspFloatType * X, DspFloatType * Y) {
            for(size_t i = 0; i < n; i++) out[i] = Tick(in[i],A[i],X[i],Y[i]);
        }
    };
}