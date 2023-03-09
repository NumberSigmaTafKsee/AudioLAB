#pragma once

#include "Core/core_cauerfilter.hpp"

namespace Filters
{
    ////////////////////////////////////////////////////////////////////////
    // IIR Filters Pole calculations
    ////////////////////////////////////////////////////////////////////////
    // -(p1 + p2)s
    DspFloatType complex_pole1(std::complex<DspFloatType> p1, std::complex<DspFloatType> p2)
    {
        return abs(p1 + p2);
    }

    // (p1*p2)
    DspFloatType complex_pole2(std::complex<DspFloatType> p1, std::complex<DspFloatType> p2)
    {
        return abs(p1 * p2);
    }

    std::complex<DspFloatType> butterworthpole(DspFloatType k, DspFloatType n)
    {
        DspFloatType p = M_PI * ((2 * k + n - 1) / (2 * n));
        return std::complex<DspFloatType>(-std::cos(p), -std::sin(p));
    }
    std::complex<DspFloatType> ButterworthPoles(DspFloatType K, DspFloatType N)
    {
        DspFloatType theta = ((2*K+N-1)/(2*N))*M_PI;
        DspFloatType sk = cos(theta);
        DspFloatType omegak = -sin(theta);
        std::complex<DspFloatType> p(sk,omegak);
        return p;
    }
    std::complex<DspFloatType> ChebyshevH0(DspFloatType N, DspFloatType r=1.0)
    {        
        DspFloatType e     = sqrt(pow(10.0,r/10.0)-1.0);
        DspFloatType theta = (M_PI/2.0)+((2*1-1.0)/(2.0*N))*M_PI;
        DspFloatType phi   = (1.0/N)*asinh(1.0/e);
        std::complex<DspFloatType> P(sinh(phi)*cos(theta),-cosh(phi)*sin(theta));        
        for(size_t K=2; K <= N; K++)
        {
            e     = sqrt(pow(10.0,r/10.0)-1.0);
            theta = (M_PI/2.0)+((2*K-1.0)/(2.0*N))*M_PI;
            phi   = (1.0/N)*asinh(1.0/e);
            std::complex<DspFloatType> p(sinh(phi)*cos(theta),-cosh(phi)*sin(theta));        
            P *= -p;        
        }
        if(fmod(N,2) == 0) return P/std::sqrt(1 + e*e);
        return P;
    }

    // Chebyshev2 = 1/pole
    std::complex<DspFloatType> ChebyshevPole(DspFloatType K, DspFloatType N, DspFloatType r=1.0)
    {      
        DspFloatType e     = sqrt(pow(10.0,r/10.0)-1.0);
        DspFloatType theta = (M_PI/2.0)+((2*K-1.0)/(2.0*N))*M_PI;
        DspFloatType phi   = (1.0/N)*asinh(1.0/e);
        std::complex<DspFloatType> p(sinh(phi)*cos(theta),-cosh(phi)*sin(theta));
        return p;
    }


    // Chebyshev2 = 1/pole
    std::complex<DspFloatType> Chebyshev2Zeros(DspFloatType K, DspFloatType N, DspFloatType r=1.0)
    {      
        DspFloatType uk    = ((2*K-1)/N)*(M_PI/2.0);
        std::complex<DspFloatType> p(0,cos(uk));
        return (DspFloatType)1.0/p;
    }

    // Chebyshev2 = 1/pole
    std::complex<DspFloatType> Chebyshev2Pole(DspFloatType K, DspFloatType N, DspFloatType r=1.0)
    {      
        DspFloatType e     = 1.0/sqrt(pow(10.0,r/10.0)-1.0);
        DspFloatType theta = (M_PI/2.0)+((2*K-1.0)/(2.0*N))*M_PI;
        DspFloatType phi   = (1.0/N)*asinh(1.0/e);
        std::complex<DspFloatType> p(sinh(phi)*cos(theta),-cosh(phi)*sin(theta));
        return (DspFloatType)1.0/p;
    }

    std::vector<std::complex<DspFloatType>> ChebyshevPoles(DspFloatType N, DspFloatType r)
    {
        std::vector<std::complex<DspFloatType>> out(N);
        for(size_t K = 1; K <= N; K++)
        {
            out.push_back(ChebyshevPole(K,N,r));
        }
        return out;
    }

    // Hn(s) = H0/d * PROD( (s^2 + ai) / (s^2 + bi*s + ci))
    // d = s+p0 n=odd
    // d = 1    n=even


    DspFloatType factorial(DspFloatType x)
    {
        if (x == 0)
            return 1;
        return x * factorial(x - 1);
    }
    DspFloatType binomial(DspFloatType n, DspFloatType k)
    {
        if(n == 0 || k == 0) return 1;
        return factorial(n) / (factorial(k) * factorial(fabs(n - k)));
    }

    DspFloatType bk(DspFloatType n, DspFloatType k) {
        DspFloatType num = factorial(2*n-k);
        DspFloatType den = pow(2.0,n-k)*factorial(k)*factorial(n-k);
        return num/den;
    }

    // calculate bessel filter coefficients qn(s) order n 
    // it will return the polynomial coefficients starting at b0
    // b0,b1*s,b2*s^2 + ... bn*s^n
    // these use the polynomial root solver to find the poles (complex)
    std::vector<DspFloatType> qn(size_t n) 
    {
        std::vector<DspFloatType> r(n+1);
        r[0] = bk(n,0);
        for(size_t k=1; k <= n; k++) {
            DspFloatType b = bk(n,k);
            r[k] = b;
        }
        return r;
    }



    // will calculate to the 8th order
    // I do not know how to calculate renorm[x]
    #define MAXORDER 9
    void bessel_coefficients(int order, std::vector<DspFloatType> & coeff, bool D=true)
    {
        int i,N,index,indexM1,indexM2;
        DspFloatType B[3][MAXORDER];
        DspFloatType A,renorm[MAXORDER];
        renorm[2] = 0.72675;
        renorm[3] = 0.57145;
        renorm[4] = 0.46946;
        renorm[5] = 0.41322;
        renorm[6] = 0.37038;
        renorm[7] = 0.33898;
        renorm[8] = 0.31546;
        index = 1;
        indexM1 = 0;
        indexM2 = 2;

        memset(B,0,3*MAXORDER*sizeof(DspFloatType));
        B[0][0] = 1.0;
        B[1][0] = 1.0;
        B[1][1] = 1.0;
        coeff.resize(order+1);
        for(N=2; N <= order; N++)
        {
            index = (index+1)%3;
            indexM1= (indexM1 + 1)%3;
            indexM2= (indexM2 + 1)%3;
            for(i=0; i < N; i++) B[index][i] = (2*N-1) * B[indexM1][i];
            for(i=2; i <= N; i++) B[index][i] = B[index][i] + B[indexM2][i-2];
            if(D) 
                for(i=0; i <= order; i++) coeff[i] = B[index][i];
            else 
                for(i=0; i <= order; i++) coeff[i] = B[index][i] * pow(A,order-i);
        }
    }

///////////////////////////////////
// Bessel filter pole/zero tables
///////////////////////////////////


    // i will make the polynormial root solver version soon anyways

    // if you want more use octave
    // [Z,P,K] = bessel(order,1,'s')
    // you only want P (poles) as bessel has no finite zeros
    // make sure it is 1 for normalized
    // make sure it is 's' for analog
    // only use lowpass the other filters can be frequency transforms from it
    // you can do this for any filters really but I alreayd had the others
    // it is much easier to make a script to generate them already


    std::complex<DspFloatType> bessel_poles_2[] = { std::complex<DspFloatType>(-0.5001,0.8660), 
                                            std::complex<DspFloatType>(-0.5001,-0.8660)};
    std::complex<DspFloatType> bessel_poles_3[] = { std::complex<DspFloatType>(-0.9416,0),
                                            std::complex<DspFloatType>(-0.7456,0.7114),
                                            std::complex<DspFloatType>(-0.7456,-0.7114)};                                          
    std::complex<DspFloatType> bessel_poles_4[] = { std::complex<DspFloatType>(-0.6572,0.0302), 
                                            std::complex<DspFloatType>(-0.6572,0.0302), 
                                            std::complex<DspFloatType>(-0.9048,0.2709),                                          
                                            std::complex<DspFloatType>(-0.9048,-0.2709)};
    std::complex<DspFloatType> bessel_poles_5[] = { std::complex<DspFloatType>(-0.9264,0),
                                            std::complex<DspFloatType>(-0.5906,0.9072),
                                            std::complex<DspFloatType>(-0.5906,-0.9072),
                                            std::complex<DspFloatType>(-0.8516,0.4427),
                                            std::complex<DspFloatType>(-0.8515,-0.4427)};
    std::complex<DspFloatType> bessel_poles_6[] = { 
                                            std::complex<DspFloatType>(-0.5386,0.9617),
                                            std::complex<DspFloatType>(-0.5386,-0.9617),
                                            std::complex<DspFloatType>(-0.7997,0.5622),
                                            std::complex<DspFloatType>(-0.7997,-0.5622),
                                            std::complex<DspFloatType>(-0.9094,0.1857),
                                            std::complex<DspFloatType>(-0.9094,-0.1857)};
    std::complex<DspFloatType> bessel_poles_7[] = { std::complex<DspFloatType>(-0.9195,0),
                                            std::complex<DspFloatType>(-0.4967,1.0025), 
                                            std::complex<DspFloatType>(-0.4967,-1.0025),
                                            std::complex<DspFloatType>(-0.7527,0.6505),
                                            std::complex<DspFloatType>(-0.7527,-0.6505),
                                            std::complex<DspFloatType>(-0.8800,0.3217),
                                            std::complex<DspFloatType>(-0.8800,-0.3217)};
    std::complex<DspFloatType> bessel_poles_8[] = { 
                                            std::complex<DspFloatType>(-0.4622,1.0344),
                                            std::complex<DspFloatType>(-0.4622,-1.0344),
                                            std::complex<DspFloatType>(-0.7111,0.7187),
                                            std::complex<DspFloatType>(-0.7111,-0.7187),
                                            std::complex<DspFloatType>(-0.8473,0.4259),
                                            std::complex<DspFloatType>(-0.8473,-0.4259),
                                            std::complex<DspFloatType>(-0.9097,0.1412),
                                            std::complex<DspFloatType>(-0.9097,-0.1412) };
    


    std::complex<DspFloatType> *bessel_poles[] = {
        bessel_poles_2,
        bessel_poles_3,
        bessel_poles_4,
        bessel_poles_5,
        bessel_poles_6,
        bessel_poles_7,
        bessel_poles_8,
    };


//////////////////////////
// Filter Coefficients
//////////////////////////

    struct FilterCoefficients
    {
        DspFloatType a[2];
        DspFloatType b[3];
    };

//////////////////////////
// BiquadSection
//////////////////////////

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
            //memcpy(z, b.z, sizeof(z));
            //memcpy(p, b.p, sizeof(p));
            z[0] = b.z[0];
            z[1] = b.z[1];
            z[2] = b.z[2];
            p[0] = b.p[0];
            p[1] = b.p[1];
            p[2] = b.p[2];
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
            //memcpy(z, n, sizeof(z));
            //memcpy(p, d, sizeof(p));
            z[0] = n[0];
            z[1] = n[1];
            z[2] = n[2];
            p[0] = d[0];
            p[1] = d[1];
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
            //memcpy(z, b.z, sizeof(z));
            //memcpy(p, b.p, sizeof(p));
            z[0] = b.z[0];
            z[1] = b.z[1];
            z[2] = b.z[2];
            p[0] = b.p[0];
            p[1] = b.p[1];
            p[2] = b.p[2];
            return *this;
        }

        void print()
        {
            std::cout << z[0] << "," << z[1] << "," << z[2] << std::endl;
            std::cout << p[0] << "," << p[1] << "," << p[2] << std::endl;
            
        }
    };


//////////////////////////
// BiquadSOS
//////////////////////////
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

        DspFloatType Tick(DspFloatType I, DspFloatType A = 1, DspFloatType X = 0, DspFloatType Y = 0)
        {
            Undenormal denormal;
            x = I;
            y = biquad.z[0] * x + d1;
            d1 = biquad.z[1] * x - biquad.p[0] * y + d2;
            d2 = biquad.z[2] * x - biquad.p[1] * y;
            return A * y;
        }
        void ProcessBlock(size_t n, DspFloatType * in, DspFloatType * out) {
            for(size_t i = 0; i < n; i++) out[i] = Tick(in[i]);
        }
    };


/////////////////////////////////
// BiquadFilterCascade
/////////////////////////////////
    struct BiquadFilterCascade
    {
        
        std::vector<BiquadTransposedTypeII> biquads;
        size_t order;
        DspFloatType fc,sr,R,q,bw,g,ripple,rolloff,stop,pass;
        bool init = false;

        BiquadFilterCascade() 
        {
        }
        BiquadFilterCascade(const BiquadSOS &s)
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

     
/////////////////////////////////
// prewarp
/////////////////////////////////
    void prewarp(DspFloatType *a0, DspFloatType *a1, DspFloatType *a2, DspFloatType fc, DspFloatType fs)
    {
        DspFloatType wp, pi;

        pi = 4.0 * std::atan(1.0);
        wp = 2.0 * fs * std::tan(pi * fc / fs);

        *a2 = (*a2) / (wp * wp);
        *a1 = (*a1) / wp;
    }
    void prewarpR(DspFloatType *a0, DspFloatType *a1, DspFloatType *a2, DspFloatType fc, DspFloatType fs, DspFloatType R)
    {
        DspFloatType wp, pi;

        pi = 4.0 * std::atan(1.0);
        wp = 2.0 * fs * std::tan(pi * fc / fs);

        *a2 = R * R * (*a2) / (wp * wp);
        *a1 = R * (*a1) / wp;
    }
    void prewarpQ(DspFloatType *a0, DspFloatType *a1, DspFloatType *a2, DspFloatType fc, DspFloatType fs, DspFloatType Q)
    {
        DspFloatType wp, pi;

        pi = 4.0 * std::atan(1.0);
        wp = 2.0 * fs * std::tan(pi * fc / fs);

        *a2 = (*a2) / (Q * Q * wp * wp);
        *a1 = (*a1) / (Q * wp);
    }
    void prewarpRQ(DspFloatType *a0, DspFloatType *a1, DspFloatType *a2, DspFloatType fc, DspFloatType fs, DspFloatType R, DspFloatType Q)
    {
        DspFloatType wp, pi;

        pi = 4.0 * std::atan(1.0);
        wp = 2.0 * fs * std::tan(pi * fc / fs);

        *a2 = R * R * (*a2) / (Q * Q * wp * wp);
        *a1 = R * (*a1) / (Q * wp);
    }

/////////////////////////////////
// Bilinear/Z
/////////////////////////////////
    void inversebilinear(
        DspFloatType z[3], DspFloatType p[3],
        DspFloatType k,    /* overall gain factor */
        DspFloatType fs,   /* sampling rate */
        DspFloatType *coef /* pointer to 4 iir coefficients */
    )
    {
        DspFloatType ad, bd;
        DspFloatType b0 = k;
        DspFloatType b1 = coef[0];
        DspFloatType b2 = coef[1];
        DspFloatType a0 = 1;
        DspFloatType a1 = coef[2];
        DspFloatType a2 = coef[3];

        ad = 1 / 4. * a2 / (fs * fs) + a1 / (2 * fs) + a0;
        bd = 1 / 4. * b2 / (fs * fs) + b1 / (2 * fs) + b0;

        z[0] = k * bd / ad;
        z[1] = b1 * bd;
        z[2] = b2 * bd;
        p[0] = 1;
        p[1] = a1 * ad;
        p[2] = a2 * ad;
    }

    void bilinear(
        DspFloatType a0, DspFloatType a1, DspFloatType a2, /* numerator coefficients */
        DspFloatType b0, DspFloatType b1, DspFloatType b2, /* denominator coefficients */
        DspFloatType *k,                       /* overall gain factor */
        DspFloatType fs,                       /* sampling rate */
        DspFloatType *coef                     /* pointer to 4 iir coefficients */
    )
    {
        DspFloatType ad, bd;

        
        /* alpha (Numerator in s-domain) */
        ad = 4. * a2 * fs * fs + 2. * a1 * fs + a0;
        /* beta (Denominator in s-domain) */
        bd = 4. * b2 * fs * fs + 2. * b1 * fs + b0;

        /* update gain constant for this section */
        *k *= ad / bd;
    
        
        /* Nominator */
        *coef++ = (2. * a0 - 8. * a2 * fs * fs) / ad;         /* alpha1 */
        *coef++ = (4. * a2 * fs * fs - 2. * a1 * fs + a0) / ad; /* alpha2 */

        *coef++ = (2. * b0 - 8. * b2 * fs * fs) / bd;           /* beta1 */
        *coef++ = (4. * b2 * fs * fs - 2. * b1 * fs + b0) / bd; /* beta2 */        
    }


/////////////////////////////////////////////////////////////////////////////////////////////
// frequency transform
/////////////////////////////////////////////////////////////////////////////////////////////
// lp(s/wc) => ( (s/wc) - p1) * (s/wc) - p2) => (s/wc)^2 - (s/wc)*p1 - (s/wc)*p2 + p1*p2
// hp(wc/s) => ( (wc/s) - p1) * (wc/s) - p2) => (wc/s)^2 - (wc/s)*p1 - (wc/s) *p2 + p1*p2
//          => (wc^2/s^2 - (wc/s)*p1 - (wc/s)*p2 +p1p2)
//          => s^2 / (wc^2 - wc*s*p1 - wc*s*p2 + p1p2*s^2)
/////////////////////////////////////////////////////////////////////////////////////////////
    void lp2lp(std::vector<std::complex<DspFloatType>> & poles,
                DspFloatType wc, DspFloatType & gain) 
    {
            
        gain *= pow(wc,poles.size());
        for(size_t i = 0; i < poles.size(); i++)
            poles[i] = wc * poles[i];
            
    }
    void lp2hp(std::vector<std::complex<DspFloatType>> & zeros,
                std::vector<std::complex<DspFloatType>> & poles,
                DspFloatType wc, DspFloatType & gain)
    {
        std::complex<DspFloatType> prodz(1.0,0.0),prodp(1.0,0.0);
        for(size_t i = 0; i < zeros.size(); i++) prodz *= -zeros[i];
        for(size_t i = 0; i < poles.size(); i++) prodp *= -poles[i];
        gain *= prodz.real() / prodp.real();
        for(size_t i = 0; i < poles.size(); i++)
            if(abs(poles[i])) poles[i] = std::complex<DspFloatType>(wc) / poles[i];
        zeros.resize(poles.size());
        for(size_t i = 0; i < zeros.size(); i++)
            zeros[i] = std::complex<DspFloatType>(0.0);
    }               
    void lp2bp(std::vector<std::complex<DspFloatType>> & zeros,
                std::vector<std::complex<DspFloatType>> & poles,
                DspFloatType wu, DspFloatType wl, DspFloatType & gain)
    {
        DspFloatType wc = sqrt(wu*wl);
        DspFloatType bw = wu-wl;
        gain      *= pow(bw,poles.size()-zeros.size());
        std::vector<std::complex<DspFloatType>> temp;
        for(size_t i = 0; i < poles.size(); i++) 
        {
            if(abs(poles[i])) {
                std::complex<DspFloatType> first = DspFloatType(0.5) * poles[i] * bw;
                std::complex<DspFloatType> second= DspFloatType(0.5) * std::sqrt(bw*bw) * (poles[i]*poles[i]-DspFloatType(4.0)*wc*wc);
                temp.push_back(first + second);
            }
        }
        for(size_t i = 0; i < poles.size(); i++) {
            if(abs(poles[i])) {
                std::complex<DspFloatType> first = DspFloatType(0.5) * poles[i] * bw;
                std::complex<DspFloatType> second= DspFloatType(0.5) * std::sqrt(bw*bw) * (poles[i]*poles[i]-DspFloatType(4.0)*wc*wc);
                temp.push_back(first - second);
            }
        }
        zeros.resize(poles.size());
        for(size_t i = 0; i < zeros.size(); i++) {
            zeros[i] = std::complex<DspFloatType>(0);
        }
        size_t index = 0;
        poles.resize(temp.size());
        for(auto i = temp.begin(); i != temp.end(); i++) {
            poles[index] = *i;
            index++;
        }        
    }       
    void lp2bs(std::vector<std::complex<DspFloatType>> & zeros,
                std::vector<std::complex<DspFloatType>> & poles,
                DspFloatType wu, DspFloatType wl, DspFloatType & gain)
    { 
        DspFloatType bw = wu-wl;
        DspFloatType Wc = sqrt(wu*wl);
        std::complex<DspFloatType> prodz(1.0,0.0);
        std::complex<DspFloatType> prodp(1.0,0.0);
        for(size_t i = 0; i < zeros.size(); i++)
            prodz *= -zeros[i];
        for(size_t i = 0; i < poles.size(); i++)
            prodp *= -poles[i];
        gain *= prodz.real() / prodp.real();
        std::vector<std::complex<DspFloatType>> ztmp;
        for(size_t i = 0; i < zeros.size(); i++) {
            ztmp.push_back(std::complex<DspFloatType>(0.0,Wc));
            ztmp.push_back(std::complex<DspFloatType>(0.0,-Wc));            
        }
        std::vector<std::complex<DspFloatType>> ptmp;
        for(size_t i = 0; i < poles.size(); i++) {
            if(abs(poles[i])) {
                std::complex<DspFloatType> term1 = DspFloatType(0.5) * bw / poles[i];
                std::complex<DspFloatType> term2 = DspFloatType(0.5) * sqrt((bw*bw) / (poles[i]*poles[i]) - (DspFloatType(4)*Wc*Wc));
                ptmp.push_back(term1+term2);
            }
        }
        size_t index = 0;
        for(auto i = ztmp.begin(); i != ztmp.end(); i++) {
            zeros[index++] = *i;
        }
        index = 0;
        for(auto i = ptmp.begin(); i != ptmp.end(); i++) {
            poles[index++] = *i;
        }
    } 



/////////////////////////////////////////////////////////////////////////////////////////////
// AnalogBiquadSection
/////////////////////////////////////////////////////////////////////////////////////////////

    // convert analog setion to biquad type I
    BiquadSection AnalogBiquadSection(const BiquadSection &section, DspFloatType fc, DspFloatType fs)
    {
        BiquadSection ns = section;
        prewarp(&ns.z[0], &ns.z[1], &ns.z[2], fc, fs);
        prewarp(&ns.p[0], &ns.p[1], &ns.p[2], fc, fs);
        
        //std::vector<DspFloatType> coeffs(4);
        DspFloatType coeffs[4] = {0,0,0,0};
        DspFloatType k =1;
        
        bilinear(ns.z[0], ns.z[1], ns.z[2], ns.p[0], ns.p[1], ns.p[2], &k, fs, coeffs);
        ns.z[0] = k;
        ns.z[1] = k*coeffs[0];
        ns.z[2] = k*coeffs[1];
        ns.p[0] = coeffs[2];
        ns.p[1] = coeffs[3];
        ns.p[2] = 0;    
        return ns;
    }

/////////////////////////////////////////////////////////////////////////////////////////////
// AnalogBiquadCascade
/////////////////////////////////////////////////////////////////////////////////////////////

    // H(s) => Bilinear/Z => H(z)
    // convert analog sos to biquad cascade type I
    BiquadSOS AnalogBiquadCascade(const BiquadSOS &sos, DspFloatType fc, DspFloatType fs)
    {
        BiquadSOS nsos = sos;
        for (size_t i = 0; i < sos.size(); i++)
        {
            BiquadSection b = AnalogBiquadSection(sos[i], fc, fs);
            nsos[i] = b;
        }
        return nsos;
    }


/////////////////////////////////////////////////////////////////////////////////////////////
// Biquad
/////////////////////////////////////////////////////////////////////////////////////////////
    struct Biquad
    {
        DspFloatType z[3];
        DspFloatType p[3];
        DspFloatType x[2];
        DspFloatType y[2];

        Biquad() {
            x[0] = x[1] = 0;
            y[0] = y[1] = 0;
        }
        void setCoefficients(DspFloatType _z[3], DspFloatType _p[3]) {
            memcpy(z,_z,3*sizeof(DspFloatType));
            memcpy(p,_p,3*sizeof(DspFloatType));
        }
        DspFloatType Tick(DspFloatType I) {
            DspFloatType out = z[0]*I + z[1]*x[0] + z[2] * z[1];
            out = out - p[0]*y[0] - p[1]*y[1];
            x[1] = x[0];
            x[0] = I;
            y[1] = y[0];
            y[1] = out;
            return out;
        }
    };


/////////////////////////////////////////////////////////////////////////////////////////////
// Butterworth
/////////////////////////////////////////////////////////////////////////////////////////////
    BiquadSOS butterlp(int order, double Q=1.0)
    {
        BiquadSOS sos;    
        size_t n = 1;
        if(order %2 != 0) {
            BiquadSection c;
            std::complex<DspFloatType> p1  = ButterworthPoles(1,order);        
            DspFloatType x1 = abs(p1);
            DspFloatType x2 = 0;
                    
            // (s-p1)        
            c.z[0] = 1.0/x1;
            c.z[1] = 0.0;
            c.z[2] = 0.0;
            c.p[0] = 1.0;;        
            c.p[1] = 1/x1;
            c.p[2] = 0.0;
            sos.push_back(c);        
            n++;
        }
            
        for(size_t i = n; i < order; i += 2)
        {
            std::complex<DspFloatType> p1  = ButterworthPoles(i,order);
            std::complex<DspFloatType> p2  = ButterworthPoles(i+1,order);
            
            DspFloatType x1 = abs(p1*p2);
            DspFloatType x2 = abs(-p1-p2);
            
            // (s-p1)*(s-p2) = s^2 -sp1 -sp2 +p1p2
            BiquadSection c;
            c.z[0] = 1.0/x1;
            c.z[1] = 0.0;
            c.z[2] = 0.0;
            c.p[0] = 1;        
            c.p[1] = (x2/Q)/x1;
            c.p[2] = 1/x1;    
            sos.push_back(c);
        }
        return sos;
    }
    BiquadSOS butterhp(int order, double Q=1.0)
    {
        BiquadSOS sos;    
        size_t n = 1;
        
        if(order %2 != 0) {
            BiquadSection c;
            std::complex<DspFloatType> p1  = ButterworthPoles(1,order);        
            DspFloatType x1 = abs(p1);
            

            c.z[0] = 0.0;
            c.z[1] = 1.0/x1;
            c.z[2] = 0.0;
            c.p[0] = 1.0;
            c.p[1] = 1/x1;
            c.p[2] = 0.0;

            sos.push_back(c);        
            n++;
        }
            
        for(size_t i = n; i < order; i += 2)
        {
            std::complex<DspFloatType> p1  = ButterworthPoles(i,order);
            std::complex<DspFloatType> p2  = ButterworthPoles(i+1,order);
            
            DspFloatType x1 = abs(p1*p2);
            DspFloatType x2 = abs(-p1-p2);
            
            // (s-p1)*(s-p2) = s^2 -sp1 -sp2 +p1p2
            BiquadSection c;
            c.z[0] = 0.0;
            c.z[1] = 0.0;
            c.z[2] = 1.0/x1;
            c.p[0] = 1;        
            c.p[1] = (1.0/Q)*x2/x1;
            c.p[2] = 1.0/x1;    
            sos.push_back(c);
        }
        
        return sos;
    }

    BiquadSOS butterlp2hp(int order, double Q=1.0)
    {
        BiquadSOS sos;    
        size_t n = 1;
        
        if(order %2 != 0) {
            BiquadSection c;
            std::complex<DspFloatType> p1  = ButterworthPoles(1,order);        
            DspFloatType x1 = abs(p1);
            

            c.z[0] = 0.0;
            c.z[1] = 1.0/x1;
            c.z[2] = 0.0;
            c.p[0] = 1.0;
            c.p[1] = 1/x1;
            c.p[2] = 0.0;

            sos.push_back(c);        
            n++;
        }
            
        for(size_t i = n; i < order; i += 2)
        {
            BiquadSection c;
            std::complex<DspFloatType> p1  = ButterworthPoles(i,order);        
            std::complex<DspFloatType> p2  = ButterworthPoles(1+1,order);        
            
            std::vector<std::complex<DspFloatType>> zeros,poles;
            poles.push_back(p1);
            poles.push_back(p2);
            DspFloatType gain;
            lp2hp(zeros,poles,1.0,gain);
            // (1-z)(1-z) = 1 -z1-z2 +z1z2
            
            DspFloatType x1 = abs(poles[0]*poles[1]);
            c.z[0] = gain*abs(zeros[0]*zeros[1])/x1;
            c.z[1] = gain*abs(-zeros[0]-zeros[1])/x1;
            c.z[2] = 1.0/x1;
            c.p[0] = 1.0;
            c.p[1] = (1.0/Q)*abs(-poles[0]-poles[1])/x1;
            c.p[2] = 1.0/x1;

            sos.push_back(c);
        }
        
        return sos;
    }
    BiquadSOS butterlp2bp(int order, double Q=1.0)
    {
        BiquadSOS sos;    
        size_t n = 1;
                    
        for(size_t i = n; i < order; i += 2)
        {
            BiquadSection c;
            std::complex<DspFloatType> p1  = ButterworthPoles(i,order);        
            std::complex<DspFloatType> p2  = ButterworthPoles(1+1,order);        
            std::vector<std::complex<DspFloatType>> zeros,poles;
            poles.push_back(p1);
            poles.push_back(p2);
            DspFloatType gain;
            lp2bp(zeros,poles,1.0,0.5,gain);
            // (1-z)(1-z) = 1 -z1-z2 +z1z2
            
            DspFloatType x1 = abs(poles[0]*poles[1]);
            c.z[0] = gain*abs(zeros[0]*zeros[1])/x1;
            c.z[1] = gain*abs(-zeros[0]-zeros[1])/x1;
            c.z[2] = 1.0/x1;
            c.p[0] = 1.0;
            c.p[1] = (1.0/Q)*abs(-poles[0]-poles[1])/x1;
            c.p[2] = 1.0/x1;

            sos.push_back(c);
        }
        
        return sos;
    }
    BiquadSOS butterlp2bs(int order, double Q=1.0)
    {
        BiquadSOS sos;    
        size_t n = 1;
                    
        for(size_t i = n; i < order; i += 2)
        {
            BiquadSection c;
            std::complex<DspFloatType> p1  = ButterworthPoles(i,order);        
            std::complex<DspFloatType> p2  = ButterworthPoles(1+1,order);        
            std::vector<std::complex<DspFloatType>> zeros,poles;
            poles.push_back(p1);
            poles.push_back(p2);
            DspFloatType gain;
            lp2bs(zeros,poles,1.0,0.5,gain);
            // (1-z)(1-z) = 1 -z1-z2 +z1z2
            
            DspFloatType x1 = abs(poles[0]*poles[1]);
            c.z[0] = gain*abs(zeros[0]*zeros[1])/x1;
            c.z[1] = gain*abs(-zeros[0]-zeros[1])/x1;
            c.z[2] = 1.0/x1;
            c.p[0] = 1.0;
            c.p[1] = (1.0/Q)*abs(-poles[0]-poles[1])/x1;
            c.p[2] = 1.0/x1;

            sos.push_back(c);
        }
        
        return sos;
    }
    
    BiquadSection butter2lp(double Q=1.0)
    {    
        std::complex<DspFloatType> p1  = ButterworthPoles(1,2);
        std::complex<DspFloatType> p2  = ButterworthPoles(2,2);
    
        DspFloatType x1 = abs(p1*p2);
        DspFloatType x2 = abs(-p1-p2);
        //std::cout << p1 << "," << p2 << "," << x1 << "," << x2 << std::endl;
        // (s-p1)*(s-p2) = s^2 -sp1 -sp2 +p1p2
        BiquadSection c;
        c.z[0] = 1.0/x1;
        c.z[1] = 0.0;
        c.z[2] = 0.0;
        c.p[0] = 1;    
        c.p[1] = (1.0/Q)*x2/x1;
        c.p[2] = 1.0/x1;        

        return c;
    }
    BiquadSection butter2hp(double Q=1.0)
    {    
        std::complex<DspFloatType> p1  = ButterworthPoles(1,2);
        std::complex<DspFloatType> p2  = ButterworthPoles(2,2);
    
        DspFloatType x1 = abs(p1*p2);
        DspFloatType x2 = abs(-p1-p2);
        //std::cout << p1 << "," << p2 << "," << x1 << "," << x2 << std::endl;
        // (s-p1)*(s-p2) = s^2 -sp1 -sp2 +p1p2
        BiquadSection c;
        c.z[0] = 0.0;
        c.z[1] = 0.0;
        c.z[2] = 1.0/x1;
        c.p[0] = 1;    
        c.p[1] = (1.0/Q)*x2/x1;
        c.p[2] = 1.0/x1;        

        return c;
    }
    BiquadSection butterlp2hp2(double Q=1.0)
    {    
        std::complex<DspFloatType> p1  = ButterworthPoles(1,2);
        std::complex<DspFloatType> p2  = ButterworthPoles(2,2);
    
        std::vector<std::complex<DspFloatType>> zeros,poles;
        poles.push_back(p1);
        poles.push_back(p2);
        DspFloatType gain;
        lp2hp(zeros,poles,1.0,gain);
        // (1-z)(1-z) = 1 -z1-z2 +z1z2
        BiquadSection c;        
        DspFloatType x1 = abs(poles[0]*poles[1]);
        c.z[0] = gain*abs(zeros[0]*zeros[1])/x1;
        c.z[1] = gain*abs(-zeros[0]-zeros[1])/x1;
        c.z[2] = 1.0/x1;
        c.p[0] = 1.0;
        c.p[1] = (1.0/Q)*abs(-poles[0]-poles[1])/x1;
        c.p[2] = 1.0/x1;

        return c;
    }
    BiquadSection butterlp2bp2(double Q=1.0)
    {    
        std::complex<DspFloatType> p1  = ButterworthPoles(1,2);
        std::complex<DspFloatType> p2  = ButterworthPoles(2,2);
    
        std::vector<std::complex<DspFloatType>> zeros,poles;
        poles.push_back(p1);
        poles.push_back(p2);
        DspFloatType gain;
        // i dont not really know what this should be normalized 1.0,0 or 1.0,0.5?
        lp2bp(zeros,poles,1.0,0.5,gain);
        // (1-z)(1-z) = 1 -z1-z2 +z1z2
        BiquadSection c;        
        DspFloatType x1 = abs(poles[0]*poles[1]);
        c.z[0] = gain*abs(zeros[0]*zeros[1])/x1;
        c.z[1] = gain*abs(-zeros[0]-zeros[1])/x1;
        c.z[2] = 1.0/x1;
        c.p[0] = 1.0;
        c.p[1] = (1.0/Q)*abs(-poles[0]-poles[1])/x1;
        c.p[2] = 1.0/x1;
        return c;
    }
    BiquadSection butterlp2bs2(double Q=1.0)
    {    
        std::complex<DspFloatType> p1  = ButterworthPoles(1,2);
        std::complex<DspFloatType> p2  = ButterworthPoles(2,2);
    
        std::vector<std::complex<DspFloatType>> zeros,poles;
        poles.push_back(p1);
        poles.push_back(p2);
        DspFloatType gain;
        lp2bs(zeros,poles,1.0,0.5,gain);
        // (1-z)(1-z) = 1 -z1-z2 +z1z2
        BiquadSection c;        
        DspFloatType x1 = abs(poles[0]*poles[1]);
        c.z[0] = gain*abs(zeros[0]*zeros[1])/x1;
        c.z[1] = gain*abs(-zeros[0]-zeros[1])/x1;
        c.z[2] = 1.0/x1;
        c.p[0] = 1.0;
        c.p[1] = (1.0/Q)*abs(-poles[0]-poles[1])/x1;
        c.p[2] = 1.0/x1;
        return c;
    }
    
    
/////////////////////////////////////////////////////////////////////////////////////////////
// Bessel Filter
/////////////////////////////////////////////////////////////////////////////////////////////    
    BiquadSOS bessellp(int order, double Q=1.0)
    {
        
        BiquadSOS sos;    
        if(order <= 0) throw std::runtime_error("Bessel order <= 0");
        if(order > 7) throw std::runtime_error("Bessel order > 7");
        
        size_t n = 1;
        std::complex<DspFloatType> *poles = bessel_poles[order-1];

        if(order %2 != 0) {
            BiquadSection c;
            std::complex<DspFloatType> p1  =poles[0];
            DspFloatType x1 = abs(p1);
            DspFloatType x2 = 0;
                    
            // (s-p1)        
            c.z[0] = 1.0/x1;
            c.z[1] = 0.0;
            c.z[2] = 0.0;
            c.p[0] = 1.0;;        
            c.p[1] = 1/x1;
            c.p[2] = 0.0;
            sos.push_back(c);        
            n++;
        }
            
        for(size_t i = n; i < order; i += 2)
        {
            std::complex<DspFloatType> p1  = poles[i];
            std::complex<DspFloatType> p2  = poles[i+1];
            
            DspFloatType x1 = abs(p1*p2);
            DspFloatType x2 = abs(-p1-p2);
            
            // (s-p1)*(s-p2) = s^2 -sp1 -sp2 +p1p2
            BiquadSection c;
            c.z[0] = 1.0/x1;
            c.z[1] = 0.0;
            c.z[2] = 0.0;
            c.p[0] = 1;        
            c.p[1] = (1.0/Q)*x2/x1;
            c.p[2] = 1/x1;    
            sos.push_back(c);
        }
        return sos;
    }

    // todo polysolver
    BiquadSOS besselhp(int order, double Q=1.0)
    {
        
        BiquadSOS sos;    
        if(order > 8) order = 8;
        size_t n = 1;
        if(order %2 != 0) {
            BiquadSection c;
            std::complex<DspFloatType> p1  = bessel_poles_2[0];
            DspFloatType x1 = abs(p1);
            DspFloatType x2 = 0;
                    
            // (s-p1)        
            c.z[0] = 0.0;
            c.z[1] = 1.0;
            c.z[2] = 0.0;
            c.p[0] = 1.0;;        
            c.p[1] = 1/x1;
            c.p[2] = 0.0;
            sos.push_back(c);        
            n++;
        }
            
        for(size_t i = n; i < order; i += 2)
        {
            std::complex<DspFloatType> p1  = bessel_poles_2[i];
            std::complex<DspFloatType> p2  = bessel_poles_2[i+1];
            
            DspFloatType x1 = abs(p1*p2);
            DspFloatType x2 = abs(-p1-p2);
            
            // (s-p1)*(s-p2) = s^2 -sp1 -sp2 +p1p2
            BiquadSection c;
            
            c.z[0] = 0.0;
            c.z[1] = 0.0;
            c.z[2] = 1.0/x1;
            c.p[0] = 1;        
            c.p[1] = (1.0/Q)*x2/x1;
            c.p[2] = 1/x1;  

            sos.push_back(c);
        }
        return sos;
    }
    
/////////////////////////////////////////////////////////////////////////////////////////////
// Chebyshev 1
/////////////////////////////////////////////////////////////////////////////////////////////
    BiquadSOS cheby1lp(int order, double Q=1.0)
    {
        
        BiquadSOS sos;            
        std::complex<DspFloatType> H0  = ChebyshevH0(order);            
        size_t n = 1;
        if(order %2 != 0) {
            BiquadSection c;
            std::complex<DspFloatType> p1  = ChebyshevPole(1,order);
            DspFloatType x1 = abs(p1);
            DspFloatType x2 = 0;
                    
            // (s-p1)        
            c.z[0] = abs(H0)/x1;
            c.z[1] = 0.0;
            c.z[2] = 0.0;
            c.p[0] = 1.0;;        
            c.p[1] = 1/x1;
            c.p[2] = 0.0;
            sos.push_back(c);        
            n++;
        }
            
        for(size_t i = n; i < order; i += 2)
        {
            std::complex<DspFloatType> p1  = ChebyshevPole(i,order);
            std::complex<DspFloatType> p2  = ChebyshevPole(i+1,order);
            
            DspFloatType x1 = abs(p1*p2);
            DspFloatType x2 = abs(-p1-p2);
            
            // (s-p1)*(s-p2) = s^2 -sp1 -sp2 +p1p2
            BiquadSection c;
            c.z[0] = abs(H0)/x1;
            c.z[1] = 0.0;
            c.z[2] = 0.0;
            c.p[0] = 1;        
            c.p[1] = (1.0/Q)*x2/x1;
            c.p[2] = 1/x1;    
            sos.push_back(c);
        }
        return sos;
    }

    BiquadSOS cheby1hp(int order, double Q=1.0)
    {
        std::complex<DspFloatType> H0  = ChebyshevH0(order);            
        BiquadSOS sos;            
        size_t n = 1;
        if(order %2 != 0) {
            BiquadSection c;
            std::complex<DspFloatType> p1  = ChebyshevPole(1,order);
            DspFloatType x1 = abs(p1);
            DspFloatType x2 = 0;
                    
            // (s-p1)        
            c.z[0] = 0.0;
            c.z[1] = abs(H0)/x1;
            c.z[2] = 0.0;
            c.p[0] = 1.0;        
            c.p[1] = 1/x1;
            c.p[2] = 0.0;
            sos.push_back(c);        
            n++;
        }
            
        for(size_t i = n; i < order; i += 2)
        {
            std::complex<DspFloatType> p1  = ChebyshevPole(i,order);
            std::complex<DspFloatType> p2  = ChebyshevPole(i+1,order);
            
            DspFloatType x1 = abs(p1*p2);
            DspFloatType x2 = abs(-p1-p2);
            
            // (s-p1)*(s-p2) = s^2 -sp1 -sp2 +p1p2
            BiquadSection c;
            c.z[0] = 0.0;
            c.z[1] = 0.0;
            c.z[2] = abs(H0)/x1;
            c.p[0] = 1;        
            c.p[1] = (1.0/Q)*x2/x1;
            c.p[2] = 1/x1;    
            sos.push_back(c);
        }
        return sos;
    }

    

/////////////////////////////////////////////////////////////////////////////////////////////
// Chebyshev 2
/////////////////////////////////////////////////////////////////////////////////////////////
    BiquadSOS cheby2lp(int order, double Q=1.0,double rips=1.0)
    {        
        BiquadSOS sos;            
        size_t n = 1;
        if(order %2 != 0) {
            BiquadSection c;
            std::complex<DspFloatType> H0  = Chebyshev2Zeros(1,order,1.0);
            std::complex<DspFloatType> p1  = Chebyshev2Pole(1,order);

            DspFloatType x1 = abs(p1);
            DspFloatType x2 = abs(H0);
                    
            // (s-p1)        
            c.z[0] = x2/x1;
            c.z[1] = 0.0;
            c.z[2] = 0.0;
            c.p[0] = 1.0;;        
            c.p[1] = 1/x1;
            c.p[2] = 0.0;
            sos.push_back(c);        
            n++;
        }
            
        for(size_t i = n; i < order; i += 2)
        {
            std::complex<DspFloatType> H0  = Chebyshev2Zeros(1,2,1.0);
            std::complex<DspFloatType> H1  = Chebyshev2Zeros(2,2,1.0);
            std::complex<DspFloatType> p1  = Chebyshev2Pole(1,2,1);
            std::complex<DspFloatType> p2  = Chebyshev2Pole(2,2,1);
    
            DspFloatType x1 = abs(p1*p2);
            DspFloatType x2 = abs(-p1-p2);
            DspFloatType z1 = abs(H0*H1);
            DspFloatType z2 = abs(-H0-H1);

            
            // (s-p1)*(s-p2) = s^2 -sp1 -sp2 +p1p2
            BiquadSection c;

            c.z[0] = z1/x1;
            c.z[1] = z2/x1;
            c.z[2] = 0;
            c.p[0] = 1;
            // radius is the same thing but goes from 0..1
            // 0 = most resonant
            // 1 = least resonant
            c.p[1] = (1.0/Q)*x2/x1;
            c.p[2] = 1/x1;    
            
            sos.push_back(c);
        }
        return sos;
    }    
    

/////////////////////////////////////////////////////////////////////////////////////////////
// Elliptic Filter
/////////////////////////////////////////////////////////////////////////////////////////////    
    BiquadSOS ellipticlp( double Ap,  //maximum pass band loss in dB,
                        double As,  //minium stop band loss in dB
                        double wp,  //pass band cutoff frequency
                        double ws   //stop band cut off frequency
                        )
    {
        Filters::CauerFilter filter(Ap,As,wp,ws);
        filter.calculate_coefficients();
        //std::cout << filter.getHo() << "\n";
        //for(size_t i = 0; i < filter.Ai.size(); i++) std::cout << filter.Ai[i] << "\n";
        //for(size_t i = 0; i < filter.Bi.size(); i++) std::cout << filter.Bi[i] << "\n";
        //for(size_t i = 0; i < filter.Ci.size(); i++) std::cout << filter.Ci[i] << "\n";

        double hZero = filter.getHo();
        BiquadSOS sos;
        for(size_t i = 0; i < filter.r; i++)
        {
            BiquadSection c;
            c.z[0] = hZero*filter.Ai[i];
            c.z[1] = 0;
            c.z[2] = hZero;    
            c.p[0] = filter.Ci[i];
            c.p[1] = filter.Bi[i];
            c.p[2] = 1;       
            c.z[0] /= c.p[0];    
            c.z[2] /= c.p[0];    
            c.p[1] /= c.p[0];
            c.p[2] /= c.p[0];
            c.p[0] /= c.p[0];
            sos.push_back(c);
        }
        return sos;
    }       
}
