#pragma once


/*
using namespace Octopus;

IIRFilters::BiquadSOS oct_butterlp(int order, DspFloatType Q=1.0)
{
    OctopusValueList values;
    values.vlist(0) = order;
    values.vlist(1) = 1.0;
    values.vlist(2) = 's';
    values = octave_butter(values,3);    
    // octave tf2sos is too slow
    VectorXcf Z = values.vlist(0).float_complex_row_vector_value();
    VectorXcf P = values.vlist(1).float_complex_row_vector_value();
    float    K = values.vlist(2).float_scalar_value();
    IIRFilters::BiquadSOS sos; 
    size_t n = 0;
    if(order %2 != 0) {
        IIRFilters::BiquadSection c;
        std::complex<DspFloatType> p1  = P(n++);
        DspFloatType x1 = abs(p1);
        DspFloatType x2 = 0;
                
        // (s-p1)        
        c.z[0] = K*1.0/x1;
        c.z[1] = 0.0;
        c.z[2] = 0.0;
        c.p[0] = 1.0;
        c.p[1] = 1.0/x1;
        c.p[2] = 0.0;
        sos.push_back(c);                        
    }         
    for(size_t i = n; i < order; i += 2)
    {
        std::complex<DspFloatType> p1  = P(n++);
        std::complex<DspFloatType> p2  = P(n++);
        
        DspFloatType x1 = abs(p1*p2);
        DspFloatType x2 = abs(-p1-p2);
        
        // (s-p1)*(s-p2) = s^2 -sp1 -sp2 +p1p2
        IIRFilters::BiquadSection c;
        c.z[0] = K*1.0/x1;
        c.z[1] = 0.0;
        c.z[2] = 0.0;
        c.p[0] = 1;        
        c.p[1] = (1.0/Q)*x2/x1;
        c.p[2] = 1/x1;    
        sos.push_back(c);
    }   
    return sos;
}
*/

////////////////////////////////////////////////////////////////////////
// IIR Filters
////////////////////////////////////////////////////////////////////////
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

    // frequency transform
    // lp(s/wc) => ( (s/wc) - p1) * (s/wc) - p2) => (s/wc)^2 - (s/wc)*p1 - (s/wc)*p2 + p1*p2
    // hp(wc/s) => ( (wc/s) - p1) * (wc/s) - p2) => (wc/s)^2 - (wc/s)*p1 - (wc/s) *p2 + p1*p2
    //          => (wc^2/s^2 - (wc/s)*p1 - (wc/s)*p2 +p1p2)
    //          => s^2 / (wc^2 - wc*s*p1 - wc*s*p2 + p1p2*s^2)
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


///////////////////////////////////
// Bessel filter pole/zero tables
///////////////////////////////////
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

    /////////////////////////////////////////////////////////////////////////////////////////////
    // Bessel Filter
    /////////////////////////////////////////////////////////////////////////////////////////////    
    BiquadSOS bessellp(int order, double Q=1.0)
    {
        
        BiquadSOS sos;    
        
        if(order <  2) order = 2;
        if(order >  8) order = 8;
        
        size_t n = 0;
        std::complex<DspFloatType> *poles = bessel_poles[order-2];

        if(order %2 != 0) {
            BiquadSection c;
            std::complex<DspFloatType> p1=poles[n++];
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
        }
            
        for(size_t i = n; i < order; i += 2)
        {
            std::complex<DspFloatType> p1  = poles[n++];
            std::complex<DspFloatType> p2  = poles[n++];
            
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
        if(order < 2) order = 2;
        if(order > 8) order = 8;
        size_t n = 1;
        std::complex<DspFloatType> *poles = bessel_poles[order-2];

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

///////////////////////////////////////
// Butterworth filter pole/zero tables
///////////////////////////////////////
    std::complex<DspFloatType> butter_2_poles[] = {
        std::complex<DspFloatType>(-0.707107,0.707107),
        std::complex<DspFloatType>(-0.707107,-0.707107)
    };
    std::complex<DspFloatType> butter_3_poles[] = {
        std::complex<DspFloatType>(-1.000000,0.000000),
        std::complex<DspFloatType>(-0.500000,0.866025),   
        std::complex<DspFloatType>(-0.500000,-0.866025),
    };
    std::complex<DspFloatType> butter_4_poles[] = {
        std::complex<DspFloatType>(-0.382683,0.923880),
        std::complex<DspFloatType>(-0.923880,0.382683),
        std::complex<DspFloatType>(-0.923880,-0.382683),
        std::complex<DspFloatType>(-0.382683,-0.923880),
    };
    std::complex<DspFloatType> butter_5_poles[] = {
        std::complex<DspFloatType>(-1.000000,0.000000),
        std::complex<DspFloatType>(-0.309017,0.951057),
        std::complex<DspFloatType>(-0.809017,0.587785),
        std::complex<DspFloatType>(-0.809017,-0.587785),
        std::complex<DspFloatType>(-0.309017,-0.951057),
    };
    std::complex<DspFloatType> butter_6_poles[] = {
        std::complex<DspFloatType>(-0.258819,0.965926),
        std::complex<DspFloatType>(-0.707107,0.707107),
        std::complex<DspFloatType>(-0.965926,0.258819),
        std::complex<DspFloatType>(-0.965926,-0.258819),
        std::complex<DspFloatType>(-0.707107,-0.707107),
        std::complex<DspFloatType>(-0.258819,-0.965926),
    };
    std::complex<DspFloatType> butter_7_poles[] = {
        std::complex<DspFloatType>(-1.000000,0.000000),    
        std::complex<DspFloatType>(-0.222521,0.974928),
        std::complex<DspFloatType>(-0.623490,0.781831),
        std::complex<DspFloatType>(-0.900969,0.433884),    
        std::complex<DspFloatType>(-0.900969,-0.433884),
        std::complex<DspFloatType>(-0.623490,-0.781831),
        std::complex<DspFloatType>(-0.222521,-0.974928),
    };
    std::complex<DspFloatType> butter_8_poles[] = {
        std::complex<DspFloatType>(-0.195090,0.980785),
        std::complex<DspFloatType>(-0.555570,0.831470),
        std::complex<DspFloatType>(-0.831470,0.555570),
        std::complex<DspFloatType>(-0.980785,0.195090),
        std::complex<DspFloatType>(-0.980785,-0.195090),
        std::complex<DspFloatType>(-0.831470,-0.555570),
        std::complex<DspFloatType>(-0.555570,-0.831470),
        std::complex<DspFloatType>(-0.195090,-0.980785),
    };
    std::complex<DspFloatType> *butter_poles[] = {
        butter_2_poles,
        butter_3_poles,
        butter_4_poles,
        butter_5_poles,
        butter_6_poles,
        butter_7_poles,
        butter_8_poles,
    };

    BiquadSOS butterlp(int order, double Q=1.0)
    {
        BiquadSOS sos; 
        if(order <= 2) order = 2;
        if(order >  8) order = 8;
        std::complex<DspFloatType> * poles = butter_poles[order-2];
        size_t n = 0;
        if(order %2 != 0) {
            BiquadSection c;
            std::complex<DspFloatType> p1  = poles[n++];
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
        }
            
        for(size_t i = n; i < order; i += 2)
        {
            std::complex<DspFloatType> p1  = poles[n++];
            std::complex<DspFloatType> p2  = poles[n++];
            
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
    BiquadSOS butterhp(int order, double Q=1.0)
    {
        BiquadSOS sos;    
        if(order <= 2) order = 2;
        if(order >  8) order = 8;
        std::complex<DspFloatType> * poles = butter_poles[order-2];
        size_t n = 0;

        if(order %2 != 0) {
            BiquadSection c;
            std::complex<DspFloatType> p1  = poles[n++];
            DspFloatType x1 = abs(p1);
            

            c.z[0] = 0.0;
            c.z[1] = 1.0/x1;
            c.z[2] = 0.0;
            c.p[0] = 1.0;
            c.p[1] = 1/x1;
            c.p[2] = 0.0;

            sos.push_back(c);                        
        }
            
        for(size_t i = n; i < order; i += 2)
        {
            std::complex<DspFloatType> p1  = poles[n++];
            std::complex<DspFloatType> p2  = poles[n++];
            
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
        if(order <= 2) order = 2;
        if(order >  8) order = 8;
        std::complex<DspFloatType> * poles = butter_poles[order-2];
        size_t n = 0;

        if(order %2 != 0) {
            BiquadSection c;
            std::complex<DspFloatType> p1  = poles[n++];
            DspFloatType x1 = abs(p1);
            

            c.z[0] = 0.0;
            c.z[1] = 1.0/x1;
            c.z[2] = 0.0;
            c.p[0] = 1.0;
            c.p[1] = 1/x1;
            c.p[2] = 0.0;

            sos.push_back(c);                        
        }
            
        for(size_t i = n; i < order; i += 2)
        {
            BiquadSection c;
            std::complex<DspFloatType> p1  = poles[n++];
            std::complex<DspFloatType> p2  = poles[n++];
            
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
        if(order <= 2) order = 2;
        if(order >  8) order = 8;
        std::complex<DspFloatType> * poles = butter_poles[order-2];
        size_t n = 0;
        for(size_t i = 1; i < order; i += 2)
        {
            BiquadSection c;
            std::complex<DspFloatType> p1  = poles[n++];        
            std::complex<DspFloatType> p2  = poles[n++];
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
        if(order <= 2) order = 2;
        if(order >  8) order = 8;
        std::complex<DspFloatType> * poles = butter_poles[order-2];
        size_t n = 0;
        for(size_t i = 1; i < order; i += 2)
        {
            BiquadSection c;
            std::complex<DspFloatType> p1  = poles[n++];
            std::complex<DspFloatType> p2  = poles[n++];
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
        std::complex<DspFloatType> * poles = butter_poles[0];

        std::complex<DspFloatType> p1  = poles[0];
        std::complex<DspFloatType> p2  = poles[1];
    
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
        std::complex<DspFloatType> * poles = butter_poles[0];

        std::complex<DspFloatType> p1  = poles[0];
        std::complex<DspFloatType> p2  = poles[1];
    
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
        std::complex<DspFloatType> * bpoles = butter_poles[0];

        std::complex<DspFloatType> p1  = bpoles[0];
        std::complex<DspFloatType> p2  = bpoles[1];
    
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
        std::complex<DspFloatType> * bpoles = butter_poles[0];

        std::complex<DspFloatType> p1  = bpoles[0];
        std::complex<DspFloatType> p2  = bpoles[1];
    
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
        std::complex<DspFloatType> * bpoles = butter_poles[0];

        std::complex<DspFloatType> p1  = bpoles[0];
        std::complex<DspFloatType> p2  = bpoles[1];
    
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

std::complex<DSPFloatType> cheby1_2_gain(0.9826,0);
    std::complex<DspFloatType> cheby1_2_poles[] = {
        std::complex<DspFloatType>(-0.548867,-0.895129),
        std::complex<DspFloatType>(-0.548867,0.895129),
    };
    std::complex<DSPFloatType> cheby1_3_gain(0.4913,0);
    std::complex<DspFloatType> cheby1_3_poles[] = {
        std::complex<DspFloatType>(-0.494171,0.000000),
        std::complex<DspFloatType>(-0.247085,-0.965999),        
        std::complex<DspFloatType>(-0.247085,0.965999),
    };
    std::complex<DSPFloatType> cheby1_4_gain(2.4565e-01,6.1843e-18);
    std::complex<DspFloatType> cheby1_4_poles[] = {
        std::complex<DspFloatType>(-0.139536,-0.983379),
        std::complex<DspFloatType>(-0.336870,-0.407329),
        std::complex<DspFloatType>(-0.336870,0.407329),
        std::complex<DspFloatType>(-0.139536,0.983379),
    };
    std::complex<DSPFloatType> cheby1_5_gain(1.2283e-01,8.6736e-18);
    std::complex<DspFloatType> cheby1_5_poles[] = {
        std::complex<DspFloatType>(-0.089458,-0.990107),
        std::complex<DspFloatType>(-0.234205,-0.611920),
        std::complex<DspFloatType>(-0.289493,0.000000),
        std::complex<DspFloatType>(-0.234205,0.611920),
        std::complex<DspFloatType>(-0.089458,0.990107),
    };
    std::complex<DSPFloatType> cheby1_6_gain(6.1413e-02,-4.6382e-18);
    std::complex<DspFloatType> cheby1_6_poles[] = {
        std::complex<DspFloatType>(-0.062181,-0.993411),
        std::complex<DspFloatType>(-0.169882,-0.727227),
        std::complex<DspFloatType>(-0.232063,-0.266184),
        std::complex<DspFloatType>(-0.232063,0.266184),
        std::complex<DspFloatType>(-0.169882,0.727227),
        std::complex<DspFloatType>(-0.062181,0.993411),
    };
    std::complex<DSPFloatType> cheby1_7_gain(3.0707e-02, 1.3010e-18);
    std::complex<DspFloatType> cheby1_7_poles[] = {
        std::complex<DspFloatType>(-0.205414,0.000000),
        std::complex<DspFloatType>(-0.045709,-0.995284),
        std::complex<DspFloatType>(-0.128074,-0.798156),
        std::complex<DspFloatType>(-0.185072,-0.442943),        
        std::complex<DspFloatType>(-0.185072,0.442943),
        std::complex<DspFloatType>(-0.128074,0.798156),
        std::complex<DspFloatType>(-0.045709,0.995284),
    };
    std::complex<DSPFloatType> cheby1_8_gain(1.5353e-02,8.6967e-19);
    std::complex<DspFloatType> cheby1_8_poles[] = {
        std::complex<DspFloatType>(-0.035008,-0.996451),
        std::complex<DspFloatType>(-0.099695,-0.844751),
        std::complex<DspFloatType>(-0.149204,-0.564444),
        std::complex<DspFloatType>(-0.175998,-0.198206),
        std::complex<DspFloatType>(-0.175998,0.198206),
        std::complex<DspFloatType>(-0.149204,0.564444),
    };
    std::complex<DSPFloatType> cheby1_gains[] = {
        cheby1_2_gain,
        cheby1_3_gain,
        cheby1_4_gain,
        cheby1_5_gain,
        cheby1_6_gain,
        cheby1_7_gain,
        cheby1_8_gain,
    };
    std::complex<DSPFloatType> *cheby1_poles[] = {
        cheby1_2_poles,
        cheby1_3_poles,
        cheby1_4_poles,
        cheby1_5_poles,
        cheby1_6_poles,
        cheby1_7_poles,
        cheby1_8_poles,
    };

/////////////////////////////////////////////////////////////////////////////////////////////
// Chebyshev 1
/////////////////////////////////////////////////////////////////////////////////////////////
    BiquadSOS cheby1lp(int order, double Q=1.0)
    {        
        BiquadSOS sos;
        if(order < 2) order = 2;            
        if(order > 7) order = 8;
        std::complex<DspFloatType> H0  = cheby1_gains[order-2];            
        std::complex<DspFloatType> * poles = cheby1_poles[order-2];
        size_t n = 0;
        if(order %2 != 0) {
            BiquadSection c;
            std::complex<DspFloatType> p1  = poles[n++];
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
        }
            
        for(size_t i = n; i < order; i += 2)
        {
            std::complex<DspFloatType> p1  = poles[n++];
            std::complex<DspFloatType> p2  = poles[n++];
            
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
        BiquadSOS sos;
        if(order < 2) order = 2;            
        if(order > 7) order = 8;
        std::complex<DspFloatType> H0  = cheby1_gains[order-2];            
        std::complex<DspFloatType> * poles = cheby1_poles[order-2];
        size_t n = 0;
        if(order %2 != 0) {
            BiquadSection c;
            std::complex<DspFloatType> p1  = poles[n++];
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
        }
            
        for(size_t i = n; i < order; i += 2)
        {
            std::complex<DspFloatType> p1  = poles[n++];
            std::complex<DspFloatType> p2  = poles[n++];            
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
    std::complex<DspFloatType> cheby2_2_gain(0.8913,0);
    std::complex<DspFloatType> cheby2_2_zeros[] = {
        std::complex<DspFloatType>(0.000000,1.414214),
        std::complex<DspFloatType>(-0.000000,-1.414214),
    };
    std::complex<DspFloatType> cheby2_2_poles[] = {
        std::complex<DspFloatType>(-0.311324,-1.298299),
        std::complex<DspFloatType>(-0.311324,1.298299),
    };
    std::complex<DspFloatType> cheby2_3_gain(5.8957,0);
    std::complex<DspFloatType> cheby2_3_zeros[] = {
        std::complex<DspFloatType>(0.000000,1.154701),
        std::complex<DspFloatType>(-0.000000,-1.154701),
    }
    std::complex<DspFloatType> cheby2_3_poles[] = {
        std::complex<DspFloatType>(-6.106489,-0.000000),
        std::complex<DspFloatType>(-0.105405,-1.129687),        
        std::complex<DspFloatType>(-0.105405,1.129687),
    };
    std::complex<DspFloatType> cheby2_4_gain(0.8913,0);
    std::complex<DspFloatType> cheby2_4_zeros[] = {
        std::complex<DspFloatType>(0.000000,1.082392),
        std::complex<DspFloatType>(0.000000,2.613126),
        std::complex<DspFloatType>(-0.000000,-2.613126),
        std::complex<DspFloatType>(-0.000000,-1.082392),
    };
    std::complex<DspFloatType> cheby2_4_poles[] = {
        std::complex<DspFloatType>(-0.054008,-1.071629),
        std::complex<DspFloatType>(-0.701365,-2.387691),
        std::complex<DspFloatType>(-0.701365,2.387691),
        std::complex<DspFloatType>(-0.054008,1.071629),
    };
    std::complex<DspFloatType> cheby2_5_gain(9.8261,0);
    std::complex<DspFloatType> cheby2_5_zeros[] = {
        std::complex<DspFloatType>(0.000000,1.051462),
        std::complex<DspFloatType>(0.000000,1.701302),
        std::complex<DspFloatType>(-0.000000,-1.701302),
        std::complex<DspFloatType>(-0.000000,-1.051462),
    };
    std::complex<DspFloatType> cheby2_5_poles[] = {
        std::complex<DspFloatType>(-10.206345,-0.000000),
        std::complex<DspFloatType>(-0.033122,-1.045402),
        std::complex<DspFloatType>(-0.223227,-1.663234),        
        std::complex<DspFloatType>(-0.223227,1.663234),
        std::complex<DspFloatType>(-0.033122,1.045402),
    };
    std::complex<DspFloatType> cheby2_6_gain(0.8913,0);
    std::complex<DspFloatType> cheby2_6_zeros[] = {
        std::complex<DspFloatType>(0.000000,1.035276),
        std::complex<DspFloatType>(0.000000,1.414214),
        std::complex<DspFloatType>(0.000000,3.863703),
        std::complex<DspFloatType>(-0.000000,-3.863703),
        std::complex<DspFloatType>(-0.000000,-1.414214),
        std::complex<DspFloatType>(-0.000000,-1.035276),
    };
    std::complex<DspFloatType> cheby2_6_poles[] = {
        std::complex<DspFloatType>(-0.022478,-1.031356),
        std::complex<DspFloatType>(-0.113895,-1.400264),
        std::complex<DspFloatType>(-1.070345,-3.525988),
        std::complex<DspFloatType>(-1.070345,3.525988),
        std::complex<DspFloatType>(-0.113895,1.400264),
        std::complex<DspFloatType>(-0.022478,1.031356),
    };
    std::complex<DspFloatType> cheby2_7_gain(13.757,0);
    std::complex<DspFloatType> cheby2_7_zeros[] = {
        std::complex<DspFloatType>(0.000000,1.025717),
        std::complex<DspFloatType>(0.000000,1.279048),
        std::complex<DspFloatType>(0.000000,2.304765),
        std::complex<DspFloatType>(-0.000000,-2.304765),
        std::complex<DspFloatType>(-0.000000,-1.279048),
        std::complex<DspFloatType>(-0.000000,-1.025717),
    };
    std::complex<DspFloatType> cheby2_7_poles[] = {
        std::complex<DspFloatType>(-14.300043,-0.000000),
        std::complex<DspFloatType>(-0.016288,-1.022959),
        std::complex<DspFloatType>(-0.070763,-1.271995),
        std::complex<DspFloatType>(-0.326203,-2.251897),        
        std::complex<DspFloatType>(-0.326203,2.251897),
        std::complex<DspFloatType>(-0.070763,1.271995),
        std::complex<DspFloatType>(-0.016288,1.022959),
    };
    std::complex<DspFloatType> cheby2_8_gain(0.8913,0);
    std::complex<DspFloatType> cheby2_8_zeros[] = {
        std::complex<DspFloatType>(0.000000,1.019591),
        std::complex<DspFloatType>(0.000000,1.202690),
        std::complex<DspFloatType>(0.000000,1.799952),
        std::complex<DspFloatType>(0.000000,5.125831),
        std::complex<DspFloatType>(-0.000000,-5.125831),
        std::complex<DspFloatType>(-0.000000,-1.799952),
        std::complex<DspFloatType>(-0.000000,-1.202690),
        std::complex<DspFloatType>(-0.000000,-1.019591),
    };
    std::complex<DspFloatType> cheby2_8_poles[] = {
        std::complex<DspFloatType>(-0.012359,-1.017538),
        std::complex<DspFloatType>(-0.048898,-1.198450),
        std::complex<DspFloatType>(-0.162825,-1.781713),
        std::complex<DspFloatType>(-1.435344,-4.675639),
        std::complex<DspFloatType>(-1.435344,4.675639),
        std::complex<DspFloatType>(-0.162825,1.781713),
        std::complex<DspFloatType>(-0.048898,1.198450),
        std::complex<DspFloatType>(-0.012359,1.017538),
    };
    std::complex<DspFloatType> cheby2_gains[] = {
        cheby2_2_gain,
        cheby2_3_gain,
        cheby2_4_gain,
        cheby2_5_gain,
        cheby2_6_gain,
        cheby2_7_gain,
        cheby2_8_gain,
    };
    std::complex<DspFloatType> *cheby2_zeros[] = {
        cheby2_2_zeros,
        cheby2_3_zeros,
        cheby2_4_zeros,
        cheby2_5_zeros,
        cheby2_6_zeros,
        cheby2_7_zeros,
        cheby2_8_zeros,
    };
    std::complex<DspFloatType> *cheby2_poles[] = {
        cheby2_2_poles,
        cheby2_3_poles,
        cheby2_4_poles,
        cheby2_5_poles,
        cheby2_6_poles,
        cheby2_7_poles,
        cheby2_8_poles,
    };

    BiquadSOS cheby2lp(int order, double Q=1.0,double rips=1.0)
    {        
        BiquadSOS sos;  
        if(order < 2) order = 2;
        if(order > 8) order = 8;
        std::complex<DspFloatType> * czeros[] = cheby2_zeros[order-1];          
        std::complex<DspFloatType> * cpoles[] = cheby2_poles[order-1];          
        size_t n = 0;
        if(order %2 != 0) {
            BiquadSection c;
            std::complex<DspFloatType> H0  = cheby2_gains[order-1]*czeros[n];
            std::complex<DspFloatType> p1  = cpoles[n];

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
        }
            
        for(size_t i = n; i < order; i += 2)
        {
            std::complex<DspFloatType> H0  = cheby2_gains[order-1]*czeros[n]; 
            std::complex<DspFloatType> H1  = cheby2_gains[order-1]*czeros[n+1]; 
            std::complex<DspFloatType> p1  = cpoles[n++];
            std::complex<DspFloatType> p2  = cpoles[n++];
    
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

    std::complex<DspFloatType> ellip_2_gain(3.1620e-03,0);
    std::complex<DspFloatType> ellip_2_zeros[] = {
        std::complex<DspFloatType>(0.000000,7.294427),
        std::complex<DspFloatType>(-0.000000,-7.294427),
    };
    std::complex<DspFloatType> ellip_2_poles[] = {
        std::complex<DspFloatType>(-0.115703,0.720179),
        std::complex<DspFloatType>(-0.115703,-0.720179),
    };
    std::complex<DspFloatType> ellip_3_gain(0.017775,0);
    std::complex<DspFloatType> ellip_3_zeros[] = {
        std::complex<DspFloatType>(0.000000,2.280738),
        std::complex<DspFloatType>(-0.000000,-2.280738),
    };
    std::complex<DspFloatType> ellip_3_poles[] = {  
        std::complex<DspFloatType>(-0.117388,0.000000),
        std::complex<DspFloatType>(-0.049807,0.886094),
        std::complex<DspFloatType>(-0.049807,-0.886094),        
    };
    std::complex<DspFloatType> ellip_4_gain(3.1620e-03,0);
    std::complex<DspFloatType> ellip_4_zeros[] = {
        std::complex<DspFloatType>(0.000000,3.035120),
        std::complex<DspFloatType>(0.000000,1.433499),
        std::complex<DspFloatType>(-0.000000,-3.035120),
        std::complex<DspFloatType>(-0.000000,-1.433499),
    };
    std::complex<DspFloatType> ellip_4_poles[] = {  
        std::complex<DspFloatType>(-0.083369,0.450273),
        std::complex<DspFloatType>(-0.022592,0.949851),
        std::complex<DspFloatType>(-0.083369,-0.450273),
        std::complex<DspFloatType>(-0.022592,-0.949851),
    };
    std::complex<DspFloatType> ellip_5_gain(0.013033,0);
    std::complex<DspFloatType> ellip_5_zeros[] = {
        std::complex<DspFloatType>(0.000000,1.594131),
        std::complex<DspFloatType>(0.000000,1.172108),
        std::complex<DspFloatType>(-0.000000,-1.594131),
        std::complex<DspFloatType>(-0.000000,-1.172108),
    };
    std::complex<DspFloatType> ellip_5_poles[] = {  
        std::complex<DspFloatType>(-0.091161,0.000000)
        std::complex<DspFloatType>(-0.049256,0.720872),
        std::complex<DspFloatType>(-0.010192,0.977717),
        std::complex<DspFloatType>(-0.049256,-0.720872),
        std::complex<DspFloatType>(-0.010192,-0.977717),
    };
    std::complex<DspFloatType> ellip_6_gain(3.1628e-03,0);
    std::complex<DspFloatType> ellip_6_zeros[] = {
        std::complex<DspFloatType>(0.000000,2.663135),
        std::complex<DspFloatType>(0.000000,1.226580),
        std::complex<DspFloatType>(0.000000,1.072668),
        std::complex<DspFloatType>(-0.000000,-2.663135),
        std::complex<DspFloatType>(-0.000000,-1.226580),
        std::complex<DspFloatType>(-0.000000,-1.072668),
    };
    std::complex<DspFloatType> ellip_6_poles[] = {  
        std::complex<DspFloatType>(-0.074612,0.401042),
        std::complex<DspFloatType>(-0.025369,0.867234),
        std::complex<DspFloatType>(-0.004554,0.990115),
        std::complex<DspFloatType>(-0.074612,-0.401042),
        std::complex<DspFloatType>(-0.025369,-0.867234),
        std::complex<DspFloatType>(-0.004554,-0.990115),
    };
    std::complex<DspFloatType> ellip_7_gain(0.012329,0);
    std::complex<DspFloatType> ellip_7_zeros[] = {
        std::complex<DspFloatType>(0.000000,1.505495),
        std::complex<DspFloatType>(0.000000,1.094314),
        std::complex<DspFloatType>(0.000000,1.031498),
        std::complex<DspFloatType>(-0.000000,-1.505495),
        std::complex<DspFloatType>(-0.000000,-1.094314),
        std::complex<DspFloatType>(-0.000000,-1.031498),
    }
    std::complex<DspFloatType> ellip_7_poles[] = {  
        std::complex<DspFloatType>(-0.086432,0.000000),
        std::complex<DspFloatType>(-0.047099,0.684688),
        std::complex<DspFloatType>(-0.012071,0.939200),
        std::complex<DspFloatType>(-0.002024,0.995621),
        std::complex<DspFloatType>(-0.047099,-0.684688),
        std::complex<DspFloatType>(-0.012071,-0.939200),
        std::complex<DspFloatType>(-0.002024,-0.995621),        
    };
    std::complex<DspFloatType> ellip_8_gain(3.1628e-03,0);
    std::complex<DspFloatType> ellip_8_zeros[] = {
        std::complex<DspFloatType>(0.000000,2.599225),
        std::complex<DspFloatType>(0.000000,1.196001),
        std::complex<DspFloatType>(0.000000,1.040615),
        std::complex<DspFloatType>(0.000000,1.013796),
        std::complex<DspFloatType>(-0.000000,-2.599225),
        std::complex<DspFloatType>(-0.000000,-1.196001),
        std::complex<DspFloatType>(-0.000000,-1.040615),
        std::complex<DspFloatType>(-0.000000,-1.013796),
    }
    std::complex<DspFloatType> ellip_8_poles[] = {  
        std::complex<DspFloatType>(-0.072877,0.391644),
        std::complex<DspFloatType>(-0.024974,0.847710),
        std::complex<DspFloatType>(-0.005518,0.972693),
        std::complex<DspFloatType>(-0.000896,0.998063),
        std::complex<DspFloatType>(-0.072877,-0.391644),
        std::complex<DspFloatType>(-0.024974,-0.847710),
        std::complex<DspFloatType>(-0.005518,-0.972693),
        std::complex<DspFloatType>(-0.000896,-0.998063),
    };
    std::complex<DspFloatType> ellip_gains[] = {
        ellip_2_gain,
        ellip_3_gain,
        ellip_4_gain,
        ellip_5_gain,
        ellip_6_gain,
        ellip_7_gain,
        ellip_8_gain,
    };
    std::complex<DspFloatType> *ellip_zeros[] = {
        ellip_2_zeros,
        ellip_3_zeros,
        ellip_4_zeros,
        ellip_5_zeros,
        ellip_6_zeros,
        ellip_7_zeros,
        ellip_8_zeros,
    };
    std::complex<DspFloatType> *ellip_poles[] = {
        ellip_2_poles,
        ellip_3_poles,
        ellip_4_poles,
        ellip_5_poles,
        ellip_6_poles,
        ellip_7_poles,
        ellip_8_poles,
    };

    BiquadSOS elliplp(int order, double Q=1.0)
    {        
        BiquadSOS sos;  
        if(order < 2) order = 2;
        if(order > 8) order = 8;
        std::complex<DspFloatType> * czeros[] = ellip_zeros[order-1];          
        std::complex<DspFloatType> * cpoles[] = ellip_poles[order-1];          
        size_t n = 0;
        if(order %2 != 0) {
            BiquadSection c;
            std::complex<DspFloatType> H0  = ellip_gains[order-1]*czeros[n];
            std::complex<DspFloatType> p1  = cpoles[n++];

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
        }
            
        for(size_t i = n; i < order; i += 2)
        {
            std::complex<DspFloatType> H0  = ellip_gains[order-1]*czeros[n]; 
            std::complex<DspFloatType> H1  = ellip_gains[order-1]*czeros[n+1]; 
            std::complex<DspFloatType> p1  = cpoles[n++];
            std::complex<DspFloatType> p2  = cpoles[n++];
    
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
           
struct IPPIIRBiquad: public Casino::IPP::IIRBiquad<DspFloatType>
{
    IPPIIRBiquad() = default;
    IPPIIRBiquad(const BiquadSection &c) {
        setCoefficients(c);
    }
    IPPIIRBiquad(const BiquadSOS & sos) {
        setCoefficients(sos);
    }
    void setCoefficients(const BiquadSection & c)
    {
        DspFloatType buf[6];        
        buf[0] = c.z[0];
        buf[1] = c.z[1];
        buf[2] = c.z[2];
        buf[3] = 1.0;
        buf[4] = c.p[0];
        buf[5] = c.p[1];
        this->initCoefficients(blockSize,1,buf);
    }
    void setCoefficients(const BiquadSOS & sos)
    {        
        DspFloatType buf[6*sos.size()];
        int x = 0;
        for(size_t i = 0; i < sos.size(); i++)
        {    
            buf[x++] = sos[i].z[0];
            buf[x++] = sos[i].z[1];
            buf[x++] = sos[i].z[2];
            buf[x++] = 1.0;
            buf[x++] = sos[i].p[0];
            buf[x++] = sos[i].p[1];
        }
        this->initCoefficients(blockSize,sos.size(),buf);
    }
    void setCoefficients(const FilterCoefficients & c)
    {
        DspFloatType buf[6];        
        buf[0] = c.b[0];
        buf[1] = c.b[1];
        buf[2] = c.b[2];
        buf[3] = 1.0;
        buf[4] = c.a[0];
        buf[5] = c.a[1];
        this->initCoefficients(blockSize,1,buf);
    }    
    void setCoefficients(const std::vector<FilterCoefficients> & c)
    {
        DspFloatType buf[6*c.size()];
        int x = 0;
        for(size_t i = 0; i < c.size(); i++)
        {    
            buf[x++] = c[i].b[0];
            buf[x++] = c[i].b[1];
            buf[x++] = c[i].b[2];
            buf[x++] = 1.0;
            buf[x++] = c[i].a[0];
            buf[x++] = c[i].a[1];
        }
        this->initCoefficients(blockSize,c.size(),buf);
    }    
    void ProcessBlock(size_t n, DspFloatType * in, DspFloatType * out)
    {
        assert(n == this->len);
        this->Execute(in,out);
    }
};
