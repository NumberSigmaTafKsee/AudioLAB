//////////////////////////////////////
// IIR Filters
//////////////////////////////////////


// coefficient generators 

// 1p designer
// 1z designer
// 1p/1z designer
// 2p designer
// 2z designer
// biquad designer
// cascader f(f(x))
// serial f1*f2
// parallel f1+f2
// virtual analog designer 
// all poles designer

// kfr 
// zpk designer
// iir designer
// fir designer
// biquad designer 

// Spuce IIR
// DSPFilters

template<typename T>
T cascade_filter(T input, std::vector<std::function<T (T)>> funcs)
{
    std::reverse(funcs.begin(),funcs.end());
    T a = input;
    for(std::function<T (T)> & f : funcs) {
        a = f(a);
    }
    return a;
}
template<typename T>
T serial_filter(T input, std::vector<std::function<T (T)>> funcs)
{
    std::reverse(funcs.begin(),funcs.end());
    T a = input;
    for(std::function<T (T)> & f : funcs) {
        a *= f(a);
    }
    return a;
}
template<typename T>
T parallel_filter(T input, std::vector<std::function<T (T)>> funcs)
{
    std::reverse(funcs.begin(),funcs.end());
    T a = input;
    for(std::function<T (T)> & f : funcs) {
        a += f(a);
    }
    return a;
}
template<typename T>
T circulator_filter(T input, std::vector<std::function<T (T)>> funcs)
{
    std::reverse(funcs.begin(),funcs.end());
    T a = input;
    for(std::function<T (T)> & f : funcs) {
        a = std::fmod(a,f(a));
    }
    return a;
}


// zero[2]
// pole[2]
// RBJ
// ResonantLP
// DAFX Transfers
// Audio Effects Transfers
// AN Coefficients
template<typename T>
struct TBiquadCoef
{
    T b0,b1,b2;
    T a1,a2;
    T fc,fs,r;
};

template<typename T>
TBiquadCoeff<T> ResonanceLPCoeffs(T fc, T r, T fs)
{
    TBiquadCoeff<T> c;
    T f = fc/fs;
    c.fc = fc;
    c.fs = fs;
    c.r  = r;    
    c.a1 = -2*std::cos(2*M_PI*f);
    c.a2 = r*r;
    c.b0 = (1-r*r)/2;
    c.b1 = 0;
    c.b2 = -c.b0;
    return c;
}
template<typename T>
TBiquadCoeff<T> ResonanceLP2Coeffs(T fc, T Q, T fs)
{
    TBiquadCoeff<T> c;
    T f = fc/fs;
    c.fc = fc;
    c.fs = fs;
    c.r  = r; 
    T Bw = (fc/Q);
    T R = std::exp(-M_PI*Bw/fs);
    c.a1 = -2*R*std::cos(2*M_PI*f);
    c.a2 = r*r;
    c.b0 = (1-r*r)/2;
    c.b1 = 0;
    c.b2 = -R;
    return c;
}
// Stk OnePole
struct TOnePoleFilter
{
    T b0,a1,y1;

    TOnePoleFilter() {
        b0=a1=y1=0;
    }
    TOnePoleFilter(T pole) {
        if(pole > 0.0)
            b0 = 1.0-pole;
        else 
            b0 = 1.0+pole;
        a1 = pole;
        y1 = 0.0;
    }
    

    T pole(T fc, T fs) { return std::exp(-2*M_PI*(fc/fs)); }

    // Low
    // High
    // b0*x + a1*y
    T Tick(T in, T A = 1, T F = 0, T R = 0)
    {
        T y = b0*in - a1*y1;
        y1  = y;
        return y;
    }
};

// Stk OneZero
struct TOneZeroFilter
{
    // Low
    // High
    // b0*x + b1*x1 
    T b0,b1,x1;

    TOneZeroFilter(T _b0, T _b1, T sr=44100) {
        b0 = _b0;
        b1 = _b1;
    }
    TOneZeroFilter(T zero, T sr=44100) {
        if(zero > 0)
            b0 = 1.0 / (1.0 + zero);
        else
            b0 = 1.0 / (1.0 - zero);        
        b1 = zero * b0;
        x1 = 0.0;
    }
    
    T zero(T fc, T fs) { return std::exp(-2*M_PI*(fc/fs)); }

    // Low
    // High
    // b0*x + a1*y
    T Tick(T in, T A = 1, T F = 0, T R = 0)
    {
        T y = b0*in  + b1*x1;
        x1  = in;
        return y;
    }
};

struct TOnePoleOneZeroFilter
{
    // Low 
    // High
    // b0*x + b1*x1 - a*y1
    T b0,b1,a1,x1,y1;

    T Tick(T in, T A = 1, T F = 0, T R = 0)
    {
        T y = b0*in  + b1*x1 - a1 * y1;
        x1  = in;
        y1  = y;
        return y;
    }
    
};

// Stk TwoPole
struct TTwoPoleFilter
{
    // b0*x -a1*y1 -a2*y2
    T b0,a1,a2,x1,x2;

    T Tick(T in, T A = 1, T F = 0, T R = 0)
    {
        T y = b0*in  - a1*y1 - a2 * y2;
        y1  = y2;
        y2 = in;
        return y;
    }
};

// Stk TwoZero
struct TTwoZeroFilter
{
    // b0*x + b1*x1 + b2*x2
    T b0,b1,b2,x1,x2;

    T Tick(T in, T A = 1, T F = 0, T R = 0)
    {
        T y = b0*in  + b1*x1 + b2 * x2;
        x1  = x2;
        x2 = in;
        return y;
    }
};

// Stk Biquad
struct TBiquad 
{
    // b0 + b1*x1 + b2 * x2 - a1*y1 - a1*y2
    T b0,b1,b2,a1,a2;
    T x0,x1,y0,y1;

    TBiquad() 
    {
        b0=b1=b2=b1=b2=0;
        x0=x1=y1=y1=0;
    }
    TBiquad(TBiquadCoeff<T> & c)
    {
        x0=x1=y1=y1=0;
        b0 = c.b0;
        b1 = c.b1;
        b2 = c.b2;
        a1 = c.a1;
        a2 = c.a2;
    }

    T Tick(T in, T A =1 , T F = 0, T R = 0)
    {
        T y = b0*in + b1*x0 + b2 * x1 - a1*y0 - a2*y1;
        x0 = x1;
        x1 = in;
        y0 = y1;
        y1 = in;
        return y;
    }
};

struct FilterCascader
{
// f(f(x))
};

struct FilterSerializer
{
// *
};

struct FilterParallelizer
{
// +
};

struct FilterCircularizer
{
// %
};

