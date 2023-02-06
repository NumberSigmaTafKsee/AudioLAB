#pragma once 
#include <cmath>


DspFloatType signum(DspFloatType x)
{
    if(x < 0) return -1.0f;
    return 1.0f;
}
DspFloatType clamp(DspFloatType x, DspFloatType min, DspFloatType max) {
    if(x < min) return min;
    if(x > max) return max;
    return x;
}

namespace FX::Distortion
{

    DspFloatType serpent_curve(DspFloatType x, DspFloatType g=1)
    {
        DspFloatType max = 2.0f / (g + 1.0f);
        if(max == 0.0f) max = 1.0f;
        return ((2*x)/(g*x*x+1))/(2*max);
    }

    // the sigmoids
    DspFloatType sigmoider(DspFloatType in, DspFloatType g=1) {    
        DspFloatType sign = in < 0? -1.0f:1.0f;
        //DspFloatType max = (2*sign*(1.0/(1.0+exp(-2*g)))-1);
        //if(max == 0.0f) max = 1.0f;
        DspFloatType r = 2*sign*(1.0 / (1.0 + exp(-2*g*in))-1);
        return r;///max;
    }
    DspFloatType erfmoider(DspFloatType in, DspFloatType g=2)
    {
        return erf(sigmoider(g*in)*g);
    }
    DspFloatType gunderballs(DspFloatType x,DspFloatType g=2) {
        return 2*(2 * std::atan(std::tanh(g*(x/2))))-1;
    }
    DspFloatType algebraballs(DspFloatType x, DspFloatType g= 2) {
        return 4 *( x/(1 + std::pow(std::pow(std::abs(x),g),1/g)));
    }
    DspFloatType algebramoider(DspFloatType x, DspFloatType g= 2) {
        return x / (std::sqrt(1 + (g*x)*(g*x)));
    }
    DspFloatType tanhify(DspFloatType x, DspFloatType g=2) {
        return std::tanh(g*x)/std::tanh(g);
    }
    DspFloatType tanhballs(DspFloatType x, DspFloatType g=2) {
        return std::tanh(g*x/std::sqrt(g+g*x*x));
    }
    DspFloatType tanhmoider(DspFloatType x, DspFloatType g=2) {
        x = sigmoider(x*g);
        return std::tanh(g*x/std::sqrt(g+g*x*x));
    }
    DspFloatType atanballs(DspFloatType x, DspFloatType g=2) {
        return std::atan(g*x/std::sqrt(g+g*x*x));
    }
    DspFloatType atanmoider(DspFloatType x, DspFloatType g=2) {
        x = sigmoider(x);
        return std::atan(g*x/std::sqrt(g+g*x*x));
    }
    DspFloatType bipolar_tanh(DspFloatType x, DspFloatType p = 2.0, DspFloatType m =2.0)
    {
        if(x < 0) return std::tanh(m*x);
        return std::tanh(p*x);
    }
    DspFloatType quadrant_tanh(DspFloatType x, DspFloatType p1 = 2.0, DspFloatType p2 = 2.0, DspFloatType p3 = 2.0, DspFloatType p4 = 2.0) 
    {
        if(x >= 0 && x < M_PI/2) return std::tanh(x*p1);
        if(x >= M_PI/2 && x < M_PI) return std::tanh(x*p2);
        if(x >= M_PI && x < 3*M_PI/4) return std::tanh(x*p3);
        return std::tanh(x*p4);
    }

    struct ClipFunction
    {
        enum {
            SERPENT_CURVE,
            SIGMOIDER,
            ERFMOIDER,
            GUNDERBALLS,
            ALGEBRABALLS,
            ALGEBRAMOIDER,
            TANHIFY,
            TANHBALLS,
            TANHMOIDER,
            ATANBALLS,
            ATANMOIDER,
        };
        int Type = SERPENT_CURVE;
        ClipFunction(int type = SERPENT_CURVE) {
            Type = type;
        }
        DspFloatType Tick(DspFloatType I, DspFloatType A = 1, DspFloatType X = 1, DspFloatType Y = 1)
        {
            switch(Type)
            {
            case SERPENT_CURVE: return serpent_curve(I,A);
            case SIGMOIDER: return sigmoider(I,A);
            case ERFMOIDER: return erfmoider(I,A);
            case GUNDERBALLS: return gunderballs(I,A);
            case ALGEBRABALLS: return algebraballs(I,A);
            case ALGEBRAMOIDER: return alegbramoider(I,A);
            case TANHIFY: return tanhify(I,A);
            case TANHBALLS: return tanhballs(I,A);
            case TANHMOIDER: return tanhmoider(I,A);
            case ATANBALLS: return atanballs(I,A);
            case ATANMOIDER: return atanmoider(I,A);
            }
            return sigmoider(I,A);
        }
    };
    struct MorphClipper
    {
        ClipFunction clip1,clip2;

        MorphClipper(int t1, int t2) : clip1(t1),clip2(t2)
        {

        }
        DspFloatType Tick(DspFloatType I, DspFloatType A =1, DspFloatType X = 1, DspFloatType Y = 1)
        {
            DspFloatType t = clip1.Tick(I,A) 
            return t + (X*Y)*(clip2.Tick(I,A)- t);
        }
    };
}