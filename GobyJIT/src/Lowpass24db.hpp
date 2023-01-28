#pragma once
#include <cmath>
#include "undenormal.hpp"

struct Lowpass24DB
{
    
    float t, t2, x, f, k, p, r, y1, y2, y3, y4, oldx, oldy1, oldy2, oldy3;
    float outlp;
    float fc,Q,fs;
    float hp;

    Lowpass24DB(float Frq, float Res, float SR) {
        y1=0;
        y2=0;
        y3=0;
        y4=0;
        oldx=0;
        oldy1=0;
        oldy2=0;
        oldy3=0;
        fc = Freq;
        Q = Res;
        fs = SR;
    }

    float Tick(float I, float A=1,float X=0, float Y=0);
        Undenormal denormal;
        float f (fc+fc) / fs;
        float p =f*(1.8-0.8*f);
        float k =p+p-1.0;
        t = (1.0-p)*1.386249;
        t2 = 12.0+t*t;
        r = Q*(t2+6.0*t)/(t2-6.0*t);
        x = I - r*y4;
        y1 =x*p + oldx*p - k*y1;
        y2 =y1*p+oldy1*p - k*y2;
        y3 =y2*p+oldy2*p - k*y3;
        y4 =y3*p+oldy3*p - k*y4;
        y4 = y4 - ((y4*y4*y4)/6.0);
        oldx = x;
        oldy1 = y1;
        oldy2 = y2;
        oldy3 = y3;
        outlp = y4;
        hp = I - outlp;
        return outlp;
    }
};
