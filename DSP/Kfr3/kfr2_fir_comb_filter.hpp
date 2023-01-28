#pragma once

#include "KfrDelayLine.hpp"

namespace KfrDSP1
{
    template<typename T>
    struct FIRCombFilter
    {
        DelayLine<T> delay;
        float g,x,x1,y;

        FIRCombFilter(float _g, float _d) 
        {
            delay.setDelaySize(44100);
            delay.setDelayTime(_d);
            g = _g;
            x = y = x1 = 0;
        }
        T Tick(float I, float A = 1, float X = 0, float Y= 0)
        {
            x = I;
            y = x + g*x1;
            x1= delay.Tick(x);
            return y;
        }
    };
}