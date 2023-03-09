#pragma once

#include "KfrDelayLine.hpp"

namespace kfr3
{
    template<typename T>
    struct MultiTapIIRCombFilter
    {
        MultiTapDelayLine<T> delay;
        float g,x,y,y1;

        MultiTapIIRCombFilter(float _g, float _d) 
        {
            delay.setDelaySize(44100);
            delay.setDelayTime(_d);
            g = _g;
            x = y = y1 = 0;
        }
        void addTap(float t) {
            delay.addTap(t);
        }
        T Tick(float I, float A = 1, float X = 0, float Y= 0)
        {
            x = I;
            y = x - g*y1;
            y1= delay.Tick(y);
            return y;
        }
    };
}
