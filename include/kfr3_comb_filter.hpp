#pragma once

#include "KfrDelayLine.hpp"

namespace kfr3
{
    template<typename T>
    struct CombFilter
    {
        DelayLine<T> delay[2];
        T x1,y,y1;
        T gain[2];
        T delayTime[2];
        
        enum {
            X_index,
            Y_index,
        };

        CombFilter(T _g1, T _g2, T _d1, T _d2)
        {
            gain[X_index] = _g1;
            gain[Y_index] = _g2;
            delayTime[0] = _d1;
            delayTime[1] = _d2;
            for(size_t i = 0; i < 2; i++)
            {
                delay[i].setDelaySize(44100);
                delay[i].setDelayTime(delayTime[i]);
            }       
            x1=y=y1=0;
        }    

        // X modulation * depth
        // Y modulation * depth
        T Tick(T I, T A=1, T X = 0, T Y=0)
        {
            T x = I;
            y = x + gain[X_index] * x1 - gain[Y_index] * y1;
            x1 = delay[X_index].Tick(x);
            y1 = delay[Y_index].Tick(y);
            return y;
        }
    };
}
