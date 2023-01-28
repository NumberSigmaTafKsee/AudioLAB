#pragma once

#include <cmath>

namespace FXDSP
{
    template<typename T>
    void linear_pan(T control, T *l_gain, T *r_gain)
    {
        *l_gain = 1.0f - control;
        *r_gain = control;        
    }
    
    template<typename T>
    equal_power_3dB_pan(T control, T *l_gain, T *r_gain)
    {
        *l_gain = std::sin((1.0 - control) * M_PI_2);
        *r_gain = std::sin(control * M_PI_2);        
    }

    template<typename T>
    equal_power_6dB_pan(T control, T *l_gain, T *r_gain)
    {
        *l_gain = std::pow(std::sin((1.0 - control) * M_PI_2), 2.0);
        *r_gain = std::pow(std::sin(control * M_PI_2), 2.0);        
    } 
}