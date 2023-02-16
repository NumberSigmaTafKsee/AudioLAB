#pragma once

namespace DSPFX
{

    struct FxPanLaw
    {
        void linear_pan(DspFloatType control, DspFloatType * l_gain, DspFloatType * r_gain) {
            #ifdef DSPFLOATDOUBLE
            ::linear_panD(control,l_gain,r_gain);
            #else
            ::linear_pan(control,l_gain,r_gain);
            #endif
        }
        void equal_power_3db_pan(DspFloatType control, DspFloatType * l_gain, DspFloatType * r_gain) 
        {
            #ifdef DSPFLOATDOUBLE
            ::equal_power_3dB_panD(control,l_gain,r_gain);
            #else
            ::equal_power_3dB_pan(control,l_gain,r_gain);
            #endif
        }
        void equal_power_6db_pan(DspFloatType control, DspFloatType * l_gain, DspFloatType * r_gain) 
        {
            #ifdef DSPFLOATDOUBLE
            ::equal_power_6dB_panD(control,l_gain,r_gain);
            #else
            ::equal_power_6dB_pan(control,l_gain,r_gain);
            #endif
        }
    };
}
