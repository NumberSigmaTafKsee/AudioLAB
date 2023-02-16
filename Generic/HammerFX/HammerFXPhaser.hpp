#pragma once

#include "HammerFX.hpp"

namespace FX::HammerFX
{
    struct Phasor
    {
        phasor_params cp;

        Phasor() {
            phasor_create(&cp);
        }
        
        void ProcessBlock(size_t n, DspFloatType * in, DspFloatType * out) {
            phasor_filter_mono(&cp,n,in,out);
        }
    };
}

