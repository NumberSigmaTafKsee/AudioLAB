#pragma once

#include "ATK.hpp"

namespace AudioTK
{  
    struct LinkwitzRileyLowPassFilter: public ATKFilter
    {
        ATK::LinkwitzRileyLowPassCoefficients<DspFloatType> * filter;

        LinkwitzRileyLowPassFilter() : ATKFilter() {
            filter = new ATK::LinkwitzRileyLowPassCoefficients<DspFloatType>(2);
        }
        ~LinkwitzRileyLowPassFilter() {
            delete filter;
        }
        enum {
            PORT_CUTOFF
        };
        void setPort(int port, DspFloatType v) {
            if(port == PORT_CUTOFF) filter->set_cut_frequency(v);
        }
    };
    struct LinkwitzRileyHighPassFilter: public ATKFilter
    {
        ATK::LinkwitzRileyHighPassCoefficients<DspFloatType> * filter;

        LinkwitzRileyHighPassFilter() : ATKFilter() {
            filter = new ATK::LinkwitzRileyHighPassCoefficients<DspFloatType>(2);
        }
        ~LinkwitzRileyHighPassFilter() {
            delete filter;
        }
        enum {
            PORT_CUTOFF
        };
        void setPort(int port, DspFloatType v) {
            if(port == PORT_CUTOFF) filter->set_cut_frequency(v);
        }
    };
    struct LinkwitzRiley4LowPassFilter: public ATKFilter
    {
        ATK::LinkwitzRiley4LowPassCoefficients<DspFloatType> * filter;

        LinkwitzRiley4LowPassFilter() : ATKFilter() {
            filter = new ATK::LinkwitzRiley4LowPassCoefficients<DspFloatType>(2);
        }
        ~LinkwitzRiley4LowPassFilter() {
            delete filter;
        }
        enum {
            PORT_CUTOFF
        };
        void setPort(int port, DspFloatType v) {
            if(port == PORT_CUTOFF) filter->set_cut_frequency(v);
        }
    };
    struct LinkwitzRiley4HighPassFilter: public ATKFilter
    {
        ATK::LinkwitzRiley4HighPassCoefficients<DspFloatType> * filter;

        LinkwitzRiley4HighPassFilter() : ATKFilter() {
            filter = new ATK::LinkwitzRiley4HighPassCoefficients<DspFloatType>(2);
        }
        ~LinkwitzRiley4HighPassFilter() {
            delete filter;
        }
        enum {
            PORT_CUTOFF
        };
        void setPort(int port, DspFloatType v) {
            if(port == PORT_CUTOFF) filter->set_cut_frequency(v);
        }
    };
}