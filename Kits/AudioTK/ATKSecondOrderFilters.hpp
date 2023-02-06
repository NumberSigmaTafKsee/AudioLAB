#pragma once

#include "ATK.hpp"

namespace AudioTK
{  
    struct SecondOrderBandPassFilter: public ATKFilter
    {
        ATK::IIRFilter<ATK::SecondOrderBandPassCoefficients<DspFloatType>>* filter;

        SecondOrderBandPassFilter() : ATKFilter()
        {
            filter = new ATK::IIRFilter<ATK::SecondOrderBandPassCoefficients<DspFloatType>>(2);
        }
        ~SecondOrderBandPassFilter() {
            delete filter;
        }
        enum {
            PORT_CUTOFF,
            PORT_Q,         
        };
        void setPort(int port, DspFloatType v) {
            switch(port)
            {
                case PORT_CUTOFF: filter->set_cut_frequency(v); break;
                case PORT_Q: filter->set_Q(v); break;
            }
        }
    };
    struct SecondOrderLowPassFilter: public ATKFilter
    {
        ATK::IIRFilter<ATK::SecondOrderLowPassCoefficients<DspFloatType>>* filter;

        SecondOrderLowPassFilter() : ATKFilter()
        {
            filter = new ATK::IIRFilter<ATK::SecondOrderLowPassCoefficients<DspFloatType>>(2);
        }
        ~SecondOrderLowPassFilter() {
            delete filter;
        }
        enum {
            PORT_CUTOFF,
            PORT_Q,         
        };
        void setPort(int port, DspFloatType v) {
            switch(port)
            {
                case PORT_CUTOFF: filter->set_cut_frequency(v); break;
                //case PORT_Q: filter->set_Q(v); break;
            }
        }
    };
    struct SecondOrderHighPassFilter: public ATKFilter
    {
        ATK::IIRFilter<ATK::SecondOrderHighPassCoefficients<DspFloatType>>* filter;

        SecondOrderHighPassFilter() : ATKFilter()
        {
            filter = new ATK::IIRFilter<ATK::SecondOrderHighPassCoefficients<DspFloatType>>(2);
        }
        ~SecondOrderHighPassFilter() {
            delete filter;
        }
        enum {
            PORT_CUTOFF,
            PORT_Q,         
        };
        void setPort(int port, DspFloatType v) {
            switch(port)
            {
                case PORT_CUTOFF: filter->set_cut_frequency(v); break;
                //case PORT_Q: filter->set_Q(v); break;
            }
        }
    };
    struct SecondOrderBandPassPeakFilter: public ATKFilter
    {
        ATK::IIRFilter<ATK::SecondOrderBandPassPeakCoefficients<DspFloatType>>* filter;

        SecondOrderBandPassPeakFilter() : ATKFilter()
        {
            filter = new ATK::IIRFilter<ATK::SecondOrderBandPassPeakCoefficients<DspFloatType>>(2);
        }
        ~SecondOrderBandPassPeakFilter() {
            delete filter;
        }
        enum {
            PORT_CUTOFF,
            PORT_Q,    
            PORT_GAIN,     
        };
        void setPort(int port, DspFloatType v) {
            switch(port)
            {
                case PORT_CUTOFF: filter->set_cut_frequency(v); break;
                case PORT_Q: filter->set_Q(v); break;
                case PORT_GAIN: filter->set_gain(v); break;
            }
        }
    };
    struct SecondOrderAllPassFilter: public ATKFilter
    {
        ATK::IIRFilter<ATK::SecondOrderAllPassCoefficients<DspFloatType>>* filter;

        SecondOrderAllPassFilter() : ATKFilter()
        {
            filter = new ATK::IIRFilter<ATK::SecondOrderAllPassCoefficients<DspFloatType>>(2);
        }
        ~SecondOrderAllPassFilter() {
            delete filter;
        }
        enum {
            PORT_CUTOFF,
            PORT_Q,         
        };
        void setPort(int port, DspFloatType v) {
            switch(port)
            {
                //case PORT_CUTOFF: filter->set_cut_frequency(v); break;
                case PORT_Q: filter->set_Q(v); break;
            }
        }
    };
    struct SecondOrderLowShelfFilter: public ATKFilter
    {
        ATK::IIRFilter<ATK::SecondOrderLowShelvingCoefficients<DspFloatType>>* filter;

        SecondOrderLowShelfFilter() : ATKFilter()
        {
            filter = new ATK::IIRFilter<ATK::SecondOrderLowShelvingCoefficients<DspFloatType>>(2);
        }
        ~SecondOrderLowShelfFilter() {
            delete filter;
        }
        enum {
            PORT_CUTOFF,
            PORT_Q,    
            PORT_GAIN,     
        };
        void setPort(int port, DspFloatType v) {
            switch(port)
            {
                case PORT_CUTOFF: filter->set_cut_frequency(v); break;
                //case PORT_Q: filter->set_Q(v); break;
                case PORT_GAIN: filter->set_gain(v); break;
            }
        }
    };
    struct SecondOrderHighShelfFilter: public ATKFilter
    {
        ATK::IIRFilter<ATK::SecondOrderHighShelvingCoefficients<DspFloatType>>* filter;

        SecondOrderHighShelfFilter() : ATKFilter()
        {
            filter = new ATK::IIRFilter<ATK::SecondOrderHighShelvingCoefficients<DspFloatType>>(2);
        }
        ~SecondOrderHighShelfFilter() {
            delete filter;
        }
        enum {
            PORT_CUTOFF,
            PORT_Q,    
            PORT_GAIN,     
        };
        void setPort(int port, DspFloatType v) {
            switch(port)
            {
                case PORT_CUTOFF: filter->set_cut_frequency(v); break;
                //case PORT_Q: filter->set_Q(v); break;
                case PORT_GAIN: filter->set_gain(v); break;
            }
        }
    };
}
    