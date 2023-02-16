#pragma once

#include "ATK.hpp"

namespace AudioTK
{  
    struct TimeVaryingSecondOrderSVFLowPassFilter: public ATKFilter
    {
        using Filter = ATK::TimeVaryingSecondOrderSVFFilter<ATK::TimeVaryingSecondOrderSVFLowPassCoefficients<DspFloatType>>;

        Filter * filter;

        TimeVaryingSecondOrderSVFLowPassFilter() : ATKFilter()
        {
            filter = new Filter(2);
        }
        ~TimeVaryingSecondOrderSVFLowPassFilter() {
            delete filter;
        }
        enum  {            
            PORT_Q,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {               
                case PORT_Q: filter->set_Q(v); break;
            }
        }
    };
    struct TimeVaryingSecondOrderSVFHighPassFilter: public ATKFilter
    {
        using Filter = ATK::TimeVaryingSecondOrderSVFFilter<ATK::TimeVaryingSecondOrderSVFHighPassCoefficients<DspFloatType>>;
        Filter * filter;

        TimeVaryingSecondOrderSVFHighPassFilter() : ATKFilter()
        {
            filter = new Filter(2);
        }
        ~TimeVaryingSecondOrderSVFHighPassFilter() {
            delete filter;
        }
        enum  {            
            PORT_Q,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {               
                case PORT_Q: filter->set_Q(v); break;
            }
        }
    };
    struct TimeVaryingSecondOrderSVFBandPassFilter: public ATKFilter
    {
        using Filter = ATK::TimeVaryingSecondOrderSVFFilter<ATK::TimeVaryingSecondOrderSVFBandPassCoefficients<DspFloatType>>;
        Filter * filter;

        TimeVaryingSecondOrderSVFBandPassFilter() : ATKFilter()
        {
            filter = new Filter(2);
        }
        ~TimeVaryingSecondOrderSVFBandPassFilter() {
            delete filter;
        }
        enum  {            
            PORT_Q,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {               
                case PORT_Q: filter->set_Q(v); break;
            }
        }
    };
    struct TimeVaryingSecondOrderSVFNotchPassFilter: public ATKFilter
    {
        using Filter = ATK::TimeVaryingSecondOrderSVFFilter<ATK::TimeVaryingSecondOrderSVFNotchCoefficients<DspFloatType>>;
        Filter * filter;

        TimeVaryingSecondOrderSVFNotchPassFilter() : ATKFilter()
        {
            filter = new Filter(2);
        }
        ~TimeVaryingSecondOrderSVFNotchPassFilter() {
            delete filter;
        }
        enum  {            
            PORT_Q,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {                
                case PORT_Q: filter->set_Q(v); break;
            }
        }
    };
    struct TimeVaryingSecondOrderSVFPeakFilter: public ATKFilter
    {
        using Filter = ATK::TimeVaryingSecondOrderSVFFilter<ATK::TimeVaryingSecondOrderSVFPeakCoefficients<DspFloatType>>;
        Filter * filter;

        TimeVaryingSecondOrderSVFPeakFilter() : ATKFilter()
        {
            filter = new Filter(2);
        }
        ~TimeVaryingSecondOrderSVFPeakFilter() {
            delete filter;
        }
        enum  {            
            PORT_Q,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {                
                case PORT_Q: filter->set_Q(v); break;
            }
        }
    };
    struct TimeVaryingSecondOrderSVFBellFilter: public ATKFilter
    {
        using Filter = ATK::TimeVaryingSecondOrderSVFFilter<ATK::TimeVaryingSecondOrderSVFBellCoefficients<DspFloatType>>;
        Filter * filter;

        TimeVaryingSecondOrderSVFBellFilter() : ATKFilter()
        {
            filter = new Filter(2);
        }
        ~TimeVaryingSecondOrderSVFBellFilter() {
            delete filter;
        }
        enum  {            
            PORT_Q,
            PORT_GAIN,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {                
                case PORT_Q: filter->set_Q(v); break;
                case PORT_GAIN: filter->set_gain(v); break;
            }
        }
    };
    struct TimeVaryingSecondOrderSVFLowShelfFilter: public ATKFilter
    {
        using Filter = ATK::TimeVaryingSecondOrderSVFFilter<ATK::TimeVaryingSecondOrderSVFLowShelfCoefficients<DspFloatType>>;
        Filter * filter;

        TimeVaryingSecondOrderSVFLowShelfFilter() : ATKFilter()
        {
            filter = new Filter(2);
        }
        ~TimeVaryingSecondOrderSVFLowShelfFilter() {
            delete filter;
        }
        enum  {            
            PORT_Q,
            PORT_GAIN,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {                
                case PORT_Q: filter->set_Q(v); break;
                case PORT_GAIN: filter->set_gain(v); break;
            }
        }
    };
    struct TimeVaryingSecondOrderSVFHighShelfFilter: public ATKFilter
    {
        using Filter = ATK::TimeVaryingSecondOrderSVFFilter<ATK::TimeVaryingSecondOrderSVFHighShelfCoefficients<DspFloatType>>;
        Filter * filter;

        TimeVaryingSecondOrderSVFHighShelfFilter() : ATKFilter()
        {
            filter = new Filter(2);
        }
        ~TimeVaryingSecondOrderSVFHighShelfFilter() {
            delete filter;
        }
        enum  {            
            PORT_Q,
            PORT_GAIN,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {                
                case PORT_Q: filter->set_Q(v); break;
                case PORT_GAIN: filter->set_gain(v); break;
            }
        }
    };
}