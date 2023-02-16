#pragma once

#include "ATK.hpp"

namespace AudioTK
{  
    struct TimeVaryingBandPassFilter: public ATKFilter
    {   
        ATK::TimeVaryingIIRFilter<ATK::TimeVaryingBandPassCoefficients<DspFloatType>> * filter;

        TimeVaryingBandPassFilter(): ATKFilter()
        {
            filter = new ATK::TimeVaryingIIRFilter<ATK::TimeVaryingBandPassCoefficients<DspFloatType>>();
        }
        ~TimeVaryingBandPassFilter() {
            if(filter) delete filter;
        }
         enum {
            PORT_MIN,
            PORT_MAX,
            PORT_Q,
            //PORT_GAIN,
            PORT_STEPS,
            PORT_MEMORY
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_MIN: filter->set_min_frequency(v); break;
                case PORT_MAX: filter->set_max_frequency(v); break;
                case PORT_Q: filter->set_Q(v); break;
                //case PORT_GAIN: filter->set_gain(v); break;
                case PORT_STEPS: filter->set_number_of_steps(v); break;
                case PORT_MEMORY: filter->set_memory(v); break;
            }
        }
    };
    struct TimeVaryingLowPassFilter: public ATKFilter
    {   
        ATK::TimeVaryingIIRFilter<ATK::TimeVaryingLowPassCoefficients<DspFloatType>> * filter;

        TimeVaryingLowPassFilter(): ATKFilter()
        {
            filter = new ATK::TimeVaryingIIRFilter<ATK::TimeVaryingLowPassCoefficients<DspFloatType>>();
        }
        ~TimeVaryingLowPassFilter() {
            if(filter) delete filter;
        }
        enum {
            PORT_MIN,
            PORT_MAX,
            //PORT_Q,
            //PORT_GAIN,
            PORT_STEPS,
            PORT_MEMORY
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_MIN: filter->set_min_frequency(v); break;
                case PORT_MAX: filter->set_max_frequency(v); break;
                //case PORT_Q: filter->set_Q(v); break;
                //case PORT_GAIN: filter->set_gain(v); break;
                case PORT_STEPS: filter->set_number_of_steps(v); break;
                case PORT_MEMORY: filter->set_memory(v); break;
            }
        }
    };
    struct TimeVaryingHighPassFilter: public ATKFilter
    {   
        ATK::TimeVaryingIIRFilter<ATK::TimeVaryingHighPassCoefficients<DspFloatType>> * filter;

        TimeVaryingHighPassFilter(): ATKFilter()
        {
            filter = new ATK::TimeVaryingIIRFilter<ATK::TimeVaryingHighPassCoefficients<DspFloatType>>();
        }
        ~TimeVaryingHighPassFilter() {
            if(filter) delete filter;
        }
         enum {
            PORT_MIN,
            PORT_MAX,
            //PORT_Q,
            //PORT_GAIN,
            PORT_STEPS,
            PORT_MEMORY
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_MIN: filter->set_min_frequency(v); break;
                case PORT_MAX: filter->set_max_frequency(v); break;
                //case PORT_Q: filter->set_Q(v); break;
                //case PORT_GAIN: filter->set_gain(v); break;
                case PORT_STEPS: filter->set_number_of_steps(v); break;
                case PORT_MEMORY: filter->set_memory(v); break;
            }
        }
    };
    struct TimeVaryingBandPassPeakFilter: public ATKFilter
    {   
        ATK::TimeVaryingIIRFilter<ATK::TimeVaryingBandPassPeakCoefficients<DspFloatType>> * filter;

        TimeVaryingBandPassPeakFilter(): ATKFilter()
        {
            filter = new ATK::TimeVaryingIIRFilter<ATK::TimeVaryingBandPassPeakCoefficients<DspFloatType>>();
        }
        ~TimeVaryingBandPassPeakFilter() {
            if(filter) delete filter;
        }
         enum {
            PORT_MIN,
            PORT_MAX,
            PORT_Q,
            PORT_GAIN,
            PORT_STEPS,
            PORT_MEMORY
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_MIN: filter->set_min_frequency(v); break;
                case PORT_MAX: filter->set_max_frequency(v); break;
                case PORT_Q: filter->set_Q(v); break;
                //case PORT_GAIN: filter->set_gain(v); break;
                case PORT_STEPS: filter->set_number_of_steps(v); break;
                case PORT_MEMORY: filter->set_memory(v); break;
            }
        }
    };
    struct TimeVaryingAllPassFilter: public ATKFilter
    {   
        ATK::TimeVaryingIIRFilter<ATK::TimeVaryingAllPassCoefficients<DspFloatType>> * filter;

        TimeVaryingAllPassFilter(): ATKFilter()
        {
            filter = new ATK::TimeVaryingIIRFilter<ATK::TimeVaryingAllPassCoefficients<DspFloatType>>();
        }
        ~TimeVaryingAllPassFilter() {
            if(filter) delete filter;
        }
         enum {
            PORT_MIN,
            PORT_MAX,
            PORT_Q,
            PORT_GAIN,
            PORT_STEPS,
            PORT_MEMORY
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_MIN: filter->set_min_frequency(v); break;
                case PORT_MAX: filter->set_max_frequency(v); break;
                case PORT_Q: filter->set_Q(v); break;
                //case PORT_GAIN: filter->set_gain(v); break;
                case PORT_STEPS: filter->set_number_of_steps(v); break;
                case PORT_MEMORY: filter->set_memory(v); break;
            }
        }
    };
    struct TimeVaryingLowShelfFilter: public ATKFilter
    {   
        ATK::TimeVaryingIIRFilter<ATK::TimeVaryingLowShelvingCoefficients<DspFloatType>> * filter;

        TimeVaryingLowShelfFilter(): ATKFilter()
        {
            filter = new ATK::TimeVaryingIIRFilter<ATK::TimeVaryingLowShelvingCoefficients<DspFloatType>>();
        }
        ~TimeVaryingLowShelfFilter() {
            if(filter) delete filter;
        }
         enum {
            PORT_MIN,
            PORT_MAX,
            //PORT_Q,
            PORT_GAIN,
            PORT_STEPS,
            PORT_MEMORY
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_MIN: filter->set_min_frequency(v); break;
                case PORT_MAX: filter->set_max_frequency(v); break;
                //case PORT_Q: filter->set_Q(v); break;
                case PORT_GAIN: filter->set_gain(v); break;
                case PORT_STEPS: filter->set_number_of_steps(v); break;
                case PORT_MEMORY: filter->set_memory(v); break;
            }
        }
    };
    struct TimeVaryingHighShelfFilter: public ATKFilter
    {   
        ATK::TimeVaryingIIRFilter<ATK::TimeVaryingHighShelvingCoefficients<DspFloatType>>* filter;

        TimeVaryingHighShelfFilter(): ATKFilter()
        {
            filter = new ATK::TimeVaryingIIRFilter<ATK::TimeVaryingHighShelvingCoefficients<DspFloatType>>();
        }
        ~TimeVaryingHighShelfFilter() {
            if(filter) delete filter;
        }
         enum {
            PORT_MIN,
            PORT_MAX,
            //PORT_Q,
            PORT_GAIN,
            PORT_STEPS,
            PORT_MEMORY
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_MIN: filter->set_min_frequency(v); break;
                case PORT_MAX: filter->set_max_frequency(v); break;
                //case PORT_Q: filter->set_Q(v); break;
                case PORT_GAIN: filter->set_gain(v); break;
                case PORT_STEPS: filter->set_number_of_steps(v); break;
                case PORT_MEMORY: filter->set_memory(v); break;
            }
        }
    };
}