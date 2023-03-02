#pragma once

#include "ATK.hpp"

namespace AudioTK
{
    struct AllPassReverb : public ATKFilter
    {
        using Filter = ATK::AllPassReverbFilter<DspFloatType>;
        Filter * filter;

        AllPassReverb(size_t max, DspFloatType fbk) : ATKFilter() {
            filter = new Filter(max);
            filter->set_feedback(fbk);
        }
        ~AllPassReverb() {
            if(filter) delete filter;
        }
        enum {
            PORT_DELAY,
            PORT_FEEDBACK
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_DELAY: filter->set_delay(v); break;
                case PORT_FEEDBACK: filter->set_feedback(v); break;
            }
        }
    };
    struct LowPassReverb : public ATKFilter
    {
        using Filter = ATK::LowPassReverbFilter<DspFloatType>;
        Filter * filter;

        LowPassReverb(size_t max, DspFloatType fbk, DspFloatType cutoff) : ATKFilter() {
            filter = new Filter(max);
            filter->set_feedback(fbk);
            filter->set_cutoff(cutoff);
        }
        ~LowPassReverb() {
            if(filter) delete filter;
        }
        enum {
            PORT_DELAY,
            PORT_FEEDBACK,
            PORT_CUTOFF,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_DELAY: filter->set_delay(v); break;
                case PORT_FEEDBACK: filter->set_feedback(v); break;
                case PORT_CUTOFF: filter->set_cutoff(v); break;
            }
        }
    };
    
}
