#pragma once

#include "ATK.hpp"

namespace AudioTK
{
    
    // Feedback Delay Network    
    struct MonoFDNDelay : public ATKMonoFilter
    {
        using Filter = ATK::FeedbackDelayNetworkFilter<ATK::HadamardMixture<DspFloatType,1>> ;
        Filter* filter;
        MonoFDNDelay(size_t max_delay) :         
        ATKMonoFilter()
        {
            filter = new Filter(max_delay);
            this->setFilter(filter);
        }
        ~MonoFDNDelay() {
            if(filter) delete filter;
        }
        enum {
            PORT_DELAY,
            PORT_INGAIN,
            PORT_FEEDBACK,
            PORT_OUTGAIN,
        };
        void setPort(size_t port, DspFloatType value)  {
            switch(port)
            {
                case PORT_DELAY: filter->set_delay(0,value); break;
                case PORT_INGAIN: filter->set_ingain(0,value); break;
                case PORT_FEEDBACK: filter->set_feedback(0,value); break;
                case PORT_OUTGAIN: filter->set_outgain(0,value); break;

            }
        }
    };

    // Feedback Delay Network    
    struct FDNDelay : public ATKFilter
    {
        using Filter = ATK::FeedbackDelayNetworkFilter<ATK::HadamardMixture<DspFloatType,2>> ;
        Filter* filter;
        FDNDelay(size_t max_delay) :         
        ATKFilter()
        {
            filter = new Filter(max_delay);
            this->setFilter(filter);
        }
        ~FDNDelay() {
            if(filter) delete filter;
        }
        enum {
            PORT_DELAY,
            PORT_INGAIN,
            PORT_FEEDBACK,
            PORT_OUTGAIN,
        };
        void setPort(size_t port, DspFloatType value)  {
            switch(port)
            {
                case PORT_DELAY: 
                filter->set_delay(0,value); 
                filter->set_delay(1,value); 
                break;
                case PORT_INGAIN: 
                filter->set_ingain(0,value); 
                filter->set_ingain(1,value); 
                break;
                case PORT_FEEDBACK: 
                filter->set_feedback(0,value); 
                filter->set_feedback(1,value); 
                break;
                case PORT_OUTGAIN: 
                filter->set_outgain(0,value); 
                filter->set_outgain(1,value); 
                break;

            }
        }
    };
}