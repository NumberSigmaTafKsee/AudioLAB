#pragma once

#include "ATK.hpp"

namespace AudioTK
{
    
    struct FollowerTransistorClassA : public ATKFilter
    {
    protected:
        using  Filter = ATK::TransistorClassAFilter<DspFloatType>;
        Filter filter;
    public:
        FollowerTransistorClassA() : ATKFilter(), filter(Filter::build_standard_filter())  {
            
        }
        ~FollowerTransistorClassA() {
            
        }
    };
}