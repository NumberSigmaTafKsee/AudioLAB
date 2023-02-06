#pragma once

#include "ATK.hpp"

namespace AudioTK
{
    struct TransistorClassA : public ATKFilter
    {
    protected:
        using Filter = ATK::TransistorClassAFilter<DspFloatType>;
        Filter *filter;
    public:
        TransistorClassA() : ATKFilter() {
            filter = new Filter(Filter::build_standard_filter());
        }
        ~TransistorClassA() {
            if(filter) delete filter;
        }
    };
}