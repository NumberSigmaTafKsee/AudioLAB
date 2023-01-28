#pragma once

#include "ATK.hpp"

namespace Filters::AudioTK::EQ
{  
    struct RemezBasedFilter: public ATKFilter
    {
        using Filter = ATK::FIRFilter<ATK::RemezBasedCoefficients<DspFloatType>>;
        Filter * filter;

        RemezBasedFilter() : ATKFilter()
        {
            filter = new Filter(2);
        }
        ~RemezBasedFilter() {
            if(filter) delete filter;
        }
        /*
        void setTemplate(std::vector<
                            std::pair<
                                std::pair<DspFloatType,DspFloatType>,
                                std::pair<DspFloatType,DspFloatType>>> & target)
        {
            filter->set_template(target);
        }
        */
    };
}