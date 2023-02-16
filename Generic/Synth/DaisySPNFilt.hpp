#pragma once

#include "DaisySP.hpp"

namespace DaisySP::Filters
{
    // how to calculate the coefficients?
    struct NFilt : public FilterProcessorPlugin<daisysp::NlFilt>
    {
        NFilt() : FilterProcessorPlugin<daisysp::NlFilt>()
        {
            this->Init();
        }

    };
}