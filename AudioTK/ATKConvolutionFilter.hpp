#pragma once
#include "ATK.hpp"
namespace AudioTK
{
	struct ConvolutionFilter : public ATKFilter
    {
	private:
        using Filter = ATK::ConvolutionFilter<DspFloatType>;
        Filter * filter;
	public:
        ConvolutionFilter(const std::vector<DspFloatType> & ir, size_t split_size = 1) : ATKFilter() {
            filter = new Filter();
            Filter::AlignedScalarVector impulse(ir.size());
            memcpy(impulse.data(),ir.data(),sizeof(DspFloatType)*ir.size());
            filter->set_impulse(std::move(impulse));
            filter->set_split_size(split_size);
        }
        ~ConvolutionFilter() {
            if(filter) delete filter;
        }
    };
 }
