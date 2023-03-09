#pragma once

#include <spuce/filters/iir_coeff.h>
#include <spuce/filters/chebyshev1_iir.h>
#include <spuce/filters/iir_df.h>


namespace Filters
{
    struct SpuceIIRFilter
    {
        spuce::iir_coeff coef;
        spuce::iir_df<DspFloatType> filter;                
    };

    struct LowpassChebyshev1
    {
        SpuceIIRFilter filter;
        DspFloatType x = 3.0;
        DspFloatType fc,ff,fs;
        size_t order;
        
        LowpassChebyshev1(size_t O,DspFloatType cutoff, DspFloatType sr)
        {
            fs = sr;            
            order = O;
            setCutoff(cutoff);
        }
        void setCutoff(DspFloatType f)
        {            
            ff = f/fs;
            filter.coef = spuce::iir_coeff(order,spuce::filter_type::low);         
            spuce::chebyshev1_iir(filter.coef,ff,x);
            filter.filter = spuce::iir_df<DspFloatType>(filter.coef);
        }
        DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=0, DspFloatType Y=0)
        {
            return A*filter.filter.clock(I);
        }
    };
}