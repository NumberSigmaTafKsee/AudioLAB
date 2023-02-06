#pragma once

#include <spuce/filters/iir_coeff.h>
#include <spuce/filters/elliptic_iir.h>
#include <spuce/filters/iir_df.h>


namespace Filters
{
    struct SpuceIIRFilter
    {
        spuce::iir_coeff coef;
        spuce::iir_df<DspFloatType> filter;                
    };

    struct BandstopElliptic
    {
        SpuceIIRFilter filter;
        DspFloatType x = 3.0;
        DspFloatType fc,fs;
        size_t order;
        
        BandstopElliptic(size_t O, DspFloatType cutoff, DspFloatType w, DspFloatType sr)
        {
            fs = sr;            
            order = O;
            x = w;
            setCutoff(cutoff);
        }
        void setCutoff(DspFloatType f)
        {
            fc = f/fs;
            filter.coef = spuce::iir_coeff(order,spuce::filter_type::bandstop);            
            spuce::elliptic_iir(filter.coef,fc,x);
            filter.filter = spuce::iir_df<DspFloatType>(filter.coef);
        }
        DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=0, DspFloatType Y=0)
        {
            return A*filter.filter.clock(I);
        }
    };
}