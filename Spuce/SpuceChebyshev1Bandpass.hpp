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

    struct BandpassChebyshev2
    {
        SpuceIIRFilter filter;
        DspFloatType a = 40.0;
        DspFloatType fc,ff,fs;
        size_t order;
        
        BandpassChebyshev1(size_t O,DspFloatType f1, DspFloatType f2, DspFloatType atten, DspFloatType sr)
        {
            fs = sr;            
            order = O;
            a = atten;
            setCutoff(f1,f2);
        }
        void setCutoff(DspFloatType f1,DspFloatType f2)
        {            
            
            f1 /= fs;
            f2 /= fs;
            fc = f2/f1;
            if(fc >= 1.1) fc = sqrt(f2*f1);
            else fc = (f1+f2)*0.5;
            filter.coef = spuce::iir_coeff(order,spuce::filter_type::bandpass);
            
            filter.coef.set_center(fc);
            spuce::chebyshev1_iir(filter.coef,(f2+f1)/2.0,a);
            filter.filter = spuce::iir_df<DspFloatType>(filter.coef);
        }
        DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=0, DspFloatType Y=0)
        {
            return A*filter.filter.clock(I);
        }
    };
}