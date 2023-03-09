#pragma once

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <spuce/filters/farrow_upsampler.h>


namespace Filters
{
    struct FarrowUpsampler
    {                
        spuce::farrow_upsampler<double> farrow;
        DspFloatType fs;
        size_t N;

        FarrowUpsampler(DspFloatType sr, size_t n)
        {
            fs = sr;            
            N = n;            
            farrow = spuce::farrow_upsampler(N);
        }        
        DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=0, DspFloatType Y=0)
        {                    
            farrow.need_sample(time_inc,I);
            return A * farrow.output(I);            
        }
    };
}