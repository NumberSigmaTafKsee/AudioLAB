#pragma once

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <spuce/filters/design_fir.h>
#include <spuce/filters/fir.h>


namespace Filters
{
    
    struct ComplexSincBandStopFIR 
    {        
        std::vector<std::complex<double>> taps;
        spuce::fir<std::complex<double>> fir;
        DspFloatType fs,fu,fl;        
        DspFloatType alpha = 0.1;
        DspFloatType weight=100;

        size_t num_taps;

        ComplexSincBandStopFIR(size_t n, DspFloatType sr, DspFloatType fl, DspFloatType fu=0, DspFloatType a = 0.1, DspFloatType w = 100)
        {
            fs = sr;            
            num_taps = n;
            filter_type = type;
            alpha = a;
            weight= w;
            setCutoff(fu,fl);
        }
        void setCutoff(DspFloatType fu, DspFloatType fl) {
            
            case BS: taps = spuce::design_complex_fir("sinc","COMPLEX_BAND_STOP", num_taps, fl, fu, alpha,weight); break;            
            fir = spuce::fir<double>(taps);
        }
        
        DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=0, DspFloatType Y=0)
        {                    
            return A * fir.clock(I);
        }
    };
    struct ComplexBandStopFIR 
    {        
        std::vector<std::complex<double>> taps;
        spuce::fir<std::complex<double>> fir;
        DspFloatType fs,fu,fl;        
        DspFloatType alpha = 0.1;
        DspFloatType weight=100;

        size_t num_taps;

        ComplexBandStopFIR(size_t n, DspFloatType sr, DspFloatType fl, DspFloatType fu=0, DspFloatType a = 0.1, DspFloatType w = 100)
        {
            fs = sr;            
            num_taps = n;
            filter_type = type;
            alpha = a;
            weight= w;
            setCutoff(fu,fl);
        }
        void setCutoff(DspFloatType fu, DspFloatType fl) {
            
            case BS: taps = spuce::design_complex_fir("remez","COMPLEX_BAND_STOP", num_taps, fl, fu, alpha,weight); break;            
            fir = spuce::fir<double>(taps);
        }
        
        DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=0, DspFloatType Y=0)
        {                    
            return A * fir.clock(I);
        }
    };
}