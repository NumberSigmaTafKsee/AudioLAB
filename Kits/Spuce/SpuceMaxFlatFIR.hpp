#pragma once

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <spuce/filters/design_fir.h>
#include <spuce/filters/fir.h>


namespace Filters
{
    
    struct MaxFlatFIR 
    {        
        std::vector<double> taps;
        spuce::fir<double> fir;
        DspFloatType fs,fu,fl;        
        DspFloatType alpha = 0.1;
        DspFloatType weight=100;

        size_t num_taps;

        enum FirType {
            LP,
            HP,
            BP,
            BS,
        } filter_type = LP;

        MaxFlatFIR(FirType type, size_t n, DspFloatType sr, DspFloatType fl, DspFloatType fu=0, DspFloatType a = 0.1, DspFloatType w = 100)
        {
            fs = sr;            
            num_taps = n;
            filter_type = type;
            alpha = a;
            weight= w;
            setCutoff(fu,fl);
        }
        void setCutoff(DspFloatType fu, DspFloatType fl) {
            switch(filter_type)
            {
                case LP: taps = spuce::design_fir("maxflat","LOW_PASS", num_taps, fl, fu, alpha,weight); break;
                case HP: taps = spuce::design_fir("maxflat","HIGH_PASS", num_taps, fl, fu, alpha,weight); break;
                case BP: taps = spuce::design_fir("maxflat","BAND_PASS", num_taps, fl, fu, alpha,weight); break;
                case BS: taps = spuce::design_fir("maxflat","BAND_STOP", num_taps, fl, fu, alpha,weight); break;
            }
            fir = spuce::fir<double>(taps);
        }
        
        DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=0, DspFloatType Y=0)
        {                    
            return A * fir.clock(I);
        }
    };
}