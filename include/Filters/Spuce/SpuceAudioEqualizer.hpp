#pragma once

#include <spuce/filters/audio_equalizer.h>


namespace Filters
{
    
    struct AudioEqualizer : public spuce::audio_equalizer
    {        
        DspFloatType fs;        
        
        AudioEqualizer(DspFloatType sr, size_t N=10) : spuce::audio_equalizer(N)
        {
            fs = sr;            
            this->set_fs(sr);
        }
        
        DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=0, DspFloatType Y=0)
        {
            std::complex<double> in(I,0),out;
            out = spuce::audio_equalizer::run(in);
            return A * abs(out);
        }
    };
}