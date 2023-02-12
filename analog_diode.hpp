#pragma once

// https://github.com/titas2001/Diode_Clipper
// https://github.com/titas2001/Diode_Ring_Modulator

#include "FX/LiquidNeuron.hpp"
#include "FX/Diode.hpp"

namespace Analog::Distortion::Diode
{
    // i do not know wtf this was
    struct Dioder
    {        
        DspFloatType sampleRate=44100.0;
        DspFloatType v = 0.0253;
        DspFloatType e = 1.68;
        DspFloatType i = .015;
        
        Dioder(DspFloatType sr=44100.0) 
        {                    
            sampleRate = sr;            
        }    
        DspFloatType Tick(DspFloatType In, DspFloatType V=1, DspFloatType E=1, DspFloatType I=1)
        {
            //return FX::Distortion::Diode::Diode(In,v*V,e*E,i*I);
            return (I*i) * (exp(0.1*x/((E*e)*(V*v)))-1);
        }
        void ProcessSIMD(size_t n, DspFloatType * in, DspFloatType * out) {
            #pragma omp simd
            for(size_t i = 0; i < n; i++) {
                out[i] = (I*i) * (exp(0.1*in[i]/((E*e)*(V*v)))-1);
            }
        }
    };    
}