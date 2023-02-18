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
        DspFloatType Iv = 0.0253;
        DspFloatType Ie = 1.68;
        DspFloatType Ii = .015;
        
        Dioder(DspFloatType sr=44100.0) 
        {                    
            sampleRate = sr;            
        }    
        DspFloatType Tick(DspFloatType In, DspFloatType V=1, DspFloatType E=1, DspFloatType I=1)
        {
            //return FX::Distortion::Diode::Diode(In,v*V,e*E,i*I);
            return (I*Ii) * (exp(0.1*In/((E*Ie)*(V*Iv)))-1);
        }
        void ProcessSIMD(size_t n, DspFloatType * in, DspFloatType * out) {
            #pragma omp simd aligned(in,out)
            for(size_t i = 0; i < n; i++) {
                out[i] = Ii * (exp(0.1*in[i]/((Ie)*(Iv)))-1);
            }
        }
        void ProcessBlock(size_t n, DspFloatType * input, DspFloatType * output) {
            ProcessSIMD(n,input,output);
        }
            
        void ProcessInplace(size_t n, DspFloatType * input) {
            ProcessBlock(n,input,input);
        }
    };    
}
