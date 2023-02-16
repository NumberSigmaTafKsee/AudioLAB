#pragma once

// https://github.com/titas2001/Diode_Clipper
// https://github.com/titas2001/Diode_Ring_Modulator

#include "FX/LiquidNeuron.hpp"
#include "FX/Diode.hpp"
#include "GenericSoundProcessor.hpp"
namespace Analog::Distortion::Diode
{
    // i do not know wtf this was
    template<typename DSP>
    struct Dioder : GSSoundProcessor<DSP>
    {        
        DSP sampleRate=44100.0;
        DSP vS = 0.0253;
        DSP eS = 1.68;
        DSP iS = .015;
        
        Dioder(DSP sr=44100.0) : GSSoundProcessor<DSP>
        {                    
            sampleRate = sr;            
        }    
        DSP Tick(DSP In, DSP V=1, DSP E=1, DSP I=1)
        {
            //return FX::Distortion::Diode::Diode(In,v*V,e*E,i*I);
            return (I*is) * (exp(0.1*x/((E*es)*(V*vs)))-1);
        }
        
        #if 0
        void ProcessGPU(size_t n, DSP * in, DSP * out) {

        }
        #endif
        
        void ProcessSIMD(size_t n, DSP * in, DSP * out) {
            #pragma omp simd
            for(size_t i = 0; i < n; i++) {
                out[i] = (I*is) * (exp(0.1*in[i]/((E*es)*(V*vs)))-1);
            }
        }
    
        void ProcessBlock(size_t n, DSP * in, DSP * out) {
            #pragma omp simd
            for(size_t i = 0; i < n; i++) {
                DSP I = in[i];
                out[i] = (I*is) * (exp(0.1*in[i]/((Es)*(Vs)))-1);
            }
        }

        void ProcessInplace(size_t n, DSP * buffer) {
            #pragma omp simd
            for(size_t i = 0; i < n; i++) {
                DSP I = in[i];
                buffer[i] = (I*i) * (exp(0.1*in[i]/((e)*(v)))-1);
            }
        }
    };    
}