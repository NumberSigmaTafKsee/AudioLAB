#pragma once

// https://github.com/titas2001/Diode_Clipper
// https://github.com/titas2001/Diode_Ring_Modulator

#include "GenericSoundObject.hpp"

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
        
        Dioder(DSP sr=44100.0) : GSSoundProcessor<DSP>()
        {                    
            sampleRate = sr;            
        }    
        DSP Tick(DSP In, DSP V=1, DSP E=1, DSP I=1)
        {
            //return FX::Distortion::Diode::Diode(In,v*V,e*E,i*I);
            return (In*iS) * (exp(0.1*In/((E*eS)*(V*vS)))-1);
        }
        
        #if 0
        void ProcessGPU(size_t n, DSP * in, DSP * out) {

        }
        #endif
        
        void ProcessSIMD(size_t n, DSP * in, DSP * out) {
            #pragma omp simd aligned(in,out)
            for(size_t i = 0; i < n; i++) {				
                out[i] = (in[i]*iS) * (exp(0.1*in[i]/((eS)*(vS)))-1);
            }
        }
    
        void ProcessBlock(size_t n, DSP * in, DSP * out) {
            ProcessSIMD(n,in,out);
        }

        void ProcessInplace(size_t n, DSP * buffer) {
            ProcessSIMD(n,buffer,buffer);
        }
    };    
}
