#pragma once

#include "MusicFunctions.hpp"
#include "Undenormal.hpp"
#include "GenericSoundObject.hpp"

namespace Analog::Filters::DiodeLadderFilter2
{
    template<typename DSP>
    struct DiodeLadderFilter : public GSSoundProcessor<DSP>
    {
        DSP eta = 1.836f;
        DSP VT = 0.0260f;
        DSP gamma = eta * VT;
        DSP C = 1.0e-7f;
        DSP Mp = 1.0e-4f;
        
        DSP VC1, VC2, VC3, VC4;
        DSP u1, u2, u3, u4, u5;
        DSP s1, s2, s3, s4;
        DSP Vin;
        DSP Vout;
        DSP VoutPrev;
        DSP Fs;
        DSP inputFs;

        int iteration = 0;
        int maxNrIterations = 100;
        DSP biasParameter = 0.5;
        DSP gainParameter = 0.5;

        DiodeLadderFilter(DSP sr = 44100.0f) : GSSoundProcessor<DSP>()
        {
            VC1 = 0.0f;
            VC2 = 0.0f;
            VC3 = 0.0f;
            VC4 = 0.0f;
            u1 = 0.0f;
            u2 = 0.0f;
            u3 = 0.0f;
            u4 = 0.0f;
            u5 = 0.0f;
            s1 = 0.0f;
            s2 = 0.0f;
            s3 = 0.0f;
            s4 = 0.0f;
            Vout = 0.0f;
            VoutPrev = 0.0f;        
            Fs = sr;
            biasParameter = 1.0;
            gainParameter  = 0.5;
        }

        void setCutoff(DSP c) {
            biasParameter = MusicFunctions::cv2freq(c);
        }
        void setResonance(DSP q) {
            gainParameter = q*20*3.459431619;
        }
        enum {
                PORT_CUTOFF,
                PORT_Q,
            };
        void setPort(int port, DSP v) {
            switch(port) {
                case PORT_CUTOFF: setCutoff(v); break;
                case PORT_Q: setResonance(v); break;
                default: printf("No port %d\n", port);
            }
        }
        DSP Tick(DSP I, DSP A=1, DSP X=1, DSP Y=1)
        {
            Undenormal noDenormals;            
                
            auto I0 = 8 * C * VT * 2 * Fs * std::tan((M_PI * biasParameter)/ Fs);
            DSP K = gainParameter;
            Vin = I;
            iteration = 0;
            while (1) {
                u1 = std::tanh((Vin - VoutPrev) / (2 * VT));
                VC1 = (I0 / (4.0 * C * Fs)) * (u2 + u1) + s1;
                u2 = std::tanh((VC2 - VC1) / (2 * gamma));
                VC2 = (I0 / (4.0 * C * Fs) * (u3 - u2)) + s2;
                u3 = std::tanh((VC3 - VC2) / (2 * gamma));
                VC3 = (I0 / (4.0 * C * Fs) * (u4 - u3)) + s3;
                u4 = std::tanh((VC4 - VC3) / (2 * gamma));
                VC4 = (I0 / (4.0 * C * Fs) * (-u5 - u4)) + s4;
                u5 = std::tanh(VC4 / (6.0f * gamma));
                Vout = (K + 0.5f) * VC4;
                if (std::fabs(Vout - VoutPrev) >= Mp * std::fabs(VoutPrev) || iteration > maxNrIterations)
                {
                    VoutPrev = Vout;
                    break;
                }
                VoutPrev = Vout;
                iteration++;
            }
            s1 = 1 / (2 * Fs) * u1 + VC1;
            s2 = 1 / (2 * Fs) * u2 + VC2;
            s3 = 1 / (2 * Fs) * u3 + VC3;
            s4 = 1 / (2 * Fs) * u4 + VC4;
            
            return Vout;
        }    
        void ProcessSIMD(size_t n, DSP * in, DSP * out)
        {            
            Undenormal noDenormals;            
            #pragma omp simd aligned(in,out)
            for(size_t i = 0; i < n; i++) {        
                auto I0 = 8 * C * VT * 2 * Fs * std::tan((M_PI * biasParameter)/ Fs);
                DSP K = gainParameter;
                const DSP I = in[i];
                Vin = I;
                iteration = 0;
                while (1) {
                    u1 = std::tanh((Vin - VoutPrev) / (2 * VT));
                    VC1 = (I0 / (4.0 * C * Fs)) * (u2 + u1) + s1;
                    u2 = std::tanh((VC2 - VC1) / (2 * gamma));
                    VC2 = (I0 / (4.0 * C * Fs) * (u3 - u2)) + s2;
                    u3 = std::tanh((VC3 - VC2) / (2 * gamma));
                    VC3 = (I0 / (4.0 * C * Fs) * (u4 - u3)) + s3;
                    u4 = std::tanh((VC4 - VC3) / (2 * gamma));
                    VC4 = (I0 / (4.0 * C * Fs) * (-u5 - u4)) + s4;
                    u5 = std::tanh(VC4 / (6.0f * gamma));
                    Vout = (K + 0.5f) * VC4;
                    if (std::fabs(Vout - VoutPrev) >= Mp * std::fabs(VoutPrev) || iteration > maxNrIterations)
                    {
                        VoutPrev = Vout;
                        break;
                    }
                    VoutPrev = Vout;
                    iteration++;
                }
                s1 = 1 / (2 * Fs) * u1 + VC1;
                s2 = 1 / (2 * Fs) * u2 + VC2;
                s3 = 1 / (2 * Fs) * u3 + VC3;
                s4 = 1 / (2 * Fs) * u4 + VC4;
                
                out[i] = Vout;
            }
        }
        void ProcessBlock(size_t n, DSP * in, DSP * out)
        {
            ProcessSIMD(n,in,out);
        }
        void ProcessInplace(size_t n, DSP * out)
        {
            ProcessSIMD(n,out,out);
        }
    };
}
