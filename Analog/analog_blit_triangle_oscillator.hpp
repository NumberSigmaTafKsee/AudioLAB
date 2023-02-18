#pragma once

#include "GenericSoundObject.hpp"
#include "FX/OnePole.hpp"
#include "FX/Filters.h"
#include "analog_blit_square_oscillator.hpp"

namespace Analog::Oscillators
{
    ///////////////////////////////////////////////////////////////////////////////////////
    // triangle integrates the square
    ///////////////////////////////////////////////////////////////////////////////////////
    template<typename DSP>
    struct BlitTriangle : public GSSoundProcessor<DSP>
    {
        FX::Filters::OnePole b1;    
        Analog::Oscillators::BlitSquare<DSP> sqr;
        DSP _out = 0;
        DSP sampleRate=44100.0;

        BlitTriangle(DSP sampleRate=44100) : GSSoundProcessor<DSP>()
        {
            this->sampleRate = sampleRate;
            b1.setFc(10.0f/sampleRate);            
            setFrequency(440);
        }
        void setFrequency(DSP f)
        {        
            sqr.setFrequency(f);                
        }
        void setDuty(DSP d)
        {    
            sqr.setDuty(d);
        }
        void reset() {
            sqr.reset();
            _out = 0;
        }
        enum {
            PORT_FREQ,
            PORT_DUTY,
            PORT_RESET,
        };
        void setPort(int port, DSP v) {
            switch(port) {
                case PORT_FREQ: setFrequency(v); break;
                case PORT_DUTY: setDuty(v); break;
                case PORT_RESET: reset(); break;
                default: printf("No port %d\n", port);
            }
        }
        DSP Tick(DSP I=1, DSP A=1, DSP X=0, DSP Y=0) {
            DSP x = sqr.Tick();
            DSP a = 1.0 - 0.01*std::fmin(1,sqr.s1.f/1000.0);
            _out = a*_out + x/sqr.s1.p_;            
            _out -= b1.process(_out);
            return 3*_out;
        }
        void ProcessSIMD(size_t n, DSP * in, DSP * out) {
            #pragma omp simd aligned(in,out)
            for(size_t i = 0; i < n; i++)
            {
                DSP x = sqr.Tick();
                DSP a = 1.0 - 0.01*std::fmin(1,sqr.s1.f/1000.0);
                _out = a*_out + x/sqr.s1.p_;            
                _out -= b1.process(_out);
                out[i] =  3*_out;
                if(in) out[i] *= in[i];
            }
        } 
        void ProcessBlock(size_t n, DSP * in, DSP * out) {
            ProcessSIMD(n,in,out);
        }
        void ProcessInplace(size_t n, DSP * out) {
            ProcessSIMD(n,nullptr,out);
        }
    };
}
