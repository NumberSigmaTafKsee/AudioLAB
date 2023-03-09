#pragma once

#include "GenericSoundObject.hpp"
#include "FX/OnePole.hpp"
#include "FX/Filters.h"

#include "analog_blit_saw_oscillator.hpp"

namespace Analog::Oscillators
{
    ///////////////////////////////////////////////////////////////////////////////////////
    // square is made from subtracting out of phase sawtooth waves
    ///////////////////////////////////////////////////////////////////////////////////////
    template<typename DSP>
    struct BlitSquare : public GSSoundProcessor<DSP>
    {
        FX::Filters::OnePole block;
        Analog::Oscillators::BlitSaw<DSP> s1,s2;
        DSP _out = 0;
        DSP _duty = 0.5;
        DSP sampleRate=44100.0;

        BlitSquare(DSP sampleRate=44100) : GSSoundProcessor<DSP>()
        {
            this->sampleRate = sampleRate;
            block.setFc(10.0f/sampleRate);
            _out = 0;
            _duty = 0.5;
            setFrequency(440.0f);
            s1.setGain(1);
            s2.setGain(1);
        }
        void setFrequency(DSP f)
        {
            s1.setFrequency(f);
            s2.setFrequency(f);        
            s2.setPhaseOffset(s1.getPhase() + _duty*M_PI);
        }
        void setDuty(DSP d)
        {
            _duty = d;
        }
        void reset() {
            s1.reset();
            s2.reset();
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
            DSP r1 = s1.Tick();                    
            DSP r2 = s2.Tick();                
            _out = r1-r2;            
            return _out;
        }
        void ProcessSIMD(size_t n, DSP * in, DSP * out) {
            #pragma omp simd aligned(in,out)
            for(size_t i = 0; i < n; i++)
            {
                DSP r1 = s1.Tick();                    
                DSP r2 = s2.Tick();                
                _out = r1-r2;            
                out[i] = _out;
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
