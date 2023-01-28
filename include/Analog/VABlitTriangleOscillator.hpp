#pragma once

#include "SoundObject.hpp"
#include "FX/OnePole.hpp"
#include "FX/Filters.h"


namespace Analog::Oscillators
{
    ///////////////////////////////////////////////////////////////////////////////////////
    // triangle integrates the square
    ///////////////////////////////////////////////////////////////////////////////////////
    struct BlitTriangle : public OscillatorProcessor
    {
        FX::Filters::OnePole b1,b2;    
        BlitSquare s1;
        DspFloatType _out = 0;
        DspFloatType sampleRate=44100.0;

        BlitTriangle(DspFloatType sampleRate=44100) : OscillatorProcessor()
        {
            this->sampleRate = sampleRate;
            b1.setFc(10.0f/sampleRate);
            b2.setFc(10.0f/sampleRate);
            setFrequency(440);
        }
        void setFrequency(DspFloatType f)
        {        
            s1.setFrequency(f);                
        }
        void setDuty(DspFloatType d)
        {    
            s1.setDuty(d);
        }
        void reset() {
            s1.reset();
            _out = 0;
        }
        enum {
            PORT_FREQ,
            PORT_DUTY,
            PORT_RESET,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_FREQ: setFrequency(v); break;
                case PORT_DUTY: setDuty(v); break;
                case PORT_RESET: reset(); break;
                default: printf("No port %d\n", port);
            }
        }
        DspFloatType Tick(DspFloatType I=1, DspFloatType A=1, DspFloatType X=0, DspFloatType Y=0) {
            DspFloatType r1 = s1.Tick();                        
            // there's a tremendous amount of weird dc noise in this thing
            r1   -= b1.process(r1);
            // not really sure why it works but it does I think m_ = harmonic * harmonic like the fourier expansion
            _out += (r1/s1.s1.m_);                
            DspFloatType x = _out;        
            return 2*(x-b2.process(x));
        }
    };
}