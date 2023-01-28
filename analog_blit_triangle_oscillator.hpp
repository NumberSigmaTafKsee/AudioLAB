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
        FX::Filters::OnePole b1;    
        BlitSquare sqr;
        DspFloatType _out = 0;
        DspFloatType sampleRate=44100.0;

        BlitTriangle(DspFloatType sampleRate=44100) : OscillatorProcessor()
        {
            this->sampleRate = sampleRate;
            b1.setFc(10.0f/sampleRate);            
            setFrequency(440);
        }
        void setFrequency(DspFloatType f)
        {        
            sqr.setFrequency(f);                
        }
        void setDuty(DspFloatType d)
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
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_FREQ: setFrequency(v); break;
                case PORT_DUTY: setDuty(v); break;
                case PORT_RESET: reset(); break;
                default: printf("No port %d\n", port);
            }
        }
        DspFloatType Tick(DspFloatType I=1, DspFloatType A=1, DspFloatType X=0, DspFloatType Y=0) {
            DspFloatType x = sqr.Tick();
            DspFloatType a = 1.0 - 0.01*std::fmin(1,sqr.s1.f/1000.0);
            _out = a*_out + x/sqr.s1.p_;            
            _out -= b1.process(_out);
            return 3*_out;
        }
    };
}