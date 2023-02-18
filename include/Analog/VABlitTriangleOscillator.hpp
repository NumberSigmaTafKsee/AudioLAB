#pragma once

#include "SoundObject.hpp"
#include "FX/OnePole.hpp"
#include "FX/Filters.h"
#include "VABlitSquareOscillator.hpp"

namespace Analog::Oscillators
{
    ///////////////////////////////////////////////////////////////////////////////////////
    // triangle integrates the square
    ///////////////////////////////////////////////////////////////////////////////////////
    struct BlitTriangle : public OscillatorProcessor
    {
        FX::Filters::OnePole b1,b2;    
        Analog::Oscillators::BlitSquare s1;
        DspFloatType _out = 0;
        DspFloatType sampleRate=44100.0;

        BlitTriangle(DspFloatType freq=220.0,DspFloatType sampleRate=44100) : OscillatorProcessor()
        {
            this->sampleRate = sampleRate;
            b1.setFc(10.0f/sampleRate);
            b2.setFc(10.0f/sampleRate);
            setFrequency(freq);
            setDuty(0.5);
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
            if(_out < 0) _out *= 0.5;        
            _out *= 2.0;         
            return 4*((_out-b2.process(_out)+1)-0.5);
        }
        void ProcessSIMD(size_t n, DspFloatType * in, DspFloatType * out)
        {
            #pragma omp simd aligned(in,out)
            for(size_t i = 0; i < n; i++)
            {
                DspFloatType r1 = s1.Tick();                        
                // there's a tremendous amount of weird dc noise in this thing
                r1   -= b1.process(r1);
                // not really sure why it works but it does I think m_ = harmonic * harmonic like the fourier expansion
                _out += (r1/s1.s1.m_);                            
                out[i] =  4*((_out-b2.process(_out)+1)-0.5);
                if(out[i] < 0) out[i] *= 0.5;
                out[i] *= 2.0;
                if(in) out[i] *= in[i];
            }
        }

        void ProcessBlock(size_t n, DspFloatType * input, DspFloatType * output) {
            ProcessSIMD(n,input,output);
        }
            
        void ProcessInplace(size_t n, DspFloatType * input) {
            ProcessBlock(n,nullptr,input);
        }
    };
}
