#pragma once

#include <cmath>

#define TWOPI 2*M_PI

namespace Oscillators::Generators
{
    inline DspFloatType function_noise() 
    {
        // remeber to seed
        DspFloatType r = std::rand() / (DspFloatType)RAND_MAX;
        return r * 2 - 1;    
    }

    struct Function : public GeneratorProcessor
    {
        enum Polarity 
        {
            BIPOLAR,
            POSITIVE,
            NEGATIVE,
        }
        polarity = BIPOLAR;

        std::function<DspFloatType (DspFloatType)> func;

        Function(DspFloatType f,DspFloatType sr,std::function<DspFloatType (DspFloatType)> function) 
        : GeneratorProcessor(),
        phase(0.0), sampleRate(sr) {
            func = function;
            frequency = f;
            inc = f/sr;
        }
        void setFrequency(DspFloatType f)
        {
            frequency = f;
            inc = f/sampleRate;
        }

        void phaseReset(DspFloatType phaseIn) {
            // This allows you to set the phase of the oscillator to anything you like.
            phase = phaseIn;
        }

        inline void phaseIncrement() {
            phase = fmod(phase+inc,1);
        }        
        DspFloatType Tick(DspFloatType I=1, DspFloatType A = 1, DspFloatType X = 0, DspFloatType Y = 0) {            
            DspFloatType r = func(phase);                    
            phaseIncrement();      
            if(polarity == POSITIVE) r = std::abs(r);
            else if(polarity == NEGATIVE) r = -std::abs(r);
            return A*r;
        }
        DspFloatType operator()() {
            return Tick();
        }
        DspFloatType sampleRate = 44100;
        DspFloatType frequency;
        DspFloatType phase,inc;    
    };


    struct SineGenerator : public Function
    {
        SineGenerator(DspFloatType freq, DspFloatType sr = 44100.0) : Function(freq,sr,[](DspFloatType phase){ return std::sin(phase * 2.0*M_PI); }) {}
    };
    struct CosGenerator : public Function
    {
        CosGenerator(DspFloatType freq, DspFloatType sr = 44100.0) : Function(freq,sr,[](DspFloatType phase){ return std::cos(phase * 2.0*M_PI); }) {}
    };
    struct TanGenerator : public Function
    {
        TanGenerator(DspFloatType freq, DspFloatType sr = 44100.0) : Function(freq,sr,[](DspFloatType phase){ return std::tan(phase * 2.0*M_PI); }) {}
    };
    struct PhasorGenerator : public Function
    {
        PhasorGenerator(DspFloatType freq, DspFloatType sr = 44100.0) : Function(freq,sr,[](DspFloatType phase){ return phase;}) {}
    };
    struct SquareGenerator : public Function
    {
        SquareGenerator(DspFloatType freq, DspFloatType sr = 44100.0) : Function(freq,sr,[](DspFloatType phase){ return phase < 0.5f ? -1.0 : 1.0;}) {}
    };
    struct SawGenerator : public Function
    {
        SawGenerator(DspFloatType freq, DspFloatType sr = 44100.0) : Function(freq,sr,[](DspFloatType phase){ return phase * 2.0 - 1.0;}) {}
    };
    struct TriangleGenerator : public Function
    {
        TriangleGenerator(DspFloatType freq, DspFloatType sr = 44100.0) : Function(freq,sr,[](DspFloatType phase){ 
            if (phase <= 0.5) {
                return (phase - 0.25) * 4;
            } else {
                return ((1.0 - phase) - 0.25) * 4;
            }}) {}
    };
}


