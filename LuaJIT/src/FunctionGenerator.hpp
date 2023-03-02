#pragma once

#include <cmath>

#define TWOPI 2*M_PI

inline float function_noise() 
{
    // remeber to seed
    float r = std::rand() / (float)RAND_MAX;
    return r * 2 - 1;    
}

struct SineFunction
{
    enum Polarity 
    {
        BIPOLAR,
        POSITIVE,
        NEGATIVE,
    }
    polarity = BIPOLAR;

    SineFunction(float f, float sr = 44100): frequency(f), startphase(0), phase(0.0), sampleRate(sr) {}

    void phaseReset(float phaseIn) {
        // This allows you to set the phase of the oscillator to anything you like.
        phase = phaseIn;
    }

    inline void phaseIncrement() {
        // can just do an fmod
        if (phase >= 1.0f)
        phase -= 1.0f;
        phase += (1.f / (sampleRate / frequency));
    } 

    float function() 
    {
        return sin(phase * 2*M_PI);            
    }

    float phasor(float frequency, float startphase=0, float endphase=1) {
           
        float output = phase;
        
        if (phase < startphase) {
            phase = startphase;
        }
        
        if (phase >= endphase)
            phase = startphase;
        
        phase += ((endphase - startphase) / (sampleRate / frequency));
        
        return (output);
    }

    // I = dc offset
    // X = frequency modulation
    // Y = phase modulation
    float Tick(float I=0, float A = 1, float X = -1, float Y = 1) {
        float r = A*function();
        float p = phase;
        phase = phasor(frequency + 0.5*X*frequency) + Y + phase;
        phase = std::fmod(phase,1);        
        if(polarity == POSITIVE) r = std::abs(r);
        else if(polarity == NEGATIVE) r = -std::abs(r);
        phase = p;
        phaseIncrement();      
        return A*(r+I);
    }

    float sampleRate = 44100;
    float frequency;
    float phase;    
    float startphase;
};

struct PhasorFunction
{
    enum Polarity 
    {
        BIPOLAR,
        POSITIVE,
        NEGATIVE,
    }
    polarity = BIPOLAR;

    PhasorFunction(float f, float sr = 44100) : frequency(f),phase(0.0), sampleRate(sr) {}

    void phaseReset(float phaseIn) {
        // This allows you to set the phase of the oscillator to anything you like.
        phase = phaseIn;
    }
    float phasor(float frequency, float startphase=0, float endphase=1) {
           
        float output = phase;
        
        if (phase < startphase) {
            phase = startphase;
        }
        
        if (phase >= endphase)
            phase = startphase;
        
        phase += ((endphase - startphase) / (sampleRate / frequency));
        
        return (output);
    }

    inline void phaseIncrement() {
        // can just do an fmod
        if (phase >= 1.0f)
        phase -= 1.0f;
        phase += (1.f / (sampleRate / frequency));
    } 

   
    float Tick(float I=0, float A = 1, float X = -1, float Y = 1) {
        float r = A*phase;
        float p = phase;
        phase = phasor(frequency + 0.5*X*frequency) + Y + phase;
        phase = std::fmod(phase,1);        
        if(polarity == POSITIVE) r = std::abs(r);
        else if(polarity == NEGATIVE) r = -std::abs(r);
        phase = p;
        phaseIncrement();      
        return A*r;
    }

    float sampleRate = 44100;
    float frequency;
    float phase;    
};


struct Function
{
    enum Polarity 
    {
        BIPOLAR,
        POSITIVE,
        NEGATIVE,
    }
    polarity = BIPOLAR;

    std::function<float (float)> func;

    Function(float f,float sr,std::function<float (float)> function) 
    : phase(0.0), sampleRate(sr) {
        func = function;
    }

    void phaseReset(float phaseIn) {
        // This allows you to set the phase of the oscillator to anything you like.
        phase = phaseIn;
    }

    inline void phaseIncrement() {
        // can just do an fmod
        if (phase >= 1.0f)
        phase -= 1.0f;
        phase += (1.f / (sampleRate / frequency));
    } 

    float phasor(float frequency, float startphase=0, float endphase=1) {
           
        float output = phase;
        
        if (phase < startphase) {
            phase = startphase;
        }
        
        if (phase >= endphase)
            phase = startphase;
        
        phase += ((endphase - startphase) / (sampleRate / frequency));
        
        return (output);
    }

    
    float Tick(float I=0, float A = 1, float X = 0, float Y = 0) {
        float p = phase;
        phase = phasor(frequency + 0.5*X*frequency) + Y + phase;
        phase = std::fmod(phase,1);
        float r = func(phase);
        phase = p;        
        phaseIncrement();      
        if(polarity == POSITIVE) r = std::abs(r);
        else if(polarity == NEGATIVE) r = -std::abs(r);
        return A*r;
    }

    float sampleRate = 44100;
    float frequency;
    float phase;    
};


struct SineGenerator : public Function
{
    SineGenerator(float freq, float sr = 44100.0f) : Function(freq,sr,[](float phase){ return std::sin(phase * 2.0f*M_PI); }) {}
};
struct CosGenerator : public Function
{
    CosGenerator(float freq, float sr = 44100.0f) : Function(freq,sr,[](float phase){ return std::cos(phase * 2.0f*M_PI); }) {}
};
struct TanGenerator : public Function
{
    TanGenerator(float freq, float sr = 44100.0f) : Function(freq,sr,[](float phase){ return std::tan(phase * 2.0f*M_PI); }) {}
};
struct PhasorGenerator : public Function
{
    PhasorGenerator(float freq, float sr = 44100.0f) : Function(freq,sr,[](float phase){ return phase;}) {}
};
struct SquareGenerator : public Function
{
    SquareGenerator(float freq, float sr = 44100.0f) : Function(freq,sr,[](float phase){ return phase < 0.5f ? -1.0f : 1.0f;}) {}
};
struct SawGenerator : public Function
{
    SawGenerator(float freq, float sr = 44100.0f) : Function(freq,sr,[](float phase){ return phase * 2.0f - 1.0f;}) {}
};
struct TriangleGenerator : public Function
{
    TriangleGenerator(float freq, float sr = 44100.0f) : Function(freq,sr,[](float phase){ 
        if (phase <= 0.5f) {
            return (phase - 0.25f) * 4;
        } else {
            return ((1.0f - phase) - 0.25f) * 4;
        }}) {}
};


struct FunctionGenerator
{
    enum Type
    {
        NOISE,
        SINEWAVE,   
        COSWAVE,
        PHASOR,
        SQUARE,
        PULSE,
        SAW,
        REVERSE_SAW,
        TRIANGLE,
    } 
    type = SINEWAVE;

    enum Polarity
    {
        BIPOLAR,
        POSITIVE,
        NEGATIVE,
    }
    polarity = BIPOLAR;

    FunctionGenerator(float f, float sr=44100) : frequency(f),sampleRate(sr),startphase(0),phase(0.0) {}
        
    void phaseReset(float phaseIn) {
        // This allows you to set the phase of the oscillator to anything you like.
        phase = phaseIn;
    }

    inline void phaseIncrement() {
        if (phase >= 1.0f)
        phase -= 1.0f;
        phase += (1.f / (sampleRate / frequency));
    } 

    float phasor(float frequency, float startphase=0, float endphase=1) {
        
        output = phase;
        
        if (phase < startphase) {
            phase = startphase;
        }
        
        if (phase >= endphase)
            phase = startphase;
        
        phase += ((endphase - startphase) / (sampleRate / frequency));
        
        return (output);
    }
        
    float noise() {
        // White Noise
        // always the same unless you seed it.
        float r = std::rand() / (float)RAND_MAX;
        output = r * 2 - 1;
        return (output);
    }

    float sinewave() {
        return sin(phase * TWOPI);            
    }

    float coswave() {
        return cos(phase * TWOPI);            
    }

    float phasor() {
        return phase;
    }

    float square() {
        return phase < 0.5f ? -1.0f : 1.0f;            
    }

    float pulse() {
        float output;
        if (phase < duty)
            output = -1.f;
        if (phase > duty)
            output = 1.f;
        return output;
    }
    float saw() 
    { 
        return phase * 2.0f - 1.0f; 
    }

    float triangle() {     
        float output;   
        if (phase <= 0.5f) {
            output = (phase - 0.25f) * 4;
        } else {
            output = ((1.0f - phase) - 0.25f) * 4;
        }
        return (output);
    }

    float function() {
        switch(type)
        {
        case SINEWAVE: return sinewave();
        case COSWAVE: return coswave();
        case PHASOR: return phasor();
        case SQUARE: return square();
        case PULSE: return pulse();
        case SAW: return saw();
        // pesuro reverse
        case REVERSE_SAW: return 1-saw();
        }
        return triangle();
    }
  
    float Tick(float I=0, float A = 1, float X = 0, float Y = 0) {        
        float p = phase;
        phase = phasor(frequency + 0.5*X*frequency) + Y + phase;
        phase = std::fmod(phase,1);
        float r = function();
        phase = p;        
        phaseIncrement();      
        if(polarity == POSITIVE) r = std::abs(r);
        else if(polarity == NEGATIVE) r = -std::abs(r);                
        return A*(r+I);
    }

    float sampleRate;
    float frequency;
    float phase;
    float duty = 0.5f;
    float startphase;
    float endphase;
    float output;
    float tri;
};


