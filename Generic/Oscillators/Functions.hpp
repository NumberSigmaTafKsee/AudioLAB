#pragma once
#include <cmath>

inline DspFloatType dsp_osc_step(DspFloatType freq, unsigned int rate)
{
    return rate ? freq / rate : 0.0;
}

/**
 * Increment the time using the step.
 *   @t: The time.
 *   @step: The step.
 */

inline DspFloatType dsp_osc_inc(DspFloatType t, DspFloatType step)
{
    return fmod(t + step, 1.0);
}

/**
 * Decrement the time using the step.
 *   @t: The time.
 *   @step: The step.
 */

inline DspFloatType dsp_osc_dec(DspFloatType t, DspFloatType step)
{
    return fmod(fmod(t - step, 1.0) + 1.0, 1.0);
}

/**
 * Generate a sine wave data point.
 *   @t: The time.
 *   &returns: The data point.
 */

inline DspFloatType dsp_osc_sine(DspFloatType t)
{
    return sin(2.0 * M_PI * t);
}

/**
 * Generate a square wave data point.
 *   @t: The time.
 *   &returns: The data point.
 */

inline DspFloatType dsp_osc_square(DspFloatType t)
{
    return (t < 0.5) ? 1.0 : -1.0;
}

/**
 * Generate a sawtooth wave data point.
 *   @t: The time.
 *   &returns: The data point.
 */

inline DspFloatType dsp_osc_saw(DspFloatType t)
{
    return fmod(1.0 + 2.0 * t, 2.0) - 1.0;
}

/**
 * Generate a reverse sawtooth wave data point.
 *   @t: The time.
 *   &returns: The data point.
 */

inline DspFloatType dsp_osc_rsaw(DspFloatType t)
{
    return fmod(3.0 - 2.0 * t, 2.0) - 1.0;
}

/**
 * Generate a triangle wave data point.
 *   @t: The time.
 *   &returns: The data point.
 */

inline DspFloatType dsp_osc_tri(DspFloatType t)
{
    return 2.0 * fabs(2.0 * t - 1.0) - 1.0;
}



struct NonBandlimitedOsc
{
    DspFloatType phase;
    DspFloatType freq;
    DspFloatType step;        
    DspFloatType sample_rate;

    DspFloatType (*function)(DspFloatType);
    
    enum Type {
            SIN,
            SAW,
            REVERSE_SAW,
            SQUARE,
            TRIANGLE,
        } _type;

    void SetFunction(Type t) {
        _type = t;
        switch(_type)
        {
            case SIN: function = dsp_osc_sine; break;
            case SAW: function = dsp_osc_saw; break;
            case SQUARE: function = dsp_osc_square; break;
            case TRIANGLE: function = dsp_osc_tri; break;
            case REVERSE_SAW: function = dsp_osc_rsaw; break;
        }
    }

    NonBandlimitedOsc(Type t, DspFloatType f, DspFloatType sr = 44100) {
        freq = f;
        sample_rate = sr;
        step = dsp_osc_step(f,sr);
        _type = t;
        SetFunction(t);
    }

    void Increment() {
        phase = dsp_osc_inc(phase,step);        
    }

    void Decrement() {
        phase = dsp_osc_dec(phase,step);        
    }

    DspFloatType Tick() {            
        Increment();        
        return function(phase);
    }

    DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=0,DspFloatType Y=0) {                
        DspFloatType tp = phase;                
        DspFloatType tf = freq;
        phase = phase + Y;
        if(phase > 1.0) phase -=1.0;
        if(phase < 0) phase += 1.0;
        DspFloatType ts = step;
        step = dsp_osc_step(freq + X,sample_rate);
        DspFloatType r = I*A*Tick();                
        phase = tp;
        freq  = tf;
        step = dsp_osc_step(freq,sample_rate);
        return r;
    }
};

/*
struct LFO : public NonBandlimitedOsc
{        
    LFO( Type type, DspFloatType freq, DspFloatType sampleRate=44100) : NonBandlimitedOsc(type,freq,sampleRate)
    {
        
    }
    ~LFO() = default;
};
*/