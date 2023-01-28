/*
  ==============================================================================

    Phasor.h
    Created: 7 Dec 2020 12:42:12pm
    Author:  sierra

  ==============================================================================
*/

#pragma once


#include <cmath>

class Phasor
{
public:
    double getSample()
    {
        double ret = phase / M_PI_z_2;
        phase = std::fmod(phase + phase_inc, TAU); //
        
        return ret;
    }
    void setSampleRate(double v) { sampleRate = v; calculateIncrement(); }
    void setFrequency (double v) { frequency = v; calculateIncrement(); }
    void setPhase (double v) { phase = v; calculateIncrement(); }
    void Reset() { phase = 0.0; }
    
    double getPhase_inc(){ return phase_inc; }
protected:
    void calculateIncrement() { phase_inc = TAU * frequency / sampleRate; }
    double sampleRate = 44100.0;
    double frequency = 1.0;
    double phase = 0.0;
    double phase_inc = 0.0;
    const double TAU = 2*M_PI;
    const double M_PI_z_2 = M_PI/2;
};
