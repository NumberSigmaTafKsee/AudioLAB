/*
  ==============================================================================
    mertIIR.cpp
    Created: 22 Apr 2019 9:22:50pm
    Author:  Mert Salahi
  ==============================================================================
*/

#include "mertIIR.h"
#include <cmath>
#include <iostream>

// Class constructor - just intializes coeffecients and filter memory
mertIIR::mertIIR()
{
    b0 = 0.0;
    b1 = 0.0;
    b2 = 0.0;
    
    a0 = 0.0;
    a1 = 0.0;
    a2 = 0.0;
    
    s1 = 0.0;
    s2 = 0.0;
}

mertIIR::~mertIIR()
{
    
}

// Sample-wise implementation of transposed direct-form II IIR
float mertIIR::updateSample(float newSample)
{
    double out = newSample * b0 + s1;
    s1 = newSample * b1 + s2 - a1 * out;
    s2 = newSample * b2 - a2 * out;
    return out;
}

// Lowpass filter - second-order resonator
void mertIIR::updateCoeffecients(double newFcutoff, double newRes)
{
    // Determine poles
    double r = 1 - newRes * std::sin(M_PI * newFcutoff);
    a0 = 1.0;
    a1 = -2.0 * r * std::cos(M_PI * newFcutoff);
    a2 = r*r;

    // Place two zeros at z = 1
    b0 = 1.0;
    b1 = 2.0;
    b2 = 1.0;
    
    // Scale for unity gain at DC.
    // A better approach might be to ensure unity gain at fcutoff
    double dcGain = (b0 + b1 + b2)/(a0 + a1 + a2);
    b0 /= dcGain;
    b1 /= dcGain;
    b2 /= dcGain;
}