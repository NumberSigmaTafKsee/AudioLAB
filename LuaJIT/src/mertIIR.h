/*
  ==============================================================================
    mertIIR.h
    Created: 22 Apr 2019 9:22:50pm
    Author:  Mert Salahi
  ==============================================================================
*/

#pragma once
#include <vector>

class mertIIR
{
public:
    // Public methods
    mertIIR();
    ~mertIIR();
    float updateSample(float newSample);
    void updateCoeffecients(double newFcutoff, double newRes);
    
    // Public member variables
    double b0, b1, b2;
    double a0, a1, a2;
    double s1, s2;

private:
};