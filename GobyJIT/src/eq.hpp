//
//  Filter.hpp
//  GraphicEQ - AU
//
//  Created by Ryan Miller on 10/9/18.
//

#ifndef Filter_hpp
#define Filter_hpp

#define _USE_MATH_DEFINES
#include <math.h>

#include <stdio.h>

using namespace std;

class Filter{
public:
    //Constructor Definition
    Filter();
    
    //Set methods for Q, f, omega, and alpha
    void setQ(float Q);
    void setFreq(int f);
    void setOmega(int f, float Fs);
    void setAlpha(float w, float Q);
    
    //Get method to use the setOmega value as an input to the setAlpha method
    float getOmega();
    
    //Define filtering methods
    float lowShelfFilter(float input, float gain, int channel);
    float peakFilter(float input, float gain, int channel);
    float highShelfFilter(float input, float gain, int channel);
    
private:
    //Define the private variables to be used by the filter class
    float Q;
    int f;
    float gain; //Pass in a pointer to the gain AudioParameterFloat for each filter
    int Fs;
    
    //Define omega, alpha, and A variables used to calculate a/b coefficients
    float w; float alpha; float A;
    
    //Define the a/b coefficients that the filter class will use
    float b0; float b1; float b2;
    float a0; float a1; float a2;
    
    //Define input and output float variables. These represent individual sample values from a buffer that get passed into the filter and processed
    float input; float output;
    
    //Delay variables
    float delay1[2] = {0.0f}; float delay2[2] = {0.0f};
};


#endif /* Filter_hpp */