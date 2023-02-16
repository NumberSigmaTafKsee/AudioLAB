//
//  Filter.cpp
//  GraphicEQ - AU
//
//  Created by Ryan Miller on 10/9/18.
//

#include "Filter.hpp"
#include "PluginProcessor.h"

using namespace std;

//Constructor Implementation. 
Filter::Filter(){}

//Set Methods for Q, f, omega, and alpha
void Filter::setQ(float Q){
    this->Q = Q;
}

void Filter::setFreq(int f){
    this->f = f;
}

void Filter::setOmega(int f, float Fs){
    w = 2*M_PI*f/Fs;
}

void Filter::setAlpha(float w, float Q){
    alpha = sin(w)/(2*Q);
}

//Get method to access omega in order to use to calculate alpha
float Filter::getOmega(){
    return w;
}

//Low Shelf Filter Method
float Filter::lowShelfFilter(float input, float gain, int channel){
    A = sqrt(pow(10.0,(gain)/20.0));
    b0 = A * ((A+1) - (A-1)*cos(w) + 2*sqrt(A)*alpha);
    b1 = 2*A * ((A-1) - (A+1)*cos(w));
    b2 = A * ((A+1) - (A-1)*cos(w) - 2*sqrt(A)*alpha);
    a0 = (A+1) + ((A-1)*cos(w)) + 2*sqrt(A)*alpha;
    a1 = -2 * ((A-1) + (A+1)*cos(w));
    a2 = (A+1) + (A-1)*cos(w) - 2*sqrt(A)*alpha;
    
    output = (b0/a0)*input + delay1[channel];
    delay1[channel] = (b1/a0)*input + (-a1/a0)*output + delay2[channel];
    delay2[channel] = (b2/a0)*input + (-a2/a0)*output;
    
    return output;
}

//Peak Filter Method
float Filter::peakFilter(float input, float gain, int channel){
    A = sqrt(pow(10.0, (gain)/20.0));
    b0 = 1 + alpha*A;
    b1 = -2 * cos(w);
    b2 = 1 - alpha*A;
    a0 = 1 + (alpha/A);
    a1 = -2 * cos(w);
    a2 = 1 - (alpha/A);
    
    output = (b0/a0)*input + delay1[channel];
    delay1[channel] = (b1/a0)*input + (-a1/a0)*output + delay2[channel];
    delay2[channel] = (b2/a0)*input + (-a2/a0)*output;
    
    return output;
}

//High Shelf Filter Method
float Filter::highShelfFilter(float input, float gain, int channel){
    A = sqrt(pow(10.0, (gain)/20.0));
    b0 = A * ((A+1) + (A-1)*cos(w) + 2*sqrt(A)*alpha);
    b1 = -2*A * ((A-1) + (A+1)*cos(w));
    b2 = A * ((A+1) + (A-1)*cos(w) - 2*sqrt(A)*alpha);
    a0 = (A+1) - ((A-1)*cos(w)) + 2*sqrt(A)*alpha;
    a1 = 2 * ((A-1) - (A+1)*cos(w));
    a2 = (A+1) - (A-1)*cos(w) - 2*sqrt(A)*alpha;
    
    output = (b0/a0)*input + delay1[channel];
    delay1[channel] = (b1/a0)*input + (-a1/a0)*output + delay2[channel];
    delay2[channel] = (b2/a0)*input + (-a2/a0)*output;
    
    return output;
}