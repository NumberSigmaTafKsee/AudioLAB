//
//  CircularBuffer.hpp
//
//  Created by Squishy on 14/02/2018.
//

#ifndef CircularBuffer_hpp
#define CircularBuffer_hpp

#include <cmath>
#include <vector>

typedef std::vector<float> Vec_Float;
enum Selector{upperBound, lowerBound};
enum InterType{cubic, linear};

class CircularBuffer{
    Vec_Float buffer;
    int bufferLength, head, tail;
	
public:
    CircularBuffer(float inValue);
    
    float read(float numElementsToRead, InterType inValue);
    
    void write(float inValue);
    
    float interpCalcAmount(float inValue, Selector inSelector);

    float cubicInterpolation(double y0, double y1, double y2, double y3, double mu);
    
    float getSample(float inValue);
        
    void setBufferLength(float inValue);
    
    int getBufferLength();
};

#endif /* CircularBuffer_hpp */

//
//  CircularBuffer.cpp
//
//  Created by Squishy on 14/02/2018.
//

#include "CircularBuffer.hpp"

CircularBuffer::CircularBuffer(float inValue){
	setBufferLength(inValue);
}

float CircularBuffer::read(float index, InterType inValue){
	float y0, y1, y2, y3, mu, upper, lower, interpAmount;
	
	switch(inValue){
		case cubic:{
				y0 = interpCalcAmount(index - 1, lowerBound);
				y1 = interpCalcAmount(index, lowerBound);
				y2 = interpCalcAmount(index, upperBound);
				y3 = interpCalcAmount(index + 1, upperBound);

				mu = index - y1;

				return cubicInterpolation(getSample(y0), getSample(y1), getSample(y2), getSample(y3), mu);
			break;
		}
		case linear:{
				upper = interpCalcAmount(index, upperBound);
				lower = interpCalcAmount(index, lowerBound);
				interpAmount = index - lower;
				return (getSample(upper) * interpAmount + (1.0 - interpAmount) * getSample(lower));
			break;
		}
		default:
			return 0.0;
		break;
	}
}

void CircularBuffer::write(float inValue){
    head++;
    buffer[head % bufferLength] = inValue;
}

float CircularBuffer::interpCalcAmount(float inValue, Selector inSelector){
    float upper, lower, isRounded = round(inValue);
    
    if(isRounded > inValue){
        upper = isRounded;
        lower = upper - 1;
    }
    else{
        lower = isRounded;
        upper = lower + 1;
    }
	
	switch(inSelector){
		case upperBound:
			return upper;
			break;
		case lowerBound:
			return lower;
			break;
		default:
			return 0.0;
	}
}

float CircularBuffer::getSample(float inValue){
    tail = head - inValue;
    return buffer[tail % bufferLength];
}

void CircularBuffer::setBufferLength(float inValue){
    bufferLength = inValue;
    buffer.resize(inValue);
}

int CircularBuffer::getBufferLength(){
    return (int)buffer.size();
}

float CircularBuffer::cubicInterpolation(double y0, double y1, double y2, double y3, double mu){
    //Cubic interp taken from: http://paulbourke.net/miscellaneous/interpolation/
    double a0, a1, a2, a3, mu2;
    
    mu2 = mu * mu;
    a0 = y3 - y2 - y0 + y1;
    a1 = y0 - y1 - a0;
    a2 = y2 - y0;
    a3 = y1;
    
    return(a0 * mu * mu2 + a1 * mu2 + a2 * mu + a3);
}



#ifndef SchroederNestedAllPass_hpp
#define SchroederNestedAllPass_hpp

#include <stdio.h>
#include "CircularBuffer.hpp"

//delay 1.7 - 5 ms

class SchroederNestedAllPass{
    CircularBuffer CB1{44100}, CB2{44100};
    
    int Fs;
    float delayLength;
    float g;
public:
    SchroederNestedAllPass(float inValue, float inG);
    
    void process(float* samples, int bufferSize);
    
    float processSingleSample(float sample);
    
    float processSingleSample2(float sample);
    
    void setFeedback(float inValue);
    
    float getFeedback();
    
    void setDelayLength(float inValue);
    
    float getDelayLength();
    
    void setFs(int inValue);
};

//
//  SchroederNestedAllPass.cpp
//
//  Created by Squishy.
//
//  If you use this, please credit me :)

#include "SchroederNestedAllPass.hpp"

SchroederNestedAllPass::SchroederNestedAllPass(float inValue, float inG){
    delayLength = inValue;
    g = inG;
}

void SchroederNestedAllPass::process(float* samples, int bufferSize){
    for(int i = 0; i < bufferSize; i++)
        samples[i] = processSingleSample(samples[i]);
}

float SchroederNestedAllPass::processSingleSample(float sample){
    float delay = processSingleSample2(sample);
    CB1.write(sample + (delay * g));
    
    delay = delay * (1 - (g * g));
    
    float feedforward = sample * -g;
    
    return (delay + feedforward);
}

float SchroederNestedAllPass::processSingleSample2(float sample){
    float delay = CB2.read(delayLength, cubic);
    CB2.write(sample + (delay * g));
    
    delay = delay * (1 - (g * g));
    
    float feedforward = sample * -g;
    
    return (delay + feedforward);
}

void SchroederNestedAllPass::setFeedback(float inValue){
    g = inValue;
}

float SchroederNestedAllPass::getFeedback(){
    return g;
}

void SchroederNestedAllPass::setDelayLength(float inValue){
    delayLength = inValue;
    if (delayLength > CB1.getBufferLength())
        CB1.setBufferLength(delayLength);
    if (delayLength > CB2.getBufferLength())
        CB2.setBufferLength(delayLength);
}

float SchroederNestedAllPass::getDelayLength(){
    return delayLength;
}

void SchroederNestedAllPass::setFs(int inValue){
    Fs = inValue;
}