
// https://github.com/Squishy47/Schroeder-Nested-All-Pass-Filter
#pragma once

#include <cmath>
#include <vector>

#include "delays_shroeder_circular_buffer"

namespace Delays::Schroeder::Nested
{
    typedef std::vector<float> Vec_Float;
    enum Selector{upperBound, lowerBound};
    enum InterType{cubic, linear};

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
}