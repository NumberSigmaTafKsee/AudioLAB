https://github.com/Squishy47/Schroeder-All-Pass-Filter
#pragma once

#include <cmath>
#include <vector>
#include <cstdio>

#include "delays_shroeder_circular_buffer"
namespace Delays::Schreoder::AllPass
{
    

    //delay 1.7 - 5 ms

    class SchroederAllPass{
        CircularBuffer CB{44100};
        
        int Fs;
        float delayLength;
        float g;
        
        float delayedSample = 0;
        float feedFordwardSample = 0;
    public:
        SchroederAllPass(float inValue, float inG);
            
        void process(float* samples, int bufferSize);
        
        float processSingleSample(float sample);
        
        void setFeedback(float inValue);
        
        float getFeedback();
        
        void setDelayLength(float inValue);
        
        float getDelayLength();
        
        void setFs(int inValue);
    };

    SchroederAllPass::SchroederAllPass(float inValue, float inG){
        delayLength = inValue;
        g = inG;
    }

    void SchroederAllPass::process(float* samples, int bufferSize){
        for(int i = 0; i < bufferSize; i++)
            samples[i] = processSingleSample(samples[i]);
    }

    float SchroederAllPass::processSingleSample(float sample){
        delayedSample = CB.readCubic(delayLength);
        
        CB.write(sample + (delayedSample * g));

        feedFordwardSample = sample * -g;

        return (delayedSample + feedFordwardSample);
    }

    void SchroederAllPass::setFeedback(float inValue){
        g = inValue;
    }

    float SchroederAllPass::getFeedback(){
        return g;
    }

    void SchroederAllPass::setDelayLength(float inValue){
        delayLength = inValue;
        if (delayLength > CB.getBufferLength())
            CB.setBufferLength(delayLength);
    }

    float SchroederAllPass::getDelayLength(){
        return delayLength;
    }

    void SchroederAllPass::setFs(int inValue){
        Fs = inValue;
    }
}