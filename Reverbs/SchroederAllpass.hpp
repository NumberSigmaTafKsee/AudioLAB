// https://github.com/Squishy47/Schroeder-All-Pass-Filter/blob/master/SchroederAllPass.cpp
#pragma once

// USER HEADERS
#include <cstdio>
#include <cmath>
#include <vector>

#include "CircularBuffer.hpp"

namespace FX::Delays
{
    //delay 1.7 - 5 ms
    class SchroederAllPass{
        CircularBuffer CB{44100};
        
        int Fs;
        DspFloatType delayLength;
        DspFloatType g;
        
        DspFloatType delayedSample = 0;
        DspFloatType feedFordwardSample = 0;
    public:
        SchroederAllPass(DspFloatType inValue, DspFloatType inG);
            
        void process(DspFloatType* samples, int bufferSize);
        
        DspFloatType processSingleSample(DspFloatType sample);
        
        void setFeedback(DspFloatType inValue);
        
        DspFloatType getFeedback();
        
        void setDelayLength(DspFloatType inValue);
        
        DspFloatType getDelayLength();
        
        void setFs(int inValue);
    };


    SchroederAllPass::SchroederAllPass(DspFloatType inValue, DspFloatType inG){
        delayLength = inValue;
        g = inG;
    }

    void SchroederAllPass::process(DspFloatType* samples, int bufferSize){
        for(int i = 0; i < bufferSize; i++)
            samples[i] = processSingleSample(samples[i]);
    }

    DspFloatType SchroederAllPass::processSingleSample(DspFloatType sample){
        delayedSample = CB.readCubic(delayLength);
        
        CB.write(sample + (delayedSample * g));

        feedFordwardSample = sample * -g;

        return (delayedSample + feedFordwardSample);
    }

    void SchroederAllPass::setFeedback(DspFloatType inValue){
        g = inValue;
    }

    DspFloatType SchroederAllPass::getFeedback(){
        return g;
    }

    void SchroederAllPass::setDelayLength(DspFloatType inValue){
        delayLength = inValue;
        if (delayLength > CB.getBufferLength())
            CB.setBufferLength(delayLength);
    }

    DspFloatType SchroederAllPass::getDelayLength(){
        return delayLength;
    }

    void SchroederAllPass::setFs(int inValue){
        Fs = inValue;
    }
}