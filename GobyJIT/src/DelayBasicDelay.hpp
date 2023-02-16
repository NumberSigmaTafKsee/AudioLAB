/*
  ==============================================================================

    BasicDelay.h
    Created: 6 Oct 2019 3:21:22pm
    Author:  Francesco Fucci

  ==============================================================================
*/

#pragma once

#include <vector>
#include <cmath>

namespace Delay
{
    struct BasicDelay : public MonoFXProcesor {
        
        BasicDelay(float sr, int num_channels, float delay, float feedback, float mix)
        : MonoFXProcessor(),sampleRate(sr), numChannels(num_channels),delayTime(delay/1000.0f),feedbackAmount(feedback),mixAmount(mix)
        {
            //size_t delayBufferLen = delayTime * sampleRate;
            //if(delayBufferLen < 1) delayBufferLen=1;
            size_t delayBufferLen = sr;
            delayBuffer[0].resize(delayBufferLen);
            delayBuffer[1].resize(delayBufferLen);
            writePosition[0] = 0;
            writePosition[1] = 0;  	
        }
    
        void setDelayTime(float time) 
        {
            delayTime = time/1000.0f;
            size_t delayBufferLen = delayTime * sampleRate;
            if(delayBufferLen < 1) delayBufferLen=1;
            delayBuffer[0].resize(delayBufferLen);
            delayBuffer[1].resize(delayBufferLen);
        }
            
        void processBlock(int channel, size_t n, float *buffer) noexcept;
        
        void ProcessBlock(size_t n, float ** in, float ** out) {
            memcpy(out[0],in[0],n*sizeof(float));
            processBlock(0,n,out[0]);
            memcpy(out[1],in[1],n*sizeof(float));
            processBlock(1,n,out[1]);
        }
        inline void updatePosition(int channel, int size) noexcept{
            writePosition[channel] += size;
            if( writePosition[channel] >= delayBuffer[channel].size()){
                writePosition[channel] -= delayBuffer[channel].size();
            }
        }
        
        size_t numChannels;
        float feedbackAmount,delayTime,mixAmount;
        float sampleRate;
        std::vector<float> delayBuffer[2];
        size_t writePosition[2];
        const float delta = 0.3f;    
    };


    inline void BasicDelay::processBlock(int channel, size_t n, float *buffer) noexcept{
        //const auto* bufferPointer = buffer.getReadPointer(channel);
        //const auto* delayPointer = mDelayBuffer.getReadPointer(channel);

        const int bufferLength 	 = n;
        const int delayBufferLength = delayBuffer[channel].size();
        
        //delaySmoothed[channel].setTargetValue(*timeAmount);
        //feedbackSmoothed[channel].setTargetValue(*feedbackAmount);
        
        auto* dryBuffer 	     = buffer;
        auto* delayWritePointer  = delayBuffer[channel].data();
        
        //float delayTime = delaySmoothed.getNextValue();
        //float feedbackCurrentAmount = feedbackSmoothed.getNextValue();
        
        //const float nextValue = feedbackSmoothed.getNextValue();
        //String feedbackAm(nextValue);
        //Logger::outputDebugString(feedbackAm);
        int localWritePosition = writePosition[channel];
        
        for(auto sample = 0; sample < bufferLength; ++sample){
            const float in = dryBuffer[numChannels*sample + channel];
            const float fullPosition = delayBufferLength + localWritePosition -(sampleRate*delayTime);
            float readPosition = fmodf(fullPosition, delayBufferLength);

            float out = 0.0f;

            int approxReadPosition = floorf(readPosition);
            int approxReadPosition_1 = (approxReadPosition + 1);
            
            if(approxReadPosition_1 >= delayBufferLength){
                approxReadPosition_1 -= delayBufferLength;
            }
            
            // Calculate the fractional delay line (approxReadPosition + 1) % delayBufferLength
            // If the read position equals the write position then the fractional delay is zero 
            if(approxReadPosition != localWritePosition){            
                float fraction = readPosition - static_cast<float>(approxReadPosition);
                float delayP1 = delayWritePointer[approxReadPosition];
                float delayP2 = delayWritePointer[approxReadPosition_1];
                out = delayP1 + fraction*(delayP2- delayP1);
                //Write back to the input buffer
                dryBuffer[numChannels*sample+channel] = in + (out - in) * mixAmount;
                //delayWritePointer[localWritePosition] = in + out*feedbackAmount;
                delayWritePointer[localWritePosition] = in + out*(feedbackAmount);
            }
            
            if(++localWritePosition >= delayBufferLength){
                localWritePosition -= delayBufferLength;
            }
            
            /**readPosition += 1 - delta;
            
            if(readPosition >= delayBufferLength){
                readPosition -= delayBufferLength;
            }**/
        }
        updatePosition(channel,n);
    }
}
