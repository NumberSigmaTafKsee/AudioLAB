// Delay

#include "StereoVector.hpp"
#include "Undenormal.hpp"

struct Delay
{
    StereoVector delayBuffer;
    int delayBufferSamples;
    int delayBufferChannels;
    int delayWritePosition;
    float paramDelay;
    float paramFeedback;
    float paramMix;
    int numChannels;
    int numInputChannels;
    int numOutputChannels;
    float currentDelayTime;

    Delay(float sampleRate, size_t num_channels=2, float delayTime=50.0f,float maxDelayTime=1000.0f,float paramDelayTime=5.0f, flat feedback=0.9f, float mix = 1.0f) 
    {
        const double smoothTime = 1e-3;
        paramDelay = paramDelayTime/1000.0f;
        paramFeedback = feedback;
        paramMix = mix;
        float maxDelayTime = paramDelayTime;
        delayBufferSamples = (int)(maxDelayTime/1000.0f * (float)sampleRate) + 1;
        if (delayBufferSamples < 1)
            delayBufferSamples = 1;

        delayBufferChannels = num_channels;
        delayBuffer.resize(delayBufferSamples);
        delayBuffer.clear();
        numChannels = num_channels;    
        delayWritePosition = 0;
        currentDelayTime = delayTime/1000.0f;
    }

    float Tick(float in) {
            Undenormal sphinxter;
            const float in = channelData[sample];
            float out = 0.0f;

            float readPosition    = fmodf ((float)localWritePosition - currentDelayTime + (float)delayBufferSamples, delayBufferSamples);
            int localReadPosition = floorf (readPosition);
            int localWritePosition = delayWritePosition;
            if (localReadPosition != localWritePosition) {
                float fraction = readPosition - (float)localReadPosition;
                float delayed1 = delayData[(localReadPosition + 0)];
                float delayed2 = delayData[(localReadPosition + 1) % delayBufferSamples];
                out = delayed1 + fraction * (delayed2 - delayed1);
                out = in + currentMix * (out - in);
                delayData[localWritePosition] = in + out * currentFeedback;
            }            
            delayWritePosition = localWritePosition;
            return out;
    }

    void Process(StereoBuffer& input, StereoBuffer & output)
    {
        Undenormal sphinxter;

        const int numInputChannels  = numChannels;
        const int numOutputChannels = numChannels;
        const int numSamples = input.size();

        //======================================

        float currentDelayTime = paramDelayTime * (float)getSampleRate();
        float currentFeedback = paramFeedback;
        float currentMix = paramMix;

        int localWritePosition;

        for (int channel = 0; channel < numChannels; ++channel) {
            float* channelData = input[channel].data();
            float* outputData  = output[channel].data();
            float* delayData = delayBuffer[channel].data();
            localWritePosition = delayWritePosition;

            for (int sample = 0; sample < numSamples; ++sample) {
                const float in = channelData[sample];
                float out = 0.0f;

                float readPosition    = fmodf ((float)localWritePosition - currentDelayTime + (float)delayBufferSamples, delayBufferSamples);
                int localReadPosition = floorf (readPosition);

                if (localReadPosition != localWritePosition) {
                    float fraction = readPosition - (float)localReadPosition;
                    float delayed1 = delayData[(localReadPosition + 0)];
                    float delayed2 = delayData[(localReadPosition + 1) % delayBufferSamples];
                    out = delayed1 + fraction * (delayed2 - delayed1);

                    outputData[sample] = in + currentMix * (out - in);
                    delayData[localWritePosition] = in + out * currentFeedback;
                }

                if (++localWritePosition >= delayBufferSamples)
                    localWritePosition -= delayBufferSamples;
            }
        }

        delayWritePosition = localWritePosition;

    }
};



