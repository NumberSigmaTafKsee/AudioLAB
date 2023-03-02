#pragma once

#include "samples.hpp"

using namespace Sample;

enum waveformIndex {
    waveformSine = 0,
    waveformTriangle,
    waveformSawtooth,
    waveformInverseSawtooth,
};

enum interpolationIndex {
    interpolationNearestNeighbour = 0,
    interpolationLinear,
    interpolationCubic,
};

template<class T>
struct Chorus {

    std::vector<T> delayBuffer;

    int delayBufferSamples;
    int delayBufferChannels;
    int delayWritePosition;

    T lfoPhase;
    T inverseSampleRate;
    T twoPi;
        
    T sample_rate;
    T current_delay;
    T current_width;
    T current_depth;
    int   num_voices;
    T current_frequency;

    waveformIndex waveform;
    interpolationIndex interpolation;
    
    Chorus(int numChannels=2, T sampleRate=44100,T maxDelayTime=50.0f, T currentDelay=10.0, T currentWidth=10.0, T currentDepth=1.0,int numVoices=2,T currentFreq = 0.2 ) {        
        delayBufferSamples = (int)(maxDelayTime * sampleRate) + 1;
        if (delayBufferSamples < 1)delayBufferSamples = 1;        
        delayBufferChannels = numChannels;
        delayWritePosition = 0;
        delayBuffer.resize(delayBufferSamples, numChannels);
        lfoPhase = 0.0f;
        inverseSampleRate = 1.0f / (T)sampleRate;
        twoPi = 2.0f * M_PI;
        waveform = waveformSine;
        interpolation = interpolationNearestNeighbour;
        sample_rate = sampleRate;        
        current_delay = currentDelay;
        current_width = currentWidth;
        current_depth = currentDepth;
        num_voices = numVoices;
        current_frequency = currentFreq;
    }

    T lfo (T phase, int waveform) {
        T out = 0.0f;
        switch (waveform) {
            case waveformSine: {
                out = 0.5f + 0.5f * sinf (twoPi * phase);
                break;
            }
            case waveformTriangle: {
                if (phase < 0.25f)
                    out = 0.5f + 2.0f * phase;
                else if (phase < 0.75f)
                    out = 1.0f - 2.0f * (phase - 0.25f);
                else
                    out = 2.0f * (phase - 0.75f);
                break;
            }
            case waveformSawtooth: {
                if (phase < 0.5f)
                    out = 0.5f + phase;
                else
                    out = phase - 0.5f;
                break;
            }
            case waveformInverseSawtooth: {
                if (phase < 0.5f)
                    out = 0.5f - phase;
                else
                    out = 1.5f - phase;
                break;
            }
        }
        return out;
    }


    void Process(size_t n, const T* input, T* output) {    
        const int numInputChannels  = delayBufferChannels;
        const int numOutputChannels = delayBufferChannels;
        const int numSamples = n/delayBufferChannels;
        bool  stereo = delayBufferChannels == 2;

        int localWritePosition;
        T phase;

        for (int channel = 0; channel < numInputChannels; ++channel) 
        {            
            std::vector<T> channelData = get_channel(channel,n,input);
            std::vector<T> delayData   = get_channel(channel, delayBuffer);

            localWritePosition = delayWritePosition;
            phase = lfoPhase;
            
            for (int sample = 0; sample < numSamples; ++sample) {
                const T in = channelData[sample];
                T out = 0.0f;
                T phaseOffset = 0.0f;
                T weight;

                for (int voice = 0; voice < num_voices - 1; ++voice) {
                    if (stereo && num_voices > 2) {
                        weight = (T)voice / (T)(num_voices - 2);
                        if (channel != 0)
                            weight = 1.0f - weight;
                    } else {
                        weight = 1.0f;
                    }

                    T localDelayTime = (current_delay + current_width * lfo (phase + phaseOffset, waveform)) * (T)sample_rate;
                    T readPosition = fmodf ((T)localWritePosition - localDelayTime + (T)delayBufferSamples, delayBufferSamples);
                    int   localReadPosition = floorf (readPosition);

                    switch (interpolation) {
                        case interpolationNearestNeighbour: {
                            T closestSample = delayData[localReadPosition % delayBufferSamples];
                            out = closestSample;
                            break;
                        }
                        case interpolationLinear: {
                            T fraction = readPosition - (T)localReadPosition;
                            T delayed0 = delayData[(localReadPosition + 0)];
                            T delayed1 = delayData[(localReadPosition + 1) % delayBufferSamples];
                            out = delayed0 + fraction * (delayed1 - delayed0);
                            break;
                        }
                        case interpolationCubic: {
                            T fraction = readPosition - (T)localReadPosition;
                            T fractionSqrt = fraction * fraction;
                            T fractionCube = fractionSqrt * fraction;

                            T sample0 = delayData[(localReadPosition - 1 + delayBufferSamples) % delayBufferSamples];
                            T sample1 = delayData[(localReadPosition + 0)];
                            T sample2 = delayData[(localReadPosition + 1) % delayBufferSamples];
                            T sample3 = delayData[(localReadPosition + 2) % delayBufferSamples];

                            T a0 = - 0.5f * sample0 + 1.5f * sample1 - 1.5f * sample2 + 0.5f * sample3;
                            T a1 = sample0 - 2.5f * sample1 + 2.0f * sample2 - 0.5f * sample3;
                            T a2 = - 0.5f * sample0 + 0.5f * sample2;
                            T a3 = sample1;
                            out = a0 * fractionCube + a1 * fractionSqrt + a2 * fraction + a3;
                            break;
                        }
                    }

                    if (stereo && num_voices == 2)
                        channelData[sample] = (channel == 0) ? in : out * current_depth;
                    else
                        channelData[sample] += out * current_depth * weight;

                    if (num_voices == 3)
                        phaseOffset += 0.25f;
                    else if (num_voices > 3)
                        phaseOffset += 1.0f / (T)(num_voices - 1);
                }

                delayData[localWritePosition] = in;

                if (++localWritePosition >= delayBufferSamples)
                    localWritePosition -= delayBufferSamples;

                phase += current_frequency * inverseSampleRate;
                if (phase >= 1.0f)
                    phase -= 1.0f;
            }
            set_channel(delayData,channel);
            set_channel(channelData,channel);
        }

        delayWritePosition = localWritePosition;
        lfoPhase = phase;
    }

    void Process(size_t n, const std::vector<T> & input, std::vector<T> & output) {
        Process(n,input.data(),output.data());
    }
};