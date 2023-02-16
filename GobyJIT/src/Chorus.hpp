#pragma once

#include <vector>
#include <cmath>

namespace FX::Chorus
{
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


    struct MonoChorus : MonoFXProcessor() {

        std::vector<float> delayBuffer[2];

        int delayBufferSamples;
        int delayBufferChannels;
        int delayWritePosition[2];

        float lfoPhase;
        float inverseSampleRate;
        float twoPi;
        
        size_t inputChannels;
        size_t outputChannels;

        float sample_rate;
        float current_delay;
        float current_width;
        float current_depth;
        int   num_voices;
        float current_frequency;

        waveformIndex waveform;
        interpolationIndex interpolation;
        
        MonoChorus( float sampleRate=44100, 
                float maxDelayTime=50.0f,             
                int numChannels=2, 
                float currentDelay=10.0, 
                float currentWidth=10.0, 
                float currentDepth=1.0,
                int numVoices=5,
                float currentFreq = 0.2 ) 
            : MonoFXProcessor()
        {                    
            delayBufferChannels = numChannels;
            inputChannels = numChannels;
            outputChannels= numChannels;

            delayWritePosition[0] = 0;        
            delayWritePosition[1] = 0;        
            lfoPhase = 0.0f;
            sample_rate   = sampleRate;        
            inverseSampleRate = 1.0f / (float)sampleRate;
            twoPi         = 2.0f * M_PI;
            waveform      = waveformSine;
            interpolation = interpolationNearestNeighbour;        
            current_delay = currentDelay/1000.0f;
            current_width = currentWidth/1000.0f;
            current_depth = currentDepth;
            num_voices    = numVoices;
            current_frequency = currentFreq;

            setMaxDelayTime(maxDelayTime);
        }
        void setMaxDelayTime(float t) {
            float maxDelayTime = t/1000.0f;
            delayBufferSamples = (int)(maxDelayTime * sample_rate) + 1;
            if (delayBufferSamples < 1) delayBufferSamples = 1;                  
            delayBuffer[0].resize(delayBufferSamples);
            delayBuffer[1].resize(delayBufferSamples);
        }
        void setCurrentDelay(float ms) {
            current_delay = ms/1000.0f;
        }
        void setWidth(float ms) {
            current_width = ms/1000.0f;
        }
        void setDepth(float d) {
            current_depth = d;
        }
        void setNumVoices(unsigned n) {
            num_voices = n;
        }
        void setLfoFrequency(float f) {
            current_frequency = f;
        }
        float lfo (float phase, int waveform) {
            float out = 0.0f;
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

        // monophonic stereo mix
        double Tick(double I, double A = 1, double X = 0, double Y = 0)
        {
            float in[2];
            float out[2];
            in[0] = in[1] = I;
            out[0] = out[1] = 0.0f;
            Process(1,in,out);
            return 0.5*(out[0] + out[1]);
        }

        void Process(size_t n, const float * input, float * output) {
        
            const int numInputChannels  = delayBufferChannels;
            const int numOutputChannels = delayBufferChannels;
            const int numSamples = n;
            
            //======================================
            
            bool  stereo = (delayBufferChannels == 2);

            int localWritePosition[2];
            float phase;

            for (int channel = 0; channel < numInputChannels; ++channel) {
                const float * channelData = input;
                float       * delayData   = delayBuffer[channel].data();
                
                localWritePosition[channel] = delayWritePosition[channel];            
                phase = lfoPhase;
                
                for (int sample = 0; sample < numSamples; ++sample) 
                {                
                    size_t input_index  = numInputChannels*sample + channel;
                    size_t output_index = numOutputChannels*sample + channel;                
                    const float in      = channelData[input_index];
                    float out = 0.0f;
                    float phaseOffset = 0.0f;
                    float weight;
                                
                    output[output_index] = in;
                    
                    for (int voice = 0; voice < num_voices - 1; ++voice) 
                    {
                        if (stereo && num_voices > 2) {
                            weight = (float)voice / (float)(num_voices - 2);
                            if (channel != 0)
                                weight = 1.0f - weight;
                        } else {
                            weight = 1.0f;
                        }

                        float localDelayTime = (current_delay + current_width * lfo (phase + phaseOffset, waveform)) * (float)sample_rate;
                        float readPosition = fmodf ((float)localWritePosition[channel] - localDelayTime + (float)delayBufferSamples, delayBufferSamples);                    
                        int   localReadPosition = floorf (readPosition);
                        
                        switch (interpolation) {
                            case interpolationNearestNeighbour: {
                                float closestSample = delayData[localReadPosition % delayBufferSamples];
                                out = closestSample;
                                break;
                            }
                            case interpolationLinear: {
                                float fraction = readPosition - (float)localReadPosition;
                                float delayed0 = delayData[(localReadPosition + 0)];
                                float delayed1 = delayData[(localReadPosition + 1) % delayBufferSamples];
                                out = delayed0 + fraction * (delayed1 - delayed0);
                                break;
                            }
                            case interpolationCubic: {
                                float fraction = readPosition - (float)localReadPosition;
                                float fractionSqrt = fraction * fraction;
                                float fractionCube = fractionSqrt * fraction;

                                float sample0 = delayData[(localReadPosition - 1 + delayBufferSamples) % delayBufferSamples];
                                float sample1 = delayData[(localReadPosition + 0)];
                                float sample2 = delayData[(localReadPosition + 1) % delayBufferSamples];
                                float sample3 = delayData[(localReadPosition + 2) % delayBufferSamples];

                                float a0 = - 0.5f * sample0 + 1.5f * sample1 - 1.5f * sample2 + 0.5f * sample3;
                                float a1 = sample0 - 2.5f * sample1 + 2.0f * sample2 - 0.5f * sample3;
                                float a2 = - 0.5f * sample0 + 0.5f * sample2;
                                float a3 = sample1;
                                out = a0 * fractionCube + a1 * fractionSqrt + a2 * fraction + a3;
                                break;
                            }
                        }

                        if (stereo && num_voices == 2)
                            output[output_index] = (channel == 0) ? in : out * current_depth;
                        else
                            output[output_index] += out * current_depth * weight;

                        if (num_voices == 3)
                            phaseOffset += 0.25f;
                        else if (num_voices > 3)
                            phaseOffset += 1.0f / (float)(num_voices - 1);
                    }                                
                    
                    delayData[localWritePosition[channel]] = in;                                
                    if (++localWritePosition[channel] >= delayBufferSamples)
                        localWritePosition[channel] -= delayBufferSamples;

                    phase += current_frequency * inverseSampleRate;
                    if (phase >= 1.0f)
                        phase -= 1.0f;
                }                                    
            }        
            delayWritePosition[0] = localWritePosition[0];
            delayWritePosition[1] = localWritePosition[1];
            lfoPhase = phase;
        }

    };
}