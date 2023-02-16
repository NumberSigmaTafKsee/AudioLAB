/*
  ==============================================================================
    BodeShifter.h
    Created: 24 May 2020 12:49:44pm
    Author:  mw
  ==============================================================================
*/
#pragma once
#include <vector>


class QuadratureOscillator 
{
public:
    QuadratureOscillator(float frequency, unsigned int sample_rate = 44100)
    {
        setFrequency(frequency, sample_rate);
    }

    void setSampleRate(unsigned int sample_rate)
    {
        setFrequency(frequency_, sample_rate_);
    }

    void setFrequency(float frequency, unsigned int sample_rate) 
    {
        frequency_ = frequency;
        sample_rate_ = sample_rate;

        cyclesPerSample_ = frequency / (float)sample_rate;
        angleDelta_ = cyclesPerSample_ * 2.f * M_PI;
    }

    void setFrequency(float frequency)
    {
        cyclesPerSample_ = frequency / (float)sample_rate_;
        angleDelta_ = cyclesPerSample_ * 2.f * M_PI;
    }

    std::pair<float, float> processSample() 
    {  
        currentAngle_ += angleDelta_;
        fmod(currentAngle_, cyclesPerSample_);
        return std::make_pair(sinf(currentAngle_), cosf(currentAngle_));
    }

private:
    float frequency_;
    unsigned int sample_rate_;
    float cyclesPerSample_;
    float angleDelta_;
    float currentAngle_ = 0.f;
};

class IIRLowPass
{
public:
    IIRLowPass(float cutoff_frequency, unsigned int sample_rate = 44100) :
    previous_(0.f),
    sample_rate_(sample_rate),
    cutoff_frequency_(cutoff_frequency)
    {
        setCutoffFrequency(cutoff_frequency);
    }

    void setSampleRate(unsigned int sample_rate)
    {
        sample_rate_ = sample_rate;
        setCutoffFrequency(cutoff_frequency_);
    }

    void setCutoffFrequency(float cutoff_frequency)
    {
        cutoff_frequency_ = cutoff_frequency;
        float delta_t = 1.f / sample_rate_;
        float numerator = 2.f * M_PI * delta_t * cutoff_frequency;
        alpha_ = numerator / (numerator + 1.f);
    }

    float processSample(float sample)
    {
        float out = previous_ + alpha_ * (sample - previous_);
        previous_ = sample;
        return out;
    }

private:
    float alpha_;
    float previous_;
    unsigned int sample_rate_;
    float cutoff_frequency_;
};

class BodeShifter
{
public:
    BodeShifter(unsigned int sample_rate = 44100) :
        filter_1_(sample_rate, sample_rate),
        filter_2_(sample_rate, sample_rate),
        oscillator_1_(sample_rate/2, sample_rate),
        oscillator_2_(sample_rate/2, sample_rate)
    {
        setFrequency1(sample_rate);
        setFrequency2(0);
    }

    float processSample(float sample) 
    {
        auto o1 = oscillator_1_.processSample();
        auto o2 = oscillator_2_.processSample();

        float m1a = o1.first * sample;
        float m1b = o1.second * sample;

        float m2a = o2.first * sample;
        float m2b = o2.second * sample;

        float f2a = filter_1_.processSample(sample);
        float f2b = filter_2_.processSample(sample);

        return f2a * o2.first + f2b * o2.second;
    }

    void setSampleRate(unsigned int sample_rate)
    {
        oscillator_1_.setSampleRate(sample_rate);
        oscillator_2_.setSampleRate(sample_rate);
        filter_1_.setSampleRate(sample_rate);
        filter_2_.setSampleRate(sample_rate);
    }

    void setFrequency1(float frequency) 
    {
        oscillator_1_.setFrequency(frequency);
        filter_1_.setCutoffFrequency(frequency);
        filter_2_.setCutoffFrequency(frequency);
    }

    void setFrequency2(float frequency)
    {
        oscillator_2_.setFrequency(frequency);
    }

private:
    QuadratureOscillator oscillator_1_;
    QuadratureOscillator oscillator_2_;
    IIRLowPass filter_1_;
    IIRLowPass filter_2_;
};