/*
The following code was taken from the following post:
https://www.musicdsp.org/en/latest/Filters/243-butterworth-optimized-c-class.html

Other than an extra function and some renaming, the code is nearly identical
*/

#pragma once

#include <JuceHeader.h>

class ButterworthFilter
{
public:
    ButterworthFilter(void);
    ~ButterworthFilter(void);
    void SetSampleRate(float fs);
    void Set(float cutoff, float q);
    float Run(float input);
    void processBlock(const AudioSourceChannelInfo& bufferToFill);

private:
    float t0, t1, t2, t3;
    float coef0, coef1, coef2, coef3;
    float history1, history2, history3, history4;
    float gain;
    float min_cutoff, max_cutoff;
};
