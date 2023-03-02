#pragma once

#ifndef FRACDELAYLINE_H_INCLUDED
#define FRACDELAYLINE_H_INCLUDED

class FracDelayLine
{
public:
    FracDelayLine();
    ~FracDelayLine();

    void clearBuffer();
    void setSize(int size);
    int currentIndex(int idx);
    void setSample(float sample);
    float getLinInterpSample(float sampleIdx);

private:
    int idx;
    int circularBufferSize;
    float* circularBuffer;
};
