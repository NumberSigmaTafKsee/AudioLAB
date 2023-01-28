#pragma once 

#include <cstring>
#include <cstddef>

class FractionalDelayBuffer
{
public:
    FractionalDelayBuffer();
    ~FractionalDelayBuffer();
    
    void clear();
    void setBufferSize(int size);
    void addSample(float sample);
    float getSample(float sampleIndex);
    int getValidIndex(int index);
    float linearInterpolate(float* buffer, int bufferSize, int bufferPosition);
    
private:
    int index;
    int bufferSize;
    float* buffer;
};