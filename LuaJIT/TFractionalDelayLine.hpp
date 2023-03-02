/*
  ==============================================================================
    FractionalDelayBuffer.h
  ==============================================================================
*/

#pragma once

#include <iostream>
#include "SoundAlchemy.hpp"

namespace SoundAlchemy::Delays
{

template<typename T>
class FractionalDelayBuffer
{
public:
    FractionalDelayBuffer();
    ~FractionalDelayBuffer();
    
    void clear();
    void setBufferSize(int size);
    void addSample(T sample);
    T    getSample();
    int  getValidIndex(int index);
    
private:    
    int read_cursor;
    int index;
    int bufferSize;
    T* buffer;
};


template<typename T>
FractionalDelayBuffer<T>::FractionalDelayBuffer()
{    
    read_cursor =0;
    index = 0;
    bufferSize = 0;
    buffer = NULL;
}


template<typename T>
FractionalDelayBuffer<T>::~FractionalDelayBuffer()
{}


template<typename T>
void FractionalDelayBuffer<T>::clear()
{
    memset(buffer, 0, bufferSize*sizeof(T));
}


template<typename T>
void FractionalDelayBuffer<T>::setBufferSize(int size)
{
    if (buffer) {
        delete [] buffer;
    }
    
    bufferSize = size;
    buffer = new T[bufferSize];
    index = bufferSize-1;
    clear();
}


template<typename T>
void FractionalDelayBuffer<T>::addSample(T sample)
{
    index = getValidIndex(index);
    buffer[index] = sample;
    index++;
}


template<typename T>
int FractionalDelayBuffer<T>::getValidIndex(int index) {
    return index % (bufferSize);
}

template<typename T>
T FractionalDelayBuffer<T>::getSample()
{    
    int lower = read_cursor++;
    int upper = lower + 1;
    if (upper == bufferSize) upper = 0;    
    T x1 = buffer[lower];
    T x2 = buffer[upper];
    T frac = x1 - floor(x1);
    return x1 + frac*(x2-x1);
}
}