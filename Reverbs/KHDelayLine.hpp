/*
  ==============================================================================

    DelayLine.h
    New version of a basic DelayLine which uses fractional delay
    Operations for reading and writing from the delay have been
    separated out into new functions so that the can be read from
    and written to in the Comb and All pass filters seperately.
    
    Created: 15 Oct 2014 5:15:53pm
    Author:  Keith Hearne

  ==============================================================================
*/

#pragma once

#include <algorithm>
#include <cmath>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <cstring>
#include "Helper.h"

//////////////////////////////////////////////////////////
//  DELAY LINE CLASS (for use with Moorer Reverb)
//  see .cpp file for comments
//////////////////////////////////////////////////////////

namespace FX::Delays
{
    class DelayLine : public FunctionProcessor
    {

    public:
        //constructor / destructor
        DelayLine(const int sr = 44100, const DspFloatType d_ms = 0.0f, const DspFloatType d_max = 0.0f);
        ~DelayLine();
        
        //getters
        DspFloatType getDelayTimeMS();
        DspFloatType getMaxDelayTimeMS();
        
        //setters
        void setDelayTimeMS(const int sr, const DspFloatType d_ms);
        void setDelay(DspFloatType d_ms);
        
        //business functions
        DspFloatType readDelay();
        void writeDelay(DspFloatType in);
        DspFloatType next(const DspFloatType in);
        void resetBuffer();
        void resetDelay();

        DspFloatType numSamplesFromMSf(const int sr, const DspFloatType d_ms);
        DspFloatType linInterp(DspFloatType x1, DspFloatType x2, DspFloatType y1, DspFloatType y2, DspFloatType x);
    
        DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1)
        {
            return A*next(I);
        }
    private:
        int writePos, readPosA, MAX_DELAY_SAMPLES;
        DspFloatType delay_ms, delay_samples, fraction, MAX_DELAY_MS;
        DspFloatType *buffer;
        
    };

    //helper functions
    //-------------------------------------------------------------------------
    // numSamplesFromMSf :
    // Determine the number of samples from the number of milliseconds delay
    // passed to function
    //-------------------------------------------------------------------------
        
    inline DspFloatType DelayLine::numSamplesFromMSf(const int sr, const DspFloatType d_ms){
        return sr * d_ms * .001;
    }

    //----------------------------------------------------------------------
    //  Linear Interpolation Function
    //  x1  :  weighting 1
    //  x2  :  weighting 2
    //  y1  :  output y(n) at read index in buffer
    //  y2  :  output y(n-1) at read index minus 1 in buffer
    //  x   :  fractional value between samples
    //----------------------------------------------------------------------
    

    inline DspFloatType DelayLine::linInterp(DspFloatType x1, DspFloatType x2, DspFloatType y1, DspFloatType y2, DspFloatType x){
        DspFloatType denom = x2 - x1;
        if(denom == 0)
            return y1; // should not ever happen

        // calculate decimal position of x
        DspFloatType dx = (x - x1)/(x2 - x1);

        // use weighted sum method of interpolating
        DspFloatType result = dx*y2 + (1-dx)*y1;

        return result; 
    }
    //////////////////////////////////////////////////////////
    //  DELAY LINE CLASS (for use with Schroeder Reverb)
    //////////////////////////////////////////////////////////

    //-------------------------------------------------------------------------
    // Constructor :
    // Predefined sample rate = 44100, delay = 0, max delay
    // set read and write index for delay line
    // set max delay set from decay factor in reverb
    // create a max delay buffer
    //-------------------------------------------------------------------------
    DelayLine::DelayLine(const int sr, const DspFloatType d_ms, const DspFloatType d_max)
    : FunctionProcessor()
    {
        //assert(d_ms <= d_ms_max);//check bound on delay time
        
        //initialisation
        buffer = NULL;
        readPosA = writePos = 0;
        delay_samples = 0.0f;
        delay_ms = d_ms;
        
        MAX_DELAY_SAMPLES = ceil(numSamplesFromMSf(sr, d_max));
        MAX_DELAY_MS = MAX_DELAY_SAMPLES * 1000.0f / sr; //make sure DspFloatType version is set with to integer-rounded buffer size
        
        DspFloatType delay_samplesf = numSamplesFromMSf(sr, d_ms);
        delay_samples = floor(delay_samplesf);
        
        //storing fractional delay time, will be interpolated
        fraction = delay_samplesf - delay_samples;
        
        buffer = new DspFloatType[MAX_DELAY_SAMPLES];
        memset(buffer, 0, MAX_DELAY_SAMPLES*sizeof(DspFloatType));

    }

    //-------------------------------------------------------------------------
    // Destructor :
    // delete delay buffer
    //-------------------------------------------------------------------------
    DelayLine::~DelayLine(){
        delete[] buffer;
    }

    //getters
    //-------------------------------------------------------------------------
    // getDelayTimeMS :
    // return the delay in milliseconds
    //-------------------------------------------------------------------------
    DspFloatType DelayLine::getDelayTimeMS(){return delay_ms;}

    //--------------------------------------------------------------------------------
    //  getMaxDelayTimeMS :
    //  return the Max delay time in Miliiseconds
    //--------------------------------------------------------------------------------
    DspFloatType DelayLine::getMaxDelayTimeMS(){return MAX_DELAY_MS;}

    //setters
    //--------------------------------------------------------------------------------
    //  Setter function that determines read position index
    //  read position is determined by subtracting the number of samples to delay
    //  from the write position index
    //
    //  readIndex = writeIndex - number of sample delay
    //
    //--------------------------------------------------------------------------------
    void DelayLine::setDelayTimeMS(const int sr, const DspFloatType d_ms){
        assert(d_ms <= MAX_DELAY_MS);//check bound on delay time
        
        DspFloatType delay_samplesf = numSamplesFromMSf(sr, d_ms);
        delay_samples = floor(delay_samplesf);
        
        //storing fractional delay time, output will be interpolate
        fraction = delay_samplesf - delay_samples;
        
        readPosA = writePos - delay_samples; 
        
        //wrap index if bounds exceeded
        if(readPosA < 0)
            readPosA += MAX_DELAY_SAMPLES;
            
            
    }

    //-------------------------------------------------------------------------
    // setDelay :
    // set the delay in milliseconds by the delay value passed to function
    //-------------------------------------------------------------------------
    void DelayLine::setDelay(DspFloatType d){
        //receiving the delay value through here in milliseconds 0 to 2000
        delay_ms = d;
        setDelayTimeMS(44100,delay_ms);
    };

    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    //  readDelay :
    //
    //  The following function reads from delay line at read index
    //  It uses terms similar to the differential equation for delay
    //
    //      x(n)                : the input at sample n
    //      y(n)                : the output at sample n after delay processing
    //      buffer              : the circular buffer used
    //      readPosA            : the read index of the delay buffer
    //      writePos            : the write index of the delay buffer
    //      MAX_DELAY_SAMPLES   : Max size of delay buffer
    //  
    //      y(n) = x(n) + x(n - D)              'delay with no feedback
    //
    //      MAX_DELAY_SAMPLES = sr * d_ms * .001;
    //---------------------------------------------------------------------------------
    //---------------------------------------------------------------------------------
    DspFloatType DelayLine::readDelay(){
        
        //output of delay at readPos
        DspFloatType yn = buffer[readPosA];
        
        //read location at n-1, one behind yn
        int readPos_minus1 = readPosA - 1;
        if(readPos_minus1 < 0)
            readPos_minus1 = MAX_DELAY_SAMPLES - 1;   //MAX_DELAY_SAMPLES - 1 is the last location of the buffer
            
        //get y(n-1)
        DspFloatType yn_minus1 = buffer[readPos_minus1];
        
        //perform linear interpolation of : (0,yn) and (1,yn_minus1) by the ammount of fractional delay(fraction)
        return linInterp(0, 1, yn, yn_minus1, fraction);
    }

    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    //  writeDelay :
    //
    //  The following function writes into the delay buffer at write index
    //---------------------------------------------------------------------------------
    //---------------------------------------------------------------------------------
    void DelayLine::writeDelay(DspFloatType in){
        
        //write the input into buffer at write position
        buffer[writePos] = in;
        
        //increment the write index, wrap index if bounds exceeded
        writePos++;
        if(writePos >= MAX_DELAY_SAMPLES)
            writePos = 0;
            
        readPosA++;
        if(readPosA >= MAX_DELAY_SAMPLES)
            readPosA = 0;
    }

    //business function
    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    //  next    :
    //  The following function is the workhorse of the delay line, but the filters
    //  in the reverb will call the individual read/write function so as to have
    //  more control over manipulating the inputs and outputs.
    //
    //      x(n)                : the input at sample n
    //      y(n)                : the output at sample n after delay processing
    //  
    //      y(n) = x(n) + x(n - D)              'delay with no feedback
    //
    //  read the interpolated value from the delay line and write input back to buffer
    //
    //---------------------------------------------------------------------------------
    //---------------------------------------------------------------------------------
    DspFloatType DelayLine::next(const DspFloatType in){
            
        //input
        DspFloatType xn = in;

        //if delay value is zero just pass input out
        DspFloatType yn;
        if(delay_samples == 0)
            yn = xn;
        else
            yn = readDelay(); 

        writeDelay(xn);
        
        return yn;
            
    }

    //-------------------------------------------------------------------------
    // resetBuffer :
    // delete contents of buffer and instantiate a new buffer
    //-------------------------------------------------------------------------
    void DelayLine::resetBuffer(){
        
        if(buffer)
            delete [] buffer;
            
        buffer = new DspFloatType[MAX_DELAY_SAMPLES];
        
        resetDelay();
        setDelayTimeMS(44100,delay_ms);
        
        return;
        
    }
    
    //-------------------------------------------------------------------------
    // resetDelay :
    // reset the delay buffer by filling with 0's, reset all indexes and bypass
    //------------------------------------------------------------------------- 
    void DelayLine::resetDelay(){
    
        if(buffer)
            memset(buffer, 0, MAX_DELAY_SAMPLES*sizeof(DspFloatType));
        
        readPosA = writePos = 0;
        
        return;
        
    }

}