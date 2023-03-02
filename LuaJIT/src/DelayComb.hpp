/*
  ==============================================================================

    Comb.h
    Created: 15 Oct 2014 8:55:11pm
    Author:  Keith Hearne
    
    A LowPass Comb Filter class that sets delay and gain while passing
    its output through a low pass filter to mimic the effect air has
    on the high frequency content lossand as refelctions in the feedback
    loop are processed 

  ==============================================================================
*/

#pragma once

#include <cstring>
#include "DelayLine.hpp"
#include "LowPass.hpp"

//////////////////////////////////////////////////////////
//  COMB FILTER CLASS
//  see .cpp file for comments
//////////////////////////////////////////////////////////

namespace Delay
{
    class Comb : public FilterProcessor
    {
        
    public:
        //constructor / destructor
        Comb(const int sr, const float d_ms, const float d_ms_max, const float g, const float cutOff);
        ~Comb();
        
        //getters
        float getGain();
        float getDelayTimeMS();
        
        //setters
        void setGain(const float g);
        void setDelayTimeMS(const float sr, const float d_ms);
        void setLPF(const float cutoff);
        
        //business methods
        float next(const float in);
        
        double Tick(double I, double A=1, double X=1, double Y=1) {
            return next(I);
        }
    private:
        float gain;
        DelayLine *delay; 
        Lowpass *lpFilter;
        
    };



    //////////////////////////////////////////////////////////
    //  BASIC COMB FILTER CLASS
    //////////////////////////////////////////////////////////

    //-------------------------------------------------------------------------
    // Constructor :
    // Predefined sample rate = 44100, delay time, max delay time and gain
    // and lowpass cutoff frequency. creates a new DelayLine and sets 
    // the input gain and initiates the low pass filter
    //-------------------------------------------------------------------------
    Comb::Comb(const int sr, const float d_ms, const float d_ms_max, const float g, const float lp_freq)
    : FilterProcessor()
    {
        gain = g;
        delay = new DelayLine(sr, d_ms, d_ms_max);
        lpFilter = new Lowpass(44100, lp_freq);
    }

    //-------------------------------------------------------------------------
    // Destructor :
    // deletes the delay and filter
    //-------------------------------------------------------------------------
    Comb::~Comb(){
        delete delay;
        delete lpFilter;
    }

    //getters
    //-------------------------------------------------------------------------
    // getGain  :
    // return the gain scalar value
    //-------------------------------------------------------------------------
    float Comb::getGain(){return gain;}

    //-------------------------------------------------------------------------
    // getDelayTimeMS  :
    // return the delay time in milliseconds
    //-------------------------------------------------------------------------
    float Comb::getDelayTimeMS(){return delay->getDelayTimeMS();}

    //setters
    //-------------------------------------------------------------------------
    // setGain  :
    // set the scalar value for the gain of the filter
    //-------------------------------------------------------------------------
    void Comb::setGain(const float g){gain = g;}

    //--------------------------------------------------------------------------------
    //  setDelayTimeMS
    //  Setter function sets delay time and from milliseconds
    //  and passes it to the delayline which converts to discrete time samples
    //
    //--------------------------------------------------------------------------------
    void Comb::setDelayTimeMS(const float sr, const float d_ms){return delay->setDelayTimeMS(sr, d_ms);}

    //--------------------------------------------------------------------------------
    //  setLPF
    //  Setter function sets the cutoff frequency of the low pass filter in the
    //  feedback loop
    //--------------------------------------------------------------------------------
    void Comb::setLPF(const float cutoff_freq){ return lpFilter->setCutoff(44100, cutoff_freq);};

    //business methods
    //------------------------------------------------------------------
    //------------------------------------------------------------------
    //  next    : Function to process next sample input in
    //          : The comb filter process involves reading the delay 
    //          : at the current readindex position (delay->readDelay)
    //          : scaling this delay value by the combs gain (gain*dL)
    //          : passing this value through the low pass filter
    //          : and writing this value back into the delay line
    //          : at the appropriate write position 
    //          : (delay->writeDelay(dLW))
    //
    //  in      :   input sample form the audio buffer
    //  
    //------------------------------------------------------------------
    //------------------------------------------------------------------
    float Comb::next(const float in){
        //read delay
        float dL = delay->readDelay();
        
        //attenuate with gain
        float dlAttn = dL * gain;
        //pass through low pass filter
        float lpOut = lpFilter->next(dlAttn);
        
        //combine output for feedback loop
        float dLW = in + lpOut;
        //write feedback loop back to delay buffer
        delay->writeDelay(dLW);
        
        //return the initially read delay
        return dL;
        
    }

}