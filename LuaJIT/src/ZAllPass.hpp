
/*
  ==============================================================================

    AllPass.h
    Created: 15 Oct 2014 8:55:30pm
    Author:  Keith Hearne
    
    Based on model that Schroeder proposed in 1962 paper presenting his
    initial reverb designs, that uses a feedback delay line with feedforward
    line.


  ==============================================================================
*/


#pragma once

#include <cstring>
#include <iostream>
#include "DelayLine.hpp"

namespace Delay
{
    //////////////////////////////////////////////////////////
    //    ALLPASS FILTER CLASS
    //  see .cpp file for comments
    //////////////////////////////////////////////////////////

    class Allpass : public FunctionProcessor
    {

    public:
        //constructor
        Allpass(const int sr, const double d_ms, const double d_ms_max, const double g);
        ~Allpass();
        
        //getters
        double getGain();
        double getDelayTimeMS();
        
        //setters
        void setGain(const double g);
        void setDelayTimeMS(const double sr, const double d_ms);
        
        //business methods
        double next(const double in);
        
        double Tick(double in, double A=1, double X=1, double Y=1) {
            double G = getGain();
            double ms = getDelayTimeMS();
            setGain(G * X);
            setDelayTimeMS(ms * Y);
            double out = next(in);
            setGain(G);
            setDelayTimeMS(sampleRate,ms);
            return A*out;
        }
    private:
        double gain;
        DelayLine *delay;
        
    };




    //////////////////////////////////////////////////////////
    //  ALLPASS FILTER CLASS
    //////////////////////////////////////////////////////////

    //constructor
    //-------------------------------------------------------------------------
    // Constructor :
    // Predefined sample rate = 44100, delay time, max delay time and gain
    // creates a new DelayLine and sets the input gain
    //-------------------------------------------------------------------------
    Allpass::Allpass(const int sr, const double d_ms, const double d_ms_max, const double g)
    : FilterObject()
    {
        gain = g;
        delay = new DelayLine(sr, d_ms, d_ms_max);
    }

    //-------------------------------------------------------------------------
    // Destructor :
    // deletes the delay
    //-------------------------------------------------------------------------
    Allpass::~Allpass(){
        delete delay;
    }

    //getters
    //-------------------------------------------------------------------------
    // getGain  :
    // return the gain scalar value
    //-------------------------------------------------------------------------
    double Allpass::getGain(){return gain;}

    //-------------------------------------------------------------------------
    // getDelayTimeMS  :
    // return the delay time in milliseconds
    //-------------------------------------------------------------------------
    double Allpass::getDelayTimeMS(){return delay->getDelayTimeMS();}

    //setters
    //-------------------------------------------------------------------------
    // setGain  :
    // set the scalar value for the gain of the filter
    //-------------------------------------------------------------------------
    void Allpass::setGain(const double g){gain = g;}

    //--------------------------------------------------------------------------------
    //  setDelayTimeMS
    //  Setter function sets delay time and from milliseconds
    //  and passes it to the delayline which converts to discrete time samples
    //
    //--------------------------------------------------------------------------------
    void Allpass::setDelayTimeMS(const double sr, const double d_ms){return delay->setDelayTimeMS(sr, d_ms);}

    //------------------------------------------------------------------
    //------------------------------------------------------------------
    //  next    : Function to process next sample input in
    //          : The all-pass filter process involves reading the delay 
    //          : at the current readindex position (delay->readDelay)
    //          : scaling this delay value by the combs gain (gain*dL)
    //          : and writing this value back into the delay line
    //          : at the appropriate write position 
    //          : (delay->writeDelay(dLW))
    //
    //  in      :   input sample form the audio buffer
    //  
    //------------------------------------------------------------------
    //------------------------------------------------------------------
    double Allpass::next(const double in){
        //read delay value from buffer
        double dL = delay->readDelay();
        // for the filter write or feedback input to buffer
        double fW = in + (gain*dL);
        //form the output of all-pass  (delay + in*-gain)
        double out = -gain*fW + dL;  
        // write feedback loop input back to delay buffer
        delay->writeDelay(fW);
        return out;
    }
}
