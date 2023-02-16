/*
  ==============================================================================

    LowPass.h
    Created: 27 Oct 2014 9:10:20pm
    Author:  Keith Hearne

    A basic Single One Pole Low Pass Filter
    Formula for coefficients taken from R. Boulanger, 2011, p.486
    Filter cutoff frequecy can be varied on the input 
    and changed via the UI of the plugin.

  ==============================================================================
*/

#pragma once

#include <cmath>

//////////////////////////////////////////////////////////
//  BASIC LOWPASS FILTER CLASS
//  see .cpp file for comments
//////////////////////////////////////////////////////////
namespace FX::Delays
{
    class Lowpass
    {

    public:
        //constructor
        Lowpass(const int sr, const float cf_hz);
        
        //getters
        float getCutoff();
        
        //setters
        void setCutoff(const int sr, const float cf_hz);
        
        //business methods
        float next(const float in);
        
    private:
        float cutoff, coef, prev;
        
    };


    //////////////////////////////////////////////////
    //  BASIC LOWPASS FILTER CLASS
    //////////////////////////////////////////////////////////

    //-------------------------------------------------------------------------
    // Constructor :
    // Predefined sample rate = 44100, default cutoff frequency passed here
    // from these values the coefficients for a and b are calculated 
    // a = 1 + b
    // b = sqrt(squared[2 - cos(2*PI*freq/sr)] - 1) -2 + cos(2*PI*freq/sr)
    //-------------------------------------------------------------------------
    Lowpass::Lowpass(const int sr, const float cf_hz){
        cutoff = coef = prev = 0;
        cutoff = cf_hz;
        float costh = 2.0 - cos(2.0 * M_PI * cutoff / sr);
        coef = sqrt(costh * costh - 1.0) - costh;
    }

    //getters
    //-------------------------------------------------------------------------
    // getCutoff :
    // return the value of the cutoff frequency
    //-------------------------------------------------------------------------
    float Lowpass::getCutoff(){return cutoff;}

    //setters
    //-------------------------------------------------------------------------
    // setCutoff :
    // function to adjust and set the filter coefficients from the cutoff
    // frequency parameter
    //-------------------------------------------------------------------------
    void Lowpass::setCutoff(const int sr, const float cf_hz){
        cutoff = cf_hz;
        float costh = 2.0 - cos(2.0 * M_PI * cutoff / sr);
        coef = sqrt(costh * costh - 1.0) - costh;
    }

    //business methods
    //-------------------------------------------------------------------------
    // next :
    // the process function which takes the input discrete time sample and 
    // applies the coefficient one sample delay formula to it returning the
    // previous value (input * (1 + b) - (delayed sample * b)
    //-------------------------------------------------------------------------
    float Lowpass::next(const float in){
        prev = in * (1 + coef) - (prev * coef);
        return prev;
    }
}