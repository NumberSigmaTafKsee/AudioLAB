/*
  ==============================================================================
    AllPass.h
    Created: 15 Oct 2014 8:55:30pm
    Author:  Keith Hearne
    
    Based on model that Schroeder proposed in 1962 paper presenting his
    initial reverb designs, that uses a feedback delay line with feedforward
    line.
    
    A Basic All-pass IIR Filter class that sets delay and gain and allows
    access to each, while providing comb filter processing function
    
  ==============================================================================
*/

#pragma once

#include <algorithm>
#include <cmath>
#include <cassert>
#include <cstdlib>
#include <iostream>

#include "DelayLine.hpp"
#include "AllPass.hpp"
#include "Comb2.hpp"

namespace FX::Delays::KHDelays
{
//////////////////////////////////////////////////////////
//  SCHROEDER REVERB
//  see .cpp file for full comments
//////////////////////////////////////////////////////////

//predefined number of comb and allpass filters for array parsing

class Schroeder{

    static const int NUM_COMBS=4;
    static const int NUM_ALLPASSES=2;
public:
    //--------------------------------------------------------------
    //constructor setting initial values for comb delays and gains
    //comb delays must be mutually prime
    //
    //  Comb 1  : 29.7 msec delay
    //  Comb 2  : 37.1 msec delay
    //  Comb 3  : 41.1 msec delay
    //  Comb 4  : 43.7 msec delay
    //  APF 1   : 5.0 msec delay, gain 0.707
    //  APF 2   : 1.7 msec delay, gain 0.5
    //  sr      : defaulted to 44100KHz
    //  rt60    : defaulted to 3.0 seconds on initialisation
    //--------------------------------------------------------------
    Schroeder(const int sr = 44100, const DspFloatType rt60 = 3.0,
              const DspFloatType cDelay1 = 29.7, const DspFloatType cDelay2 = 37.1,
              const DspFloatType cDelay3 = 41.1, const DspFloatType cDelay4 = 43.7,
              const DspFloatType aDelay1 = 5.0, const DspFloatType aDelay2 = 1.7,
              const DspFloatType aGain1 = 0.707, const DspFloatType aGain2 = 0.5);
    ~Schroeder();
    
    //getters
    DspFloatType getDecayFactor();
    DspFloatType getCombDelay(const int id);
    DspFloatType getAllpassDelay(const int id);
    DspFloatType getAllpassGain(const int id);
    //DspFloatType getLowpassCutoff(const int id);
    bool getBypass();
    
    //setters
    void setDecayFactor(const DspFloatType df);
    void setCombDelay(const int id, const DspFloatType sr, const DspFloatType d_ms);
    void setAllpassGain(const int id, const DspFloatType g);
    void setAllpassDelay(const int id, const int sr, const DspFloatType d_ms);
    void setBypass(bool bp);
    void reset();
    
    //business methods
    DspFloatType next(const DspFloatType in);
    
    DspFloatType calcCombGain(const DspFloatType d_ms, const DspFloatType rt60);
    DspFloatType linInterp(DspFloatType x1, DspFloatType x2, DspFloatType y1, DspFloatType y2, DspFloatType x);
    DspFloatType numSamplesFromMSf(const int sr, const DspFloatType d_ms);

    private:
    DspFloatType decayFactor, ALLPASS_GAIN_LIMIT = 0.707f;//to keep the allpasses from exploding
    bool bypass;
    Comb2 *combs[NUM_COMBS];
    Allpass *allpasses[NUM_ALLPASSES];
};


//////////////////////////////////////////////////////////
//  SCHROEDER REVERB
//////////////////////////////////////////////////////////


//helper functions
//------------------------------------------------------------------
//------------------------------------------------------------------
//  calcCombGain : Function to calculate gain from decay/RT60
//
//  RT60    :   value from plugin decay parameter
//  d_ms    :   Delay value of the comb filter
//------------------------------------------------------------------
//------------------------------------------------------------------


inline DspFloatType Schroeder::calcCombGain(const DspFloatType d_ms, const DspFloatType rt60){
    return pow(10.0, ((-3.0 * d_ms) / (rt60 * 1000.0)));
}


//--------------------------------------------------------------
//--------------------------------------------------------------
//constructor setting initial values for comb delays and gains
//comb delays must be mutually prime
//
//  Comb 1  : 29.7 msec delay
//  Comb 2  : 37.1 msec delay
//  Comb 3  : 41.1 msec delay
//  Comb 4  : 43.7 msec delay
//  APF 1   : 5.0 msec delay, gain 0.707
//  APF 2   : 1.7 msec delay, gain 0.5
//  sr      : defaulted to 44100KHz
//  rt60    : defaulted to 3.0 seconds on initialisation
//--------------------------------------------------------------
//--------------------------------------------------------------
Schroeder::Schroeder(const int sr, const DspFloatType rt60,
          const DspFloatType cDelay1, const DspFloatType cDelay2, const DspFloatType cDelay3, const DspFloatType cDelay4,
          const DspFloatType aDelay1, const DspFloatType aDelay2, const DspFloatType aGain1, const DspFloatType aGain2){
    
    decayFactor = rt60;
    DspFloatType d_ms, d_ms_max = 100.0f, gain;
    bypass = false;
    
    //Comb Filter 1 setup
    d_ms = cDelay1;
    gain = calcCombGain(d_ms, decayFactor);
    combs[0] = new Comb2(sr, d_ms, d_ms_max, gain);
    setCombDelay(0,sr,d_ms);
    
    //Comb Filter 2 setup
    d_ms = cDelay2;
    gain = calcCombGain(d_ms, decayFactor);
    combs[1] = new Comb2(sr, d_ms, d_ms_max, gain);
    setCombDelay(1,sr,d_ms);
    
    //Comb Filter 3 setup
    d_ms = cDelay3;
    gain = calcCombGain(d_ms, decayFactor);
    combs[2] = new Comb2(sr, d_ms, d_ms_max, gain);
    setCombDelay(2,sr,d_ms);
    
    //Comb Filter 4 setup
    d_ms = cDelay4;
    gain = calcCombGain(d_ms, decayFactor);
    combs[3] = new Comb2(sr, d_ms, d_ms_max, gain);
    setCombDelay(3,sr,d_ms);

    d_ms_max = 20.0f;
    
    //All-pass filter setup
    allpasses[0] = new Allpass(sr, aDelay1, d_ms_max, aGain1);
    allpasses[1] = new Allpass(sr, aDelay2, d_ms_max, aGain2);


}

//-------------------------------------------------------------------------
// Destructor :
// delete all combs and allpasses
//-------------------------------------------------------------------------
Schroeder::~Schroeder(){
    for(int i = 0; i < NUM_COMBS; i++){
        delete combs[i];
    }
    for(int i = 0; i < NUM_ALLPASSES; i++){
        delete allpasses[i];
    }
}

//getters
//-------------------------------------------------------------------------
// getDecayFactor :
// return the decay factor used for determining RT60 and filter gain
//-------------------------------------------------------------------------
DspFloatType Schroeder::getDecayFactor(){return decayFactor;}

//-------------------------------------------------------------------------
// getCombDelay : comb id
// get the specified delay time in milliseconds of the indexed comb filter
//-------------------------------------------------------------------------
DspFloatType Schroeder::getCombDelay(const int id){return combs[id]->getDelayTimeMS();}

//-------------------------------------------------------------------------
// getAllpassDelay : allpass id
// get the specified delay time in milliseconds of the indexed allpass filter
//-------------------------------------------------------------------------
DspFloatType Schroeder::getAllpassDelay(const int id){return allpasses[id]->getDelayTimeMS();}

//-------------------------------------------------------------------------
// getAllpassGain : comb id
// get the specified gain scalar value of the indexed comb filter
//-------------------------------------------------------------------------
DspFloatType Schroeder::getAllpassGain(const int id){return allpasses[id]->getGain();}

//-------------------------------------------------------------------------
// getBypass : 
// return the status of the boolean for bypassing the plugin processing
//-------------------------------------------------------------------------
bool Schroeder::getBypass(){return bypass;}

//setters
//-------------------------------------------------------------------------
// setDecayFactor : decayfactor value in seconds
// decay time/desired RT60 is passed from UI to this function
// and the required comb filter gain values that will adhere to that RT60
// are calculated based on this factor
//-------------------------------------------------------------------------
void Schroeder::setDecayFactor(const DspFloatType df){
    decayFactor = df;
    
    //cycle through each comb and reset the comb gain value according to
    //the new decayFactor
    for(int i = 0; i < NUM_COMBS; i++){
        combs[i]->setGain(calcCombGain(combs[i]->getDelayTimeMS(), decayFactor));
    }
};

//-------------------------------------------------------------------------
// setCombDelay : id of comb, sample rate, delay time in milliseconds
// sets the gain and the delaytime in milliseconds of the Comb filters
// delay buffer when a value is changed through the UI
//-------------------------------------------------------------------------
void Schroeder::setCombDelay(const int id, const DspFloatType sr, const DspFloatType d_ms){
    combs[id]->setGain(calcCombGain(d_ms, decayFactor));
    combs[id]->setDelayTimeMS(sr, d_ms);
}

//-------------------------------------------------------------------------
// setAllpassGain : id of comb, gain
// sets the gain for the allpass filter, scaling by the GAIN_LIMIT so as
// not to blow the filter up due to the unstable nature of IIR filters
//-------------------------------------------------------------------------
void Schroeder::setAllpassGain(const int id, const DspFloatType g){allpasses[id]->setGain(g * ALLPASS_GAIN_LIMIT);}

//-------------------------------------------------------------------------
// setAllpassDelay : id of comb, sample rate, delay in milliseconds
// sets the delay time in milliseconds of the all pass delay line
//-------------------------------------------------------------------------
void Schroeder::setAllpassDelay(const int id, const int sr, const DspFloatType d_ms){allpasses[id]->setDelayTimeMS(sr, d_ms);}

//-------------------------------------------------------------------------
// setBypass : boolean bypass value
// sets a boolean which determines if processing should be bypassed in the
// worker next function
//-------------------------------------------------------------------------
void Schroeder::setBypass(bool bp){bypass = bp;}

//-------------------------------------------------------------------------
// reset : 
// resets the delay lines in the combs and allpass filters
//-------------------------------------------------------------------------
void Schroeder::reset(){
    for(int i = 0; i < NUM_COMBS; i++){
        combs[i]->resetDelay();
    }
    
    for(int i = 0; i < NUM_ALLPASSES; i++){
        allpasses[i]->resetDelay();
    }
    
}

//------------------------------------------------------------------
//------------------------------------------------------------------
//  next    : Function to process next sample input in
//          : each input sample is passed to each of the comb 
//          : filters in turn (scaling to prevent clipping)
//          : the output value is accumulated/summed
//          : the result is then passed in series through the 
//          : all-pass filters
//
//  in      :   input sample form the audio buffer
//  
//------------------------------------------------------------------
//------------------------------------------------------------------
DspFloatType Schroeder::next(const DspFloatType in){
    
    // if bypass is enabled on plugin then just pass back
    // the unprocessed input sample
    if(bypass)
        return in;
        
    DspFloatType out = 0.0f;
    
    //step through each of the 4 comb filters and pass a scaled input
    for(int i = 0; i < NUM_COMBS; i++){
        out += combs[i]->next(in * 0.50f); //scale down to avoid clipping
    }
    
    DspFloatType passOut = 0.0f;
    DspFloatType passOut2 = 0.0f;
    
    passOut = allpasses[0]->next(out);          //1 stage all-pass filtering
    passOut2 = allpasses[1]->next(passOut);     //2 stage all-pass filtering
  
    return passOut2;
    
}


}
