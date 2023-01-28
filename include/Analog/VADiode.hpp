#pragma once

// https://github.com/titas2001/Diode_Clipper
// https://github.com/titas2001/Diode_Ring_Modulator

#include "FX/LiquidNeuron.hpp"
#include "FX/Diode.hpp"

namespace Analog::Distortion::Diode
{
    // i do not know wtf this was
    struct Dioder
    {
        Liquid::LiquidPole vt,eta,is;
        DspFloatType Vt,Eta,Is,sampleRate=44100.0;
        int ctr=0;
        Dioder(DspFloatType sr=44100.0)  :vt(1/0.02,1000.0),eta(1/0.02,1000),is(1/0.02,1000)
        {        
            Random noise;
            sampleRate = sr;
            Vt  = noise.frand()*0.01;
            Eta = noise.frand()*0.1;
            Is  = noise.frand()*0.001;
            ctr=1000;
        }    
        DspFloatType Tick(DspFloatType In, DspFloatType V=1, DspFloatType E=1, DspFloatType I=1)
        {
            ctr++;
            DspFloatType v = 0.0253+vt.Tick(Vt);        
            DspFloatType e = 1.68+eta.Tick(Eta);
            DspFloatType i = .015+is.Tick(Is);
            
            if(ctr >= sampleRate/1000) {
                Random noise;
                Vt  = noise.frand()*0.01;
                Eta = noise.frand()*0.1;
                Is  = noise.frand()*0.001;
                ctr = 0;
            }                    
            return FX::Distortion::Diode::Diode(In,v,e,i);
        }
    };    
}