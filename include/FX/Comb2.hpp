#pragma once

namespace FX::Delays
{
    //////////////////////////////////////////////////////////
    //  COMB FILTER CLASS
    //////////////////////////////////////////////////////////

    class Comb2{
        
    public:
        //constructor / destructor
        Comb2(const int sr, const DspFloatType d_ms, const DspFloatType d_ms_max, const DspFloatType g);
        ~Comb2();
        
        //getters
        DspFloatType getGain();
        DspFloatType getDelayTimeMS();
        
        //setters
        void setGain(const DspFloatType g);
        void setDelayTimeMS(const DspFloatType sr, const DspFloatType d_ms);
        
        //business methods
        DspFloatType next(const DspFloatType in);
        
        void resetDelay() { delay->resetDelay(); }

    private:
        DspFloatType gain;
        DelayLine *delay;  
        Lowpass *lpFilter;
        
    };

    //////////////////////////////////////////////////////////
    //  BASIC COMB FILTER CLASS
    //////////////////////////////////////////////////////////

    //constructor / destructor
    Comb2::Comb2(const int sr, const DspFloatType d_ms, const DspFloatType d_ms_max, const DspFloatType g){
        gain = g;
        delay = new DelayLine(sr, d_ms, d_ms_max);
        lpFilter = new Lowpass(44100, 10000.0f);
    }

    Comb2::~Comb2(){
        delete delay;
        delete lpFilter;
    }

    //getters
    DspFloatType Comb2::getGain(){return gain;}
    DspFloatType Comb2::getDelayTimeMS(){return delay->getDelayTimeMS();}

    //setters
    void Comb2::setGain(const DspFloatType g){gain = g;}
    void Comb2::setDelayTimeMS(const DspFloatType sr, const DspFloatType d_ms){return delay->setDelayTimeMS(sr, d_ms);}

    //business methods
    DspFloatType Comb2::next(const DspFloatType in){
        //KH***return delay->next(in, gain);
        DspFloatType dL = delay->readDelay();
        
        //DspFloatType lpOut = dL * gain;
        //DspFloatType lpRetVal = lpFilter->next(lpOut);
        DspFloatType lpRetVal = lpFilter->next(dL);
        
        //DspFloatType dLW = in + lpRetVal;
        DspFloatType dLW = in + lpRetVal*gain;
        delay->writeDelay(dLW);
        return dL;
        
    }

}