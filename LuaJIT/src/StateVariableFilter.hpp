#pragma once 
#include "Undenormal.hpp"
#include <cmath>

namespace Filters::SVF
{
    struct StateVariableFilter
    {
        /*
        cutoff = cutoff freq in Hz
        fs = sampling frequency //(e.g. 44100Hz)
        f = 2 sin (pi * cutoff / fs) //[approximately]
        q = resonance/bandwidth [0 < q <= 1]  most res: q=1, less: q=0
        low = lowpass output
        high = highpass output
        band = bandpass output
        notch = notch output
        */
        float cutoff,scale,fs,low,high,band,notch;
            
        StateVariableFilter(float Fc, float Q, float Fs) {
            scale = Q;
            cutoff= Fc;
            fs    = Fs;
            low=high=band=notch=0;
        }
        void setCutoff(float F) { cutoff = F; }
        void setResonance(float R) { scale = 1.25*(1.0-R); }
        float Tick(float I, float A = 1, float X=1, float Y=1)
        {
            Undenormal denormal;
            float f     = 2 * std::sin(2 * M_PI * cutoff/fs);        
            //--beginloop
            //I = tanhify(I);
            float c = cutoff;
            float q = Q;
            setCutoff(c * fabs(X));
            setResonance(q * fabs(Y));
            low = low + f * band;
            high = scale * I - low - scale*band;
            band = f * high + band;
            notch = high + low;
            setCutoff(c);
            setResonance(q);
            return A*low;
        }
    };

    struct StateVariableLowPassFilter : public FilterProcessor
    {
        StateVariableFilter * filter;

        StateVariableHighPassFilter(float Fc, float Q, float Fs) : FilterProcessor() {
            filter = new StateVariableFilter(Fc,Q,Fs);
        }
        ~StateVariableLowPassFilter() {
            if(filter) delete filter;
        }
        float Tick(float I, float A = 1, float X=1, float Y=1) {
            double o = filter->Tick(I,A,X,Y);            
            return o;
        }
    };

    struct StateVariableHighPassFilter : public FilterProcessor
    {
        StateVariableFilter * filter;

        StateVariableHighPassFilter(float Fc, float Q, float Fs) : FilterProcessor() {
            filter = new StateVariableFilter(Fc,Q,Fs);
        }
        ~StateVariableHighPassFilter() {
            if(filter) delete filter;
        }
        float Tick(float I, float A = 1, float X=1, float Y=1) {
            double o = filter->Tick(I,A,X,Y);
            o = filter->high;
            return o;
        }
    };

    struct StateVariableBandPassFilter : public FilterProcessor
    {
        StateVariableFilter * filter;

        StateVariableHighPassFilter(float Fc, float Q, float Fs) : FilterProcessor() {
            filter = new StateVariableFilter(Fc,Q,Fs);
        }
        ~StateVariableHighPassFilter() {
            if(filter) delete filter;
        }
        float Tick(float I, float A = 1, float X=1, float Y=1) {
            double o = filter->Tick(I,A,X,Y);
            o = filter->band;
            return o;
        }
    };

    struct StateVariableNotchFilter : public FilterProcessor
    {
        StateVariableFilter * filter;

        StateVariableNotchFilter(float Fc, float Q, float Fs) : FilterProcessor() {
            filter = new StateVariableFilter(Fc,Q,Fs);
        }
        ~StateVariableNotchFilter() {
            if(filter) delete filter;
        }
        float Tick(float I, float A = 1, float X=1, float Y=1) {
            double o = filter->Tick(I,A,X,Y);
            o = filter->notch;
            return o;
        }
    };
}
