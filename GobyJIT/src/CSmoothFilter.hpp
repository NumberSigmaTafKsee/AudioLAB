#pragma once

#include <cmath>

namespace FX::Filters::Smooth
{
    // is just the one pole again
    class CSmoothFilter : public FilterProcessor
     {
    public:
                
        CSmoothFilter(double Fs=44100.0f, double Fc=0.1f) : FilterProcessor() {
            a0 = 1.0; b1 = 0.0; z1 = 0.0;
            fs = Fs;            
            setCutoff(Fc);        
        };
        
        void setCutoff(double Fc)
        {            
            f = Fc/fs;                    
            b1 = exp(-2.0 * M_PI *f);
            a0 = 1.0 - b1;
        }
        void setTarget(double t)
        {
            target = t;
        }
        bool isAtTarget() {
            return fabs(z1 - target) < 1e-3;
        }	    	
        double process()
        {
            return z1 = target * a0 + z1 * b1;
        }        
        double process(double in) {            
            return z1 = in * a0 + z1 * b1;
        }
        double Tick(double I, double A=1, double X=0, double Y=0) {
            return process(I);
        }
    protected:    
        double a0, b1, z1, f, fs;
        double target = 1.0;
    };

    class CParamSmooth : public FilterProcessor
    {
    public:
        CParamSmooth() : FilterProcessor() { a = 0.99f; b = 1.f - a; z = 0; };
        ~CParamSmooth();
        inline double Process(double in) { z = (in * b) + (z * a); return z; }

        double Tick(double I, double A=1, double X=0, double Y=0) {
            return Process(I);
        }
    private:
        double a, b, z;
    };

    class CSmoothTime : public FilterProcessor
    {
    public:

        CSmoothTime(double smoothingTimeInMs, double samplingRate)
        : FilterProcessor()
        {
            const double c_twoPi = 6.283185307179586476925286766559f;

            a = exp(-c_twoPi / (smoothingTimeInMs * 0.001f * samplingRate));
            b = 1.0f - a;
            z = 0.0f;
        }

        ~CSmoothTime()
        {

        }

        inline double process(double in)
        {
            z = (in * b) + (z * a);
            return z;
        }

        double Tick(double I, double A=1, double X=0, double Y=0) {
            return process(I);
        }
    private:
        double a;
        double b;
        double z;
    };
}