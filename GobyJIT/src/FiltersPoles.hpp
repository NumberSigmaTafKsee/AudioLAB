#pragma once

namespace Filters::Poles
{

    ////////////////////////////////////////////////////////////////////////////////////////////
    // One Pole/One Zero
    // Mainly used for smoothing or blocking DC or Nyquuist
    ////////////////////////////////////////////////////////////////////////////////////////////
    struct OnePole
    {
        double b0,a1;
        double y1;
        double x,y;
        double fc,fs;

        enum Type
        {
            Lowpass,
            Highpass
        }
        filterType = Lowpass;

        OnePole(double f, double sr)
        {
            y1 = 0.0f;        
            fc = f/sr;
            fs = sr;
            a1 = -std::exp(-2*M_PI*fc);
            b0 = 1.0 + a1;
            
        }
        void setCutoff(double f) {
            fc = f/fs;
            if(filterType == Lowpass)
                setLowpass(f);
            else
                setHighpass(f);

        }   
        inline void setLowpass(double Fc) {
            a1 = std::exp(-2*M_PI*fc);
            b0 = 1.0 - a1;
        }

        inline void setHighpass(double Fc) {
            a1 = exp(-2.0 * M_PI * (0.5 - Fc));
            b0 = 1.0 - a1;
        }

        double Tick(double in) {
            x = in;
            y = b0*x - a1*y1;
            y1= y;
            return y;

        }
    };

    ////////////////////////////////////////////////////////////////////////////////////////////
    // slightly resonant low/high pass
    ////////////////////////////////////////////////////////////////////////////////////////////
    struct TwoPoles
    {
        double b0,a1,a2;
        double y1,y2;
        double x,y;
        double fc,fs;
        float  R;
        enum Type
        {
            Lowpass,
            Highpass,        
        }
        filterType = Lowpass;

        TwoPoles(double sr,double f)
        {
            y1 = 0.0f;   
            y2 = 0.0f;     
            fc = f/sr;
            fs = sr;        
            R  = 0.0f;
            setCutoff(f);        
        }
        // this makes the filter
        void setCutoff(double f) {
            fc = f/fs;
            if(filterType == Lowpass)
                setLowpass();
            else if(filterType == Highpass)
                setHighpass();
        }           
        void setQ(double Q) {
            if(R == 0.0f) R = 0.005f;            
            R = Q;
        }
        inline void setLowpass() {
            /* experiment
            double K = tan(M_PI*fc);
            double bottom = (K*K*q+K+q);
            a1 = (2*q*(K*K-1))/bottom;
            a2 = (K*K*q-K+q)/bottom;
            b0 = (K*K*q)/bottom;
            */
            b0 = (1+a1)*(1+a1);
            a1 = -2*R*cos(2*M_PI*fc);
            a2 = R*R;
            
        }
        inline void setHighpass() {
            /* experiment
            double K = tan(M_PI*fc);
            double bottom = (K*K*q+K+q);
            a1 = (2*q*(K*K-1))/bottom;
            a2 = (K*K*q-K+q)/bottom;
            b0 = (K*K*q)/bottom;
            */
            b0 = (1-a1)*(1-a1);
            a1 = -2*R*cos(2*M_PI*(0.5-fc));
            a2 = R*R;
            
        }

        

        double Tick(double in) {                
            x = in;
            y = b0*x - a1*y1 - a2*y2;
            y2 = y1;
            y1= y;                
            return y;
        }
    };
}
