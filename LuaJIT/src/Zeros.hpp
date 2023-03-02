#pragma once

namespace Filters
{
    ////////////////////////////////////////////////////////////////////////////////////////////
    // a kind of low/high pass effect
    // It can be calculaated with windowed sinc as it is a FIR filter with only 1 coefficients
    ////////////////////////////////////////////////////////////////////////////////////////////
    struct OneZero
    {
        double b0,b1;
        double x1;
        double x,y;
        double fc,fs;
        double hp;

        enum Type
        {
            Lowpass,
            Highpass
        }
        filterType = Lowpass;

        OneZero(double f, double sr)
        {
            x1 = 0.0f;        
            fc = f/sr;
            fs = sr;
            setLowpass(f);
        }
        void setCutoff(double f) {
            fc      = f/fs;
            double wc= 2*M_PI*fc;
            double K = tan(wc/2);        
            b0 = K/(K+1);
            b1 = K/(K+1);

        }
        // low pass
        inline void setLowpass(double Fc) {
            double wc      = M_PI*fc;
            double K = tan(wc);        
            b0 = K/(K+1);
            b1 = K/(K+1);
        }
    
        double Tick(double in) {       
            double wc= M_PI*fc;
            double K = tan(wc);        
            b0 =  K/(K+1);
            b1 =  K/(K+1);       
            x = in;
            y = b0*x + b1*x1;
            x1= x;        
            hp = 1 - y;
            return y;

        }
    };

    ////////////////////////////////////////////////////////////////////////////////////////////
    // it's used to create a peculiar notch
    // It can be calculaated with windowed sinc as it is a FIR filter with only 2 coefficients
    ////////////////////////////////////////////////////////////////////////////////////////////
    struct TwoZeros
    {
        double b0,b1,b2;
        double x1,x2;
        double x,y;
        double fc,fs,R;
        double hp;

        
        TwoZeros(double sr, double f, float r)
        {
            x1 = 0.0f;        
            fc = f/sr;
            fs = sr;
            R  = r;
            setCutoff(f);
        }
        // this doesn't really do much
        void setCutoff(double f) {
            fc = f/fs;
            b0 = 1;
            b1 = -2*R*cos(2*M_PI*fc);
            b2 = R*R;
        }
        // this creates the notch
        void setQ(double q) 
        {
            R  = q;
            if(R == 0.0f) R = 0.005f;
            if(R == 1.0f) R = 0.995f;
            b0 = 1;
            b1 = -2*R*cos(2*M_PI*fc);
            b2 = R*R;                
        }
        double Tick(double in) {               
            x = in;
            y = b0*x + b1*x1 + b2 * x2;
            x2=x1;
            x1= x;                
            return y;

        }
    };
}
