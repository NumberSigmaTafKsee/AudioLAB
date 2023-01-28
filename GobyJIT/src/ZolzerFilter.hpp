#pragma once
#include <cmath>
namespace Filters::Zolzer
{
    ////////////////////////////////////////////////////////////////////////////////////////////
    // Zolzer
    ////////////////////////////////////////////////////////////////////////////////////////////
    struct ZolzerBiquad
    {
        double a[2];
        double b[3];
        double fs,fc,q,g;
        double x1,x2,y1,y2;
        double x,y; 
        double res;   

        enum FilterType filter_type;

        ZolzerBiquad(FilterType type, double Fs, double Fc, double G = 1, double Q=0.707)
        {
            fs = Fs;
            fc = Fc;
            q  = Q;
            g = G;
            res = 0;
            x1=x2=y1=y2=0;
            filter_type = type;
            init_filter(Fc,Q);        
        }

        void init_filter(double Fc, double Q, double gain=1)
        {
            fc = Fc/fs*0.5;
            
            q = Q;
            g = gain;

            switch(filter_type)
            {
                case Lowpass: lowpass(fc,q); break;
                case Highpass: highpass(fc,q); break;
                case Bandpass: bandpass(fc,q); break;
                case Notch: notch(fc,q); break;
                // fc/q dont matter q must be 0
                case Allpass: allpass(fc,0); break;
                //have to find it
                //case Peak: peak(fc,q,gain); break;
                //case Lowshelf: lowshelf(fc,q); break;
                //case Highshelf: highshelf(fc,q); break;
                default: assert(1==0);
            }
        }
        
        void setCutoff(double f) {
            fc = f;
            init_filter(fc,q,g);
        }
        void setQ(double Q) {
            q  = Q;
            init_filter(fc,q,g);
        }
        void setResonance(double R) 
        {
            res = R;
        }
        void setGain(double G) {
            g = G;
            init_filter(fc,q,g);
        }

        void notch(double f, double Q) {
            fc = f;
            q  = Q;
            double K = tan(M_PI*fc);
            double Kq = Q*(1+K*K) ;
            double Kk = (K*K*Q+K+Q);        
            b[0] = Kq/Kk;
            b[1] = (2*Kq)/Kk;
            b[2] = Kq/Kk;
            a[0] = (2*Q*(K*K-1))/Kk;
            a[1] = (K*K*Q-K+Q)/Kk;
        }
        void lowpass1p(double f)
        {
            fc = f;
            q  = 0;
            double K = tan(M_PI*fc);
            b[0] = K/(K+1);
            b[1] = K/(K+1);
            b[2] = 0;
            a[0] = (K-1)/(K+1);
            a[1] = 0;
        }
        void highpass1p(double f)
        {
            fc = f;
            q  = 0;
            double K = tan(M_PI*fc);
            b[0] = 1/(K+1);
            b[1] = -1/(K+1);
            b[2] = 0;
            a[0] = (K-1)/(K+1);
            a[1] = 0;
        }
        void allpass1p(double f)
        {
            fc = f;
            q  = 0;
            double K = tan(M_PI*fc);
            b[0] = (K-1)/(K+1);
            b[1] = 1;
            b[2] = 0;
            a[0] = (K-1)/(K+1);
            a[1] = 0;
        }
        void lowpass(double f, double Q) {
            fc = f;
            q  = Q;
            double K = tan(M_PI*fc);
            double Kk = (K*K*Q+K+Q);        
            double Kq = (K*K*Q);
            b[0] = Kq/Kk;
            b[1] = (2*Kq) /Kk;
            b[2] =  Kq / Kk;
            a[0] = (2*Q*(K*K-1))/Kk;
            a[1] = (K*K*Q-K+Q)/Kk;
        }
        void allpass(double f, double Q) {
            fc = f;                
            q  = Q;
            double K = tan(M_PI*fc);
            double Kk = (K*K*Q+K+Q);        
            double Km = (K*K*Q-K+Q);
            double Kq = 2*Q*(K*K-1);
            b[0] = Km/Kk;
            b[1] = Kq/Kk;
            b[2] = 1.0f;
            a[0] = Kq/Kk;
            a[1] = Km/Kk;
        }
        void highpass(double f, double Q) {
            fc = f;
            q  = Q;
            double K = tan(M_PI*fc);
            double Kk = (K*K*Q+K+Q); 
            double Kq = 2*Q*(K*K-1);
            double Km = (K*K*Q-K+Q);
            b[0] = Q / Kk;
            b[1] = -(2*Q)/Kk;
            b[2] = Q / Kk;
            a[1] = Kq/Kk;
            a[2] = Km/Kk;
        }    
        void bandpass(double f, double Q) {
            fc = f;
            q  = Q;
            double K = tan(M_PI*fc);
            double Kk = (K*K*Q+K+Q); 
            b[0] = K / Kk;
            b[1] = 0;
            b[2] = -b[0];
            a[0] = (2*Q*(K*K-1))/Kk;
            a[1] = (K*K*Q-K+Q)/Kk;
        }
        // lowshelf
        void lfboost(double f, double G)
        {
            fc = f;
            g  = G;
            double K = tan(M_PI*fc);
            double V0= pow(10,G/20.0);
            double Kaka1 = sqrt(2*V0) * K + V0*K*K;
            double Kaka2 = 1 + sqrt(2)*K + K*K;
            b[0] = (1+Kaka1)/Kaka2;
            b[1] = (2*(V0*K*K-1))/ Kaka2;
            b[2] = (1 - Kaka1)/Kaka2;
            a[0] = (2*(K*K-1))/Kaka2;
            a[1] = (1-sqrt(2)*K+K*K)/Kaka2;
        }
        // lowshelf
        void lfcut(double f, double G)
        {
            fc = f;
            g  = G;
            double K = tan(M_PI*fc);
            double V0= pow(10,G/20.0);
            double Kaka = V0 + sqrt(2*V0)*K + K*K;
            b[0] = (V0*(1+sqrt(2)*K+K*K))/Kaka;
            b[1] = (2*V0*(K*K-1))/ Kaka;
            b[2] = (V0*(1-sqrt(2)*K+K*K))/Kaka;
            a[0] = (2*(K*K-V0))/Kaka;
            a[1] = (V0-sqrt(2*V0)*K+K*K)/Kaka;
        }
        // hishelf
        void hfboost(double f, double G)
        {
            fc = f;
            g  = G;
            double K = tan(M_PI*fc);
            double V0= pow(10,G/20.0);            
            double Kaka = 1 + sqrt(2)*K + K*K;
            b[0] = (V0 + sqrt(2*V0)*K + K*K)/Kaka;
            b[1] = (2*(K*K-V0))/Kaka;
            b[2] = (V0 - sqrt(2*V0)*K + K*K)/Kaka;
            a[0] = (2*(K*K-1))/Kaka;
            a[1] = (1-sqrt(2*K)+K*K)/Kaka;
        }
        // hishelf
        void hfcut(double f, double G)
        {
            fc = f;
            g  = G;
            double K = tan(M_PI*fc);
            double V0= pow(10,G/20.0);            
            double Kaka = 1 + sqrt(2*V0)*K + V0*K*K;
            b[0] = (V0*(1 + sqrt(2)*K + K*K))/Kaka;
            b[1] = (2*V0*(K*K-1))/Kaka;
            b[2] = (V0*(1 - sqrt(2)*K + K*K))/Kaka;
            a[0] = (2*(V0*K*K-1))/Kaka;
            a[1] = (1-sqrt(2*V0)*K + V0*K*K)/Kaka;
        }
        // peak
        void boost(double f, double Q, double G)
        {
            fc = f;
            g  = G;
            q  = Q;
            double K = tan(M_PI*fc);
            double V0= pow(10,G/20.0);            
            double Kaka = 1 + (1/Q)*K + K*K;
            b[0] = (1+(V0/Q)*K + K*K)/Kaka;
            b[1] = (2*(K*K-1))/Kaka;
            b[2] = (1- (V0/Q)*K + K*K)/Kaka;
            a[0] = (2*(K*K-1))/Kaka;
            a[1] = (1 - (1/Q)*K + K*K)/Kaka;
        }
        //peak
        void cut(double f, double Q, double G)
        {
            fc = f;
            g  = G;
            q  = Q;
            double K = tan(M_PI*fc);
            double V0= pow(10,G/20.0);            
            double Kaka = 1 + (1/(V0*Q)*K + K*K);
            b[0] = (1 + (1/Q)*K + K*K)/Kaka;
            b[1] = (2*(K*K-1))/Kaka;
            b[2] = (1 - (1/Q)*Kaka *K*K)/Kaka;
            a[0] = (2*(K*K-1))/Kaka;
            a[1] = (1 - (1/(V0*Q)*K +K*K))/Kaka;
        }
        
        double Tick(double I, double A = 1, double X = 0, double Y = 0)
        {
            Undenormal denormal;
            x = I; 
            y = b[0]*x + b[1]*x1 + b[2]*x2 - a[0]*y1 - a[1]*y2;        
            y2 = y1;
            y1 = y;
            x2 = x1;
            x1 = x;        
            return y;
        }
    };
}        
