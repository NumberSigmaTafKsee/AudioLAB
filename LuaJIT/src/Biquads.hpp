#pragma once

////////////////////////////////////////////////////////////////////////////////////////////
// -6db One Pole/One Zero
////////////////////////////////////////////////////////////////////////////////////////////

#include "Filters.hpp"
#include "SoundObject.hpp"

namespace Filters::Biquads
{
    enum FilterType
    {
        Lowpass,
        Highpass,
        Bandpass,  
        Bandpass2, // this is used in RBJ for the cszap whatevr
        Notch,
        Bandstop,
        Allpass,
        Peak,
        Lowshelf,
        Highshelf,        
    };

    struct Biquad6DB : public FilterProcessor
    {
        double a[2];
        double b[3];
        double fs,fc;
        double x1,x2,y1,y2;
        double x,y;

        FilterType filterType = Lowpass;

        Biquad6DB(FilterType type, double Fs, double Fc) : FilterProcessor() {
            fs = Fs;
            fc = Fc/Fs;
            setFilter(type);
        }
        void setFilter(FilterType type) {
            filterType = type;
            switch(type) {
                case Lowpass: lowpass(fc); break;
                case Highpass: highpass(fc); break;
                case Allpass: allpass(fc); break;
            }
        }
        void setCutoff(float f) {
            fc = f/fs;
            setFilter(filterType);
        }
        void setQ(float q) {
            // does nothing right now
        }

        void lowpass(double fc)
        {
            double K = tan(M_PI*fc);
            b[0] = K/(K+1);
            b[1] = K/(K+1);
            b[2] = 0.0;
            a[0] = (K-1)/(K+1);
            a[1] = 0.0;
        }
        void highpass(double fc)
        {
            double K = tan(M_PI*fc);
            b[0] = 1/(K+1);
            b[1] = -1/(K+1);
            b[2] = 0.0;
            a[0] = (K-1)/(K+1);
            a[1] = 0.0;
        }
        void allpass(double fc)
        {
            double K = tan(M_PI*fc);
            b[0] = (K-1)/(K+1);
            b[1] = 1;
            b[2] = 0.0;
            a[0] = (K-1)/(K+1);
            a[1] = 0.0;
        }
        
        double Tick(double I, double A=1, double X=0, double Y=0)
        {
            Undenormal denormal;        
            x = I;
            y = b[0]*x + b[1]*x1 - a[0] * y1;
            x1 = x;
            y1 = y;
            return y;
        }
        
    };

    ////////////////////////////////////////////////////////////////////////////////////////////
    // -12db Two Pole/Two Zero 1 section
    ////////////////////////////////////////////////////////////////////////////////////////////
    struct Biquad12DB : public FilterProcessor
    {
        double a[2];
        double b[3];
        double fs,fc,q,g;
        double x1,x2,y1,y2;
        double x,y;
        double bQ;

        FilterType filterType = Lowpass;

        Biquad12DB(FilterType type, double Fs, double Fc, double G = 1, double Q=0.707)
        : FilterProcessor()
        {
            fs = Fs;
            fc = Fc;
            q  = Q;
            bQ = Q;
            g = G;
            x1=x2=y1=y2=0;
            filterType = type;
            init_filter(Fc,Q);        
        }
        Biquad12DB(Filters::FilterCoefficients& bq, double Fs, double Fc, double G = 1, double Q=0.707)
        {
            fs = Fs;
            fc = Fc;
            q  = Q;
            g = G;
            x1=x2=y1=y2=0;        
            setCoefficients(bq);
        }
        

        void init_filter(double Fc, double Q, double gain=1)
        {
            fc = Fc;
            q  = Q;
            g  = gain;

            switch(filterType)
            {
                case Lowpass: lowpass(fc,q); break;
                case Highpass: highpass(fc,q); break;
                case Bandpass: bandpass(fc,q); break;
                case Notch: notch(fc,q); break;
                // fc/q dont matter q must be 0
                case Allpass: allpass(fc,0); break;
                case Peak: peak(fc,q,gain); break;
                case Lowshelf: lowshelf(fc,q); break;
                case Highshelf: highshelf(fc,q); break;
                default: assert(1==0);
            }
        }

        void setCoefficients(Filters::FilterCoefficients p)
        {
            a[0] = p.a[0];
            a[1] = p.a[1];
            b[0] = p.b[0];
            b[1] = p.b[1];
            b[2] = p.b[2];
        }
        void setCutoff(double f) {
            fc = f;
            init_filter(fc,q+bQ,g);
        }
        void setQ(double Q) {
            q  = Q;
            init_filter(fc,q+bQ,g);
        }
        void setGain(double G) {
            g = G;
            init_filter(fc,q,g);
        }


        void notch(double f, double Q) {
            fc = f;
            q  = Q;
	    Filters::FilterCoefficients c = Filters::NotchBiquad(fc,fs,q);
            setCoefficients(c);
        }
        void lowpass(double f, double Q) {
            fc = f;
            q  = Q;
	    Filters::FilterCoefficients c = Filters::LowpassBiquad(fc,fs,q);
            setCoefficients(c);
        }
        void allpass(double f, double Q) {
            fc = f;
            q  = Q;
	    Filters::FilterCoefficients c = Filters::AllpassBiquad(fc,fs,q);
            setCoefficients(c);
        }
        void highpass(double f, double Q) {
            fc = f;
            q  = Q;
 	    Filters::FilterCoefficients c = Filters::HighpassBiquad(fc,fs,q);
            setCoefficients(c);
        }
        void peak(double f, double Q, double gain) {
            fc = f;
            q  = Q;
            g  = gain;
	    Filters::FilterCoefficients c = Filters::PeakBiquad(fc,fs,q,gain);
            setCoefficients(c);

        }
        void lowshelf(double f, double Q) {
            fc = f;
            q  = Q;
	    Filters::FilterCoefficients c = Filters::LowshelfBiquad(fc,fs,q);
            setCoefficients(c);
        }
        void highshelf(double f, double Q) {
            fc = f;
            q  = Q;
	    Filters::FilterCoefficients c = Filters::HighshelfBiquad(fc,fs,q);
            setCoefficients(c);
        }
        void bandpass(double f, double Q) {
            fc = f;
            q  = Q;
	    Filters::FilterCoefficients c = Filters::BandpassBiquad(fc,fs,q);
            setCoefficients(c);
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




    ////////////////////////////////////////////////////////////////////////////////////////////
    // -24db
    ////////////////////////////////////////////////////////////////////////////////////////////

    struct Biquad24DB : public FilterProcessor
    {
        Biquad12DB *first,*second;
        double fs,fc;
        double Q1,Q2;
        double bQ1,bQ2;

        Biquad24DB(FilterType type, double Fs, double Fc, double G = 1, double Q1 = 0.54119610f, double Q2=1.3065460f)
        : FilterProcessor()
        {
            fs = Fs;
            fc = Fc;
            bQ1 = Q1;
            bQ2 = Q2;        
            first = new Biquad12DB(type,Fs,Fc,G,Q1);
            second = new Biquad12DB(type,Fs,Fc,G,Q2);
            this->Q1 = Q1;
            this->Q2 = Q2;
        }
        ~Biquad24DB()
        {
            if(first) delete first;
            if(second) delete second;
        }
        void setCutoff(double f) {
            fc = f;
            first->lowpass(fc,Q1+bQ1);
            second->lowpass(fc,Q2+bQ2);
        }
        void setQ(double q1,double q2) {
            Q1 = q1;
            Q2 = q2;        
            first->lowpass(fc,Q1+bQ1);
            second->lowpass(fc,Q2+bQ2);
        }
        double Tick(double I, double A=1, double X=0, double Y=0)   {
            return second->Tick(first->Tick(I));
        }
    };


    ////////////////////////////////////////////////////////////////////////////////////////////
    // -36db
    ////////////////////////////////////////////////////////////////////////////////////////////

    struct Biquad36DB : public FilterProcessor
    {
        Biquad12DB *first,*second,*third;
        double fs,fc;
        double Q1,Q2,Q3;
        double bQ1,bQ2,bQ3;        

        Biquad36DB(FilterType type, double Fs, double Fc, double G = 1, double Q2=0.70710678, double Q3=1.9318517)
        : FilterProcessor()
        {
            fs = Fs;
            fc = Fc;
            this->Q1 = Q1;
            this->Q2 = Q2;
            this->Q3 = Q3;        
            bQ1 = Q1;
            bQ2 = Q2;
            bQ3 = Q3;        
            first  = new Biquad12DB(type,Fs,Fc,G,Q1);
            second = new Biquad12DB(type,Fs,Fc,G,Q2);
            third  = new Biquad12DB(type,Fs,Fc,G,Q3);
        }
        ~Biquad36DB()
        {
            if(first) delete first;
            if(second) delete second;
            if(third) delete third;
        }
        void setCutoff(double f) {
            fc = f;
            first->lowpass(fc,Q1+bQ1);
            second->lowpass(fc,Q2+bQ2);
            third->lowpass(fc,Q3+bQ3);
        }
        void setQ(double q1,double q2, double q3) {
            Q1 = q1;
            Q2 = q2;
            Q3 = q3;        
            first->lowpass(fc,Q1+bQ1);
            second->lowpass(fc,Q2+bQ2);
            third->lowpass(fc,Q3+bQ3);
        }
        double Tick(double I, double A=1, double X=0, double Y=0)   {
            return third->Tick(second->Tick(first->Tick(I)));
        }
    };


    ////////////////////////////////////////////////////////////////////////////////////////////
    // -48db
    ////////////////////////////////////////////////////////////////////////////////////////////

    struct Biquad48DB : public FilterProcessor
    {
        Biquad12DB *first,*second,*third,*fourth;
        double fs,fc;
        double Q1,Q2,Q3,Q4;
        double bQ1,bQ2,bQ3,bQ4;

        Biquad48DB(FilterType type, double Fs, double Fc, double G = 1, double Q1 = 0.50979558, double Q2=0.60134489, double Q3=0.89997622, double Q4=2.5629154)
        : FilterProcessor()
        {
            fs = Fs;
            fc = Fc;        
            this->Q1 = Q1;
            this->Q2 = Q2;
            this->Q3 = Q3;        
            this->Q4 = Q4;        
            bQ1 = Q1;
            bQ2 = Q2;
            bQ3 = Q3;        
            bQ4 = Q4;        
            first = new Biquad12DB(type,Fs,Fc,G,Q1);
            second = new Biquad12DB(type,Fs,Fc,G,Q2);
            third = new Biquad12DB(type,Fs,Fc,G,Q3);
            fourth = new Biquad12DB(type,Fs,Fc,G,Q4);
        }
        ~Biquad48DB()
        {
            if(first) delete first;
            if(second) delete second;
            if(third) delete third;
            if(fourth) delete fourth;
        }
        void setCutoff(double f) {
            fc = f;
            first->lowpass(fc,Q1+bQ1);
            second->lowpass(fc,Q2+bQ2);
            third->lowpass(fc,Q3+bQ3);
            fourth->lowpass(fc,Q4+bQ4);
        }
        void setQ(double q1,double q2, double q3, double q4) {
            Q1 = q1;
            Q2 = q2;
            Q3 = q3;
            Q4 = q4;        
            first->lowpass(fc,Q1+bQ1);
            second->lowpass(fc,Q2+bQ2);
            third->lowpass(fc,Q3+bQ3);
            fourth->lowpass(fc,Q4+bQ4);
        }
        double Tick(double I, double A=1, double X=0, double Y=0)   {
            return fourth->Tick(third->Tick(second->Tick(first->Tick(I))));
        }
    };
}

