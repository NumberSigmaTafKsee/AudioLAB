#pragma once

////////////////////////////////////////////////////////////////////////////////////////////
// -6db One Pole/One Zero
////////////////////////////////////////////////////////////////////////////////////////////

namespace KfrDSP1
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

    struct FilterBase
    {
        virtual double Tick(double I, double A=1, double X=1, double Y=1) = 0;

        void ProcessBlock(size_t n, float * in, float * out) {
            for(size_t i = 0; i < n; i++) out[i] = Tick(in[i]);
        }
        void ProcessBlock(size_t n, float * in, float * out, float * A, float *X=NULL, float *Y=NULL) {
            for(size_t i = 0; i < n; i++) out[i] = Tick(in[i],A[i], X? X[i]:1.0,Y? Y[i]:1.0);
        }
        void ProcessBlock(size_t n, double * in, double * out) {
            for(size_t i = 0; i < n; i++) out[i] = Tick(in[i]);
        }
        void ProcessBlock(size_t n, double * in, double * out, double * A, double *X=NULL, double *Y=NULL) {
            for(size_t i = 0; i < n; i++) out[i] = Tick(in[i],A[i], X? X[i]:1.0,Y? Y[i]:1.0);
        }
    };
    struct Biquad6DB : public FilterBase
    {
        double a[2];
        double b[3];
        double fs,fc;
        double x1,x2,y1,y2;
        double x,y;

        FilterType filterType = Lowpass;

        Biquad6DB(FilterType type, double Fs, double Fc) {
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
            double K = std::tan(M_PI*fc);
            b[0] = K/(K+1);
            b[1] = K/(K+1);
            b[2] = 0.0;
            a[0] = (K-1)/(K+1);
            a[1] = 0.0;
        }
        void highpass(double fc)
        {
            double K = std::tan(M_PI*fc);
            b[0] = 1/(K+1);
            b[1] = -1/(K+1);
            b[2] = 0.0;
            a[0] = (K-1)/(K+1);
            a[1] = 0.0;
        }
        void allpass(double fc)
        {
            double K = std::tan(M_PI*fc);
            b[0] = (K-1)/(K+1);
            b[1] = 1;
            b[2] = 0.0;
            a[0] = (K-1)/(K+1);
            a[1] = 0.0;
        }
        
        double Tick(double I)
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
    struct Biquad12DB : public FilterBase
    {
        double a[2];
        double b[3];
        double fs,fc,q,g;
        double x1,x2,y1,y2;
        double x,y;

        FilterType filterType = Lowpass;

        Biquad12DB() = default;

        Biquad12DB(FilterType type, double Fs, double Fc, double G = 1, double Q=0.707) 
        {
            fs = Fs;
            fc = Fc;
            q  = Q;
            g = G;
            x1=x2=y1=y2=0;            
            init_filter(type,Fc,Q,G);        
        }
        Biquad12DB(const kfr::biquad_params<double>& bq, double Fs, double Fc, double G = 1, double Q=0.707)
        {
            fs = Fs;
            fc = Fc;
            q  = Q;
            g = G;
            x1=x2=y1=y2=0;        
            setCoefficients(bq);
        }
        

        void init_filter(FilterType type, double Fc, double Q=0.707, double gain=1, double Fs=44100.0)
        {
            fs = Fs;
            filterType = type;
            fc = Fc/fs*0.99;        
            q = Q;
            g = gain;

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

        void setCoefficients(kfr::biquad_params<double> p)
        {
            if(p.a0 == 0) p.a0=1.0f;        
            a[0] = p.a1/p.a0;
            a[1] = p.a2/p.a0;
            b[0] = p.b0/p.a0;
            b[1] = p.b1/p.a0;
            b[2] = p.b2/p.a0;
        }
        void setCutoff(double f) {
            fc = f;
            init_filter(filterType,fc,q,g,fs);
        }
        void setQ(double Q) {
            q  = Q;
            init_filter(filterType,fc,q,g,fs);
        }
        void setGain(double G) {
            g = G;
            init_filter(filterType,fc,q,g,fs);
        }
        void setSampleRate(double sr) {
            fs = sr;
            init_filter(filterType,fc,q,g,sr);
        }
        void setType(FilterType type) {
            filterType = type;
            init_filter(filterType,fc,q,g,fs);
        }
        void notch(double f, double Q) {
            fc = f;
            q  = Q;
            kfr::biquad_params<double> p  = kfr::biquad_notch(fc,q);  
            setCoefficients(p);
        }
        void lowpass(double f, double Q) {
            fc = f;
            q  = Q;
            kfr::biquad_params<double> p  = kfr::biquad_lowpass(fc,q);        
            setCoefficients(p);
        }
        void allpass(double f, double Q) {
            fc = f;
            q  = Q;
            kfr::biquad_params<double> p  = kfr::biquad_allpass(fc,q);        
            setCoefficients(p);
        }
        void highpass(double f, double Q) {
            fc = f;
            q  = Q;
            kfr::biquad_params<double> p  = kfr::biquad_highpass(fc,q);        
            setCoefficients(p);
        }
        void peak(double f, double Q, double gain) {
            fc = f;
            q  = Q;
            kfr::biquad_params<double> p  = kfr::biquad_peak(fc,q, gain);        
            setCoefficients(p);
        }
        void lowshelf(double f, double Q) {
            fc = f;
            q  = Q;
            kfr::biquad_params<double> p  = kfr::biquad_lowshelf(fc,q);        
            setCoefficients(p);
        }
        void highshelf(double f, double Q) {
            fc = f;
            q  = Q;
            kfr::biquad_params<double> p  = kfr::biquad_highshelf(fc,q);        
            setCoefficients(p);
        }
        void bandpass(double f, double Q) {
            fc = f;
            q  = Q;
            kfr::biquad_params<double> p  = kfr::biquad_bandpass(fc,q);        
            setCoefficients(p);
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

    struct BiquadCascade : FilterBase
    {
        std::vector<Biquad12DB> biquads;

        BiquadCascade(size_t num, FilterType type = Lowpass, double fc = 1000.0, double q = 0.707, double gain = 1.0, double sr=44100.0 ) {
            biquads.resize(num);
            for(size_t i = 0; i < num; i++)
                biquads[i].init_filter(type,fc,pow(q,num),gain,sr);
        }
        void setCutoff(double f) {
            int num = biquads.size();
            for(size_t i = 0; i < num; i++)
                biquads[i].setCutoff(f);            
        }
        void setQ(double q) {
            int num = biquads.size();
            for(size_t i = 0; i < num; i++)
                biquads[i].setQ(pow(q,num));
        }
        void setGain(double g) {
            int num = biquads.size();
            for(size_t i = 0; i < num; i++)
                biquads[i].setGain(g);
        }
        void setType(FilterType type) {
            int num = biquads.size();
            for(size_t i = 0; i < num; i++)
                biquads[i].setType(type);            
        }
        double Tick(double I, double A=1, double X=1, double Y=1) {
            double o = biquads[0].Tick(I,A,X,Y);
            for(size_t i = 1; i < biquads.size(); i++)
                o = biquads[i].Tick(o,A,X,Y);
            return o;
        }        
    };
}