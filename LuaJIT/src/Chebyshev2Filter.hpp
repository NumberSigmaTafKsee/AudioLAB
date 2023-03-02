#pragma once

namespace FX::Filters::Chebyshev2
{
    struct Chebyshev2Filter : public FilterProcessor
    {
        double          a0,
                        a1,
                        a2,
                        b1,
                        b2;		/* coefficients */
        double          x1,x2,y1,y2,x,y;
        double          Fs,Fc,ripple;
        int             lowpass;

        Chebyshev2Filter(double sr, double fc, double _ripple, int low)	: FilterProcessor()       
        {

            x1=x2=y1=y2=x=y=0;
            Fs = sr;
            Fc = fc;
            ripple = _ripple*100.0;
            lowpass = low;
            setCoefficients(Fc);
            
        }
        void setCoefficients(double fc)
        {

            double          x,
                            y,
                            z,
                            c,
                            v,
                            t,
                            r,
                            om,
                            m,
                            x0,
                            x1,
                            x2,
                            y1p,
                            y2,
                            k,
                            d,
                            tt,
                            tt2;
            Fc = fc;
            c = -cos(M_PI / 4);
            v = sin(M_PI / 4);
            if (ripple > 0) {
                t = 100.0 / (100.0 - ripple);
                x = sqrt(t * t - 1);
                t = 1 / x;
                r = t + sqrt(t / x);
                y = 0.5 * log(r + 1);
                z = 0.5 * log(r - 1);
                t = exp(z);
                z = (t + 1 / t) / 2;
                t = exp(y);
                c *= (t - 1 / t) / 2 / z;
                v *= (t + 1 / t) / 2 / z;
            }
            tt = 2 * tan(0.5);
            tt2 = tt * tt;
            om = 2 * M_PI * Fc / Fs;
            m = c * c + v * v;
            d = 4 - 4 * c * tt + m * tt2;
            x0 = tt2 / d;
            x1 = x0 * 2;
            x2 = x0;
            y1p = (8 - 2 * m * tt2) / d;
            y2 = (-4 - 4 * c * tt - m * tt2) / d;
            if (lowpass)
                k = sin(0.5 - om / 2) / sin(0.5 + om / 2);
            else
                k = -cos(om / 2 + 0.5) / cos(om / 2 - 0.5);
            d = 1 + k * (y1p - y2 * k);
            a0 = (x0 - k * (x1 - x2 * k)) / d;        
            a1 = 2 * a0;
            a2 = a0;        
            b1 = (k * (2 + y1p * k - 2 * y2) + y1p) / d;
            b2 = (-k * (k + y1p) + y2) / d;
            if (!lowpass) {
                a1 = -a1;
                b1 = -b1;        
            }
        }
        void setCutoff(double f) {        
            Fc = f;
        }
        void setRipple(double r) {
            ripple = r;
        }
        double Tick(double I, double A=1, double X=1, double Y=1)
        {
            Undenormal denormals;
            double f = Fc;
            double r = ripple;
            Fc *= X;
            ripple *= Y;
            setCoefficients(Fc);
            x = I;
            y = x * a0 + x1 * a1 + x2 * a2 + b1 * y1 + b2 * y2;
            x2 = x1;
            x1 = x;
            y2 = y1;
            y1 = y;
            ripple = r;
            setCoefficients(f);
            return A*y;
        }
    };
}