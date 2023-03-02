#pragma once

struct TwoZeroNotch
{
    double b0,b1,b2;
    double fc,sr,R;
    double y,x0,x1,x2;

    TwoZeroNotch()
    {
        fc = 440.0f;
        sr = sampleRate;
        // is complex sqrt(xp^2 + yp^2) but only real here
        // it must be < 1
        R  = 0.9;
        y = x1 = x2 = 0;
        setCutoff(fc);
    }
    void setRadius(double r) {
        if(r < 0 || r >= 1.0) return;
        R = r;
    }
    void setCutoff(double f)
    {
        if(f < 0 || f > sr/2) return;
        fc = f;
        double theta = 2*M_PI*fc/sr;
        b0 = 1.0;
        b1 = -2.0 * R * cos(theta);
        b2 = R*R;
        if(b1 > 0) b0 = 1 / ( 1 + b1 + b2);
        else b0 = 1 / (1 - b1 + b2);
        b1 *= b0;
        b2 *= b0;
    }
    double Tick(double I, double A=1, double X=1, double Y=1)
    {
        double f = fc;
        double r = R;
        R *= fabs(Y);
        setCutoff(f*fabs(X));
        x0 = I;
        y = b0*x0 + b1*x1 + b2*x2;
        x2 = x1;
        x1 = x0;
        R = r;
        setCutoff(f);
        return A*y;
    }
};
