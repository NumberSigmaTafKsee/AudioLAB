/*
# Ideal impulse response
Low pass	2*fC*sin(nwC)/nwC
High pass	−2*fC*sin(nwC)/nwC
Band pass	2*f2*sin(nw2)/nw2−2*f1*sin(nw1)/nw1
Band stop	2*f1*sin(nw1)/nw1−2*f2*sin(nw2)/nw2

High pass = -LP
Bandpass  = HP - LP
Bandstop  = -BP
Lowshelf  = GLP + HP
HighShelf = LP + GHP
Peak/Notch= GBP + BS
*/


#pragma once

#include <stdlib.h>
#include <cmath>
#include <string.h>
#include <string.h>
#include <vector>

namespace Filters::FIR
{
    class FIR_filter
    {
    public:
        FIR_filter( int taps=0, DspFloatType f1=0, DspFloatType f2=0, char* type="", 
                    char* window="");
        ~FIR_filter();

        std::vector<DspFloatType> getCoefficients();

        DspFloatType filter(DspFloatType new_sample);

    private:
        std::vector<DspFloatType> lowPass_coefficient( int taps, DspFloatType f);
        std::vector<DspFloatType> highPass_coefficient(int taps, DspFloatType f);
        std::vector<DspFloatType> bandPass_coefficient(int taps, DspFloatType f1, DspFloatType f2);
        std::vector<DspFloatType> bandStop_coefficient(int taps, DspFloatType f1, DspFloatType f2);

        std::vector<DspFloatType> window_hamming(int taps);
        std::vector<DspFloatType> window_triangle(int taps);
        std::vector<DspFloatType> window_hanning(int taps);
        std::vector<DspFloatType> window_blackman(int taps);

        std::vector<DspFloatType> h;       // FIR coefficients
        std::vector<DspFloatType> samples; // FIR delay

        int idx;        // Round robin index
        int taps;       // Number of taps of the filter
    };

    static DspFloatType sinc(const DspFloatType x)
    {
        if (x == 0)
            return 1;

        return sin(M_PI * x) / (M_PI * x);
    }

    FIR_filter::FIR_filter( int taps, DspFloatType f1, DspFloatType f2, char* type,
                            char* window): h(taps, 0), samples(taps, 0)
    {
        this->idx     = 0;
        this->taps    = taps;

        std::vector<DspFloatType> h;  // Buffer FIR coefficients
        std::vector<DspFloatType> w;  // Buffer window coefficients

        // Calculate the coefficient corresponding to the filter type
        if (!strcmp(type, "lp")) {
            h = lowPass_coefficient(taps, f1);
        }
        else if (!strcmp(type, "hp")) {
            h = highPass_coefficient(taps, f1);
        }
        else if (!strcmp(type, "bp")) {
            h = bandPass_coefficient(taps, f1, f2);
        }
        else if (!strcmp(type, "sb")) {
            h = bandStop_coefficient(taps, f1, f2);
        }

        // Calculate the window to improve the FIR filter
        if (!strcmp(window, "hamming")) {
            w = window_hamming(taps);
        }
        else if (!strcmp(window, "triangle")) {
            w = window_triangle(taps);
        }
        else if (!strcmp(window, "hanning")) {
            w = window_hanning(taps);
        }
        else if (!strcmp(window, "blackman")) {
            w = window_blackman(taps);
        }

        if (!strcmp(window, "")) {
            this->h = h;
        }
        else
        {
            for(int n = 0; n < taps; n++) {
                this->h[n] = h[n] * w[n];
            }
        }
    }

    FIR_filter::~FIR_filter()
    {

    }

    std::vector<DspFloatType> FIR_filter::getCoefficients()
    {
        return this->h;
    }

    std::vector<DspFloatType> FIR_filter::lowPass_coefficient(int taps, DspFloatType f)
    {
        std::vector<int>    n(taps, 0);
        std::vector<DspFloatType> h(taps, 0);

        for(int i = 0; i < taps; i++) {
            n[i] = i - int(taps/2);
        }

        for(int i = 0; i < taps; i++) {
            h[i] = 2.0*f*sinc(2.0*f*n[i]);
        }

        return h;
    }

    std::vector<DspFloatType> FIR_filter::highPass_coefficient(int taps, DspFloatType f)
    {
        std::vector<int>    n(taps, 0);
        std::vector<DspFloatType> h(taps, 0);

        for(int i = 0; i < taps; i++) {
            n[i] = i - int(taps/2);
        }

        for(int i = 0; i < taps; i++) {
            h[i] = sinc(n[i]) - 2.0*f*sinc(2.0*f*n[i]);
        }

        return h;
    }

    std::vector<DspFloatType> FIR_filter::bandPass_coefficient(int taps, DspFloatType f1, DspFloatType f2)
    {
        std::vector<int>    n(taps, 0);
        std::vector<DspFloatType> h(taps, 0);

        for(int i = 0; i < taps; i++) {
            n[i] = i - int(taps/2);
        }

        for(int i = 0; i < taps; i++) {
            h[i] = 2.0*f1*sinc(2.0*f1*n[i]) - 2.0*f2*sinc(2.0*f2*n[i]);
        }

        return h;
    }

    std::vector<DspFloatType> FIR_filter::bandStop_coefficient(int taps, DspFloatType f1, DspFloatType f2)
    {
        std::vector<int>    n(taps, 0);
        std::vector<DspFloatType> h(taps, 0);

        for(int i = 0; i < taps; i++) {
            n[i] = i - int(taps/2);
        }

        for(int i = 0; i < taps; i++) {
            h[i] = 2.0*f1*sinc(2.0*f1*n[i]) - 2.0*f2*sinc(2.0*f2*n[i]) + sinc(n[i]);
        }

        return h;
    }

    std::vector<DspFloatType> FIR_filter::window_hamming(int taps)
    {
        std::vector<int>    n(taps, 0);
        std::vector<DspFloatType> w(taps, 0);

        DspFloatType alpha   = 0.54;
        DspFloatType beta    = 0.46;

        for(int i = 0; i < taps; i++) {
            w[i] = alpha - beta * cos(2.0 * M_PI * i / (taps - 1));
        }

        return w;
    }

    std::vector<DspFloatType> FIR_filter::window_hanning(int taps)
    {
        std::vector<DspFloatType> w(taps, 0);

        for(int i = 0; i < taps; i++) {
            w[i] =  sin(((DspFloatType) M_PI * i) / (taps - 1)) *
                    sin(((DspFloatType) M_PI * i) / (taps - 1));
        }

        return w;
    }

    std::vector<DspFloatType> FIR_filter::window_triangle(int taps)
    {
        std::vector<DspFloatType> w(taps, 0);

        DspFloatType l = taps;

        for(int i = 0; i < taps; i++) {
            w[i] = 1 - abs((i - (((DspFloatType)(taps-1)) / 2.0)) / (((DspFloatType)l) / 2.0));
        }

        return w;
    }

    std::vector<DspFloatType> FIR_filter::window_blackman(int taps)
    {
        std::vector<DspFloatType> w(taps, 0);

        DspFloatType alpha0 = 0.42;
        DspFloatType alpha1 = 0.5;
        DspFloatType alpha2 = 0.08;

        for(int i = 0; i < taps; i++) {
            w[i] = alpha0 - alpha1 * cos(2.0 * M_PI * i / (taps - 1))
                        - alpha2 * cos(4.0 * M_PI * i / (taps - 1));
        }

        return w;
    }

    DspFloatType FIR_filter::filter(DspFloatType new_sample)
    {
        DspFloatType result = 0;

        // Save the new sample
        this->samples[this->idx] = new_sample;

        // Calculate the output
        for(int n = 0; n < this->taps; n++)
            result += this->samples[(this->idx + n) % this->taps] * this->h[n];

        // Increase the round robin index
        this->idx = (this->idx + 1) % this->taps;

        return result;
    }
}
