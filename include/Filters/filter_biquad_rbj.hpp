
#pragma once

// todo: invididual filter classes for each type(lp,hp,etc..)

namespace Filters::IIR::RBJFilters
{
//////////////////////////////////////////////////////////////////////////////////////////
// RBJ
//////////////////////////////////////////////////////////////////////////////////////////
  
    struct RBJBiquadFilter : public BiquadTransposedTypeII
    {

        FilterType filter_type;
        DspFloatType Fc, Fs, Q, G, R;

        RBJBiquadFilter(FilterType type = LOWPASS, DspFloatType sampleRate = 44100.0) : BiquadTransposedTypeII()
        {
            Fc = 1000.0;
            Fs = sampleRate;
            Q = 0.5;
            G = 1.0;
            filter_type = type;
            setCoefficients(Fc, Q, G);
        }        
        RBJBiquadFilter(FilterType type, DspFloatType fc, DspFloatType sr, DspFloatType q = 0.5, DspFloatType g = 1) : BiquadTransposedTypeII()
        {
            filter_type = type;
            Fs = sr;
            Q = q;
            G = g;
            setCoefficients(Fc, Q, G);
        }
        void setCoefficients(DspFloatType fc, DspFloatType q, DspFloatType g = 1)
        {            
            FilterCoefficients c;
            Fc = fc/4.0;
            Q = q;
            G = g;            
            switch (filter_type)
            {
            case LOWPASS:
                c = RBJLowpassBiquad(Fc, Fs, Q);
                break;
            case HIGHPASS:
                c = RBJHighpassBiquad(Fc, Fs, Q);
                break;
            case BANDPASS:
            case ZERODBBANDPASS:
                c = RBJBandpassConstant0dbBiquad(Fc, Fs, Q);
                break;
            case SKIRTBANDPASS:            
                c = RBJBandpassConstantSkirtBiquad(Fc, Fs, Q);
                break;
            case NOTCH:
                c = RBJNotchBiquad(Fc, Fs, Q);
                break;
            case PEAK:
                c = RBJPeakBiquad(Fc, Fs, Q, G);
                break;
            case LOWSHELF:
                c = RBJLowshelfBiquad(Fc, Fs, Q, G);
                break;
            case HIGHSHELF:
                c = RBJHighshelfBiquad(Fc, Fs, Q, G);
                break;
            case ALLPASS:
                c = RBJAllpassBiquad(Fc, Fs, Q);
                break;
            }
            biquad.setCoefficients(c);
        }
        void setCutoff(DspFloatType fc)
        {
            if(fc < 0 || fc >= Fs/2.0) return;
            setCoefficients(fc, Q, G);
        }
        void setQ(DspFloatType q)
        {
            if(q < 0 || q > 1000.0) return;
            setCoefficients(Fc, q, G);
        }
        void setRadius(DspFloatType r)
        {
            if(r < 0 || r >= 1.0) return;
            //setCoefficientsRadius(Fc, Q, (1 - r), G);
        }
        void setGain(DspFloatType g)
        {
            if(g < 0 || g > 10.0) return;
            setCoefficients(Fc, Q, g);
        }
    };

//////////////////////////////////////////////////////////////////////////////////////////
// lowpass
//////////////////////////////////////////////////////////////////////////////////////////

    struct RBJLowPassFilter : public RBJBiquadFilter
    {        
        DspFloatType Fc, Fs, Q, G, R;

        RBJLowPassFilter(DspFloatType sampleRate = 44100.0) : RBJBiquadFilter(LOWPASS,sampleRate)
        {
        }
        RBJLowPassFilter(DspFloatType fc, DspFloatType sr, DspFloatType q = 0.5, DspFloatType g = 1) 
        : RBJBiquadFilter(LOWPASS,fc,sr,q,g)
        {
        }
    };   

//////////////////////////////////////////////////////////////////////////////////////////
// highpass
//////////////////////////////////////////////////////////////////////////////////////////

    struct RBJHighPassFilter : public RBJBiquadFilter
    {        
        DspFloatType Fc, Fs, Q, G, R;

        RBJHighPassFilter(DspFloatType sampleRate = 44100.0) : RBJBiquadFilter(HIGHPASS,sampleRate)
        {
        }
        RBJHighPassFilter(DspFloatType fc, DspFloatType sr, DspFloatType q = 0.5, DspFloatType g = 1) 
        : RBJBiquadFilter(HIGHPASS,fc,sr,q,g)
        {
        }        
    };   

//////////////////////////////////////////////////////////////////////////////////////////
// allpass
//////////////////////////////////////////////////////////////////////////////////////////

    struct RBJAllPassFilter : public RBJBiquadFilter
    {        
        DspFloatType Fc, Fs, Q, G, R;

        RBJAllPassFilter() : RBJBiquadFilter(LOWPASS)
        {
        }
        RBJAllPassFilter(DspFloatType fc, DspFloatType sr, DspFloatType q = 0.5, DspFloatType g = 1) 
        : RBJBiquadFilter(ALLPASS,fc,sr,q,g)
        {
        }
    };   

//////////////////////////////////////////////////////////////////////////////////////////
// bandpass
//////////////////////////////////////////////////////////////////////////////////////////

    struct RBJBandPassFilter : public RBJBiquadFilter
    {        
        DspFloatType Fc, Fs, Q, G, R;

        RBJBandPassFilter() : RBJBiquadFilter(BANDPASS)
        {
        }
        RBJBandPassFilter(DspFloatType fc, DspFloatType sr, DspFloatType q = 0.5, DspFloatType g = 1) 
        : RBJBiquadFilter(BANDPASS,fc,sr,q,g)
        {
        }
    };   

//////////////////////////////////////////////////////////////////////////////////////////
// skirtband
//////////////////////////////////////////////////////////////////////////////////////////

    struct RBJSkirtBandPassFilter : public RBJBiquadFilter
    {        
        DspFloatType Fc, Fs, Q, G, R;

        RBJSkirtBandPassFilter(DspFloatType sampleRate = 44100.0) : RBJBiquadFilter(SKIRTBANDPASS,sampleRate)
        {
        }
        RBJSkirtBandPassFilter(DspFloatType fc, DspFloatType sr, DspFloatType q = 0.5, DspFloatType g = 1) 
        : RBJBiquadFilter(SKIRTBANDPASS,fc,sr,q,g)
        {
        }
    };   

//////////////////////////////////////////////////////////////////////////////////////////
// bandstop
//////////////////////////////////////////////////////////////////////////////////////////

    struct RBJBandStopFilter : public RBJBiquadFilter
    {        
        DspFloatType Fc, Fs, Q, G, R;

        RBJBandStopFilter() : RBJBiquadFilter(NOTCH)
        {
        }
        RBJBandStopFilter(DspFloatType fc, DspFloatType sr, DspFloatType q = 0.5, DspFloatType g = 1) 
        : RBJBiquadFilter(NOTCH,fc,sr,q,g)
        {
        }
    };   

//////////////////////////////////////////////////////////////////////////////////////////
// peak
//////////////////////////////////////////////////////////////////////////////////////////

    struct RBJPeakFilter : public RBJBiquadFilter
    {        
        DspFloatType Fc, Fs, Q, G, R;

        RBJPeakFilter(DspFloatType sampleRate = 44100.0) : RBJBiquadFilter(PEAK,sampleRate)
        {
        }
        RBJPeakFilter(DspFloatType fc, DspFloatType sr, DspFloatType q = 0.5, DspFloatType g = 1) 
        : RBJBiquadFilter(PEAK,fc,sr,q,g)
        {
        }
    }; 

//////////////////////////////////////////////////////////////////////////////////////////
// lowshelf
//////////////////////////////////////////////////////////////////////////////////////////

    struct RBJLowShelfFilter : public RBJBiquadFilter
    {        
        DspFloatType Fc, Fs, Q, G, R;

        RBJLowShelfFilter(DspFloatType sampleRate = 44100.0) : RBJBiquadFilter(LOWSHELF,sampleRate)
        {
        }
        RBJLowShelfFilter(DspFloatType fc, DspFloatType sr, DspFloatType q = 0.5, DspFloatType g = 1) 
        : RBJBiquadFilter(LOWSHELF,fc,sr,q,g)
        {
        }
    };   

//////////////////////////////////////////////////////////////////////////////////////////
// highshelf
//////////////////////////////////////////////////////////////////////////////////////////

    struct RBJHighShelfFilter : public RBJBiquadFilter
    {        
        DspFloatType Fc, Fs, Q, G, R;

        RBJHighShelfFilter() : RBJBiquadFilter(HIGHSHELF)
        {
        }
        RBJHighShelfFilter(DspFloatType fc, DspFloatType sr, DspFloatType q = 0.5, DspFloatType g = 1) 
        : RBJBiquadFilter(HIGHSHELF,fc,sr,q,g)
        {
        }
    };   
}