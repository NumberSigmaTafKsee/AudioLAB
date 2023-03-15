#pragma once
#include "Filters/IIRFilters.hpp"

namespace Filters::IIR::ZolzerFilters
{

//////////////////////////////////////////////////////////////////////////////////////////
// Zolzer Biquad
//////////////////////////////////////////////////////////////////////////////////////////

    struct ZolzerBiquadFilter : public BiquadTransposedTypeII
    {
		enum ZolzerType {
			LOWPASS,
			LOWPASS1P,
			HIGHPASS,
			HIGHPASS1P,
			BANDPASS,
			NOTCH,
			PEAKBOOST,
			PEAKCUT,
			LOWSHELFBOOST,
			LOWSHELFCUT,
			HIGHSHELFBOOST,
			HIGHSHELFCUT,
			ALLPASS,
			ALLPASS1P
		};
		

        ZolzerType filter_type;
        DspFloatType Fc, Fs, Q, G, R;

        ZolzerBiquadFilter(ZolzerType type = LOWPASS,DspFloatType sampleRate = 44100.0) : BiquadTransposedTypeII()
        {
            Fc = 1000.0;
            Fs = sampleRate;
            Q = 0.5;
            G = 1.0;
            filter_type = type;
            setCoefficients(Fc, Q, G);
        }
        ZolzerBiquadFilter(ZolzerType type, DspFloatType fc, DspFloatType sr, DspFloatType q = 0.5, DspFloatType g = 1) : BiquadTransposedTypeII()
        {
            filter_type = type;
            Fs = sr;
            Q = q;
            G = g;
            setCoefficients(Fc, Q, G);
        }        
        
        void setCoefficients(DspFloatType fc, DspFloatType q, DspFloatType g=1.0)
        {                    
            FilterCoefficients c;   
            Fc =fc;         
            Q = q;
            G = g;
            
            switch (filter_type)
            {
            case LOWPASS:
                c = ZolzerLowpass(Fc, Fs, Q);
                break;
            case LOWPASS1P:
                c = ZolzerLowpass1p(Fc, Fs);
                break;
            case HIGHPASS:
                c = ZolzerHighpass(Fc, Fs, Q);
                break;
            case HIGHPASS1P:
                c = ZolzerHighpass1p(Fc, Fs);
                break;
            case BANDPASS:
                c = ZolzerBandpass(Fc, Fs, Q);
                break;
            case NOTCH:
                c = ZolzerNotch(Fc, Fs, Q);
                break;
            case PEAKBOOST:
                c = ZolzerBoost(Fc, Fs, Q, G);
                break;
            case PEAKCUT:
                c = ZolzerCut(Fc, Fs, Q, G);
                break;
            case LOWSHELFBOOST:
                c = ZolzerLFBoost(Fc, Fs, G);
                break;
            case LOWSHELFCUT:
                c = ZolzerLFCut(Fc, Fs, G);
                break;
            case HIGHSHELFBOOST:
                c = ZolzerHFBoost(Fc, Fs, G);
                break;
            case HIGHSHELFCUT:
                c = ZolzerHFCut(Fc, Fs, G);
                break;
            case ALLPASS:
                c = ZolzerAllpass(Fc, Fs, Q);
                break;
            case ALLPASS1P:
                c = ZolzerAllpass1p(Fc, Fs);
                break;
            }
            biquad.setCoefficients(c);
        }
        void setType(ZolzerType t) {
			filter_type = t;
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
        }
        void setGain(DspFloatType g)
        {
            if(g < 0 || g > 10.0) return;
            setCoefficients(Fc, Q, g);
        }
    };


//////////////////////////////////////////////////////////////////////////////////////////
// Lowpass
//////////////////////////////////////////////////////////////////////////////////////////

    struct ZolzerLowPassFilter : public ZolzerBiquadFilter
    {

        
        DspFloatType Fc, Fs, Q, G, R;

        ZolzerLowPassFilter(DspFloatType sampleRate = 44100.0) : ZolzerBiquadFilter(LOWPASS,sampleRate)
        {
        
        }

        ZolzerLowPassFilter(DspFloatType fc, DspFloatType sr, DspFloatType q = 0.5, DspFloatType g = 1) 
        : ZolzerBiquadFilter(LOWPASS,fc,sr,q,g)
        {

        }
    };

//////////////////////////////////////////////////////////////////////////////////////////
// Highpass
//////////////////////////////////////////////////////////////////////////////////////////

    struct ZolzerHighPassFilter : public ZolzerBiquadFilter
    {        
        DspFloatType Fc, Fs, Q, G, R;

        ZolzerHighPassFilter(DspFloatType sampleRate = 44100.0) : ZolzerBiquadFilter(HIGHPASS,sampleRate)
        {
        
        }

        ZolzerHighPassFilter(DspFloatType fc, DspFloatType sr, DspFloatType q = 0.5, DspFloatType g = 1) 
        : ZolzerBiquadFilter(HIGHPASS,fc,sr,q,g)
        {

        }
    };

//////////////////////////////////////////////////////////////////////////////////////////
// Allpass
//////////////////////////////////////////////////////////////////////////////////////////

    struct ZolzerAllPassFilter : public ZolzerBiquadFilter
    {        
        DspFloatType Fc, Fs, Q, G, R;

        ZolzerAllPassFilter(DspFloatType sampleRate = 44100.0) : ZolzerBiquadFilter(ALLPASS,sampleRate)
        {
        
        }

        ZolzerAllPassFilter(DspFloatType fc, DspFloatType sr, DspFloatType q = 0.5, DspFloatType g = 1) 
        : ZolzerBiquadFilter(ALLPASS,fc,sr,q,g)
        {

        }
    };

//////////////////////////////////////////////////////////////////////////////////////////
// Lowpass 1p
//////////////////////////////////////////////////////////////////////////////////////////

    struct ZolzerLowPass1pFilter : public ZolzerBiquadFilter
    {

        
        DspFloatType Fc, Fs, Q, G, R;

        ZolzerLowPass1pFilter(DspFloatType sampleRate = 44100.0) : ZolzerBiquadFilter(LOWPASS1P,sampleRate)
        {
        
        }

        ZolzerLowPass1pFilter(DspFloatType fc, DspFloatType sr, DspFloatType q = 0.5, DspFloatType g = 1) 
        : ZolzerBiquadFilter(LOWPASS1P,fc,sr,q,g)
        {

        }
    };

//////////////////////////////////////////////////////////////////////////////////////////
// Highpass 1p
//////////////////////////////////////////////////////////////////////////////////////////

    struct ZolzerHighPass1pFilter : public ZolzerBiquadFilter
    {        
        DspFloatType Fc, Fs, Q, G, R;

        ZolzerHighPass1pFilter(DspFloatType sampleRate = 44100.0) : ZolzerBiquadFilter(HIGHPASS1P,sampleRate)
        {
        
        }

        ZolzerHighPass1pFilter(DspFloatType fc, DspFloatType sr, DspFloatType q = 0.5, DspFloatType g = 1) 
        : ZolzerBiquadFilter(HIGHPASS1P,fc,sr,q,g)
        {

        }
    };

//////////////////////////////////////////////////////////////////////////////////////////
// Allpass 1p
//////////////////////////////////////////////////////////////////////////////////////////

    struct ZolzerAllPass1pFilter : public ZolzerBiquadFilter
    {        
        DspFloatType Fc, Fs, Q, G, R;

        ZolzerAllPass1pFilter(DspFloatType sampleRate = 44100.0) : ZolzerBiquadFilter(ALLPASS1P,sampleRate)
        {
        
        }

        ZolzerAllPass1pFilter(DspFloatType fc, DspFloatType sr, DspFloatType q = 0.5, DspFloatType g = 1) 
        : ZolzerBiquadFilter(ALLPASS1P,fc,sr,q,g)
        {

        }
    };

//////////////////////////////////////////////////////////////////////////////////////////
// Bandpass
//////////////////////////////////////////////////////////////////////////////////////////

    struct ZolzerBandPassFilter : public ZolzerBiquadFilter
    {

        
        DspFloatType Fc, Fs, Q, G, R;

        ZolzerBandPassFilter(DspFloatType sampleRate = 44100.0) : ZolzerBiquadFilter(BANDPASS,sampleRate)
        {
        
        }

        ZolzerBandPassFilter(DspFloatType fc, DspFloatType sr, DspFloatType q = 0.5, DspFloatType g = 1) 
        : ZolzerBiquadFilter(BANDPASS,fc,sr,q,g)
        {

        }
    };

//////////////////////////////////////////////////////////////////////////////////////////
// Notch
//////////////////////////////////////////////////////////////////////////////////////////

    struct ZolzerNotchFilter : public ZolzerBiquadFilter
    {

        
        DspFloatType Fc, Fs, Q, G, R;

        ZolzerNotchFilter(DspFloatType sampleRate = 44100.0) : ZolzerBiquadFilter(NOTCH,sampleRate)
        {
        
        }

        ZolzerNotchFilter(DspFloatType fc, DspFloatType sr, DspFloatType q = 0.5, DspFloatType g = 1) 
        : ZolzerBiquadFilter(NOTCH,fc,sr,q,g)
        {

        }
    };


//////////////////////////////////////////////////////////////////////////////////////////
// peak boost
//////////////////////////////////////////////////////////////////////////////////////////

    struct ZolzerPeakBoostFilter : public ZolzerBiquadFilter
    {

        
        DspFloatType Fc, Fs, Q, G, R;

        ZolzerPeakBoostFilter(DspFloatType sampleRate = 44100.0) : ZolzerBiquadFilter(PEAKBOOST,sampleRate)
        {
        
        }

        ZolzerPeakBoostFilter(DspFloatType fc, DspFloatType sr, DspFloatType q = 0.5, DspFloatType g = 1) 
        : ZolzerBiquadFilter(PEAKBOOST,fc,sr,q,g)
        {

        }
    };

//////////////////////////////////////////////////////////////////////////////////////////
// peak cut
//////////////////////////////////////////////////////////////////////////////////////////

    struct ZolzerPeakCutFilter : public ZolzerBiquadFilter
    {

        
        DspFloatType Fc, Fs, Q, G, R;

        ZolzerPeakCutFilter(DspFloatType sampleRate = 44100.0) : ZolzerBiquadFilter(PEAKCUT,sampleRate)
        {
        
        }

        ZolzerPeakCutFilter(DspFloatType fc, DspFloatType sr, DspFloatType q = 0.5, DspFloatType g = 1) 
        : ZolzerBiquadFilter(PEAKCUT,fc,sr,q,g)
        {

        }
    };

//////////////////////////////////////////////////////////////////////////////////////////
// low shelf boost
//////////////////////////////////////////////////////////////////////////////////////////

    struct ZolzerLowShelfBoostFilter : public ZolzerBiquadFilter
    {

        
        DspFloatType Fc, Fs, Q, G, R;

        ZolzerLowShelfBoostFilter(DspFloatType sampleRate = 44100.0) : ZolzerBiquadFilter(LOWSHELFBOOST,sampleRate)
        {
        
        }

        ZolzerLowShelfBoostFilter(DspFloatType fc, DspFloatType sr, DspFloatType q = 0.5, DspFloatType g = 1) 
        : ZolzerBiquadFilter(LOWSHELFBOOST,fc,sr,q,g)
        {

        }
    };

//////////////////////////////////////////////////////////////////////////////////////////
// low shelf cut
//////////////////////////////////////////////////////////////////////////////////////////

    struct ZolzerLowShelfCutFilter : public ZolzerBiquadFilter
    {        
        DspFloatType Fc, Fs, Q, G, R;

        ZolzerLowShelfCutFilter(DspFloatType sampleRate = 44100.0) : ZolzerBiquadFilter(LOWSHELFCUT,sampleRate)
        {
        
        }

        ZolzerLowShelfCutFilter(DspFloatType fc, DspFloatType sr, DspFloatType q = 0.5, DspFloatType g = 1) 
        : ZolzerBiquadFilter(LOWSHELFCUT,fc,sr,q,g)
        {

        }
    };

//////////////////////////////////////////////////////////////////////////////////////////
// highshelf boost
//////////////////////////////////////////////////////////////////////////////////////////

    struct ZolzerHighShelfBoostFilter : public ZolzerBiquadFilter
    {

        
        DspFloatType Fc, Fs, Q, G, R;

        ZolzerHighShelfBoostFilter(DspFloatType sampleRate = 44100.0) : ZolzerBiquadFilter(HIGHSHELFBOOST,sampleRate)
        {
        
        }

        ZolzerHighShelfBoostFilter(DspFloatType fc, DspFloatType sr, DspFloatType q = 0.5, DspFloatType g = 1) 
        : ZolzerBiquadFilter(HIGHSHELFBOOST,fc,sr,q,g)
        {

        }
    };

//////////////////////////////////////////////////////////////////////////////////////////
// highshelf cut
//////////////////////////////////////////////////////////////////////////////////////////

    struct ZolzerHighShelfCutFilter : public ZolzerBiquadFilter
    {        
        DspFloatType Fc, Fs, Q, G, R;

        ZolzerHighShelfCutFilter(DspFloatType sampleRate = 44100.0) : ZolzerBiquadFilter(HIGHSHELFCUT,sampleRate)
        {
        
        }

        ZolzerHighShelfCutFilter(DspFloatType fc, DspFloatType sr, DspFloatType q = 0.5, DspFloatType g = 1) 
        : ZolzerBiquadFilter(HIGHSHELFCUT,fc,sr,q,g)
        {

        }
    };

}
