
#pragma once

#include "Filters/IIRFilters.hpp"

namespace Filters::IIR::RBJFilters
{
//////////////////////////////////////////////////////////////////////////////////////////
// RBJ
//////////////////////////////////////////////////////////////////////////////////////////
  
    struct RBJBiquadFilter : public BiquadTransposedTypeII
    {

		enum RBJType {
			LOWPASS,
			HIGHPASS,
			BANDPASS,
			ZERODB,
			SKIRT,
			NOTCH,
			PEAK,
			LOWSHELF,
			HIGHSHELF,
			ALLPASS,
			LPBW,
			HPBW,
			LPR,
			HPR,
			SKIRTBW,
			ZERODBBW,
			NOTCHBW,			
			APFBW,
			PEAKBW,
			LOWSHELFSLOPE,
			HIGHSHELFSLOPE,			
		};
		
        RBJType filter_type;
        DspFloatType Fc, Fs, Q, G, R,S,BW;
		
        RBJBiquadFilter(RBJType type = LOWPASS, DspFloatType sampleRate = 44100.0) : BiquadTransposedTypeII()
        {
            Fc = 1000.0;
            Fs = sampleRate;
            Q = 0.5;
            G = 1.0;
            S = 0.0;
            BW=0.0;
            filter_type = type;
            setCoefficients(Fc);
        }        
        RBJBiquadFilter(RBJType type, DspFloatType fc, DspFloatType sr, DspFloatType q = 0.5, DspFloatType g = 1) : BiquadTransposedTypeII()
        {
            filter_type = type;
            Fs = sr/2.0;
            Q = q;
            G = g;
            setCoefficients(Fc);
        }
        void setCoefficients(DspFloatType fc)
        {            
            FilterCoefficients c;            
            Fc = fc;
            
            switch (filter_type)
            {
            case LOWPASS:
                c = RBJLowpassBiquad(Fc, Fs, Q);
                break;
            case HIGHPASS:
                c = RBJHighpassBiquad(Fc, Fs, Q);
                break;
            case BANDPASS:
				c = RBJBandpassConstant0dbBiquadBW(Fc, Fs, Q);
				break;
            case ZERODB:
                c = RBJBandpassConstant0dbBiquad(Fc, Fs, Q);
                break;
            case SKIRT:            
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
            case LPBW:
				c = RBJLowpassBiquadBW(Fc,Fs,BW);
				break;
			case HPBW:
				c = RBJHighpassBiquadBW(Fc,Fs,BW);
				break;
			case LPR:
				c = RBJLowpassBiquadR(Fc,Fs,Q,R);
				break;
			case HPR:
				c = RBJHighpassBiquadR(Fc,Fs,Q,R);
				break;
			case SKIRTBW:
				c = RBJBandpassConstantSkirtBiquadBW(Fc,Fs,BW);
				break;
			case ZERODBBW:
				c = RBJBandpassConstant0dbBiquadBW(Fc,Fs,BW);
				break;
			case NOTCHBW:
				c = RBJNotchBiquadBW(Fc,Fs,BW);
				break;
			case APFBW:
				c = RBJAllpassBiquadBW(Fc,Fs,BW);
				break;
			case PEAKBW:
				c = RBJPeakBiquadBW(Fc,Fs,BW,G);
				break;
			case LOWSHELFSLOPE:
				c = RBJLowshelfBiquadSlope(Fc,Fs,S,G);
				break;
			case HIGHSHELFSLOPE:
				c = RBJHighshelfBiquadSlope(Fc,Fs,S,G);
				break;
            }
            biquad.setCoefficients(c);
        }
        void setType(RBJType type) {
			filter_type = type;
		}
        void setCutoff(DspFloatType fc)
        {
            if(fc < 0 || fc >= Fs/2.0) return;
            Fc = fc;
            setCoefficients(fc);
        }
        void setQ(DspFloatType q)
        {
			Q = q;            
            setCoefficients(Fc);
        }        
        void setGain(DspFloatType g)
        {            
            G = g;
            setCoefficients(Fc);
        }
        void setBandWidth(DspFloatType bw) {
			BW = bw;
			setCoefficients(Fc);
		}
        void setSlope(DspFloatType s) {
			S = s;
			setCoefficients(Fc);
		}
		
        float getCutoff() const { return Fc; }
        float getQ() const { return Q; }
        float getGain() const { return G; }
        float getSlope() const { return S; }
        
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

        RBJSkirtBandPassFilter(DspFloatType sampleRate = 44100.0) : RBJBiquadFilter(SKIRT,sampleRate)
        {
        }
        RBJSkirtBandPassFilter(DspFloatType fc, DspFloatType sr, DspFloatType q = 0.5, DspFloatType g = 1) 
        : RBJBiquadFilter(SKIRT,fc,sr,q,g)
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
