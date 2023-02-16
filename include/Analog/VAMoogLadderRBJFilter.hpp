#pragma once

#include "VAMoogLadderFilters.hpp"

namespace Analog::Moog
{
    class RBJFilter : public BiQuadBase
    {
    public:
        
        enum FilterType
        {
            LOWPASS,
            HIGHPASS,
            BANDPASS,
            ALLPASS,
            NOTCH,
            PEAK,
            LOW_SHELF,
            HIGH_SHELF
        };
        
        
        RBJFilter(FilterType type = FilterType::LOWPASS, DspFloatType cutoff = 1, DspFloatType sampleRate = 44100) : sampleRate(sampleRate), t(type)
        {
            Q = 1;
            A = 1;

            a = {{0.0f, 0.0f, 0.0f}};
            b = {{0.0f, 0.0f, 0.0f}};

            SetCutoff(cutoff);
        }
        
        ~RBJFilter()
        {
            
        }
        
        void UpdateCoefficients()
        {
            cosOmega = std::cos(omega);
            sinOmega = std::sin(omega);
            
            switch (t)
            {
                case LOWPASS:
                {
                    alpha = sinOmega / (2.0 * Q);
                    b[0] = (1 - cosOmega) / 2;
                    b[1] = 1 - cosOmega;
                    b[2] = b[0];
                    a[0] = 1 + alpha;
                    a[1] = -2 * cosOmega;
                    a[2] = 1 - alpha;
                } break;
                    
                case HIGHPASS:
                {
                    alpha = sinOmega / (2.0 * Q);
                    b[0] = (1 + cosOmega) / 2;
                    b[1] = -(1 + cosOmega);
                    b[2] = b[0];
                    a[0] = 1 + alpha;
                    a[1] = -2 * cosOmega;
                    a[2] = 1 - alpha;
                } break;
                    
                case BANDPASS:
                {
                    alpha = sinOmega * sinhf(logf(2.0) / 2.0 * Q * omega/sinOmega);
                    b[0] = sinOmega / 2;
                    b[1] = 0;
                    b[2] = -b[0];
                    a[0] = 1 + alpha;
                    a[1] = -2 * cosOmega;
                    a[2] = 1 - alpha;
                } break;
                    
                case ALLPASS:
                {
                    alpha = sinOmega / (2.0 * Q);
                    b[0] = 1 - alpha;
                    b[1] = -2 * cosOmega;
                    b[2] = 1 + alpha;
                    a[0] = b[2];
                    a[1] = b[1];
                    a[2] = b[0];
                } break;
                    
                case NOTCH:
                {
                    alpha = sinOmega * sinhf(logf(2.0) / 2.0 * Q * omega/sinOmega);
                    b[0] = 1;
                    b[1] = -2 * cosOmega;
                    b[2] = 1;
                    a[0] = 1 + alpha;
                    a[1] = b[1];
                    a[2] = 1 - alpha;
                } break;
                    
                case PEAK:
                {
                    alpha = sinOmega * sinhf(logf(2.0) / 2.0 * Q * omega/sinOmega);
                    b[0] = 1 + (alpha * A);
                    b[1] = -2 * cosOmega;
                    b[2] = 1 - (alpha * A);
                    a[0] = 1 + (alpha / A);
                    a[1] = b[1];
                    a[2] = 1 - (alpha / A);
                } break;
                    
                case LOW_SHELF:
                {
                    alpha = sinOmega / 2.0 * sqrt( (A + 1.0 / A) * (1.0 / Q - 1.0) + 2.0);
                    b[0] = A * ((A + 1) - ((A - 1) * cosOmega) + (2 * sqrtf(A) * alpha));
                    b[1] = 2 * A * ((A - 1) - ((A + 1) * cosOmega));
                    b[2] = A * ((A + 1) - ((A - 1) * cosOmega) - (2 * sqrtf(A) * alpha));
                    a[0] = ((A + 1) + ((A - 1) * cosOmega) + (2 * sqrtf(A) * alpha));
                    a[1] = -2 * ((A - 1) + ((A + 1) * cosOmega));
                    a[2] = ((A + 1) + ((A - 1) * cosOmega) - (2 * sqrtf(A) * alpha));
                } break;
                    
                case HIGH_SHELF:
                {
                    alpha = sinOmega / 2.0 * sqrt( (A + 1.0 / A) * (1.0 / Q - 1.0) + 2.0);
                    b[0] = A * ((A + 1) + ((A - 1) * cosOmega) + (2 * sqrtf(A) * alpha));
                    b[1] = -2 * A * ((A - 1) + ((A + 1) * cosOmega));
                    b[2] = A * ((A + 1) + ((A - 1) * cosOmega) - (2 * sqrtf(A) * alpha));
                    a[0] = ((A + 1) - ((A - 1) * cosOmega) + (2 * sqrtf(A) * alpha));
                    a[1] = 2 * ((A - 1) - ((A + 1) * cosOmega));
                    a[2] = ((A + 1) - ((A - 1) * cosOmega) - (2 * sqrtf(A) * alpha));
                } break;
            }
            
            // Normalize filter coefficients
            DspFloatType factor = 1.0f / a[0];
            
            std::array<DspFloatType, 2> aNorm;
            std::array<DspFloatType, 3> bNorm;

            aNorm[0] = a[1] * factor;
            aNorm[1] = a[2] * factor;

            bNorm[0] = b[0] * factor;
            bNorm[1] = b[1] * factor;
            bNorm[2] = b[2] * factor;
            
            SetBiquadCoefs(bNorm, aNorm);
        }
        
        void SetSampleRate(DspFloatType sr)
        {
            sampleRate = sr;
        }

        // In Hertz, 0 to Nyquist
        void SetCutoff(DspFloatType c)
        {
            omega = HZ_TO_RAD(c) / sampleRate;
            UpdateCoefficients();
        }
        
        DspFloatType GetCutoff()
        {
            return omega;
        }
        
        // Arbitrary, from 0.01f to ~20
        void SetQValue(DspFloatType q)
        {
            Q = q;
            UpdateCoefficients();
        }
        
        DspFloatType GetQValue()
        {
            return Q;
        }

        void SetType(FilterType newType)
        {
            t = newType;
            UpdateCoefficients();
        }
        
        FilterType GetType()
        {
            return t;
        }
        
    private:

        DspFloatType sampleRate;
        
        DspFloatType omega;
        DspFloatType cosOmega;
        DspFloatType sinOmega;
        
        DspFloatType Q;
        DspFloatType alpha;
        DspFloatType A;

        std::array<DspFloatType, 3> a;
        std::array<DspFloatType, 3> b;
        
        FilterType t;
    };


    struct Filter : public FilterProcessor
    {
        enum FilterType {
            LOWPASS = RBJFilter::FilterType::LOWPASS,
            HIGHPASS = RBJFilter::FilterType::HIGHPASS,
            BANDPASS = RBJFilter::FilterType::BANDPASS,
            ALLPASS = RBJFilter::FilterType::ALLPASS,
            NOTCH = RBJFilter::FilterType::NOTCH,
            PEAK = RBJFilter::FilterType::PEAK,
            LOW_SHELF = RBJFilter::FilterType::LOW_SHELF,
            HIGH_SHELF = RBJFilter::FilterType::HIGH_SHELF,
        };
        
        RBJFilter * filter;
        Filter(FilterType type, DspFloatType cutoff, DspFloatType sample_rate) : FilterProcessor() {
            filter = new RBJFilter((RBJFilter::FilterType)type,cutoff,sample_rate);
            assert(filter != nullptr);
        }
        ~Filter() {
            if(filter) delete filter;
        }

        void UpdateCoefficients() { filter->UpdateCoefficients(); }
        void SetCutoff(DspFloatType c) { filter->SetCutoff(c); }
        void SetQ(DspFloatType q) { filter->SetQValue(q); }
        DspFloatType GetCutoff() { return filter->GetCutoff(); }
        DspFloatType GetQ() { return filter->GetQValue(); }

        enum {
            PORT_CUTOFF,
            PORT_RESONANCE,
        };
        void setPort(int port, DspFloatType v) {
            switch(port)
            {
                case PORT_CUTOFF: SetCutoff(v); break;
                case PORT_RESONANCE: SetQ(v); break;
            }
        }
        void ProcessBlock(size_t n, DspFloatType * in, DspFloatType * out) { filter->ProcessBlock(n,in,out); }
        void ProcessInplace(size_t n, DspFloatType * p) { filter->ProcessInplace(n,p); }
        
        DspFloatType Tick(DspFloatType s, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) { return filter->Tick(s,A,X,Y); }
    };

    struct RBJLowPassFilter : public Filter
    {
        RBJLowPassFilter(DspFloatType cutoff, DspFloatType sample_rate=44100) : Filter(LOWPASS,cutoff,sample_rate) {

        }
    };

    struct RBJHighPassFilter : public Filter
    {
        RBJHighPassFilter(DspFloatType cutoff, DspFloatType sample_rate=44100) : Filter(HIGHPASS,cutoff,sample_rate) {

        }
    };

    struct RBJBandPassFilter : public Filter
    {
        RBJBandPassFilter(DspFloatType cutoff, DspFloatType sample_rate=44100) : Filter(BANDPASS,cutoff,sample_rate) {

        }
    };
    struct RBJAllPassFilter : public Filter
    {
        RBJAllPassFilter(DspFloatType cutoff, DspFloatType sample_rate=44100) : Filter(ALLPASS,cutoff,sample_rate) {

        }
    };

    struct RBJNotchFilter : public Filter
    {
        RBJNotchFilter(DspFloatType cutoff, DspFloatType sample_rate=44100) : Filter(NOTCH,cutoff,sample_rate) {

        }
    };

    struct RBJPeakFilter : public Filter
    {
        RBJPeakFilter(DspFloatType cutoff, DspFloatType sample_rate=44100) : Filter(PEAK,cutoff,sample_rate) {

        }
    };

    struct RBJLowShelfFilter : public Filter
    {
        RBJLowShelfFilter(DspFloatType cutoff, DspFloatType sample_rate=44100) : Filter(LOW_SHELF,cutoff,sample_rate) {

        }
    };

    struct RBJHighShelfFilter : public Filter
    {
        RBJHighShelfFilter(DspFloatType cutoff, DspFloatType sample_rate=44100) : Filter(HIGH_SHELF,cutoff,sample_rate) {

        }
    };

}
