#pragma once

#include <cmath>
#include <cstdlib>

namespace Analog::Filters::DinkFilter
{
    class DinkyFilter : public FilterProcessor {
    public:
        enum FilterMode {
            FILTER_MODE_LOWPASS = 0,
            FILTER_MODE_HIGHPASS,
            FILTER_MODE_BANDPASS,
            kNumFilterModes
        };
        DinkyFilter() :
        cutoff(0.99),
        resonance(0.0),
        mode(FILTER_MODE_LOWPASS),
        buf0(0.0),
        buf1(0.0),
        buf2(0.0),
        buf3(0.0)
        {
            calculateFeedbackAmount();
        };

        DspFloatType process(DspFloatType inputValue);
        void setCutoff(DspFloatType newCutoff) { cutoff = newCutoff; calculateFeedbackAmount(); };
        void setResonance(DspFloatType newResonance) { resonance = newResonance; calculateFeedbackAmount(); };
        void setFilterMode(FilterMode newMode) { mode = newMode; }

        enum {
          PORT_CUTOFF,
          PORT_RESONANCE,
          PORT_FILTERMODE,          
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_CUTOFF: setCutoff(v); break;
                case PORT_RESONANCE: setResonance(v); break;
                case PORT_FILTERMODE: setFilterMode((FilterMode)v); break;  
                default: printf("No port %d\n", port);              
            }
        }
        DspFloatType Tick(DspFloatType in, DspFloatType A=1, DspFloatType X=0, DspFloatType Y=0) { 
            return A*process(in); 
        }
        
        void ProcessSIMD(size_t n, DspFloatType * in, DspFloatType * out) {
            #pragma omp simd
            for(size_t i = 0; i < n; i++)
            {
                DspFloatType inputValue = in[i];
                buf0 += cutoff * (inputValue - buf0 + (buf0-buf1)*feedbackAmount);
                buf1 += cutoff * (buf0 - buf1);
                buf2 += cutoff * (buf1 - buf2);
                buf3 += cutoff * (buf2 - buf3);
                switch (mode) {
                    case FILTER_MODE_LOWPASS:
                        out[i] =  buf3;
                        break;
                    case FILTER_MODE_HIGHPASS:
                        out[i] =  inputValue - buf0;
                        break;
                    case FILTER_MODE_BANDPASS:
                        out[i] =  buf0 - buf3;
                        break;
                    default:
                        out[i] = buf3;
                }
            }
        }
        void ProcessBlock(size_t n, DspFloatType * input, DspFloatType * output) {
            ProcessSIMD(n,input,output);
        }
            
        void ProcessInplace(size_t n, DspFloatType * input) {
            ProcessBlock(n,input,input);
        }

    private:
        DspFloatType cutoff;
        DspFloatType resonance;
        FilterMode mode;
        DspFloatType feedbackAmount;
        void calculateFeedbackAmount() { feedbackAmount = resonance + resonance/(1.0 - cutoff); }
        DspFloatType buf0;
        DspFloatType buf1;
        DspFloatType buf2;
        DspFloatType buf3;
    };

    inline DspFloatType DinkyFilter::process(DspFloatType inputValue) {
        buf0 += cutoff * (inputValue - buf0 + (buf0-buf1)*feedbackAmount);
        buf1 += cutoff * (buf0 - buf1);
        buf2 += cutoff * (buf1 - buf2);
        buf3 += cutoff * (buf2 - buf3);
        switch (mode) {
            case FILTER_MODE_LOWPASS:
                return buf3;
            case FILTER_MODE_HIGHPASS:
                return inputValue - buf0;
            case FILTER_MODE_BANDPASS:
                return buf0 - buf3;
            default:
                return 0.0;
        }
    }
}
