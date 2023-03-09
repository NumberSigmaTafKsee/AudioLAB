#pragma once

#include <stdio.h>
#include <math.h>
#include "GenericSoundObject.hpp"


namespace Analog::MoogFilters::MoogFilter4 
{
	
	template<typename DSP>
    class MoogFilter : public GSSoundProcessor<DSP> {
        DSP frequency;
        DSP g;
        DSP resonance;
        DSP drive;
        int sampleRate;
        
        DSP y_a;
        DSP y_b;
        DSP y_c;
        DSP y_d;
        DSP y_d_1;
        
    public:
        MoogFilter();
        ~MoogFilter();
        
        
        DSP getFrequency();
        DSP getResonance();
        DSP getDrive();
        

        void setFrequency(DSP f);
        void setResonance(DSP r);
        void setSampleRate(int s);
        void setDrive(DSP d);

        enum
        {
            PORT_CUTOFF,
            PORT_RESONANCE,			
            PORT_DRIVE,
        };
        void setPort(int port, DSP v)
        {
            switch (port)
            {
            case PORT_CUTOFF:
                setFrequency(v);
                break;
            case PORT_RESONANCE:
                setResonance(v);
                break;
            case PORT_DRIVE:
                setDrive(v);
                break;
            }
        }
        DSP Tick(DSP I, DSP A=1, DSP X=1, DSP Y=1) {
            DSP samples = tanhf(I * drive);
            y_a = y_a + g * (std::tanh(samples - resonance * ((y_d_1 + y_d)/2) - std::tanh(y_a)));
            y_b = y_b + g * (std::tanh(y_a) - std::tanh(y_b));
            y_c = y_c + g * (std::tanh(y_b) - std::tanh(y_c));
            y_d_1 = y_d;
            y_d = y_d + g * (std::tanh(y_c) - std::tanh(y_d));
            return A*y_d;
        }
        void ProcessSIMD(size_t n, DSP *in, DSP * out);
        
        void ProcessBlock(size_t n, DSP *in, DSP * out) {
			ProcessSIMD(n,in,out);
		}
		void ProcessInplace(size_t n, DSP *in) {
			ProcessSIMD(n,in,in);
		}        
    };

    MoogFilter::MoogFilter() : GSSoundProcessor<DSP> {
        y_a = 0;
        y_b = 0;
        y_c = 0;
        y_d = 0;
        y_d_1 = 0;
        
        frequency = 2000;
        resonance = 1;
        drive = 1;
    }
    MoogFilter::~MoogFilter() {
        
    }
    void MoogFilter::ProcessSIMD(size_t n, DSP * in, DSP * out) {
		Undenormal denormals;
		#pragma omp simd aligned(in,out)
        for (int i = 0; i < 2 * numSamples; i++) {
            samples[i/2] = tanhf(samples[i/2] * drive);
            y_a = y_a + g * (tanhf(samples[i/2] - resonance * ((y_d_1 + y_d)/2) - tanhf(y_a)));
            y_b = y_b + g * (tanhf(y_a) - tanhf(y_b));
            y_c = y_c + g * (tanhf(y_b) - tanhf(y_c));
            y_d_1 = y_d;
            y_d = y_d + g * (tanhf(y_c) - tanhf(y_d));
            samples[i/2] = y_d;
        }
    }
    
    DSP MoogFilter::getFrequency() {
        return frequency;
    }
    DSP MoogFilter::getResonance() {
        return resonance;
    }
    DSP MoogFilter::getDrive() {
        return drive;
    }

    void MoogFilter::setFrequency(DSP f) {
        if (f > 12000.0f) f = 12000.0f;
        if (f < 0.0f) f = 0.0f;
        frequency = f;
        g = 1 - expf(-2 * tanf(2 * M_PI * frequency/(2 * sampleRate)));
    }
    void MoogFilter::setResonance(DSP r) {
        if (r > 5.0f) r = 5.0f;
        if (r < 0.0f) r = 0.0f;
        resonance = r;
    }
    void MoogFilter::setSampleRate(int s) {
        sampleRate = s;
    }
    void MoogFilter::setDrive(DSP d) {
        if (d > 10.0f) d = 10.0f;
        if (d < 0.1f) d = 0.1f;
        drive = d;
    }
}
