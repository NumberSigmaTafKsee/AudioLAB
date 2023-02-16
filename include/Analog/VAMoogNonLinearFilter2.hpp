#pragma once

namespace Analog::Filters::Moog::NonLinear2
{
    class NonLinearMoogFilter2 {
        DspFloatType frequency;
        DspFloatType g;
        DspFloatType resonance;
        DspFloatType drive;
        int sampleRate;
        
        DspFloatType y_a;
        DspFloatType y_b;
        DspFloatType y_c;
        DspFloatType y_d;
        DspFloatType y_d_1;
        
    public:
        NonLinearMoogFilter2();
        ~NonLinearMoogFilter2();
        
        
        DspFloatType getFrequency();
        DspFloatType getResonance();
        DspFloatType getDrive();
        
        void setFrequency(DspFloatType f);
        void setResonance(DspFloatType r);
        void setSampleRate(int s);
        void setDrive(DspFloatType d);

        enum {
            PORT_CUTOFF,
            PORT_RESONANCE,
            PORT_DRIVE,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_CUTOFF: setFrequency(v); break;
                case PORT_RESONANCE: setResonance(v); break;
                case PORT_DRIVE: setDrive(v); break;
            }
        }

		void ProcessSIMD(size_t n, DspFloatType *input, DspFloatType * output);
		void ProcessBlock(size_t n, DspFloatType * in, DspFloatType * out) {
				ProcessSIMD(n,in,out);
		}
		void ProcessInplace(size_t n, DspFloatType * out) {
			ProcessSIMD(n,out,out);
		}
    };

    NonLinearMoogFilter2::NonLinearMoogFilter2() {
        y_a = 0;
        y_b = 0;
        y_c = 0;
        y_d = 0;
        y_d_1 = 0;
        
        frequency = 2000;
        resonance = 1;
        drive = 1;
    }
    NonLinearMoogFilter2::~NonLinearMoogFilter2() {
        
    }
    void NonLinearMoogFilter2::ProcessSIMD(size_t n, DspFloatType *input, DspFloatType * samples) {
        #pragma omp simd
        for (int i = 0; i < 2 * n; i++) {
            samples[i/2] = std::tanh(input[i/2] * drive);
            y_a = y_a + g * (std::tanh(samples[i/2] - resonance * ((y_d_1 + y_d)/2) - std::tanh(y_a)));
            y_b = y_b + g * (std::tanh(y_a) - std::tanh(y_b));
            y_c = y_c + g * (std::tanh(y_b) - std::tanh(y_c));
            y_d_1 = y_d;
            y_d = y_d + g * (std::tanh(y_c) - std::tanh(y_d));
            samples[i/2] = y_d;
        }
    }
    
    DspFloatType NonLinearMoogFilter2::getFrequency() {
        return frequency;
    }
    DspFloatType NonLinearMoogFilter2::getResonance() {
        return resonance;
    }
    DspFloatType NonLinearMoogFilter2::getDrive() {
        return drive;
    }

    void NonLinearMoogFilter2::setFrequency(DspFloatType f) {
        if (f > 12000.0f) f = 12000.0f;
        if (f < 0.0f) f = 0.0f;
        frequency = f;
        g = 1 - expf(-2 * tanf(2 * M_PI * frequency/(2 * sampleRate)));
    }
    void NonLinearMoogFilter2::setResonance(DspFloatType r) {
        if (r > 5.0f) r = 5.0f;
        if (r < 0.0f) r = 0.0f;
        resonance = r;
    }
    void NonLinearMoogFilter2::setSampleRate(int s) {
        sampleRate = s;
    }
    void NonLinearMoogFilter2::setDrive(DspFloatType d) {
        if (d > 10.0f) d = 10.0f;
        if (d < 0.1f) d = 0.1f;
        drive = d;
    }
    
};
