#pragma once

namespace Filters::Agave
{
    // todo: add ports
    class RCFilter : public FilterProcessor
    {
    // THIS CLASS IMPLEMENTES AN LTI FIRST-ORDER LOWPASS FILTER DERIVED FROM THE TRANSFER FUNCTION 
    // OF A SIMPLE PASSIVE RC FILTER. THE FILTER IS PARAMETRIZED BY SETTING THE CUTOFF FREQUENCY IN HZ
    // I.E. fc = 1 / (2*pi*R*C)
    // 
    // Usage example:
    // 	RCFilter filter(100.0f,44100.0f);
    // 	filter.process(x);
    // 	filter.getLowpassOutput();
    // 
    private:

        float sampleRate  = 44.1e3f;
        float fc = 1.0e3f; 	// Cutoff frequency (in Hz)
        float wc;			// Cutoff frequency (in rad/sec)

        float previousInput = 0.0f;
        float lowpassOutput = 0.0f;
        float highpassOutput = 0.0f;

    public:

        RCFilter() : FilterProcessor() { setCutoff(); }
        RCFilter(float cutoffFrequency, float SR) : FilterProcessor() { 
            fc = cutoffFrequency;
            sampleRate = SR;
            setCutoff();
        }
        ~RCFilter() {}

        void setSampleRate(float SR) {
            sampleRate = SR;
            setCutoff();
        }

        enum {
            PORT_CUTOFF,
            PORT_LOWPASS,
            PORT_HIGHPASS,
        };
        void setPort(int port, double v) {
            switch(port)
            {
                case PORT_CUTOFF: fc = v; setCutoff(); break;
            }
        }

        void setCutoff() {
            float wa = 2.0f*M_PI*fc; // analog cutoff freq
            wc = 2.0f*std::atan(0.5f*wa/sampleRate)*sampleRate;	// digital cutoff freq
        }

        void process(float input) {

            float alpha = 2.0f*sampleRate/wc;

            // Compute filter output
            lowpassOutput = ( (alpha - 1.0f)*lowpassOutput + input + previousInput ) / (1.0f + alpha);
            highpassOutput = input - lowpassOutput;

            // Update State
            previousInput = input;

        }
        double Tick(double I, double A=1, double X=0, double Y=0)
        {
            return A*process(I);
        }
        float getLowpassOutput() {
            return lowpassOutput;
        }

        float getHighpassOutput() {
            return highpassOutput;
        }
        double getPort(int port) {
            switch(port) {
                case PORT_LOWPASS: return getLowpassOutput();
                case PORT_HIGHPASS: return getHighpassOutput();
            }
            return 0;
        }
    };

    // todo: add ports
    class DCBlocker : public FilterProcessor {

    // THIS CLASS IMPLEMENTES AN LTI IIR DC BLOCKER BASED ON J. PEKONEN'S DESIGN, DESCRIBED IN
    // "FILTER-BASED ALIAS REDUCTION FOR DIGITAL CLASSICAL WAVEFORM SYNTHESIS" (ICASSP 2008)
    // 
    private: 

        // Default parameters. Use constructor to overwrite.
        float sampleRate = 44.1e3f;
        float fc = 1.0e3f;

        float xState = 0.0f;
        float yState = 0.0f;
        float p = 0.0f;
        float output = 0.0f;

    public:

        DCBlocker() : FilterProcessor() { setPole(); }
        DCBlocker(float cutoffFrequency, float SR) : FilterProcessor() { 
            fc = cutoffFrequency;
            sampleRate = SR;
            setPole();
        }
        ~DCBlocker() {}

        void setSampleRate(float SR) {
            sampleRate = SR;
            setPole();
        }

        enum {
            PORT_CUTOFF
        };
        void setPort(int port, double v) {
            switch(port)
            {
                case PORT_CUTOFF: fc = v; setPole(); break;
            }
        }
        void setPole() {
            p = std::tan(0.25f*M_PI - M_PI*fc/sampleRate); // Filter pole
        }

        void process(float input) {

            output = 0.5f*(1.0 + p) * ( input - xState + p*yState ); 

            // Update State
            xState = input;
            yState = output;

        }
        double Tick(double I, double A=1, double X=0, double Y=0) {
            return A*process(I);
        }
        float getFilteredOutput() {
            return output;
        }
    };

}