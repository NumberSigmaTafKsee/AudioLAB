class AnalogFilter {
public:
    AnalogFilter() {}
    
    AnalogFilter(const unsigned int p) : poles(p) {}
    
    AnalogFilter(const unsigned int p, const double cf) : poles(p), cornerFrequency(cf) {}

    virtual ~AnalogFilter() {}

    void bindSampleRate (double* sampleRate) { this->sampleRate = sampleRate; }
    
    virtual void setPoles (unsigned int p) noexcept = 0;
    
    virtual void setCornerFrequency (double frequency) noexcept = 0;
    
    virtual void process (double& buffer, unsigned int channel) noexcept = 0;
    
    template <typename SampleType>
    SampleType getProcessed (SampleType sample, unsigned int channel) noexcept;

protected:
    double* sampleRate = nullptr;
    double tau = 0.0;
    unsigned int poles = 1.0;
    double cornerFrequency = 1.0;
    double caps[32][32] = {0.0};
    double* cap_ptr;
    double sampleDelta = 0.0;
    double pi = 3.1415926535;
};

emplate <typename SampleType>
inline SampleType AnalogFilter::getProcessed (SampleType sample, unsigned int channel) noexcept
{
    double sample2 = sample;
    process (sample2, channel);
    return sample2;
}



class LowPassFilter : public AnalogFilter {
    public:
        LowPassFilter();
        
        LowPassFilter(const unsigned int poles);
    
        LowPassFilter(const unsigned int poles, const double cornerFrequency);
    
        void setPoles (unsigned int p) noexcept override;
    
        void setCornerFrequency(double frequency) noexcept override;
    
        void process (double& sample, unsigned int channel) noexcept override;
    };


#include <cmath>

#include "LowPassFilter.h"

LowPassFilter::LowPassFilter ()
{
    poles = 1;
    setCornerFrequency(1.0e5);
}

LowPassFilter::LowPassFilter(const unsigned int poles) : AnalogFilter ()
{
    this->poles = poles;
    setCornerFrequency(1.0e5);
}

LowPassFilter::LowPassFilter (const unsigned int poles, const double cornerFrequency) : AnalogFilter ()
{
    this->poles = poles;
    this->cornerFrequency = cornerFrequency;
    setCornerFrequency(cornerFrequency);
}

void LowPassFilter::setPoles(unsigned int p) noexcept
{
    if (p > 0 and p <= 32 and p % 2 == 0) poles = p;
    setCornerFrequency (cornerFrequency);
}

inline void LowPassFilter::setCornerFrequency (double frequency) noexcept
{
    cornerFrequency = frequency;
    // tau = 1.0 / (1.0 + frequency / (0.5 * poles)) * (44100.0 / (sampleRate ? *sampleRate : 44100.0));
    
    if (poles > 0)
    {
        tau = 1.0 / (2.0 * pi * cornerFrequency)
                   * (sampleRate ? *sampleRate : 44100.0)
                   * std::sqrt (std::pow (2.0, 1.0 / poles) - 1.0);
    }
}

inline void LowPassFilter::process (double& sample, unsigned int channel) noexcept
{
    cap_ptr = caps[channel];
    for (int p = 0; p < poles; p++, cap_ptr++)
    {
        *cap_ptr += (sample - *cap_ptr) / (1.0 + tau);
        sample = *cap_ptr;
    }
}

#include "LowPassFilter.h"

class ResonantLowPassFilter : public LowPassFilter {
public:
    ResonantLowPassFilter();

    ResonantLowPassFilter(const unsigned int poles) : LowPassFilter(poles), resonance(0.0) {}
    
    ResonantLowPassFilter(const unsigned int poles, const double cornerFrequency)
        : LowPassFilter(poles),
          resonance(0.0)
    {
        setCornerFrequency (cornerFrequency);
    }

    void process (double& sample, unsigned int channel) noexcept override;
    
    void setResonance (double res) noexcept;
    
protected:
    double resonance;
    double sampleFeedbackBuffer[32] = {0.0};
};

#include "ResonantLowPassFilter.h"
#include <cmath>

ResonantLowPassFilter::ResonantLowPassFilter () : LowPassFilter()
{
    resonance = 0.0;
    tau = 1.0;
}

inline void ResonantLowPassFilter::process (double& sample, unsigned int channel) noexcept {

    sample -= resonance * sampleFeedbackBuffer[channel];
    LowPassFilter::process(sample, channel);
    sampleFeedbackBuffer[channel] = sample;
    sample *= 1.0 + resonance;

}

void ResonantLowPassFilter::setResonance (double res) noexcept {
    resonance = res * 0.99;
}