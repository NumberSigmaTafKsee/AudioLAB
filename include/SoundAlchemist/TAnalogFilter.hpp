#pragma once
#include <cmath>

namespace SoundAlchemy::Filters
{
    template<typename T>
    class TAnalogFilter {
    public:
        TAnalogFilter() {}
        
        TAnalogFilter(const unsigned int p) : poles(p) {}
        
        TAnalogFilter(const unsigned int p, const double cf) : poles(p), cornerFrequency(cf) {}

        virtual ~TAnalogFilter() {}

        void bindSampleRate (double* sampleRate) { this->sampleRate = sampleRate; }
        
        virtual void setPoles (unsigned int p) noexcept = 0;
        
        virtual void setCornerFrequency (double frequency) noexcept = 0;
        
        virtual void process (double& buffer, unsigned int channel) noexcept = 0;
        
        T getProcessed (T sample, unsigned int channel) noexcept;

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

    template <typename SampleType>
    inline SampleType TAnalogFilter<SampleType>::getProcessed (SampleType sample, unsigned int channel) noexcept
    {
        double sample2 = sample;
        process (sample2, channel);
        return sample2;
    }


    template<typename T>
    class TLowPassFilter : public TAnalogFilter<T> {
        public:
            TLowPassFilter();
            
            TLowPassFilter(const unsigned int poles);
        
            TLowPassFilter(const unsigned int poles, const double cornerFrequency);
        
            void setPoles (unsigned int p) noexcept override;
        
            void setCornerFrequency(double frequency) noexcept override;
        
            void process (double& sample, unsigned int channel) noexcept override;
        };


    template<typename T>
    inline TLowPassFilter<T>::TLowPassFilter ()
    {
        this->poles = 1;
        setCornerFrequency(1.0e5);
    }

    template<typename T>
    inline TLowPassFilter<T>::TLowPassFilter(const unsigned int poles) : TAnalogFilter<T> ()
    {
        this->poles = poles;
        setCornerFrequency(1.0e5);
    }

    template<typename T>
    inline TLowPassFilter<T>::TLowPassFilter (const unsigned int poles, const double cornerFrequency) : TAnalogFilter<T> ()
    {
        this->poles = poles;
        this->cornerFrequency = cornerFrequency;
        setCornerFrequency(cornerFrequency);
    }

    template<typename T>
    inline void TLowPassFilter<T>::setPoles(unsigned int p) noexcept
    {
        if (p > 0 and p <= 32 and p % 2 == 0) this->poles = p;
        setCornerFrequency (this->cornerFrequency);
    }

    template<typename T>
    inline void TLowPassFilter<T>::setCornerFrequency (double frequency) noexcept
    {
        this->cornerFrequency = frequency;
        // tau = 1.0 / (1.0 + frequency / (0.5 * poles)) * (44100.0 / (sampleRate ? *sampleRate : 44100.0));
        
        if (this->poles > 0)
        {
            this->tau = 1.0 / (2.0 * M_PI * this->cornerFrequency)
                    * (this->sampleRate ? *this->sampleRate : 44100.0)
                    * std::sqrt (std::pow (2.0, 1.0 / this->poles) - 1.0);
        }
    }

    template<typename T>
    inline void TLowPassFilter<T>::process (double& sample, unsigned int channel) noexcept
    {
        this->cap_ptr = this->caps[channel];
        for (int p = 0; p < this->poles; p++, this->cap_ptr++)
        {
            *this->cap_ptr += (sample - *this->cap_ptr) / (1.0 + this->tau);
            sample = *this->cap_ptr;
        }
    }

    template<typename T>
    class TResonantLowPassFilter : public TLowPassFilter<T> {
    public:
        TResonantLowPassFilter();

        TResonantLowPassFilter(const unsigned int poles) : TLowPassFilter<T>(poles), resonance(0.0) {}
        
        TResonantLowPassFilter(const unsigned int poles, const double cornerFrequency)
            : TLowPassFilter<T>(poles),
            resonance(0.0)
        {
            this->setCornerFrequency (cornerFrequency);
        }

        void process (double& sample, unsigned int channel) noexcept override;
        
        void setResonance (double res) noexcept;
        
    protected:
        double resonance;
        double sampleFeedbackBuffer[32] = {0.0};
    };


    template<typename T>
    inline TResonantLowPassFilter<T>::TResonantLowPassFilter () : TLowPassFilter<T>()
    {
        resonance = 0.0;
        this->tau = 1.0;
    }

    template<typename T>
    inline void TResonantLowPassFilter<T>::process (double& sample, unsigned int channel) noexcept {

        sample -= resonance * sampleFeedbackBuffer[channel];
        TLowPassFilter<T>::process(sample, channel);
        sampleFeedbackBuffer[channel] = sample;
        sample *= 1.0 + resonance;

    }

    template<typename T>
    inline void TResonantLowPassFilter<T>::setResonance (double res) noexcept {
        resonance = res * 0.99;
    }
}

using LowpassFilter32 = SoundAlchemy::Filters::TLowPassFilter<float>;
using LowpassFilter64 = SoundAlchemy::Filters::TLowPassFilter<double>;
using ResonantLowpassFilter32 = SoundAlchemy::Filters::TResonantLowPassFilter<float>;
using ResonantLowpassFilter64 = SoundAlchemy::Filters::TResonantLowPassFilter<double>;