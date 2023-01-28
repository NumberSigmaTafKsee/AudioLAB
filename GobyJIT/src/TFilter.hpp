#pragma once

// One Pole/Zero Zolzer
// Biquad + Resonant Biquad + RBJ + Zolzer 
// IIR + Butterworth + Chebyshev + Bessel + Elliptical + Legnedre + ZPK
// FIR + Remez Exchange + Windowed Sinc FIR + Raised Cosine
// microDSPFilters 
// microIIR1
// microFastFilter
// microFastFir
// microSpuce
// microCppFilters
// microCFilters

#include "TModulator.hpp"
#include "SoundAlchemy.hpp"

namespace SoundAlchemy
{

enum FilterType
    {
        LOWPASS,
        HIGHPASS,
        BANDPASS,
        NOTCH,
        PEAK,
        LOWSHELF,
        HIGHSHELF,
        ALLPASS,
    };
    
template<typename T>
class Filter {
    
public:
    
    FilterType type = LOWPASS;
    T Fc,Fs,R;

    virtual ~Filter() = default;
    
    virtual void coefficients(T sampleRate,T frequency, T resonance) = 0;
    
};

template<typename T>
struct OnePoleFilter : public Filter<T>
{
    T a1;
    T b0;
    T y1;
    T x,y;
    T hp;

    
    OnePoleFilter() = default;

    void coefficients(T sampleRate,T frequency, T resonance) 
    {

    }

    /*
    LP:
    recursion: tmp = (1-p)*in + p*tmp with output = tmp
    coefficient: p = (2-cos(x)) - sqrt((2-cos(x))^2 - 1) with x = 2*pi*cutoff/samplerate
    coeficient approximation: p = (1 - 2*cutoff/samplerate)^2

    HP:
    recursion: tmp = (p-1)*in - p*tmp with output = tmp
    coefficient: p = (2+cos(x)) - sqrt((2+cos(x))^2 - 1) with x = 2*pi*cutoff/samplerate
    coeficient approximation: p = (2*cutoff/samplerate)^2
    */

    void LowPass(T sampleRate,T frequency) {
        T x = 2.0*M_PI*frequency/sampleRate;
        T c = std::cos(x);
        this->Fc = frequency;
        this->Fs = sampleRate;
        this->R  = 0;
        a1 = (2.0 - c) - std::sqrt((2.0+c*c)-1);
        b0 = 1.0-a1;
    }
    void HighPass(T sampleRate,T frequency) {
        T x = 2.0*M_PI*frequency/sampleRate;
        T c = std::cos(x);
        this->Fc = frequency;
        this->Fs = sampleRate;
        this->R  = 0;
        a1 = (2.0 - c) - std::sqrt((2.0-c*c)-1);
        b0 = a1-1.0;
    }

    /*
    Process loop (lowpass):
    out = a0*in - b1*tmp;
    tmp = out;

    Simple HP version: subtract lowpass output from the input (has strange behaviour towards nyquist):
    out = a0*in - b1*tmp;
    tmp = out;
    hp = in-out;

    Coefficient calculation:
    x = exp(-2.0*pi*freq/samplerate);
    a0 = 1.0-x;
    b1 = -x;
    */

    void LPF(T sampleRate,T frequency) {
        T x = 2.0*M_PI*frequency/sampleRate;
        T p = std::exp(-x);
        this->Fc = frequency;
        this->Fs = sampleRate;
        this->R  = 0;
        a1 = -x;
        b0 = 1.0-x;
    }
    
    T Tick(T I, T A = 1, T X = 0, T Y = 0) {
        x = I;
        y = b0 * x - a1 * y1;
        y1 = y;
        hp = y - in;
        return y;
    }
};

template<typename T>
struct OneZeroFilter : public Filter<T>
{
    T b1;
    T b0;
    T x1;
    T x,y;

    OneZeroFilter() = default;

    void coefficients(T sampleRate,T frequency, T resonance) 
    {

    }

    T Tick(T I, T A = 1, T X = 0, T Y = 0) {
        x = I;
        y = b0 * x + b1 * x1;
        x1 = x;
        return y;
    }
};


template<typename T>
struct OnePoleOneZeroFilter : public Filter<T>
{
    T a1;
    T b0,b1;
    T x1;
    T y1;
    T x;
    T y;

    OnePoleOneZeroFilter() = default;

    void coefficients(T sampleRate,T frequency, T resonance) 
    {

    }

    void SetLPF(float fCut, float fSampling)
    {
        float w = 2.0 * fSampling;
        float Norm;

        this->Fc = fCut;
        this->Fs = fSampling;
        this->R  = 0;
        fCut *= 2.0F * PI;
        Norm = 1.0 / (fCut + w);
        a1 = (w - fCut) * Norm;
        b0 = b1 = fCut * Norm;
    }

    void SetHPF(float fCut, float fSampling)
    {
        float w = 2.0 * fSampling;
        float Norm;

        this->Fc = fCut;
        this->Fs = fSampling;
        this->R  = 0;

        fCut *= 2.0F * PI;
        Norm = 1.0 / (fCut + w);
        b0 = w * Norm;
        b1 = -a0;
        a1 = (w - fCut) * Norm;
    }

    T Tick(T in, T A = 1, T X = 0, T Y = 0) {
        x = in;
        y = b0*x + b1*x1 - a1*y1;        
        y1 = y;        
        x1 = x;
        return y;
    }    
};


template<typename T>
struct TwoPoleFilter : public Filter<T>
{
    T a1,a2;
    T b0;
    T y1,y2;
    T x,y;

    TwoPoleFilter() = default;

    void coefficients(T sampleRate,T frequency, T resonance) 
    {

    }

    T Tick(T I, T A = 1, T X = 0, T Y = 0) {
        x = I;
        y = b0 * x - a1 * y1 - a2 * y2;
        y2 = y1;
        y1 = y;
        return y;
    }
};

template<typename T>
struct TwoZeroFilter : public Filter<T>
{
    T b0;
    T b1;
    T b2;
    T x1,x2;
    T x,y;

    TwoZeroFilter() = default;

    void coefficients(T sampleRate,T frequency, T resonance) 
    {

    }

    T Tick(T I, T A = 1, T X = 0, T Y = 0) {
        x = I;
        y = b0 * x + b1 * x1 + b2 * x2;
        x2 = x1;
        x1 = x;
        return y;
    }
};

template<typename T>
struct BiquadFilter : public Filter<T>
{
    T a1,a2;
    T b0,b1,b2;
    T x2,x1;
    T y2,y1;
    T x;
    T y;

    BiquadFilter() = default;

    void coefficients(T sampleRate,T frequency, T resonance) 
    {

    }
    /*
    a and b are revrse (a should be b and b should be a)
    Lowpass:
      c = 1.0 / tan(pi * f / sample_rate);

      a1 = 1.0 / ( 1.0 + r * c + c * c);
      a2 = 2* a1;
      a3 = a1;
      b1 = 2.0 * ( 1.0 - c*c) * a1;
      b2 = ( 1.0 - r * c + c * c) * a1;

    Hipass:
      c = tan(pi * f / sample_rate);

      a1 = 1.0 / ( 1.0 + r * c + c * c);
      a2 = -2*a1;
      a3 = a1;
      b1 = 2.0 * ( c*c - 1.0) * a1;
      b2 = ( 1.0 - r * c + c * c) * a1;
    */
    T Tick(T in, T A = 1, T X = 0, T Y = 0) {
        x = in;
        y = b0*x + b1*x1 + b2*x2 - a1*y1 - a2*y2;
        y2 = y1;
        y1 = y;
        x2 = x2;
        x1 = x;
        return y;
    }
};


template<typename T>
struct FilterCascade
{
    std::vector<Biquad*> filters;
    T x,y;

    FilterCascade() = default;

    T Tick(T in, T A = 1, T X = 0, T Y = 0) {
        std::vector<Biquad*> t = filts;
        std::reverse(filts.begin(),filter.end());
        x = in;
        y = x;
        for(size_t i = 0; i < filts.size(); i++)
            y = filters[i].Tick(y,A,X,Y);
        return y;
    }
};

template<typename T>
struct FilterSerializer
{
    std::vector<Biquad*> filters;
    T x,y;

    FilterSerializer() = default;

    T Tick(T in, T A = 1, T X = 0, T Y = 0) {
        std::vector<Biquad*> t = filts;
        std::reverse(filts.begin(),filter.end());
        x = in;        
        BiquadFilter* first = filt[0];
        y = first(x);
        for(size_t i = 1; i < filts.size(); i++)
            y *= filters[i].Tick(x,A,X,Y);
        return y;
    }
};


template<typename T>
struct FilterParallelizer
{
    std::vector<Biquad*> filters;
    T x,y;

    FilterParallelizer() = default;

    T Tick(T in, T A = 1, T X = 0, T Y = 0) {
        std::vector<Biquad*> t = filts;
        std::reverse(filts.begin(),filter.end());
        x = in;        
        BiquadFilter* first = filt[0];
        y = first(x);
        for(size_t i = 1; i < filts.size(); i++)
            y += filters[i].Tick(x,A,X,Y);
         y /= filters.size();
         return y;
    }
};


template<typename T>
class OnePoleHighPassFilter : public Filter<T>, public ModTarget<T>
{   
public:
    virtual void coefficients(float sampleRate, float frequency, float resonance) override;
    virtual void process(float *in, float *out,int numSamples);
    virtual void setModulator(Modulator* mod) override;
    virtual void setModAmount(float amount) override;
    
    OnePoleHighPassFilter();
    virtual ~OnePoleHighPassFilter();
    
private:
    std::shared_ptr<OnePoleFilter<T>> filter1,filter2;
    std::shared_ptr<Modulator> modulator;
    float frequency;
    float modAmount;
    float resonance;    
};

template<typename T>
class OnePoleOneZeroHighPassFilter : public Filter<T>, public ModTarget<T>
{   
public:
    virtual void coefficients(float sampleRate, float frequency, float resonance) override;
    virtual void process(float *in, float *out,int numSamples);
    virtual void setModulator(Modulator* mod) override;
    virtual void setModAmount(float amount) override;
    
    OnePoleOneZeroHighPassFilter();
    virtual ~OnePoleOneZeroHighPassFilter();
    
private:
    std::shared_ptr<OnePoleOneZeroFilter<T>> filter1,filter2;
    std::shared_ptr<Modulator> modulator;
    float frequency;
    float modAmount;
    float resonance;    
};


template<typename T>
class HighPassFilter : public Filter<T>, public ModTarget<T>
{   
public:
    virtual void coefficients(float sampleRate, float frequency, float resonance) override;
    virtual void process(float *in, float *out,int numSamples);
    virtual void setModulator(Modulator* mod) override;
    virtual void setModAmount(float amount) override;
    
    HighPassFilter();
    virtual ~HighPassFilter();
    
private:
    std::shared_ptr<BiquadFilter<T>> filter1,filter2;
    std::shared_ptr<Modulator> modulator;
    float frequency;
    float modAmount;
    float resonance;    
};

HighPassFilter::HighPassFilter() {
    this->filter1 = new IIRFilter();
    this->filter2 = new IIRFilter();
    this->modulator = 0;
	this->modAmount = 0;
}

HighPassFilter::~HighPassFilter() {
    this->filter1 = nullptr;
    this->filter2 = nullptr;
}


void HighPassFilter::coefficients(float sampleRate, float frequency, float resonance) {
    
    this->frequency = frequency / 10;
    this->resonance = resonance;
    
    if (frequency >= sampleRate / 2) {
        frequency = sampleRate / 2;
    }
    
    IIRCoefficients ic1  = IIRCoefficients::makeHighPass (sampleRate, frequency / 10, resonance);
    filter1->setCoefficients(ic1);
    filter2->setCoefficients(ic1);
}

void HighPassFilter::process(float *in, float *out, int numSamples) {
    
    float f = frequency;
    
    if (this->modulator != 0) {
        
        
        if (SynthLab::ADSR* env = dynamic_cast<SynthLab::ADSR*>(this->modulator)) {
            
            if(env->getState() != SynthLab::ADSR::env_idle) {
                f =  this->frequency + (modulator->getOutput() * this->modAmount * (22000 - this->frequency));
            }
            else {
                env->reset();
            }
            
        }
        else {
            f =  this->frequency + (modulator->getOutput() * this->modAmount * 1000);
            // modulator->process();
        }
        
        if (f <= 0) {
            f = 0.1;
        }
        if (f > 22000) {
            f = 22000;
        }
        
        IIRCoefficients ic1  = IIRCoefficients::makeHighPass (44100, f, this->resonance);
        
        filter1->setCoefficients(ic1);
        filter2->setCoefficients(ic1);
    }
    
    this->filter1->processSamples(in,numSamples);
    // in -= numSamples;
    // this->filter2->processSamples(in,numSamples);
}

void HighPassFilter::setModulator(Modulator* mod) {
    this->modulator = mod;
}

void HighPassFilter::setModAmount(float amount) {
    this->modAmount = amount;
}

template<typename T>
struct MoogFilter1 : public Filter<T>
{
//Init
// cutoff = cutoff freq in Hz
//fs = sampling frequency //(e.g. 44100Hz)
//res = resonance [0 - 1] //(minimum - maximum)

    float f,fs,k,p,scale,r,y1,y2,y3,y4,oldx,oldy1,oldy2,oldy3;

    MoogFilter1(float cutoff, float resonance, float sampleRate) {
        
        coefficients(sampleRate,cutoff,resonance);
        y1=y2=y3=y4=oldx=oldy1=oldy2=oldy3=0;
    }

    void coefficients(T sampleRate,T frequency, T resonance) 
    {
        f = 2 * cutoff / fs; //[0 - 1]
        k = 3.6*f - 1.6*f*f -1; //(Empirical tuning)
        p = (k+1)*0.5;
        scale = std::exp((1-p)*1.386249);
        r = resonance*scale;
        
    }
    T Tick(T input, T A=1, T X=0, T Y=0)
    {
        //Loop
        //--Inverted feed back for corner peaking
        x = input - r*y4;

        //Four cascaded onepole filters (bilinear transform)
        y1=x*p + oldx*p - k*y1;
        y2=y1*p+oldy1*p - k*y2;
        y3=y2*p+oldy2*p - k*y3;
        y4=y3*p+oldy3*p - k*y4;

        //Clipper band limited sigmoid
        y4 = y4 - (y4^3)/6;

        oldx = x;
        oldy1 = y1;
        oldy2 = y2;
        oldy3 = y3;

        return y4;
    }
};


class MultimodeFilter : public Filter, public StereoEffect, public ModTarget {

public:
    enum Mode {
        HIGHPASS,
        LOWPASS
    };
    
    MultimodeFilter();
    virtual ~MultimodeFilter();
    
    virtual void coefficients(float sampleRate, float frequency, float resonance) override;
    virtual void processStereo(float *const left, float *const right, const int numSamples) override;
    virtual void setModulator(Modulator* mod) override;
    virtual void setModAmount(float amount) override;
    void setMode(Mode mode);

private:
    
    juce::ScopedPointer<LowPassFilter> lowPassLeft;
    juce::ScopedPointer<HighPassFilter> highPassLeft;
    
    juce::ScopedPointer<LowPassFilter> lowPassRight;
    juce::ScopedPointer<HighPassFilter> highPassRight;
    
    Mode mode;
    JUCE_LEAK_DETECTOR(MultimodeFilter);
    
};

MultimodeFilter::MultimodeFilter() {

    this->lowPassLeft = new LowPassFilter();
    this->highPassLeft = new HighPassFilter();

    this->lowPassRight = new LowPassFilter();
    this->highPassRight = new HighPassFilter();
    
    this->mode = LOWPASS;
}

MultimodeFilter::~MultimodeFilter() {
    this->lowPassLeft = nullptr;
    this->lowPassRight = nullptr;
    this->highPassLeft = nullptr;
    this->highPassRight = nullptr;
}


void MultimodeFilter::setMode(Mode mode) {
    this->mode = mode;
}

void MultimodeFilter::coefficients(float sampleRate, float frequency, float resonance) {
	if (frequency == 0) {
		frequency = 0.1;
	}
    this->lowPassLeft->coefficients(sampleRate, frequency, resonance);
    this->lowPassRight->coefficients(sampleRate, frequency, resonance);
    this->highPassLeft->coefficients(sampleRate, frequency, resonance);
    this->highPassRight->coefficients(sampleRate, frequency, resonance);
}

void MultimodeFilter::processStereo(float *const left, float *const right, const int numSamples) {
    if (this->enabled) {
        if (this->mode == Mode::LOWPASS) {
            this->lowPassLeft->process(left, 0, numSamples);
            this->lowPassRight->process(right, 0, numSamples);
        }
        else {
            this->highPassLeft->process(left, 0, numSamples);
            this->highPassRight->process(right, 0, numSamples);
        }
   
    }
}

void MultimodeFilter::setModulator(Modulator* mod) {
    this->lowPassLeft->setModulator(mod);
    this->lowPassRight->setModulator(mod);
    this->highPassLeft->setModulator(mod);
    this->highPassRight->setModulator(mod);
}

void MultimodeFilter::setModAmount(float amount) {
    this->lowPassLeft->setModAmount(amount);
    this->lowPassRight->setModAmount(amount);
    this->highPassLeft->setModAmount(amount);
    this->highPassRight->setModAmount(amount);
}
};
