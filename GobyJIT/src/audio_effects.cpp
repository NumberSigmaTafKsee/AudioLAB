///////////////////////////////////////////////////////////////////
// Audio Effects Rack
///////////////////////////////////////////////////////////////////

#include <cassert>
#include <random>
#include <chrono>
#include <complex>
#include <cstdio>
#include <cstring>
#include <cstdlib>

#include "Vc/vectorclass.h"
#include "Vc/vectormath_lib.h"

#include "Std/StdObject.h"
#include "Std/StdRandom.h"

#include "AudioMidi/audiosystem.h" 



///////////////////////////////////////
// Samples/DSP
///////////////////////////////////////

#include "samples/sample.hpp"
#include "samples/sample_dsp.hpp"
#include "Threads.hpp"

#include "MusicFunctions.hpp"

//#include "DSP/KfrDSP1.hpp"
#include "DSP/DSPResamplers.hpp"

#include "Delay/DelayLines.hpp"
#include "Delay/Delays.hpp"


#include "Distortion/12ax7_table.h"
#include "Distortion/Amplifiers.hpp"
#include "Distortion/Amplifier.hpp"
#include "Distortion/ClipFunctions.hpp"
#include "Distortion/DistortionFunctions.hpp"
#include "Distortion/Waveshapers.hpp"


#include "Waveshapers/waveshapers.h"
#include "Distortion/SstWaveshaper.hpp"

#include "Envelopes/Modulation.hpp"
#include "Envelopes/ADSR.hpp"
#include "Envelopes/GammaEnvelope.hpp"
#include "Envelopes/OnePoleEnvelope.hpp"

/* Soundtouch snips
    soundtouch::SoundTouch *cSoundTouch1;
    cSoundTouch1 = new soundtouch::SoundTouch();
    //===============================================================
    cSoundTouch1->setSampleRate((uint)this->sampleRate);
    cSoundTouch1->setPitch(1.0);
    cSoundTouch1->setTempo(1);
    cSoundTouch1->setRate(1);
    cSoundTouch1->setChannels(2);
    cSoundTouch1->setPitchSemiTones(0);
    
Process Loop
    cSoundTouch1->setChannels(2);
    cSoundTouch1->putSamples(pBuffer, sampleFrames);
    cSoundTouch1->receiveSamples(pBuffer, sampleFrames);
*/


#include "Filters/CSmoothFilter.hpp"
#include "Filters/biquad.hpp"
#include "Filters/VecSVF.hpp"
#include "Filters/Filters.h"
#include "Filters/Filters.hpp"
#include "Filters/cppfilters.hpp"
#include "Filters/StateVariableFilters.hpp"
#include "Filters/MoogFilters.hpp"
#include "Filters/FaustFilters.h"
#include "Filters/DCBlock.hpp"


#include "FX/AudioEffects.hpp"
#include "FX/FreeVerb3.hpp"
#include "FX/Limiter.hpp"
#include "FX/AudioDSP_Delay.h"
#include "FX/AudioDSP_VibraFlange.hpp"
#include "FX/AudioDSP_Chorus.hpp"
#include "FX/AudioDSP_Phaser.h"
#include "FX/FXChorusDelay.hpp"
#include "FX/FXYKChorus.hpp"
#include "FX/FXCE2Chorus.hpp"
#include "FX/Stereofiyer.hpp"

#include "Reverbs/MVerb.h"

#include "Oscillators/BandLimitedOscillators.hpp"
#include "Oscillators/FunctionGenerator.hpp"
#include "Oscillators/FunctionLFO.hpp"
#include "Oscillators/WaveTable.hpp"
#include "Oscillators/FourierWave.hpp"
#include "Oscillators/WaveGenerators.hpp"
#include "Oscillators/WaveFile.hpp"

#include "SndFile.hpp"

#include "Filters/CSmoothFilter.hpp"
#include "Filters/OBXFilter.hpp"

#include "Analog/VCO.hpp"
#include "Analog/VCF.hpp"
#include "Analog/SmoothSVF.hpp"
#include "Analog/VCA.hpp"

#include "Analog/RandomSmoother.hpp"
#include "Analog/SlewLimiter.hpp"
#include "Analog/Diode.hpp"

#include "SndFile.hpp"

#include "Analog/LiquidMoog.hpp"
#include "Distortion/SstWaveshaper.hpp"
#include "Distortion/Waveshaping.hpp"
#include "Distortion/ChebyDistortion.hpp"
#include "Analog/MDFM-1000.hpp"

#include "Analog/VirtualAnalogDiodeLadderFilter.hpp"
#include "Analog/VirtualAnalogStateVariableFilter.hpp"

using namespace std;

const int BufferSize = 256;







Std::RandomMersenne noise;
double  sampleRate        = 44100.0f;
double  inverseSampleRate = 1/44100.0f;

/// ENvelopes = Filters, Curves, Splines, Beziers, Interpolator, Stretch, Compress
Envelopes::ADSR adsr(0.3,0.2,0.9,0.2,sampleRate);
Envelopes::ADSR adsr2(0.25,0.9,0.4,0.5,sampleRate);

// ANALOG-9000 
Analog::VCO vco1(sampleRate,Analog::VCO::SAWTOOTH);
Analog::VCF vcf1(sampleRate,1,0.5);
Analog::VCA vca1(3.0);

OBXFilter   obxfilter;
Moog::MoogFilter1 filter(sampleRate,100,0.5);

// WAVEGENERATORS = Functions, Waves, Wave Tables, Wave Cycles, FFT Wave Table, Fractilizer, Neural Network
SinewaveGenerator lfo1(0.005,sampleRate);
SinewaveGenerator lfo2(0.01,sampleRate);
SinewaveGenerator lfo3(0.1,sampleRate);
SinewaveGenerator lfo4(0.2,sampleRate);

// IR-9000 = Convolution, Overlap, Convolution Filter Overlap, Correlation, IR processor, Impulses, Dirac, LaPlace
DSP::Convolution::ImpulseConvolver convolver("irs/01 Halls 1 Large Hall.wav",BufferSize);
//PffftConvolver convolver("irs/01 Halls 1 Large Hall.wav");

// FX-9000 - IIR/FIR, DelayLine, CombFilter, Modulation, Delays, Reverbs, Filters, Convolution
FX::Chorus chorus;
FX::RingModulator ring;
FX::Tremolo trem;
FX::DistortionEffect dist(FX::DistortionEffect::_dcSigmoid,12);
FX::AnalogSVF svf(sampleRate,1000,0.5);
FX::AutoWah awah;
FX::StereoCompressor compressor;
FX::DelayEffect delay;
FX::Flanger flanger;
FX::Phaser phaser;
FX::PingPongDelay pingpong;
FX::PVPitchShifter pitch;
FX::Vibrato vibrato;
FX::ChorusDelay cdelay;

FX::EarlyReflectionReverb rev1(sampleRate);
FX::HallReverb rev2(sampleRate);
FX::PlateReverb rev3(sampleRate);

FX::EarlyDelayReverb rev4(sampleRate);

JoonasFX::VibraFlange vflanger;
JoonasFX::Chorus jchorus;
JoonasFX::Phaser jphaser;

FX::YKChorus::ChorusEngine ykchorus(sampleRate);
FX::Ce2Chorus ce2L(sampleRate),ce2R(sampleRate);

// FX::
FX::Stereofyier stereofy;
//Limiter<float> limiter(sampleRate,.01,.01,.1,60.0,-20);


// FILTERS-9000 = IIR/FIR/Convolution,Designers,ZPK,SOS,Biquad, Difference, Laplace, Inverse Z, Transfer Function
Filters::FilterBase::FilterType type = Filters::FilterBase::FilterType::LOWPASS ;
Filters::BiquadFilter     biquad(type,1000,sampleRate,0.5);
Filters::Butterworth24db  butter(sampleRate,800.0,0.5);
Filters::OnePoleFilter    onepole(type,sampleRate,1000.0,0.5);
Filters::OneZeroFilter    onezero(type,sampleRate,1000.0,0.5);
Filters::OnePoleOneZeroFilter polezero(type,sampleRate,1000.0,0.5);
Filters::TwoZeroFilter    twozero(type,sampleRate,1000.0,0.5);
Filters::TwoPoleFilter    twopole(type,sampleRate,1000.0,0.5);



float **temp1;
float **temp2;





float Freq=1.0f;
float Vel=1.0f;
double Fcutoff = 1;
double Q = 0.5;
double Qn = 0.5;
oscwave_t wave = OT_SAW;
bool hardSync = true;
float osc1_tune = 0;
float osc2_tune = 7.0/12.0;
float osc_mix = 1.0;

WaveFile * file;

DelayLines::delayline   delayline(0.5);
DelayLines::biquaddelay bqdelay(0.2,0.5);

Distortion::AmplifierN<5> amplifier;
Distortion::QuadAmplifier quadamp;
Distortion::BipolarAmplifier biamp;
Distortion::BipolarAmplifier2 biamp2;
Distortion::TwinAmplifier tamp;
Distortion::RangeAmplifier ramp;
Distortion::Range4Amplifier ramp2;

//Analog::BezierDistortion bdist;

Liquid::LiquidNeuron neuron;
Liquid::LiquidNeuronDelay ndelay(0.5*sampleRate);
Liquid::LiquidMoog moog;
Analog::Slew glide1(1/0.01,sampleRate);
Analog::Slew glide2(1/0.01,sampleRate);
// this always sticks why?
Liquid::LiquidPole cutoff(1/0.01,sampleRate);
Liquid::LiquidPole resonance(1/0.01,sampleRate);

// SstWaveShaper waveshaper(SstWaveShaper::WaveshaperType::wst_soft);


Distortion::ChebyDistortion<4> chebyd;

FFTProcessor fft(256);


// FM5150 FM/Phase Audio Distortion
double fmdistortion(double x, double A=1, double B=1, double X=0, double Y=0)
{
    return sin(2*M_PI*(A*sin(2*M_PI*(B*x+Y))+X));
}
double pmdistortion(double x, double X=0)
{
    return sin(2*M_PI*(x+X));
}



VirtualAnalogDiodeLadderFilter dfilter;
VirtualAnalogStateVariableFilterProcessor svf_plugin;

Distortion::Diode::DiodeClipperNR diode_nr;
Distortion::Diode::DiodeClipperFP diode_fp;

// butterworth
// test all filters

// Filter Designer Handbook
// https://github.com/ruohoruotsi/Butterworth-Filter-Design
// https://github.com/adis300/filter-c
// https://github.com/DeadlyRedCube/ButterworthIIR

// RBU
// https://www.musicdsp.org/en/latest/_downloads/550af030275b020f8f85b422e699ce5f/EQ-Coefficients.pdf

// Zolzer

// AN6
// http://www.apogeebio.com/ddx/PDFs/AN-06.pdf

// equations
// https://vicanek.de/articles/BiquadFits.pdf
// https://vicanek.de/articles/ShelvingFits.pdf

// Designers
// Kfr
// Spuce => just need coefficient and math
// DspFilters => just need coefficient nad math
// Coefficient Designer

// polynomial
// roots
// any H(s) bilinear/z
// biquad section bilinear/z
// all poles bilinear/z
// all zeros bilinear/z

// H(z)
// inverse bilinear/z => H(s)
// prewarp
// bilinear/z => H(z)


namespace Filters
{
    // it's an all pole buterworth no zeros
    // 1 / (1+s) = exp(-2*PI*F/Fs)
    // 2 pole = 1 - 2 * R * cos + R*R
    struct OnePoleFilter : public FilterBase
    {
        double a1;
        double b0;
        double y1;
        double x,y;
        double hp;

        
        OnePoleFilter(FilterType type, double freq, double sample_rate, double resonance) 
        : FilterBase(type,freq,sample_rate,resonance)
        {
            setCoefficients();
        }

        void setCutoff(double f) {
            Fc = f;
        }
        // doesn't do anything
        void setQ(double q) {
            Q = q;
        }

        void setCoefficients() {
            switch(filter_type) {
                case LOWPASS:  setLowPass(Fs,Fc); break;
                case HIGHPASS: setHighPass(Fs,Fc); break;
                default: assert(1==0);
            }
        }
    
       

        void setLowPass(double sampleRate,double frequency) {            
            a1 = std::exp(-2 * M_PI * frequency/sampleRate);
            b0 = 1.0-a1;
        }
        void setHighPass(double sampleRate,double frequency) {
            a1 = -std::exp(-2 * M_PI * frequency/sampleRate);
            b0 = 1.0+a1;
        }
        
        double Tick(double I, double A = 1, double X = 0, double Y = 0) {
            Undenormal denormal;
            setCoefficients();
            x = I;
            y = b0 * x - a1 * y1;
            y1 = y;
            hp = y - I;
            return y;
        }
    };

    // it's a 1 tap fir or a zero butterworth
    // all zero butterworht no poles
    // (1+q) = exp(-2*PI*F/Fs)
    // 2 zero = 1 - 2*R*Cos + R*R
    // 3 zero = (1+q)(q^2 - 2*q +1)
    struct OneZeroFilter : public FilterBase
    {
        double b1;
        double b0;
        double x1;
        double x,y;

        OneZeroFilter(FilterType type, double sample_rate, double freq, double resonance) 
        : FilterBase(type,freq,sample_rate,resonance)
        {            
            setCoefficients();
        }

        void setCoefficients()
        {            
            switch(filter_type) {
                case LOWPASS:  setLowPass(Fc); break;
                case HIGHPASS: setHighPass(Fc); break;
                default: assert(1==0);
            }
        }
        void setCutoff(double f) {
            Fc = f;
        }
        // doesn't do anything
        void setQ(double q) {
            Q = q;
        }
        // I dont know why it's backwards this ought to be lowpass but it's highpass
        void setHighPass(double fc) {
            Fc = fc;
            fc /= Fs;    
            fc  = 0.5 - 0.5*fc;                                
            b0 = (0.54 - 0.46*std::cos(2*M_PI/2)) * std::sin(2*M_PI*fc)/M_PI;
            b1 = (0.54 - 0.46*std::cos(4*M_PI/2)) * std::sin(4*M_PI*fc)/(2*M_PI);            
        }
        // this should be highpass but it's lowpass?
        void setLowPass(double fc) {            
            Fc = fc;            
            fc /= Fs;          
            fc  = 0.5*fc + 0.5;  
            b0 = -(0.54 - 0.46*std::cos(M_PI))   * std::sin(2*M_PI*fc)/M_PI;
            b1 = -(0.54 - 0.46*std::cos(2*M_PI)) * std::sin(4*M_PI*fc)/(2*M_PI);        
        }

        double Tick(double I, double A = 1, double X = 0, double Y = 0) {
            Undenormal denormal;
            setCoefficients();
            x = I;
            y = b0 * x + b1 * x1;
            x1 = x;
            return y;
        }
    };

    // 1 pole butterworth with zero=1
    // 1 / (1+s)
    // 2 pole 1 z = 1 / (s^2 -2*s + 1)    
    // 3 pole 1 z = 1 / (1+s)(s^2 -2*s + 1)
    struct OnePoleOneZeroFilter : public FilterBase
    {
        double a1;
        double b0,b1;
        double x1;
        double y1;
        double x;
        double y;

        OnePoleOneZeroFilter(FilterType type, double freq, double sample_rate, double resonance) 
        : FilterBase(type,freq,sample_rate,resonance)
        {
            setCoefficients();
        }

        void setCoefficients()
        { 
            switch(filter_type) {
                case LOWPASS:  setLowPass(Fs,Fc); break;
                case HIGHPASS: setHighPass(Fs,Fc); break;
                default: assert(1==0);
            }
        }
        void setCutoff(float f) {
            Fc = f;
        }
        void setLowPass(double fSampling, double fCut)
        {
            double w = 2.0 * fSampling;
            double Norm;            
            fCut = fCut/fSampling * 2.0F * M_PI;
            Norm = 1.0 / (fCut + w);
            a1 = (w - fCut) * Norm;
            b0 = 
            b1 = fCut * Norm;
        }
        // it's not very good it breaks apart at lower cutoff
        void setHighPass(double fSampling,double fCut)
        {
            if(fCut < 10) fCut = 10;   
            b0 = fSampling / (fSampling+fCut);
            b1 = -b0;
            a1 = (fSampling-fCut)/(fSampling+fCut);
        }

        double Tick(double in, double A = 1, double X = 0, double Y = 0) {
            Undenormal denormal;
            setCoefficients();
            x = in;
            y = b0*x + b1*x1 - a1*y1;        
            y1 = y;        
            x1 = x;           
            return y;
        }    
    };
    // 2 pole butterworth etc
    // 2 tap sinc FIR
    struct TwoPoleFilter : public FilterBase
    {
        double a1,a2;
        double b0;
        double y1,y2;
        double x,y;

        TwoPoleFilter(FilterType type, double sample_rate, double freq, double resonance) 
        : FilterBase(type,freq,sample_rate,resonance)
        {
            setCoefficients();
        }
        void setCoefficients()
        {
            switch(filter_type) {
                case LOWPASS:  setLowPass(Fs,Fc); break;
                case HIGHPASS: setHighPass(Fs,Fc); break;
                default: assert(1==0);
            }
        }
        void setCutoff(double f) {
            Fc = f;
        }
        // it's resonant
        void setQ(double q) {
            if(q < 0.01) q = 0.01;
            if(q > 0.98) q = 0.98;
            Q = q;
        }
        void setLowPass(double fs, double fc) {            
            fc = 0.5*fc/fs;
            a1 = -2*Q*std::cos(2*M_PI*fc);
            a2 = Q*Q;
            b0 = (1+a1)*(1+a1);
        }
        // it goes the wrong direction of the knob but can't figure it out right 
        void setHighPass(double fs, double fc) {            
            float R = Q;
            fc = (0.5+0.5*fc/fs);            
            a1 = 2*R*std::cos(2*M_PI*(fc));
            a2 = R*R;
            b0 = (1-a1)*(1-a1);
        }
        
        double Tick(double I, double A = 1, double X = 0, double Y = 0) {
            Undenormal denormal;
            setCoefficients();
            x = I;
            y = b0 * x - a1 * y1 - a2 * y2;
            y2 = y1;
            y1 = y;
            return y;
        }
    };
    // Sinc FIR Filter n-taps
    // Resonator Solver (poles)
    // Antiresonator Solver (zeros)

    // it's just a two tap FIR with windowed sinc
    // or 2 zero butterworth no poles
    struct TwoZeroFilter : public FilterBase
    {
        double b0;
        double b1;
        double b2;
        double x1,x2;
        double x,y;
        
        TwoZeroFilter(FilterType type, double sample_rate, double freq, double resonance) 
        : FilterBase(type,freq,sample_rate,resonance)
        {            
            setCoefficients();
        }
        void setCoefficients()
        {            
            switch(filter_type) {
                case LOWPASS:  setLowPass(Fc); break;
                case HIGHPASS: setHighPass(Fc); break;
                default: assert(1==0);
            }
        }
        void setCutoff(double f) {
            Fc = f;
        }
        // doesn't do anything
        void setQ(double q) {
            Q = q;
        }
        // I dont know why it's backwards this ought to be lowpass but it's highpass
        void setHighPass(double fc) {
            Fc = fc;
            fc /= Fs;    
            fc  = 0.5 - 0.5*fc;                                
            b0 = (0.54 - 0.46*std::cos(2*M_PI/3))   * std::sin(2*M_PI*fc)/M_PI;
            b1 = (0.54 - 0.46*std::cos(4*M_PI/3)) * std::sin(4*M_PI*fc)/(2*M_PI);
            b2 = (0.54 - 0.46*std::cos(6*M_PI/3)) * std::sin(6*M_PI*fc)/(3*M_PI);
        }
        // this should be highpass but it's lowpass?
        void setLowPass(double fc) {            
            Fc = fc;            
            fc /= Fs;          
            fc  = 0.5*fc + 0.5;  
            b0 = -(0.54 - 0.46*std::cos(M_PI))   * std::sin(2*M_PI*fc)/M_PI;
            b1 = -(0.54 - 0.46*std::cos(2*M_PI)) * std::sin(4*M_PI*fc)/(2*M_PI);
            b2 = -(0.54 - 0.46*std::cos(3*M_PI)) * std::sin(6*M_PI*fc)/(3*M_PI);
        }
        double Tick(double I, double A = 1, double X = 0, double Y = 0) {
            Undenormal denormal;
            setCoefficients();
            x = I;
            y = b0 * x + b1 * x1 + b2 * x2;
            x2 = x1;
            x1 = x;
            return y;
        }
    };

    const double BUDDA_Q_SCALE = 6.0;
    // puts Q into butterworth
    // Butterworth12db
    // cascade
    class Butterworth24db
    {
    public:
        Butterworth24db(double Fs, double Fc, double Q);
        ~Butterworth24db(void);

        void SetSampleRate(double fs);
        void Set(double cutoff, double q);
        double Run(double input);

        void setCutoff(double f) {
            fc = f/fs;
            if(fc < 0.01) fc = 0.01;
            if(fc > 0.95) fc = 0.95;
            this->Set(fc,q);
        }
        void setQ(double Q) {     
            if(Q < 0.01) Q = 0.01;                 
            if(Q > 999) Q = 999;
            this->Set(fc,Q);
        }
        void InplaceProcess(size_t n, float * inputs) {
            for(size_t i = 0; i < n; i++)
                inputs[i] = Run(inputs[i]);
        }
        void ProcessBlock(size_t n, float * inputs, float * outputs) {
            for(size_t i = 0; i < n; i++)
                outputs[i] = Run(inputs[i]);
        }
    private:
        double t0, t1, t2, t3;
        double coef0, coef1, coef2, coef3;
        double history1, history2, history3, history4;
        double gain,fs,fc,q;
        double min_cutoff, max_cutoff;
    };
    inline Butterworth24db::Butterworth24db(double Fs, double Fc, double Q)
    {
        this->history1 = 0.f;
        this->history2 = 0.f;
        this->history3 = 0.f;
        this->history4 = 0.f;        
        fs = Fs;
        fc = Fc/Fs;
        q  = Q;
        this->SetSampleRate(Fs);
        this->Set(fc,q);
    }    
    inline Butterworth24db::~Butterworth24db(void)
    {
    }
    inline void Butterworth24db::SetSampleRate(double fs)
    {
        double pi = 4.f * atanf(1.f);

        this->t0 = 4.f * fs * fs;
        this->t1 = 8.f * fs * fs;
        this->t2 = 2.f * fs;
        this->t3 = pi / fs;

        this->min_cutoff = fs * 0.01f;
        this->max_cutoff = fs * 0.45f;
    }
    inline void Butterworth24db::Set(double cutoff, double q)
    {
        if (cutoff < this->min_cutoff)
                cutoff = this->min_cutoff;
        else if(cutoff > this->max_cutoff)
                cutoff = this->max_cutoff;

        if(q < 0.f)
                q = 0.f;
        else if(q > 1.f)
                q = 1.f;

        double wp = this->t2 * tanf(this->t3 * cutoff);
        double bd, bd_tmp, b1, b2;

        q *= BUDDA_Q_SCALE;
        q += 1.f;

        b1 = (0.765367f / q) / wp;
        b2 = 1.f / (wp * wp);

        bd_tmp = this->t0 * b2 + 1.f;

        bd = 1.f / (bd_tmp + this->t2 * b1);

        this->gain = bd * 0.5f;

        this->coef2 = (2.f - this->t1 * b2);

        this->coef0 = this->coef2 * bd;
        this->coef1 = (bd_tmp - this->t2 * b1) * bd;

        b1 = (1.847759f / q) / wp;

        bd = 1.f / (bd_tmp + this->t2 * b1);

        this->gain *= bd;
        this->coef2 *= bd;
        this->coef3 = (bd_tmp - this->t2 * b1) * bd;
    }
    inline double Butterworth24db::Run(double input)
    {
        double output = input * this->gain;
        double new_hist;

        output -= this->history1 * this->coef0;
        new_hist = output - this->history2 * this->coef1;

        output = new_hist + this->history1 * 2.f;
        output += this->history2;

        this->history2 = this->history1;
        this->history1 = new_hist;

        output -= this->history3 * this->coef2;
        new_hist = output - this->history4 * this->coef3;

        output = new_hist + this->history3 * 2.f;
        output += this->history4;

        this->history4 = this->history3;
        this->history3 = new_hist;

        return output;
    }

    struct BiquadSection
    {
        double z[3];
        double p[2];

        BiquadSection() {
            memset(z,0,sizeof(z));
            memset(p,0,sizeof(p));
        }
        BiquadSection(double z1, double z2, double z3, double p1, double p2) {
            z[0] = z1;
            z[1] = z2;
            z[2] = z3;
            p[0] = p1;
            p[2] = p2;
        }
        BiquadSection(const BiquadSection & b) {
            memcpy(z,b.z,sizeof(z));
            memcpy(p,b.p,sizeof(p));
        }
        void setCoefficients(double z1, double z2, double z3, double p1, double p2) {
            z[0] = z1;
            z[1] = z2;
            z[2] = z3;
            p[0] = p1;
            p[2] = p2;
        }
        void setCoefficients(double n[3], double d[2]) {
            memcpy(z,n,sizeof(z));
            memcpy(p,d,sizeof(z));
        }
        BiquadSection& operator = (const BiquadSection& b) {
            memcpy(z,b.z,sizeof(z));
            memcpy(p,b.p,sizeof(p));
        }
    };

    using BiquadSOS = std::vector<BiquadSection>;


    struct BiquadTypeI 
    {
        BiquadSection biquad;
        double x,y,x1,x2,y1,y2;

        BiquadTypeI() {
            x=y=0;
            x1=x2=y1=y2=0;
        }
        BiquadTypeI(const BiquadSection &b) : biquad(b)
        {
            x=y=0;
            x1=x2=y1=y2=0;
        }
        BiquadTypeI& operator = (const BiquadTypeI & b) {
            biquad = b.biquad;
            x = b.x;
            y = b.y;
            x1 = b.x1;
            x2 = b.x2;
            y1 = b.y1;
            y2 = b.y2;
            return *this;
        }
        void setBiquad(const BiquadSection & b) {
            biquad = b;
        }

        double Tick(double I, double A=1, double X=0, double Y=0)
        {
            x = I;
            y = biquad.z[0]*x + biquad.z[1]*x1 + biquad.z[2]*x2;
            y = y - biquad.p[0]*y1 - biquad.p[1]*y2;
            x2 = x1;
            x1 = x2;
            y2 = y1;
            y1 = y;
            return A*y;
        }
    };


    struct BiquadTypeICascade
    {  
        BiquadSOS sos;
        std::vector<BiquadTypeI> biquads;
    };

    
    // first 10 butterworth polynomials up to order 10
    struct ButterWorth10
    {
        static BiquadSection poly[10];    
    };

    // butterworth solver  any order possible
    // product 2*k+n-1/2n
    struct ButterworthSolver
    {

    };

    // Polynomial roots
    // H(s) transfer function bilinear/z
    // H(z) transfer function

    // Q equations
    // https://www.earlevel.com/main/2016/09/29/cascading-filters/
    // http://www.sengpielaudio.com/calculator-bandwidth.htm
    
    double Bandwidth(double f1, double f2)
    {
        return f2 - f1;
    }
    double BandwidthQ(double f1, double Q)
    {
        return f1/Q;
    }
    double Q(double f0, double BW) {
        return f0/BW;
    }
    double f0(double BW, double Q) {
        return BW * q;
    }
    double f0ff(double f1, double f2) {
        return sqrt(f1*f2);
    }
    double f1(double f0, fouble f2) {
        return (f0*f0)/f2;
    }
    double f2bw(double f2, double BW) {
        return f2 - BW;
    }
    double f2(double f0, double f1)
    {
        return (f0*f0)/f1;
    }
    double f2bw(double f1, double BW) {
        return f1 + BW;
    }
    double OctaveBWToQ(double N)
    {
        return sqrt(pow(2,N))/(pow(2,N)-1.0);
    }
    double QtoOctaveBandwidth1(double y) {
        return log(y)/log(2);
    }
    double QToOctaveBandwidth2(double Q) {
        double a = log(1 + (1 / (2*Q*Q)));
        double x = 2 + (1/(Q*Q));
        double x1 = x*x;
        double b = x1/4.0 - 1;
        b = sqrt(b);
        return (a+b)/log(2);
    }
    double QToOctaveNSinh(double Q) {
        double a = 2.0/ln(2.0);
        return a * sinh(1/(2*Q));
    }
    double QToOctaveBandwidth4(double Q)
    {
        double a = 2*Q*Q+1/(2*Q*Q);
        double x = (2*Q*Q+1)/Q*Q;
        double x1 = x*x;
        double b = x1/4 - 1;
        b = sqrt(b);
        return log(a+b)/log(2);
    }
    double OctaveRatio(double N) {
        return pow(2.0,N);
    }
    double OctaveRatioF(double f1, double f2) {
        return f2/f1;
    }

    // FIR
    // Window Sinc all bands
    /*
    # Ideal impulse response
    Low pass	2*fC*sin(nwC)/nwC
    High pass	−2*fC*sin(nwC)/nwC
    Band pass	2*f2*sin(nw2)/nw2−2*f1*sin(nw1)/nw1
    Band stop	2*f1*sin(nw1)/nw1−2*f2*sin(nw2)/nw2

    High pass = -LP
    Bandpass  = HP - LP
    Bandstop  = -BP
    Lowshelf  = GLP + HP
    HighShelf = LP + GHP
    Peak/Notch= GBP + BS
    */

    /*
    Digital Filter
    1  = pass
    0  = stop
    -1 = cut
    lowpass
    highpass
    bandpass
    bandreject
    peak boost
    peak cut
    shelf boost
    shelf cut
    boost
    cut
    Digital Allpass = all 1's
    Digital Lowshelf Boost = low(1) + high(0)
    Digital Lowshelf Cut   = low(-1) + high(0)
    */
}

// ADSR
// I = 
// A = amp
// X = base * X
// Y = rate * Y
int audio_callback( const void *inputBuffer, void *outputBuffer,
                            unsigned long framesPerBuffer,
                            const PaStreamCallbackTimeInfo* timeInfo,
                            PaStreamCallbackFlags statusFlags,
                            void *userData )
{    
    float *inputs = (float*)inputBuffer;
    float *output = (float*)outputBuffer;
        
    memset(temp2[0],0,BufferSize*sizeof(float));    
    memset(temp2[1],0,BufferSize*sizeof(float));
    
    float *p;
    float pan  = 0.5;    
    float gain = pow(10.0,3.0/20.0f);    
    float dist = pow(10.0,3/20.0);
    vca1.Randomize();

    //svf.setCutoff(fcq);
    //svf.setQ(Q);

    //ssvf.setCutoff(Fcutoff);
    //ssvf.setResonance(Qn);

    //vsvfL.setCutoff(fcq);
    //vsvfL.setQ(Q);

    //vsvfR.setCutoff(fcq);
    //vsvfR.setQ(Q);

    //osc.setWaveform(wave);
    //oscS.setWaveform(wave);
    
    //if(hardSync) osc.lpS = oscS.lpO;
    //else osc.lpS = NULL;

    // it should be oversampled but I'm not doing it yet.
    
    if(Fcutoff < 0.005) Fcutoff=0.005;
    
    for(size_t i = 0; i < framesPerBuffer; i++) 
    {
        float fcv1  = glide1.Tick(Freq+osc1_tune);
        float fcv2  = glide2.Tick(Freq+osc2_tune);
        float frq1  = cv2freq(fcv1);
        float frq2  = cv2freq(fcv2);    
        float fcutq   = cutoff.Tick(Fcutoff);
        float fresq   = resonance.Tick(Qn);
        float fq      = neuron.Tick(Q);

        bqdelay.setCutoff(frq1);
        bqdelay.setQ(fq);

        vco1.setFrequency(frq1);

        vcf1.setCutoff(fcutq);
        vcf1.setResonance(fresq);

        //moog.setCutoff(fcutq);
        //moog.setResonance(fresq);        
        //svf_plugin.setFilterType(0);        
        //svf_plugin.setCutoffFreq(cv2freq(fcutq));        
        //svf_plugin.setResonance(cv2freq(fcutq));
        //lpf.setCutoff(cv2freq(Fcutoff));
        //lpf.setResonance(fresq);
        //cheby2.setCutoff(cv2freq(fcutq));
        //rcf.setCutoff(cv2freq(Fcutoff));
        //svcf.setCutoff(cv2freq(Fcutoff));        
        //svcf.setResonance(fresq);
        //dfilter.setCutoff(fcutq);
        //dfilter.setResonance(fresq);
        //obxfilter.setResonance(fresq);
        
        float e1  = adsr.Tick();
        float e2  = adsr2.Tick();                
        float l1  = lfo1.Tick();
        float l2  = lfo2.Tick();        
        float l3  = lfo3.Tick();        
        float l4  = lfo4.Tick();        

        float x2  = Vel*vco1.Tick();        
        x2 = vcf1.Tick(x2);//,e1,l1,l2);
                       
        float tick = vca1.Tick(x2,gain,-0.9+0.1*l3,0.9+0.1);
        
        //tick = ndelay.Tick(tick);
        temp1[0][i] = tick * sin((1-pan)*M_PI/2);
        temp1[1][i] = tick * cos(pan*M_PI/2);
        
        temp1[0][i] = ce2L.Tick(temp1[0][i]);
        temp1[1][i] = ce2R.Tick(temp1[1][i]);
        
    }
    //fft.ProcessBuffer(framesPerBuffer,temp1[0]);
    //fft.ProcessBuffer(framesPerBuffer,temp1[1]);
    stereofy.InplaceProcess(framesPerBuffer,temp1);
    
    p = output;
    for(size_t i = 0; i < framesPerBuffer;i++)
    {
        *p++ = temp1[0][i];
        *p++ = temp1[1][i];
    }    
    
    return 0;
}            

#include "effects_midi"
#include "effects_gui"


//#include "effects_repl"
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Lua REPL
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "LuaJIT.hpp"
LuaJIT * lua;

int test(lua_State *L)
{
    printf("test\n");
    return 0;
}

bool is_cmd(const char * cmd, const char * key)
{
    return !strcmp(cmd,key);
}

void strupr(char * s)
{
    for(size_t i = 0; i < strlen(s); i++)
        s[i] = toupper(s[i]);
}

int RandomDist(lua_State *L)
{
    amplifier.RandomClip();
    quadamp.RandomClip();
    biamp.RandomClip();
    biamp2.RandomClip();
    tamp.RandomClip();
    ramp.RandomClip();
    ramp2.RandomClip();
    //bdist.init();
    return 0;
}
void connectLua()
{
    lua = new LuaJIT("main.lua");
    lua->CreateCFunction("randd",RandomDist);
}

void repl() {
    std::string cmd;
    std::cin >> cmd;
    lua->DoCmd(cmd);
}



void testpoly()
{
    float polynomial[] = {0,0,0,0};
    std::vector<float> points(1024);
    for(size_t i = 0; i < 4; i++)
        polynomial[i] = noise.randint(-5,5);
    
    for(size_t i = 0; i < 1024; i++)
    {
        float f = 2*((float)i / 1024.0f)-1;
        points[i] = polynomial[1]*f + polynomial[2]*f*f + polynomial[3]*f*f*f;
    }
    float max = -9999;
    for(size_t i = 0; i < 1024; i++)
        if(fabs(points[i]) > max) max =fabs(points[i]);
    printf("max=%f\n",max);
    
    for(size_t i = 0; i < 1024; i++)
    {
        float f = 2*((float)i / 1024.0f)-1;
        float out = polynomial[1]*f + polynomial[2]*f*f + polynomial[3]*f*f*f;
        printf("%f,",out/max);
    }        
    
        
}

int main()
{    
 
    Init();
    srand(time(NULL));
    int channel=1;
    
    connectLua();
    
    stk::Stk::setSampleRate(sampleRate);
    stk::Stk::setRawwavePath( "rawwaves/" );

    //mdelay.addTap(0.1);
    //delay.addTap(0.25);      
    //pitch.setPitchShift(-12);
    //delay.addTap(0.3);

    // hammerzhang
    //init_sin_lookup_table();    
    //autowah_init(&autowah);
    //chorus_create(&chorus);
    //delay_create(&delay);
    
    temp1 = (float**)calloc(2,sizeof(float*));
    temp2 = (float**)calloc(2,sizeof(float*));
    for(size_t i = 0; i < 2; i++)
    {
        temp1[i] = (float*)calloc(BufferSize,sizeof(float));
        temp2[i] = (float*)calloc(BufferSize,sizeof(float));
    }


    //plugins = new LV2Plugins;
    /*
    lv2flanger =plugin->LoadPlugin("http://polyeffects.com/lv2/flanger");

    lv2flanger->connections[0][0] = 0.5;
    lv2flanger->connections[1][0] = 0.255;
    lv2flanger->connections[2][0] = 0.9;
    lv2flanger->connections[3][0] = 1;
    lv2flanger->connections[4][0] = 0.9;
    lv2flanger->connections[5][1] = 0;                            
    */

    //plugin->infoPlugin(&plugin->plugin);
    //host = new Lv2Host(0,sampleRate,256,"http://polyeffects.com/lv2/flanger",0);
 

    int num_midi = GetNumMidiDevices();

    for(size_t i=0; i < num_midi; i++)
    {
        printf("midi device #%lu: %s\n", i, GetMidiDeviceName(i));
    }
    
    int num_audio = GetNumAudioDevices();
    int pulse = 7;

    for(size_t i=0; i < num_audio; i++)
    {
        if(!strcmp(GetAudioDeviceName(i),"jack")) pulse = i;
        printf("audio device #%lu: %s\n", i, GetAudioDeviceName(i));
    }


    //file = new WaveFile("BabyElephantWalk60.wav");     
    //osc.setSlave(oscS.lpO);
 
    set_note_on_func(note_on);
    set_note_off_func(note_off);
    set_audio_func(audio_callback);
    set_repl_func(repl);
    set_control_change_func(control_change);
    //Thread ui(runfltk,NULL);
    InitMidiDevice(1,3,3);
    InitAudioDevice(pulse,pulse,2,sampleRate,BufferSize);
    //Thread audio(run,NULL);
    //runfltk();
    RunAudio();
    StopAudio();
 
}