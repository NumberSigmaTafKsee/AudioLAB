//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Audio Effects Rack
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Analog   = VAxxx
// Delay    = Zxxx
// Dynamics = Gxxx
// Distortion = Dxxx
// Waveshaper = WSxxx
// FX         = FXxxx
// AudioEffects = Reiss
// DAFX         
// ACA          = libACA/Audio Analysis
// HammerFX     = GNuitar
// RackFX       = RackaRack/3
// Digital Control = Waveshapers + Fourier + Math + Z-domain filters
// DCOxxx
// DCAxxx
// DCFxxx
// Voltage Control = Neural Neworks + Spice + Wave Digital Filters
// VCOxxx
// VCAxxx
// VCFxxx
// IIR Filter = IIRxxx
// FIR filter = FIRxxx
// Convolution = CNVxxx
// Correlation = CORxxx
// FFT         = FFTxxx
// FIR Convolution Filter = CFIRxxx
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
#include "SndFile.hpp"
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


#include "Synthesizer/Envelopes/Modulation.hpp"
#include "Synthesizer/Envelopes/ADSR.hpp"
#include "Synthesizer/Envelopes/GammaEnvelope.hpp"
#include "Synthesizer/Envelopes/OnePoleEnvelope.hpp"


#include "Synthesizer/Oscillators/FunctionGenerator.hpp"
#include "Synthesizer/Oscillators/FunctionLFO.hpp"
#include "Synthesizer/Oscillators/WaveTable.hpp"
#include "Synthesizer/Oscillators/FourierWave.hpp"
#include "Synthesizer/Oscillators/WaveGenerators.hpp"
#include "Synthesizer/Oscillators/WaveFile.hpp"

#include "Filters/IIRFilters.hpp"
#include "Filters/IIRRBJFilters.hpp"
#include "Filters/IIRButterworth.hpp"

#include "FX/CSmoothFilter.hpp"
#include "FX/DelayLines.hpp"
#include "FX/Delays.hpp"
#include "FX/Amplifiers.hpp"
#include "FX/Amplifier.hpp"
#include "FX/ClipFunctions.hpp"
#include "FX/DistortionFunctions.hpp"
#include "FX/Waveshapers.hpp"
#include "FX/Waveshaping.hpp"
#include "FX/ChebyDistortion.hpp"



#include "FX/AudioEffects.hpp"
#include "FX/AudioDSP_Delay.hpp"
#include "FX/AudioDSP_VibraFlange.hpp"
#include "FX/AudioDSP_Chorus.hpp"
#include "FX/AudioDSP_Phaser.hpp"

/*
#include "FX/12ax7_table.h"
#include "FX/Limiter.hpp"
#include "FX/FXChorusDelay.hpp"
#include "FX/FXYKChorus.hpp"
#include "FX/FXCE2Chorus.hpp"
#include "FX/Stereofiyer.hpp"
*/
#include "Reverbs/MVerb.h"
#include "Reverbs/FreeVerb3.hpp"


#include "Analog/VAMoogLadderFilters.hpp"
#include "Analog/VirtualAnalogDiodeLadderFilter.hpp"
#include "Analog/VirtualAnalogStateVariableFilter.hpp"
#include "Analog/VCO.hpp"
#include "Analog/VCF.hpp"
#include "Analog/VCA.hpp"

#include "FX/RandomSmoother.hpp"
#include "FX/SlewLimiter.hpp"
#include "FX/Diode.hpp"
#include "FX/LiquidMoog.hpp"
#include "FX/MDFM-1000.hpp"
#include "FX/SstWaveshaper.hpp"
#include "FX/SstFilters.hpp"
#include "FX/BandLimitedOscillators.hpp"
#include "FX/StateVariableFilters.hpp"

#include "StkHeaders.hpp"
#include "Gamma.hpp"

#include "Filters/ATK.hpp"
#include "Filters/ATKAdaptiveFilters.hpp"

using namespace std;

const int BufferSize = 256;

Std::RandomMersenne noise;
double sampleRate = 44100.0f;
double inverseSampleRate = 1 / 44100.0f;

/// ENvelopes = Filters, Curves, Splines, Beziers, Interpolator, Stretch, Compress
Envelopes::ADSR adsr(0.3, 0.2, 0.9, 0.2, sampleRate);
Envelopes::ADSR adsr2(0.25, 0.9, 0.4, 0.5, sampleRate);

// ANALOG-9000
Analog::VCO vco1(sampleRate, Analog::VCO::SAWTOOTH);
Analog::VCF vcf1(sampleRate, 1, 0.5);
Analog::VCA vca1(3.0);

//Moog::MoogFilter1 filter(sampleRate, 100, 0.5);

SinewaveGenerator lfo1(0.005, sampleRate);
SinewaveGenerator lfo2(0.01, sampleRate);
SinewaveGenerator lfo3(0.1, sampleRate);
SinewaveGenerator lfo4(0.2, sampleRate);


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


JoonasFX::VibraFlange vflanger;
JoonasFX::Chorus jchorus;
JoonasFX::Phaser jphaser;

/*
FX::ChorusDelay cdelay;
FX::EarlyReflectionReverb rev1(sampleRate);
FX::HallReverb rev2(sampleRate);
FX::PlateReverb rev3(sampleRate);

FX::EarlyDelayReverb rev4(sampleRate);


FX::YKChorus::ChorusEngine ykchorus(sampleRate);
FX::Ce2Chorus ce2L(sampleRate),ce2R(sampleRate);

// FX::
FX::Stereofyier stereofy;
//Limiter<float> limiter(sampleRate,.01,.01,.1,60.0,-20);
*/

float **temp1;
float **temp2;

float Freq = 1.0f;
float Vel = 1.0f;
double Fcutoff = 1;
double Q = 0.5;
double Qn = 0.5;
bool hardSync = true;
float osc1_tune = 0;
float osc2_tune = 7.0 / 12.0;
float osc_mix = 1.0;

WaveFile *file;

DelayLines::delayline delayline(0.5);
DelayLines::biquaddelay bqdelay(0.2, 0.5);

/*
Distortion::AmplifierN<5> amplifier;
Distortion::QuadAmplifier quadamp;
Distortion::BipolarAmplifier biamp;
Distortion::BipolarAmplifier2 biamp2;
Distortion::TwinAmplifier tamp;
Distortion::RangeAmplifier ramp;
Distortion::Range4Amplifier ramp2;
*/
// Analog::BezierDistortion bdist;

Liquid::LiquidNeuron neuron;
Liquid::LiquidNeuronDelay ndelay(0.5 * sampleRate);
Liquid::LiquidMoog moog;
Analog::RateLimiters::Slew glide1(1 / 0.01, sampleRate);
Analog::RateLimiters::Slew glide2(1 / 0.01, sampleRate);
Liquid::LiquidPole cutoff(1 / 0.01, sampleRate);
Liquid::LiquidPole resonance(1 / 0.01, sampleRate);

// SstFilters
// SstWaveShaper waveshaper(SstWaveShaper::WaveshaperType::wst_soft);

Distortion::ChebyDistortion<4> chebyd;

VirtualAnalogDiodeLadderFilter dfilter;
VirtualAnalogStateVariableFilterProcessor svf_plugin;

Distortion::Diode::DiodeClipperNR diode_nr;
Distortion::Diode::DiodeClipperFP diode_fp;

// FFTProcessor fft(256);

// FM5150 FM/Phase Audio Distortion
double fmdistortion(double x, double A = 1, double B = 1, double X = 0, double Y = 0)
{
    return std::sin(2 * M_PI * (A * std::sin(2 * M_PI * (B * x + Y)) + X));
}
double pmdistortion(double x, double X = 0)
{
    return std::sin(2 * M_PI * (x + X));
}

template <typename SIMD>
SIMD vector_function(SIMD input)
{
    // computer function
    SIMD out = 1 + exp(input);
    return out;
}


double ampTodB(double power)
{
    return pow(10.0, power / 20.0);
}
double dBToAmp(double db)
{
    return 20.0 * log10(db);
}



// FFT
// Convolution
// Correlation
// Biquads
// IIR
// FIR
// Convolution Filter
// Resample
// Mixer
// Signal Morphing/Blending
// Modulation *,%


namespace Filters
{   

    // H(s) transfer function bilinear/z
    struct TransferFunctionLaplace
    {
        std::vector<double> num, den;
        
        TransferFunctionLaplace(size_t order) {
            num.resize(order);
            den.resize(order);
        }
        void setNumS(int power, double val) {
            num[power] = val;
        }
        void setDenS(int power, double val) {
            den[power] = val;
        }
    };

    // H(z) transfer function
    struct TransferFunctionDigital
    {
        std::vector<double> num, den;
    };


    struct Integrator
    {
        BiquadTypeI b;

        Integrator()
        {
            BiquadSection s;
            s.z[0] = 1;
            s.z[1] = 0;
            s.z[2] = 0;
            s.p[0] = 1;
            s.p[1] = 1;
            s.p[2] = 0;
            //            BiquadSection section = AnalogBiquadSection(s[0], fc, sampleRate);
            b = Filters::IIR::AnalogBiquadSection(s,sampleRate/2.0,sampleRate);
        }
        double Tick(double I, double A=1, double X=0, double Y=0)
        {
            return b.Tick(I,A,X,Y);
        }
    };
    //https://en.wikipedia.org/wiki/Differentiator
    struct Differentiator
    {
        BiquadTypeI b;

        Differentiator()
        {
            BiquadSection s;
            s.z[0] = 1;
            s.z[1] = 1;
            s.z[2] = 0;
            s.p[0] = 1;
            s.p[1] = 0;
            s.p[2] = 0;
            //            BiquadSection section = AnalogBiquadSection(s[0], fc, sampleRate);
            b = Filters::IIR::AnalogBiquadSection(s,sampleRate/2.0,sampleRate);
        }
        double Tick(double I, double A=1, double X=0, double Y=0)
        {
            return b.Tick(I,A,X,Y);
        }
    };

}





double Fc = 0;
double Kc = 0;

// RBJ Filters
Filters::IIR::RBJFilters::RBJBiquadFilter rbj;
Filters::IIR::RBJFilters::RBJLowPassFilter rbjlp;
Filters::IIR::RBJFilters::RBJHighPassFilter rbjhp;
Filters::IIR::RBJFilters::RBJAllPassFilter rbjap;
Filters::IIR::RBJFilters::RBJBandPassFilter rbjbp;
Filters::IIR::RBJFilters::RBJSkirtBandPassFilter rbjsbp;
Filters::IIR::RBJFilters::RBJPeakFilter rbjpeak;
Filters::IIR::RBJFilters::RBJLowShelfFilter rbjlsh;
Filters::IIR::RBJFilters::RBJHighShelfFilter rbjhsh;

#include "Filters/IIRZolzerFilter.hpp"

// Zolzer Filters
Filters::IIR::ZolzerFilters::ZolzerBiquadFilter zolz;
Filters::IIR::ZolzerFilters::ZolzerLowPassFilter zolzlp;
Filters::IIR::ZolzerFilters::ZolzerHighPassFilter zolzhp;
Filters::IIR::ZolzerFilters::ZolzerAllPassFilter zolzap;
Filters::IIR::ZolzerFilters::ZolzerLowPass1pFilter zolzlp1;
Filters::IIR::ZolzerFilters::ZolzerHighPass1pFilter zolzhp1;
Filters::IIR::ZolzerFilters::ZolzerAllPass1pFilter zolzap1;
Filters::IIR::ZolzerFilters::ZolzerBandPassFilter zolzbp;
Filters::IIR::ZolzerFilters::ZolzerNotchFilter zolzbs;
Filters::IIR::ZolzerFilters::ZolzerPeakBoostFilter zolpb;
Filters::IIR::ZolzerFilters::ZolzerPeakCutFilter zolpc;
Filters::IIR::ZolzerFilters::ZolzerLowShelfBoostFilter zolshb;
Filters::IIR::ZolzerFilters::ZolzerLowShelfCutFilter zolshc;
Filters::IIR::ZolzerFilters::ZolzerHighShelfBoostFilter zohshb;
Filters::IIR::ZolzerFilters::ZolzerHighShelfCutFilter zohshc;

// Butterworth Filters
Filters::IIR::ButterworthFilters::ButterworthLowPassCascadeFilter           butter(3);
Filters::IIR::ButterworthFilters::ButterworthResonantLowPassCascadeFilter   rbutter(3);
Filters::IIR::ButterworthFilters::ButterworthDampedLowPassCascadeFilter     dbutter(3);

Filters::IIR::ButterworthFilters::ButterworthLowPassFilter12db blp;
Filters::IIR::ButterworthFilters::ButterworthResonantLowPassFilter12db brlp;
Filters::IIR::ButterworthFilters::ButterworthDampedLowPassFilter12db bdlp;

// Resonant = Radius
// Damped   = Q
Filters::IIR::ButterworthFilters::ButterworthHighPassFilter12db bhp;
Filters::IIR::ButterworthFilters::ButterworthBandPassFilter12db bbp;
Filters::IIR::ButterworthFilters::ButterworthBandStopFilter12db bsp;

#include "Filters/IIRChebyshevFilters.hpp"
// Chebyshev 1
Filters::IIR::ChebyshevFilters::ChebyshevILowPassFilter12db c1lp;

// Chebyshev 2
Filters::IIR::ChebyshevFilters::ChebyshevIILowPassFilter12db c2lp;



// Elliptical 
// Filter Designer Handbook

// Bessel
// Filter Designer Handbook

// ZPK
// zeros [z1,z2...]
// poles [p1,p2...]
// k = gain

// Polynomial = 1 + s + s^2 (least to most)
// Laguerre Root Finder
// RootFinder
// TransferFunctionAnalog
// TransferFunctionDigital


// Octave
// butter
// cheby1
// cheby2
// bessel = analog
// elliptical

// DSP Filters
// todo get analog prototype of dsp filters for R and Q
#include "Filters/IIRBesselFilterProcessor.hpp"
#include "Filters/IIRChebyshevFilterProcessors.hpp"
Filters::IIR::Bessel::LowPassFilter bessellp(4,1000.0f,sampleRate);
//Filters::IIR::ChebyshevI::LowPassFilter bessellp(4,1000.0f,sampleRate);

// Spuce
// todo - useful has good stuff


int audio_callback(const void *inputBuffer, void *outputBuffer,
                   unsigned long framesPerBuffer,
                   const PaStreamCallbackTimeInfo *timeInfo,
                   PaStreamCallbackFlags statusFlags,
                   void *userData)
{
    float *inputs = (float *)inputBuffer;
    float *output = (float *)outputBuffer;

    memset(temp2[0], 0, BufferSize * sizeof(float));
    memset(temp2[1], 0, BufferSize * sizeof(float));

    float *p;
    float pan = 0.5;
    float gain = pow(10.0, 3.0 / 20.0f);
    float dist = pow(10.0, 3 / 20.0);

    vca1.Randomize();

    // svf.setCutoff(fcq);
    // svf.setQ(Q);

    // ssvf.setCutoff(Fcutoff);
    // ssvf.setResonance(Qn);

    // vsvfL.setCutoff(fcq);
    // vsvfL.setQ(Q);

    // vsvfR.setCutoff(fcq);
    // vsvfR.setQ(Q);

    // osc.setWaveform(wave);
    // oscS.setWaveform(wave);

    // if(hardSync) osc.lpS = oscS.lpO;
    // else osc.lpS = NULL;

    // it should be oversampled but I'm not doing it yet.

    for (size_t i = 0; i < framesPerBuffer; i++)
    {
        float fcv1 = glide1.Tick(Freq + osc1_tune);
        float fcv2 = glide2.Tick(Freq + osc2_tune);
        float frq1 = cv2freq(fcv1);
        float frq2 = cv2freq(fcv2);
        float fcutq = cutoff.Tick(Fcutoff);
        float fresq = resonance.Tick(Qn);
        float fq = neuron.Tick(Q);

        //bqdelay.setCutoff(frq1);
        //bqdelay.setQ(fq);

        vco1.setFrequency(frq1);

        vcf1.setCutoff(fcutq);
        vcf1.setResonance(fresq);

        //rbj.setCutoff(Fc+Kc);
        //rbj.setQ(fq);

        //c1lp.setCutoff(Fc);
        //c1lp.setQ(fq);
        c2lp.setCutoff(Fc);
        c2lp.setQ(fq);

        
        

        // moog.setCutoff(fcutq);
        // moog.setResonance(fresq);
        // svf_plugin.setFilterType(0);
        // svf_plugin.setCutoffFreq(cv2freq(fcutq));
        // svf_plugin.setResonance(cv2freq(fcutq));
        // lpf.setCutoff(cv2freq(Fcutoff));
        // lpf.setResonance(fresq);
        // cheby2.setCutoff(cv2freq(fcutq));
        // rcf.setCutoff(cv2freq(Fcutoff));
        // svcf.setCutoff(cv2freq(Fcutoff));
        // svcf.setResonance(fresq);
        // dfilter.setCutoff(fcutq);
        // dfilter.setResonance(fresq);
        // obxfilter.setResonance(fresq);

        float e1 = adsr.Tick();
        float e2 = adsr2.Tick();
        float l1 = lfo1.Tick();
        float l2 = lfo2.Tick();
        float l3 = lfo3.Tick();
        float l4 = lfo4.Tick();

        float x2 = Vel * vco1.Tick();
        x2 = vcf1.Tick(x2); //,e1,l1,l2);

        float tick = vca1.Tick(x2, gain, -0.9 + 0.1 * l3, 0.9 + 0.1);

        // tick = ndelay.Tick(tick);
        temp1[0][i] = tick * std::sin((1 - pan) * M_PI / 2);
        temp1[1][i] = tick * std::cos(pan * M_PI / 2);

        // temp1[0][i] = ce2L.Tick(temp1[0][i]);
        // temp1[1][i] = ce2R.Tick(temp1[1][i]);
    }
    // fft.ProcessBuffer(framesPerBuffer,temp1[0]);
    // fft.ProcessBuffer(framesPerBuffer,temp1[1]);
    // stereofy.InplaceProcess(framesPerBuffer,temp1);

    p = output;
    for (size_t i = 0; i < framesPerBuffer; i++)
    {
        *p++ = temp1[0][i];
        *p++ = temp1[1][i];
    }

    return 0;
}

#include "effects_midi"
//#include "effects_gui"
//#include "effects_repl"


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Lua REPL
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "LuaJIT.hpp"
LuaJIT *lua;

int test(lua_State *L)
{
    printf("test\n");
    return 0;
}

bool is_cmd(const char *cmd, const char *key)
{
    return !strcmp(cmd, key);
}

void strupr(char *s)
{
    for (size_t i = 0; i < strlen(s); i++)
        s[i] = toupper(s[i]);
}

int RandomDist(lua_State *L)
{
    /*
    amplifier.RandomClip();
    quadamp.RandomClip();
    biamp.RandomClip();
    biamp2.RandomClip();
    tamp.RandomClip();
    ramp.RandomClip();
    ramp2.RandomClip();
    // bdist.init();
    */
    return 0;
}
void connectLua()
{
    lua = new LuaJIT("main.lua");
    lua->CreateCFunction("randd", RandomDist);
}

void repl()
{
    std::string cmd;
    std::cin >> cmd;
    lua->DoCmd(cmd);
}

void testpoly()
{
    float polynomial[] = {0, 0, 0, 0};
    std::vector<float> points(1024);
    for (size_t i = 0; i < 4; i++)
        polynomial[i] = noise.randint(-5, 5);

    for (size_t i = 0; i < 1024; i++)
    {
        float f = 2 * ((float)i / 1024.0f) - 1;
        points[i] = polynomial[1] * f + polynomial[2] * f * f + polynomial[3] * f * f * f;
    }
    float max = -9999;
    for (size_t i = 0; i < 1024; i++)
        if (fabs(points[i]) > max)
            max = fabs(points[i]);
    printf("max=%f\n", max);

    for (size_t i = 0; i < 1024; i++)
    {
        float f = 2 * ((float)i / 1024.0f) - 1;
        float out = polynomial[1] * f + polynomial[2] * f * f + polynomial[3] * f * f * f;
        printf("%f,", out / max);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Main
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main()
{   
    Init();
    srand(time(NULL));
    int channel = 1;

    connectLua();

    stk::Stk::setSampleRate(sampleRate);
    stk::Stk::setRawwavePath("rawwaves/");

    temp1 = (float **)calloc(2, sizeof(float *));
    temp2 = (float **)calloc(2, sizeof(float *));
    for (size_t i = 0; i < 2; i++)
    {
        temp1[i] = (float *)calloc(BufferSize, sizeof(float));
        temp2[i] = (float *)calloc(BufferSize, sizeof(float));
    }

    // plugins = new LV2Plugins;
    /*
    lv2flanger =plugin->LoadPlugin("http://polyeffects.com/lv2/flanger");

    lv2flanger->connections[0][0] = 0.5;
    lv2flanger->connections[1][0] = 0.255;
    lv2flanger->connections[2][0] = 0.9;
    lv2flanger->connections[3][0] = 1;
    lv2flanger->connections[4][0] = 0.9;
    lv2flanger->connections[5][1] = 0;
    */

    // plugin->infoPlugin(&plugin->plugin);
    // host = new Lv2Host(0,sampleRate,256,"http://polyeffects.com/lv2/flanger",0);

    int num_midi = GetNumMidiDevices();

    for (size_t i = 0; i < num_midi; i++)
    {
        printf("midi device #%lu: %s\n", i, GetMidiDeviceName(i));
    }

    int num_audio = GetNumAudioDevices();
    int pulse = 7;

    for (size_t i = 0; i < num_audio; i++)
    {
        if (!strcmp(GetAudioDeviceName(i), "jack"))
            pulse = i;
        printf("audio device #%lu: %s\n", i, GetAudioDeviceName(i));
    }

    // file = new WaveFile("BabyElephantWalk60.wav");
    // osc.setSlave(oscS.lpO);

    set_note_on_func(note_on);
    set_note_off_func(note_off);
    set_audio_func(audio_callback);
    set_repl_func(repl);
    set_control_change_func(control_change);
    // Thread ui(runfltk,NULL);
    InitMidiDevice(1, 3, 3);
    InitAudioDevice(pulse, pulse, 2, sampleRate, BufferSize);
    // Thread audio(run,NULL);
    // runfltk();
    RunAudio();
    StopAudio();
}