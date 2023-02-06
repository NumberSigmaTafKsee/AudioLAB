#include <cassert>
#include <random>
#include <chrono>
#include <complex>
#include <cstdio>
#include <cstring>
#include <cstdlib>

#include "SoundObject.hpp"
#include "SndFile.hpp"
#include "AudioMidi/audiosystem.h"
#include "samples/sample.hpp"
#include "samples/sample_dsp.hpp"
#include "Threads.hpp"
#include "MusicFunctions.hpp"

///////////////////////////////////////
// Samples/DSP
///////////////////////////////////////


//#include "DSP/KfrDSP1.hpp"
#include "DSP/DSPResamplers.hpp"

#include "FX/Modulation.hpp"
#include "FX/ADSR.hpp"
#include "FX/GammaEnvelope.hpp"
#include "FX/FunctionGenerator.hpp"
#include "FX/FunctionLFO.hpp"
#include "FX/WaveTable.hpp"
#include "FX/FourierWave.hpp"
#include "FX/WaveGenerators.hpp"
#include "FX/WaveFile.hpp"

#include "Filters/IIRFilters.hpp"
#include "Filters/IIRRBJFilters.hpp"
#include "Filters/IIRButterworth.hpp"

#include "FX/CSmoothFilters.hpp"
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

#include "FX/RandomSmoother.hpp"
#include "FX/SlewLimiter.hpp"
#include "FX/Diode.hpp"
#include "FX/LiquidMoog.hpp"
#include "FX/MDFM-1000.hpp"
#include "FX/SstWaveshaper.hpp"
#include "FX/SstFilters.hpp"
#include "FX/BandLimitedOscillators.hpp"
#include "FX/StateVariableFilters.hpp"

//#include "FX/ATK.hpp"
//#include "FX/ATKAdaptiveFilters.hpp"

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
#include "Reverbs/FV3Processor.hpp"

#include "Analog/VAMoogLadderFilters.hpp"
#include "Analog/VirtualAnalogDiodeLadderFilter.hpp"
#include "Analog/VirtualAnalogStateVariableFilter.hpp"
#include "Analog/VCO.hpp"
#include "Analog/VCF.hpp"
#include "Analog/VCA.hpp"


#include "Synthesizer/StkHeaders.hpp"
#include "Synthesizer/Gamma.hpp"


using namespace std;

const int BufferSize = 256;

Default noise;
DspFloatType sampleRate = 44100.0f;
DspFloatType inverseSampleRate = 1 / 44100.0f;

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

//////////////////////////////////////////////////////
// FX
// AudioEffects
// DAFX
// HammerFX
// RackFX
// LV2
// Faust
//////////////////////////////////////////////////////
FX::Chorus chorus;
FX::RingModulator ring;
FX::Tremolo trem;
FX::DistortionEffect dist(FX::DistortionEffect::_dcSigmoid,12);
FX::AnalogSVF svf(sampleRate,1000,0.5);
FX::AutoWah awah(44100.0);
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
*/

/////////////////////////////////////////////
// Dynamics
/////////////////////////////////////////////
// Compressors
// Limiters
// Gates
// Expanders
/////////////////////////////////////////////
//Limiter<DspFloatType> limiter(sampleRate,.01,.01,.1,60.0,-20);

DspFloatType **temp1;
DspFloatType **temp2;

DspFloatType Freq = 1.0f;
DspFloatType Vel = 1.0f;
DspFloatType Fcutoff = 1;
DspFloatType Q = 0.5;
DspFloatType Qn = 0.5;
bool hardSync = true;
DspFloatType osc1_tune = 0;
DspFloatType osc2_tune = 7.0 / 12.0;
DspFloatType osc_mix = 1.0;

WaveFile *file;

// Gamma RingBuffer

///////////////////////////////////////////////////////////////////
// Delays.hpp
// DelayLines.hpp
///////////////////////////////////////////////////////////////////
FX::Delays::delayline delayline(0.5);
FX::Delays::biquaddelay bqdelay(0.2, 0.5);


/////////////////////////////////////////////
// FX and Analog
/////////////////////////////////////////////

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

FX::Distortion::Chebyshev::ChebyDistortion<4> chebyd;

Analog::Filters::VirtualAnalogDiodeLadderFilter dfilter;
Analog::Filters::VirtualAnalogStateVariableFilterProcessor svf_plugin;

FX::Distortion::Diode::DiodeClipperNR diode_nr;
FX::Distortion::Diode::DiodeClipperFP diode_fp;

// FFTProcessor fft(256);

// FM5150 FM/Phase Audio Distortion
DspFloatType fmdistortion(DspFloatType x, DspFloatType A = 1, DspFloatType B = 1, DspFloatType X = 0, DspFloatType Y = 0)
{
    return std::sin(2 * M_PI * (A * std::sin(2 * M_PI * (B * x + Y)) + X));
}
DspFloatType pmdistortion(DspFloatType x, DspFloatType X = 0)
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


DspFloatType ampTodB(DspFloatType power)
{
    return pow(10.0, power / 20.0);
}
DspFloatType dBToAmp(DspFloatType db)
{
    return 20.0 * log10(db);
}

////////////////////////////////////////////////////////////////////////////////////////////////
// Synthesizer
////////////////////////////////////////////////////////////////////////////////////////////////
// FM
// Granular
// Physical Models
// Oscillators
// Noise
// Bristol
// VCVRack
// Surge/SST
// SstFilters
// SstWaveShaper waveshaper(SstWaveShaper::WaveshaperType::wst_soft);
// libSurge/vcvsurge
// libZynAddSubFX
////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////
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
////////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////////////
// Filters
// IIR
// FIR
// FIR Convolution
// Transfer Functions
////////////////////////////////////////////////////////////////////////////////////////////////

DspFloatType Fc = 0;
DspFloatType Kc = 0;

///////////////////////////////////////////////////////////////
// RBJ Filters
///////////////////////////////////////////////////////////////
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

///////////////////////////////////////////////////////////////
// Zolzer Filters
///////////////////////////////////////////////////////////////
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

// these are experimental and not stable

///////////////////////////////////////////////////////////////
// Butterworth Filters
///////////////////////////////////////////////////////////////
Filters::IIR::ButterworthFilters::ButterworthLowPassCascadeFilter           butter(4);
Filters::IIR::ButterworthFilters::ButterworthResonantLowPassCascadeFilter   rbutter(3);
Filters::IIR::ButterworthFilters::ButterworthDampedLowPassCascadeFilter     dbutter(3);

Filters::IIR::ButterworthFilters::ButterworthLowPassFilter12db blp;
Filters::IIR::ButterworthFilters::ButterworthResonantLowPassFilter12db brlp;
Filters::IIR::ButterworthFilters::ButterworthDampedLowPassFilter12db bdlp;

Filters::IIR::ButterworthFilters::ButterworthHighPassFilter12db bhp;

///////////////////////////////////////////////////////////////
// Chebyshev 1
///////////////////////////////////////////////////////////////
#include "Filters/IIRChebyshevFilters.hpp"
Filters::IIR::ChebyshevFilters::ChebyshevILowPassFilter12db c1lp;


///////////////////////////////////////////////////////////////
// Chebyshev 2
///////////////////////////////////////////////////////////////
#include "Filters/IIRChebyshevFilterProcessors.hpp"
#include "Filters/IIRChebyshev2FilterProcessors.hpp"
Filters::IIR::ChebyshevFilters::ChebyshevIILowPassFilter12db c2lp;


///////////////////////////////////////////////////////////////
// Elliptical
///////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////
// Bessel
///////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////
// KFR 
///////////////////////////////////////////////////////////////
// IIR
// FIR
// Biquad
// Convolution Filter
// Digital Convolution Filter
///////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////
// OctaveFilters
// butter
// cheby1
// cheby2
// bessel
// elliptical
///////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////
// Poles
// Analog IIR
// Filter Designer Handbook
// C++ DSP
// Analog Filters
// VA Filter Design
// IIR Filter Design
// DSPFilters
// Spuce
// CppFilters
///////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////
// FIR
// FIR.hpp
// frequency sampling = https://github.com/devtibo/FilterDesign
// firpm = parks-mclellans
// Fir1  = https://github.com/berndporr/fir1
// kiss_fastfir
// Stk::Fir
// AudioTK::FIRFilter
// fast filters = https://github.com/jpcima/fast-filters
// https://github.com/jontio/FastFIR
///////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Spectrum
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Gamma = DFT,STFT,SDFT,CFFT,RFFT
// audiofft+fftconvolver
// pffft+pffastconv
// vcvfft = pffft
// rackdsp = convolution + dsp
// fftwpp = fftw3 + convolution
// Matx   = GPU
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// https://github.com/Flixor/IRBaboon
// https://github.com/grahman/RTConvolve
// https://github.com/GuillaumeLeNost/ConvolutionFilter
// https://github.com/deddy11/test_overlap_add
// https://github.com/fbdp1202/DSP_Overlap-Save_Overlap_Add
// https://github.com/ybdarrenwang/PhaseVocoder
// https://github.com/ljtrevor/PhaseVocoder
// https://github.com/stekyne/PhaseVocoder
// https://github.com/sevagh/pitch-detection
// https://github.com/xstreck1/LibPyin
// https://github.com/fftune/fftune
// https://github.com/karlmess/pitch-detection
// #include "OverlapAddBuffer.hpp"
// #include "SimpleOverlapAddProcessor.hpp"
// #include "ConvolutionProcessors.hpp"
// #include "ConvolutionFilter.hpp"
// #include "CxxConvolver.hpp"
// #include "OverlappingFFTProcessor.hpp"
// #include "STFTProcessors.hpp"
// #include "FIRConvolutionFilters.hpp"
// ConvolutionMatrix.hpp
// FFTConvolver.hpp
// ImpulseResponse.hpp = SpectrumConvolutionFilters.hpp
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include "Spectrum/Spectrum.hpp"
#include "Spectrum/FFTProcessors.hpp"


////////////////////////////////////////////////////////
// Convolution
////////////////////////////////////////////////////////
#include "Spectrum/SpectrumConvolutionFilters.hpp"

////////////////////////////////////////////////////////
// Correlation
////////////////////////////////////////////////////////

Spectrum::FFT::FFTShifter fftL,fftR;


struct Resampler2x
{
    DSP::Decimator9 dec;
    std::vector<DspFloatType> interpolated;
    std::vector<DspFloatType> decimated;

    Resampler2x(size_t n = BufferSize)
    {
        interpolated.resize(n*2);
        decimated.resize(n);
    }
    void Interp(size_t n, DspFloatType * in) {
        size_t x = 0;
        for(size_t i = 0; i < n; i++)
        {            
            DspFloatType xm1 = in[(i-1) < 0? 0 : i-1];
            DspFloatType x   = in[i];
            DspFloatType x1  = in[(i+1) > n? n : i+1];
            DspFloatType x2  = in[(i+2) > n? n : i+2];
            interpolated[x++] = in[i];
            interpolated[x++] = cubic_interpolate<DspFloatType>(0.5,xm1,x,x1,x2);
        }
    }
    void Decimate(size_t n, DspFloatType * in) {
        size_t x = 0;
        for(size_t i = 0; i < n; i+=2)
        {
            DspFloatType x = dec.Calc(in[i],in[i+1]);
            decimated[x++] = x;
        }
    }
};

Resampler2x resampler;


/////////////////////////////////////////////
// AMSynth
/////////////////////////////////////////////
#include "Synthesizer/AMSynthEffects.hpp"
#include "Synthesizer/AMSynth.hpp"

Synth::AMSynth::Distortion distortion;
Synth::AMSynth::SoftLimiter softlimit;
    



int audio_callback(const void *inputBuffer, void *outputBuffer,
                   unsigned long framesPerBuffer,
                   const PaStreamCallbackTimeInfo *timeInfo,
                   PaStreamCallbackFlags statusFlags,
                   void *userData)
{
    float *inputs = (float*)inputBuffer;
    float *output = (float*)outputBuffer;
    float *p;

    memset(temp2[0], 0, BufferSize * sizeof(DspFloatType));
    memset(temp2[1], 0, BufferSize * sizeof(DspFloatType));

    
    DspFloatType pan = 0.5;
    DspFloatType gain = pow(10.0, 3.0 / 20.0f);
    DspFloatType dist = pow(10.0, 3 / 20.0);

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
        DspFloatType fcv1 = glide1.Tick(Freq + osc1_tune);
        DspFloatType fcv2 = glide2.Tick(Freq + osc2_tune);
        DspFloatType frq1 = MusicFunctions::cv2freq(fcv1);
        DspFloatType frq2 = MusicFunctions::cv2freq(fcv2);
        DspFloatType fcutq = cutoff.Tick(Fcutoff);
        DspFloatType fresq = resonance.Tick(Qn);
        DspFloatType fq = neuron.Tick(Q);

        //bqdelay.setCutoff(frq1);
        //bqdelay.setQ(fq);

        vco1.setFrequency(frq1);
        vcf1.setCutoff(fcutq);
        vcf1.setResonance(fresq);

        blp.setCutoff(Fc+Kc);
        blp.setQ(fq);

        //c1lp.setCutoff(Fc);
        //c1lp.setQ(fq);
        //c2lp.setCutoff(Fc);
        //c2lp.setQ(fq);

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

        DspFloatType e1 = adsr.Tick();
        DspFloatType e2 = adsr2.Tick();
        DspFloatType l1 = lfo1.Tick();
        DspFloatType l2 = lfo2.Tick();
        DspFloatType l3 = lfo3.Tick();
        DspFloatType l4 = lfo4.Tick();

        DspFloatType x2 = Vel * vco1.Tick();
        x2 = vcf1.Tick(x2,e1,l1,l2);

        DspFloatType tick = vca1.Tick(x2, gain, -0.9 + 0.1 * l3, 0.9 + 0.1);

        // tick = ndelay.Tick(tick);
        temp1[0][i] = tick * std::sin((1 - pan) * M_PI / 2);
        temp1[1][i] = tick * std::cos(pan * M_PI / 2);

        // temp1[0][i] = ce2L.Tick(temp1[0][i]);
        // temp1[1][i] = ce2R.Tick(temp1[1][i]);
    }
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
Lua::LuaJIT *lua;


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
    lua = new Lua::LuaJIT("main.lua");
    lua->CreateCFunction("randd", RandomDist);
}

void repl()
{
    std::string cmd;
    std::cin >> cmd;
    lua->DoCmd(cmd);
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
    stk::Stk::setRawwavePath("Data/rawwaves/");

    distortion.SetCrunch(1.25);
    softlimit.SetSampleRate(sampleRate);

    temp1 = (DspFloatType **)calloc(2, sizeof(DspFloatType *));
    temp2 = (DspFloatType **)calloc(2, sizeof(DspFloatType *));
    for (size_t i = 0; i < 2; i++)
    {
        temp1[i] = (DspFloatType *)calloc(BufferSize, sizeof(DspFloatType));
        temp2[i] = (DspFloatType *)calloc(BufferSize, sizeof(DspFloatType));
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
    int pulse = 6;    
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