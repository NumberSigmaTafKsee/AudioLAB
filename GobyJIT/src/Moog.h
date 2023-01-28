#pragma once

#include <cassert>


#include "Filters.h"
#include "NoiseGenerator.h"
#include "LadderFilterBase.h"
#include "HuovilainenModel.h"
#include "ImprovedModel.h"
#include "KrajeskiModel.h"
#include "MicrotrackerModel.h"
#include "MusicDSPModel.h"
#include "OberheimVariationModel.h"
#include "RKSimulationModel.h"
#include "SimplifiedModel.h"
#include "StilsonModel.h"
#include "RingBuffer.h"
#include "std.h"

namespace SoundWave {
float clamp(float x, float a, float b) {
    return x < a? a: x > b? b: x;
}
enum FilterType {
    LOWPASS = ::RBJFilter::FilterType::LOWPASS,
    HIGHPASS = ::RBJFilter::FilterType::HIGHPASS,
    BANDPASS = ::RBJFilter::FilterType::BANDPASS,
    ALLPASS = ::RBJFilter::FilterType::ALLPASS,
    NOTCH = ::RBJFilter::FilterType::NOTCH,
    PEAK = ::RBJFilter::FilterType::PEAK,
    LOW_SHELF = ::RBJFilter::FilterType::LOW_SHELF,
    HIGH_SHELF = ::RBJFilter::FilterType::HIGH_SHELF,
};

struct Filter
{
    ::RBJFilter * filter;

    Filter(FilterType type, float cutoff, float sample_rate) {
        filter = new ::RBJFilter((::RBJFilter::FilterType)type,cutoff,sample_rate);
        assert(filter != nullptr);
    }
    ~Filter() {
        if(filter) delete filter;
    }

    void UpdateCoefficients() { filter->UpdateCoefficients(); }
    void SetCutoff(float c) { filter->SetCutoff(c); }
    void SetQ(float q) { filter->SetQValue(q); }
    float GetCutoff() { return filter->GetCutoff(); }
    float GetQ() { return filter->GetQValue(); }

    void Process(float * samples, const uint32_t n) { filter->Process(samples,n); }
    float Tick(float s) { return filter->Tick(s); }
};

struct RBJLowPassFilter : public Filter
{
    RBJLowPassFilter(float cutoff, float sample_rate=44100) : Filter(LOWPASS,cutoff,sample_rate) {

    }
};

struct RBJHighPassFilter : public Filter
{
    RBJHighPassFilter(float cutoff, float sample_rate=44100) : Filter(HIGHPASS,cutoff,sample_rate) {

    }
};

struct RBJBandPassFilter : public Filter
{
    RBJBandPassFilter(float cutoff, float sample_rate=44100) : Filter(BANDPASS,cutoff,sample_rate) {

    }
};
struct RBJAllPassFilter : public Filter
{
    RBJAllPassFilter(float cutoff, float sample_rate=44100) : Filter(ALLPASS,cutoff,sample_rate) {

    }
};

struct RBJNotchFilter : public Filter
{
    RBJNotchFilter(float cutoff, float sample_rate=44100) : Filter(NOTCH,cutoff,sample_rate) {

    }
};

struct RBJPeakFilter : public Filter
{
    RBJPeakFilter(float cutoff, float sample_rate=44100) : Filter(PEAK,cutoff,sample_rate) {

    }
};

struct RBJLowShelfFilter : public Filter
{
    RBJLowShelfFilter(float cutoff, float sample_rate=44100) : Filter(LOW_SHELF,cutoff,sample_rate) {

    }
};

struct RBJHighShelfFilter : public Filter
{
    RBJHighShelfFilter(float cutoff, float sample_rate=44100) : Filter(HIGH_SHELF,cutoff,sample_rate) {

    }
};


struct CRingBuffer
{
    RingBuffer *ring;

    CRingBuffer(size_t n) {
        ring = new RingBuffer(n);
        assert(ring != nullptr);
    }
    ~CRingBuffer() { if(ring != nullptr) delete ring; }

    void resize(size_t n) { ring->resize(n); }
    void clear() { ring->clear(); }
    size_t getSize() const { return ring->getSize(); }
    size_t getAvailableWrite() const { return ring->getAvailableWrite(); }
    size_t getAvailableRead() const { return ring->getAvailableRead(); }

    bool write(const float * array, size_t count) { return ring->write(array,count); }
    bool read(float * array, size_t count) { return ring->read(array,count); }

};

enum NoiseType
	{
		WHITE,
		PINK,
		BROWN,
	};

struct NoiseSamples
{
    ::NoiseGenerator * noise;
    NoiseType type;
    int sampleRate;
    int channels;

    NoiseSamples(NoiseType type, int sampleRate, int channels) {
        noise = new ::NoiseGenerator();
        assert(noise != nullptr);
    }
    ~NoiseSamples() {
        if(noise) delete noise;
    }

    Std::StdVector<float> produce(float seconds) {
        Std::StdVector<float> r;
        r.vector = noise->produce((::NoiseGenerator::NoiseType)type,sampleRate,channels,seconds);
        return r;
    }
};

enum MoogModelType
{
    FINN_MOOG,
    IMPROVED_MOOG,
    POLISH_MOOG,
    MICROTRACKER_MOOG,
    MUSICDSP_MOOG,
    OBERHEIM_MOOG,
    RK_MOOG,
    SIMPLIFIED_MOOG,
    STILSON_MOOG
};

struct MoogLadderFilter
{
    LadderFilterBase * moog;
    float sr;

    MoogLadderFilter(MoogModelType type, float sample_rate=44100) {
        sr = sample_rate;
        switch(type) {
        case FINN_MOOG: moog = new HuovilainenMoog(sample_rate); break;
        case IMPROVED_MOOG: moog = new ImprovedMoog(sample_rate); break;
        case POLISH_MOOG: moog = new KrajeskiMoog(sample_rate); break;
        case MICROTRACKER_MOOG: moog = new MicrotrackerMoog(sample_rate); break;
        case MUSICDSP_MOOG: moog = new MusicDSPMoog(sample_rate); break;
        case OBERHEIM_MOOG: moog = new OberheimVariationMoog(sample_rate); break;
        case RK_MOOG: moog = new RKSimulationMoog(sample_rate); break;
        case STILSON_MOOG: moog = new StilsonMoog(sample_rate); break;
        case SIMPLIFIED_MOOG: moog = new SimplifiedMoog(sample_rate); break;
        }
        assert(moog != nullptr);
    }
    ~MoogLadderFilter() {
        if(moog) delete moog;
    }
    void Process(float * samples, uint32_t n) {
        moog->Process(samples,n);
    }
    void Process(float * samples, float * modulation, uint32_t n) {
        moog->Process(samples,modulation, n);
    }
    void SetResonance(float r) {
        moog->SetResonance(r);
    }
    void SetCutoff(float c) {
        moog->SetCutoff(c);
    }
    float Tick(float I,float A=1, float X=0, float Y=0) {
        float o = I;
        float c = moog->GetCutoff();
        float r = moog->GetResonance();
        SetCutoff(clamp(c + c*(X*0.5),0,sr/2));
        SetResonance(clamp(r + Y,0,1));
        Process(&o,1);
        SetCutoff(c);
        SetResonance(r);
        return o;
    }
};

};
