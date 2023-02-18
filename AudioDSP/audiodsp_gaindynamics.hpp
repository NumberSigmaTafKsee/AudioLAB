// https://github.com/buosseph/JuceCompressor
#include <math.h>

class GainDynamics {
public:
    GainDynamics(float sampleRate, float attackTime, float releaseTime);
    ~GainDynamics();
    
    float tick(float inputSample);
    void setDetector(float sampleRate);
    void setAttack(float attackTime);
    void setRelease(float releaseTime);
    
private:
    float fs, outputGain;
    double b0Attack, b0Release, b0;
    float attackTime, releaseTime;     // in seconds
};

class PeakLevelDetector {
public:
    PeakLevelDetector(float sampleRate);
    ~PeakLevelDetector();
    
    float tick(float inputSample);
    void setDetector(float sampleRate);
    
private:
    float fs, inputAbs, peakOutput;
    float b0Attack, b0Release, b0, a1;
    float releaseTime = 0.100f;     // seconds
};

// Times are in seconds (e.g. 100ms = 0.1f, 1.2s = 1.2f)
GainDynamics::GainDynamics(float sampleRate, float newAttackTime, float newReleaseTime) {
    attackTime = newAttackTime;
    releaseTime = newReleaseTime;
    setDetector(sampleRate);
}

GainDynamics::~GainDynamics() {}

float GainDynamics::tick(float inputGain) {
    if (inputGain < outputGain) {   // Isn't this suppose to be (input > lastOutput)?
        b0 = b0Attack;
    }
    else {
        b0 = b0Release;
    }
    
    // Simplified filter equation (out = b0 * input + a1 * lastOut)
    outputGain += b0 * (inputGain - outputGain);
    
    return outputGain;
}

void GainDynamics::setDetector(float sampleRate) {
    fs = sampleRate;
    outputGain = 0.f;
    setAttack(attackTime);
    setRelease(releaseTime);
}

void GainDynamics::setAttack(float newAttackTime) {
    attackTime = newAttackTime;
    b0Attack = 1. - expf(-1. / (attackTime * fs));;
}

void GainDynamics::setRelease(float newReleaseTime) {
    releaseTime = newReleaseTime;
    b0Release = 1. - expf(-1. / (releaseTime * fs));;
}

PeakLevelDetector::PeakLevelDetector(float sampleRate) {
    setDetector(sampleRate);
}

PeakLevelDetector::~PeakLevelDetector() {}

float PeakLevelDetector::tick(float inputSample) {
    inputAbs = fabs(inputSample);
    
    if (inputAbs > peakOutput) {
        b0 = b0Attack;
    }
    else {
        b0 = b0Release;
    }
    
    // Simplified filter equation (out = b0 * input + a1 * lastOut)
    peakOutput += b0 * (inputAbs - peakOutput);
    
    return peakOutput;
}

void PeakLevelDetector::setDetector(float sampleRate) {
    fs = sampleRate;
    peakOutput = 0.f;
    
    // set coefficients for leaky integrator
    b0Attack = 1.f;
    a1 = expf(-1 / (releaseTime * fs));
    b0Release = 1.f - a1;
}
