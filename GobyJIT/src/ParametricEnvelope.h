//
// Created by eddie on 2020-10-29.
//
#pragma once 

namespace SoundWave {

enum state_enum {
    idle, attack, decay, sustain, release
};

struct parameter_struct {
    double attackTime = 1.0;
    double attackSlope = 0.0;
    double decayTime = 1.0;
    double decaySlope = 0.0;
    double sustainLevel = 0.5;
    double releaseTime = 1.0;
    double releaseSlope = 0.0;
};


class ParametricEnvelope {
private:
    double currentValue;
    double releaseInitialValue;
    unsigned long currentStep;
    state_enum currentState;
    state_enum previousState;

    // Amplitude scale
    double maxLevel;
    // Time scale. Essentially, how often is step() called?
    unsigned int stepFrequency;
    parameter_struct parameters;

    void gotoState(state_enum newState);

    double calculateAttackValue(unsigned long currentStep, double time, double slope);

    double calculateDecayValue(unsigned long currentStep, double time, double slope, double targetLevel);

    double calculateReleaseValue(unsigned long currentStep, double time, double slope, double originLevel);

public:
    ParametricEnvelope(double maxLevel = 10.0, unsigned int stepFrequency = 10000);

    //
    // Setters for all envelope parameters, with reasonable defaults
    //
    void setAttackTime(double time);

    void setAttackSlope(double slope);

    void setDecayTime(double time);

    void setDecaySlope(double slope);

    void setSustainLevel(double level);

    void setReleaseTime(double time);

    void setReleaseSlope(double slope);

    //
    // Get current time _within the current phase_
    // NB This gets reset at the beginning of the A, D and R phases and isn't updated during S
    double getCurrentTime();

    //
    // Even handlers for gate transitions
    // onGateOn will transition the envelop into STATE_ATTACK
    //
    void onGateOn();

    //
    // onGateOff will transition the envelope into STATE_RELEASE
    //
    void onGateOff();

    //
    // step is called at regular intervals to get the current signal value
    //
    double step();

    //
    // Reset the envelope
    //
    void reset();

    float Tick() {
        return step();
    }
    float Tick(float in) {
        return in*step();
    }

    void Process(size_t n, float * input, float * output) {
        for(size_t i = 0; i < n; i++) output[i] = input[i]*step();
    }
    void Process(float * input, size_t n) {
        for(size_t i = 0; i < n; i++) input[i] = input[i]*step();
    }
    //
    // Get the current state.
    //
    state_enum getState();

    //
    // Return current parameters
    parameter_struct getParameters();
};


}