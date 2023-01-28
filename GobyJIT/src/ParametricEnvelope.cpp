//
// Created by eddie on 2020-10-29.
//
#include <cstdlib>
#include "ParametricEnvelope.h"
#include <math.h>

using namespace SoundWave;

ParametricEnvelope::ParametricEnvelope(double maxLevel, unsigned int stepFrequency) {
    this->maxLevel = maxLevel;
    this->stepFrequency = stepFrequency;
    reset();
}

void ParametricEnvelope::reset() {
    currentState = idle;
    previousState = idle;
    currentValue = 0.0;
    currentStep = 0;
}

void ParametricEnvelope::setAttackTime(double time) {
    parameters.attackTime = time;
}

void ParametricEnvelope::setAttackSlope(double slope) {
    parameters.attackSlope = slope;
}

void ParametricEnvelope::setDecayTime(double time) {
    parameters.decayTime = time;
}

void ParametricEnvelope::setDecaySlope(double slope) {
    parameters.decaySlope = slope;
}

void ParametricEnvelope::setReleaseTime(double time) {
    parameters.releaseTime = time;
}

void ParametricEnvelope::setReleaseSlope(double slope) {
    parameters.releaseSlope = slope;
}

void ParametricEnvelope::setSustainLevel(double level) {
    if (level>maxLevel) {
        parameters.sustainLevel = 1.0;
    } else {
        parameters.sustainLevel = level/maxLevel;
    }
}

parameter_struct ParametricEnvelope::getParameters() {
    return parameters;
}

state_enum ParametricEnvelope::getState() {
    return currentState;
}

//
// step is called periodically to update the inner state and the currentValue of this envelope. It returns the current value
double ParametricEnvelope::step() {
    //TODO Implement value update
    switch (currentState) {
        case idle:
            break;
        case attack:
            if (previousState == idle) {
                currentStep = 0;
                previousState = attack;
            }
            currentValue = calculateAttackValue(currentStep, parameters.attackTime, parameters.attackSlope);
            if (currentValue >= 1.0) {
                gotoState(decay);
            } else {
                currentStep++;
            }
            break;
        case decay:
            if (previousState != decay) {
                currentStep = 0;
                previousState = decay;
            }
            currentValue = calculateDecayValue(currentStep, parameters.decayTime, parameters.decaySlope,
                                               parameters.sustainLevel);
            if (currentValue <= parameters.sustainLevel) {
                gotoState(sustain);
            } else {
                currentStep++;
            }
            break;
        case sustain:
            if (previousState != sustain) {
                previousState = sustain;
                currentValue = parameters.sustainLevel;
            }
            break;
        case release:
            if (previousState != release) {
                currentStep = 0;
                previousState = release;
                releaseInitialValue = currentValue;
            }
            currentValue = calculateReleaseValue(currentStep, parameters.releaseTime, parameters.releaseSlope,
                                                 releaseInitialValue);
            if (currentValue < 0.0) {
                currentValue = 0.0;
                gotoState(idle);
            } else {
                currentStep++;
            }
    }
    return currentValue * maxLevel;
}

double ParametricEnvelope::getCurrentTime() {
    return currentStep / (double)stepFrequency;
}

void ParametricEnvelope::onGateOn() {
    gotoState(attack);
}

void ParametricEnvelope::onGateOff() {
    gotoState(release);
}

void ParametricEnvelope::gotoState(state_enum newState) {
    previousState = currentState;
    currentState = newState;
}

double ParametricEnvelope::calculateAttackValue(unsigned long currentStep, double time, double slope) {
    return std::pow(getCurrentTime() / time, std::pow(2.0, -slope));
}

double ParametricEnvelope::calculateDecayValue(unsigned long currentStep, double time, double slope, double targetLevel) {
    return std::pow(getCurrentTime() / time, std::pow(2.0, -slope)) * (targetLevel - 1.0) + 1.0;
}

double ParametricEnvelope::calculateReleaseValue(unsigned long currentStep, double time, double slope, double originLevel) {
    return originLevel * ( 1- std::pow(getCurrentTime() / time, std::pow(2.0, -slope))) ;
}
