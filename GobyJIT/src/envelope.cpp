#include "SoundWave.h


ADSR :: ADSR( void )
{
  target_ = 0.0;
  value_ = 0.0;
  attackRate_ = 0.001;
  decayRate_ = 0.001;
  releaseRate_ = 0.005;
  releaseTime_ = -1.0;
  sustainLevel_ = 0.5;
  state_ = IDLE;
  Stk::addSampleRateAlert( this );
}

ADSR :: ~ADSR( void )
{
  Stk::removeSampleRateAlert( this );
}

void ADSR :: sampleRateChanged( float newRate, float oldRate )
{
  if ( !ignoreSampleRateChange_ ) {
    attackRate_ = oldRate * attackRate_ / newRate;
    decayRate_ = oldRate * decayRate_ / newRate;
    releaseRate_ = oldRate * releaseRate_ / newRate;
  }
}

void ADSR :: keyOn()
{
  if ( target_ <= 0.0 ) target_ = 1.0;
  state_ = ATTACK;
}

void ADSR :: keyOff()
{
  target_ = 0.0;
  state_ = RELEASE;

  // FIXED October 2010 - Nick Donaldson
  // Need to make release rate relative to current value!!
  // Only update if we have set a TIME rather than a RATE,
  // in which case releaseTime_ will be -1
  if ( releaseTime_ > 0.0 )
	  releaseRate_ = value_ / ( releaseTime_ * Stk::sampleRate() );
}

void ADSR :: setAttackRate( float rate )
{
  if ( rate < 0.0 ) {
    oStream_ << "ADSR::setAttackRate: argument must be >= 0.0!";
    handleError( StkError::WARNING ); return;
  }

  attackRate_ = rate;
}

void ADSR :: setAttackTarget( float target )
{
  if ( target < 0.0 ) {
    oStream_ << "ADSR::setAttackTarget: negative target not allowed!";
    handleError( StkError::WARNING ); return;
  }

  target_ = target;
}

void ADSR :: setDecayRate( float rate )
{
  if ( rate < 0.0 ) {
    oStream_ << "ADSR::setDecayRate: negative rates not allowed!";
    handleError( StkError::WARNING ); return;
  }

  decayRate_ = rate;
}

void ADSR :: setSustainLevel( float level )
{
  if ( level < 0.0 ) {
    oStream_ << "ADSR::setSustainLevel: negative level not allowed!";
    handleError( StkError::WARNING ); return;
  }

  sustainLevel_ = level;
}

void ADSR :: setReleaseRate( float rate )
{
  if ( rate < 0.0 ) {
    oStream_ << "ADSR::setReleaseRate: negative rates not allowed!";
    handleError( StkError::WARNING ); return;
  }

  releaseRate_ = rate;

  // Set to negative value so we don't update the release rate on keyOff()
  releaseTime_ = -1.0;
}

void ADSR :: setAttackTime( float time )
{
  if ( time <= 0.0 ) {
    oStream_ << "ADSR::setAttackTime: negative or zero times not allowed!";
    handleError( StkError::WARNING ); return;
  }

  attackRate_ = 1.0 / ( time * Stk::sampleRate() );
}

void ADSR :: setDecayTime( float time )
{
  if ( time <= 0.0 ) {
    oStream_ << "ADSR::setDecayTime: negative or zero times not allowed!";
    handleError( StkError::WARNING ); return;
  }

  decayRate_ = (1.0 - sustainLevel_) / ( time * Stk::sampleRate() );
}

void ADSR :: setReleaseTime( float time )
{
  if ( time <= 0.0 ) {
    oStream_ << "ADSR::setReleaseTime: negative or zero times not allowed!";
    handleError( StkError::WARNING ); return;
  }

  releaseRate_ = sustainLevel_ / ( time * Stk::sampleRate() );
  releaseTime_ = time;
}

void ADSR :: setAllTimes( float aTime, float dTime, float sLevel, float rTime )
{
  this->setAttackTime( aTime );
  this->setSustainLevel( sLevel );
  this->setDecayTime( dTime );
  this->setReleaseTime( rTime );
}

void ADSR :: setTarget( float target )
{
  if ( target < 0.0 ) {
    oStream_ << "ADSR::setTarget: negative target not allowed!";
    handleError( StkError::WARNING ); return;
  }

  target_ = target;

  this->setSustainLevel( target_ );
  if ( value_ < target_ ) state_ = ATTACK;
  if ( value_ > target_ ) state_ = DECAY;
}

void ADSR :: setValue( float value )
{
  state_ = SUSTAIN;
  target_ = value;
  value_ = value;
  this->setSustainLevel( value );
  lastFrame_[0] = value;
}


float ADSR::process() {
    switch (state) {
        case env_idle:
            break;
        case env_attack:
            output = attackBase + output * attackCoef;
            if (output >= 1.0) {
                output = 1.0;
                state = env_decay;
            }
            break;
        case env_decay:
            output = decayBase + output * decayCoef;
            if (output <= sustainLevel) {
                output = sustainLevel;
                state = env_sustain;
            }
            break;
        case env_sustain:
            break;
        case env_release:
            output = releaseBase + output * releaseCoef;
            if (output <= 0.0) {
                output = 0.0;
                state = env_idle;
            }
    }
    return output;
}

void ADSR::gate(int gate) {
    if (gate)
        state = env_attack;
    else if (state != env_idle)
        state = env_release;
}

int ADSR::getState() {
    return state;
}

void ADSR::reset() {
    state = env_idle;
    output = 0.0;
}

float ADSR::getOutput() {
    return output;
}


}

//
// Created by eddie on 2020-10-29.
//

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
