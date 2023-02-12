//
// Created by eddie on 2020-10-29.
// The Truth by Macho Charli 5
//
#pragma once 

#include "SoundObject.hpp"

namespace Envelopes
{
    enum state_enum {
        idle, attack, decay, sustain, release
    };


    struct parameter_struct {
        DspFloatType attackTime = 1.0;
        DspFloatType attackSlope = 0.0;
        DspFloatType decayTime = 1.0;
        DspFloatType decaySlope = 0.0;
        DspFloatType sustainLevel = 0.5;
        DspFloatType releaseTime = 1.0;
        DspFloatType releaseSlope = 0.0;
    };


    // todo: add ports
    class ADSR2 : public GeneratorProcessor {
    private:
        DspFloatType currentValue;
        DspFloatType releaseInitialValue;
        unsigned long currentStep;
        state_enum currentState;
        state_enum previousState;

        // Amplitude scale
        DspFloatType maxLevel;
        // Time scale. Essentially, how often is step() called?
        unsigned int stepFrequency;
        parameter_struct parameters;

        void gotoState(state_enum newState);

        DspFloatType calculateAttackValue(unsigned long currentStep, DspFloatType time, DspFloatType slope);

        DspFloatType calculateDecayValue(unsigned long currentStep, DspFloatType time, DspFloatType slope, DspFloatType targetLevel);

        DspFloatType calculateReleaseValue(unsigned long currentStep, DspFloatType time, DspFloatType slope, DspFloatType originLevel);

    public:
        ADSR2(DspFloatType maxLevel = 10.0, unsigned int stepFrequency = 10000);

        //
        // Setters for all envelope parameters, with reasonable defaults
        //
        void setAttackTime(DspFloatType time);

        void setAttackSlope(DspFloatType slope);

        void setDecayTime(DspFloatType time);

        void setDecaySlope(DspFloatType slope);

        void setSustainLevel(DspFloatType level);

        void setReleaseTime(DspFloatType time);

        void setReleaseSlope(DspFloatType slope);

        //
        // Get current time _within the current phase_
        // NB This gets reset at the beginning of the A, D and R phases and isn't updated during S
        DspFloatType getCurrentTime();

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
        DspFloatType step();

        //
        // Reset the envelope
        //
        void reset();

        enum {
            PORT_ATK,
            PORT_ATKSLOPE,
            PORT_DCY,
            PORT_DCYSLOPE,
            PORT_SUS,
            PORT_REL,
            PORT_RELSLOPE,
            PORT_RESET,
            PORT_GATEON,
            PORT_GATEOFF,
        };
        void printPorts() {
            printf("PORTS\n");
            printf("PORT_ATK=0\n");
            printf("PORT_ATKSLOPE=1\n");
            printf("PORT_DCY=2\n");
            printf("PORT_DCYSLOPE=3\n");
            printf("PORT_SUS=4\n");
            printf("PORT_REL=5\n");
            printf("PORT_RELSLOPE=6\n");
            printf("PORT_RESET=7\n");
            printf("PORT_GATEON=8\n");
            printf("PORT_GATEOFF=9\n");
        }
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_ATK: setAttackTime(v); break;
                case PORT_ATKSLOPE: setAttackSlope(v); break;
                case PORT_DCY: setDecayTime(v); break;
                case PORT_DCYSLOPE: setDecaySlope(v); break;
                case PORT_SUS: setSustainLevel(v); break;
                case PORT_REL: setReleaseTime(v); break;
                case PORT_RELSLOPE: setReleaseSlope(v); break;
                case PORT_RESET: if(v != 0.0) reset(); break;
                case PORT_GATEON: if(v != 0.0) onGateOn(); break;
                case PORT_GATEOFF: if(v != 0.0) onGateOff(); break;
            }
        }

        DspFloatType Tick(DspFloatType I=1, DspFloatType A=1, DspFloatType X=0, DspFloatType Y=0)
        {
            return I*step();
        }
        void ProcessBlock(size_t n, DspFloatType * input, DspFloatType * output) {
            #pragma omp simd
            for(size_t i = 0; i < n; i++) output[i] = input[i]*step();
        }
        void ProcessInplace(DspFloatType * input, size_t n) {
            #pragma omp simd
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


    inline ADSR2::ADSR2(DspFloatType maxLevel, unsigned int stepFrequency) 
    : GeneratorProcessor()
    {
        this->maxLevel = maxLevel;
        this->stepFrequency = stepFrequency;
        reset();
    }

    inline void ADSR2::reset() {
        currentState = idle;
        previousState = idle;
        currentValue = 0.0;
        currentStep = 0;
    }

    inline void ADSR2::setAttackTime(DspFloatType time) {
        parameters.attackTime = time;
    }

    inline void ADSR2::setAttackSlope(DspFloatType slope) {
        parameters.attackSlope = slope;
    }

    inline void ADSR2::setDecayTime(DspFloatType time) {
        parameters.decayTime = time;
    }

    inline void ADSR2::setDecaySlope(DspFloatType slope) {
        parameters.decaySlope = slope;
    }

    inline void ADSR2::setReleaseTime(DspFloatType time) {
        parameters.releaseTime = time;
    }

    inline void ADSR2::setReleaseSlope(DspFloatType slope) {
        parameters.releaseSlope = slope;
    }

    inline void ADSR2::setSustainLevel(DspFloatType level) {
        if (level>maxLevel) {
            parameters.sustainLevel = 1.0;
        } else {
            parameters.sustainLevel = level/maxLevel;
        }
    }

    parameter_struct ADSR2::getParameters() {
        return parameters;
    }

    state_enum ADSR2::getState() {
        return currentState;
    }

    //
    // step is called periodically to update the inner state and the currentValue of this envelope. It returns the current value
    inline DspFloatType ADSR2::step() {
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

    inline DspFloatType ADSR2::getCurrentTime() {
        return currentStep / (DspFloatType)stepFrequency;
    }

    inline void ADSR2::onGateOn() {
        gotoState(attack);
    }

    inline void ADSR2::onGateOff() {
        gotoState(release);
    }

    inline void ADSR2::gotoState(state_enum newState) {
        previousState = currentState;
        currentState = newState;
    }

    inline DspFloatType ADSR2::calculateAttackValue(unsigned long currentStep, DspFloatType time, DspFloatType slope) {
        return std::pow(getCurrentTime() / time, std::pow(2.0, -slope));
    }

    inline DspFloatType ADSR2::calculateDecayValue(unsigned long currentStep, DspFloatType time, DspFloatType slope, DspFloatType targetLevel) {
        return std::pow(getCurrentTime() / time, std::pow(2.0, -slope)) * (targetLevel - 1.0) + 1.0;
    }

    inline DspFloatType ADSR2::calculateReleaseValue(unsigned long currentStep, DspFloatType time, DspFloatType slope, DspFloatType originLevel) {
        return originLevel * ( 1- std::pow(getCurrentTime() / time, std::pow(2.0, -slope))) ;
    }
}