#pragma once

#include <cstdlib>
#include <math.h>

#include "SABase.hpp"

namespace SoundAlchemy
{
    // todo: exponential one pole filter
    template<typename T>
    class TADSR : public TObject<T>  {
    public:
        TADSR<T>(void);
        ~TADSR<T>(void);
        T process(void);
        T getOutput(void);
        int getState(void);
        void gate(int on);

        void setSustainLevel(T level);

        void setAttackTime(T ms) { attackRate = ms/this->sample_rate; }
        void setDecayTime(T ms) { decayRate = ms/this->sample_rate; }
        void setReleaseTime(T ms) { releaseRate = ms/this->sample_rate; }

        void setAllTimes(T attack, T decay, T sustain, T release) {
            setAttackTime(attack);
            setDecayTime(decay);
            setSustainLevel(sustain);
            setReleaseTime(release);
        }

        void setAttackRate(T rate);
        void setDecayRate(T rate);
        void setReleaseRate(T rate);

        void setAllRates(T a, T d, T s, T r) {
            setAttackRate(a);
            setDecayRate(d);
            setSustainLevel(s);
            setReleaseRate(r);
        }
        
        void setTargetRatioA(T targetRatio);
        void setTargetRatioDR(T targetRatio);
        
        void reset(void);

        T Tick() {
            return this->output_gain*process();
        }
        T Tick(T I, T A = 1, T X = 1, T Y = 1) {
            I *= this->input_gain;
            T x = std::abs(X);
            T y = -std::abs(Y);
            T r = this->output_gain*(A*I*X*Y*Tick());
            return clamp(output_gain*r,this->minA,this->maxA);
        }
                
        enum envState {
            env_idle = 0,
            env_attack,
            env_decay,
            env_sustain,
            env_release
        };

    protected:
        int state;
        T output;
        T attackRate;
        T decayRate;
        T releaseRate;
        T attackCoef;
        T decayCoef;
        T releaseCoef;
        T sustainLevel;
        T targetRatioA;
        T targetRatioDR;
        T attackBase;
        T decayBase;
        T releaseBase;

        T calcCoef(T rate, T targetRatio);
    };

    template<typename T>
    inline T TADSR<T>::process() {
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

    template<typename T>
    inline void TADSR<T>::gate(int gate) {
        if (gate)
            state = env_attack;
        else if (state != env_idle)
            state = env_release;
    }

    template<typename T>
    inline int TADSR<T>::getState() {
        return state;
    }

    template<typename T>
    inline void TADSR<T>::reset() {
        state = env_idle;
        output = 0.0;
    }

    template<typename T>
    inline T TADSR<T>::getOutput() {
        return output;
    }


    template<typename T>
    TADSR<T>::TADSR(void) {
        reset();
        setAttackRate(0);
        setDecayRate(0);
        setReleaseRate(0);
        setSustainLevel(1.0);
        setTargetRatioA(0.3);
        setTargetRatioDR(0.0001);
    }

    template<typename T>
    TADSR<T>::~TADSR(void) {
    }

    template<typename T>
    void TADSR<T>::setAttackRate(T rate) {
        attackRate = rate;
        attackCoef = calcCoef(rate, targetRatioA);
        attackBase = (1.0 + targetRatioA) * (1.0 - attackCoef);
    }

    template<typename T>
    void TADSR<T>::setDecayRate(T rate) {
        decayRate = rate;
        decayCoef = calcCoef(rate, targetRatioDR);
        decayBase = (sustainLevel - targetRatioDR) * (1.0 - decayCoef);
    }

    template<typename T>
    void TADSR<T>::setReleaseRate(T rate) {
        releaseRate = rate;
        releaseCoef = calcCoef(rate, targetRatioDR);
        releaseBase = -targetRatioDR * (1.0 - releaseCoef);
    }

    template<typename T>
    T TADSR<T>::calcCoef(T rate, T targetRatio) {
        return (rate <= 0) ? 0.0 : std::exp(-std::log((1.0 + targetRatio) / targetRatio) / rate);
    }

    template<typename T>
    void TADSR<T>::setSustainLevel(T level) {
        sustainLevel = level;
        decayBase = (sustainLevel - targetRatioDR) * (1.0 - decayCoef);
    }

    template<typename T>
    void TADSR<T>::setTargetRatioA(T targetRatio) {
        if (targetRatio < 0.000000001)
            targetRatio = 0.000000001;  // -180 dB
        targetRatioA = targetRatio;
        attackCoef = calcCoef(attackRate, targetRatioA);
        attackBase = (1.0 + targetRatioA) * (1.0 - attackCoef);
    }

    template<typename T>
    void TADSR<T>::setTargetRatioDR(T targetRatio) {
        if (targetRatio < 0.000000001)
            targetRatio = 0.000000001;  // -180 dB
        targetRatioDR = targetRatio;
        decayCoef = calcCoef(decayRate, targetRatioDR);
        releaseCoef = calcCoef(releaseRate, targetRatioDR);
        decayBase = (sustainLevel - targetRatioDR) * (1.0 - decayCoef);
        releaseBase = -targetRatioDR * (1.0 - releaseCoef);
    }

}
