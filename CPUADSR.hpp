#pragma once

#include "CPUSoundObject.hpp"


namespace CPU::Envelopes
{
    // todo: add ports
    template<typename DSPType=Vector4f, Float=float>
    class ADSR : public CPUGeneratorProcessor<T> {
    public:
        
        Float min=0.0;
        Float max=1.0;
        
        ADSR(Float sample_rate=44100.0f) : CPUGeneratorProcessor<T>()
        {        
            sr = sample_rate;
            reset();
            setAttackRate(0);
            setDecayRate(0);
            setReleaseRate(0);
            setSustainLevel(1.0);
            setTargetRatioA(0.3);
            setTargetRatioDR(0.0001);        
        }
        ~ADSR(void) {

        }

        void setVoices(SIMD v) {
            voices = v;
        }

        ADSR(Float a, Float d, Float s, Float r, Float sample_rate=44100.0f, Float cutoff=10.0f) 
        : GeneratorProcessor()
        {     
            sr = sample_rate;
            reset();
            setAttackRate(a * sample_rate);
            setDecayRate(d * sample_rate);
            setSustainLevel(s);
            setReleaseRate(r * sample_rate);                
            setTargetRatioA(0.3);
            setTargetRatioDR(0.0001);
        }
        
        SIMD process(void);
        SIMD getOutput(void);
        int getState(void);
        void gate(int on);
        
        void setAllTimes(Float a, Float d, Float s, Float r) {
            reset();
            setAttackTime(a);
            setDecayTime(d);
            setSustainLevel(s);
            setReleaseTime(r);
        }

        enum {
            PORT_ATKTIME,
            PORT_DCYTIME,
            PORT_RELTIME,
            PORT_ATKRATE,
            PORT_DCYRATE,
            PORT_RELRATE,
            PORT_SUSTAIN,
            PORT_RATIOA,
            PORT_RATIODR,
            PORT_RESET,
            PORT_NOTEON,
            PORT_NOTEOFF,
        };
        void printPorts() {
            printf("PORTS\n");
            printf("PORT_ATKTIME=0\n");
            printf("PORT_DCYTIME=1\n");
            printf("PORT_RELTIME=2\n");
            printf("PORT_ATKRATE=3\n");
            printf("PORT_DCYRATE=4\n");
            printf("PORT_RELRATE=5\n");
            printf("PORT_SUSTAIN=6\n");
            printf("PORT_RATIOA=7\n");
            printf("PORT_RADIODR=8\n");
            printf("PORT_RESET=9\n");
            printf("PORT_NOTEON=10\n");
            printf("PORT_NOTEOFF=11\n");
        }
        void setPort(int port, Float v) {
            switch(port)
            {
                case PORT_ATKTIME: setAttackTime(v); break;
                case PORT_DCYTIME: setDecayTime(v); break;
                case PORT_RELTIME: setReleaseTime(v); break;
                case PORT_ATKRATE: setAttackRate(v); break;
                case PORT_DCYRATE: setDecayRate(v); break;
                case PORT_RELRATE: setReleaseRate(v); break;
                case PORT_SUSTAIN: setSustainLevel(v); break;
                case PORT_RATIOA: setTargetRatioA(v); break;
                case PORT_RATIODR: setTargetRatioDR(v); break;
                case PORT_NOTEON: if(v != 0.0) noteOn(); break;
                case PORT_NOTEOFF: if(v != 0.0) noteOff(); break;
                case PORT_RESET: if(v != 0.0) reset(); break;                
            }
        }

        void setAttackTime(Float rate)  { setAttackRate(rate*sr);}
        void setDecayTime(Float rate)   { setDecayRate(rate*sr); }
        void setReleaseTime(Float rate) { setReleaseRate(rate*sr);}

        void setAttackRate(Float rate);
        void setDecayRate(Float rate);
        void setReleaseRate(Float rate) ;
        void setSustainLevel(Float level);

        void setTargetRatioA(Float targetRatio);
        void setTargetRatioDR(Float targetRatio);

        void noteOn() { gate(true); }
        void noteOff() { gate(false); }
        void reset(void);

        
        SIMD Tick(SIM I=1, Float A=1,Float X=1, Float Y=1) {            
            T r = process();            
            return A*r;
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
        SIMD output;
        Float attackRate;
        Float decayRate;
        Float releaseRate;
        Float attackCoef;
        Float decayCoef;
        Float releaseCoef;
        Float sustainLevel;
        Float targetRatioA;
        Float targetRatioDR;
        Float attackBase;
        Float decayBase;
        Float releaseBase;
        Float sr;
        // this controls gain on each simd channel
        SIMD voices;
        Float calcCoef(Float rate, Float targetRatio);
    };

    inline SIMD ADSR::process() {
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
        return voices*(output*(max-min) + min);    
    }

    inline void ADSR::gate(int gate) {
        if (gate)
            state = env_attack;
        else if (state != env_idle)
            state = env_release;
    }

    inline int ADSR::getState() {
        return state;
    }

    inline void ADSR::reset() {
        state = env_idle;
        output = 0.0;
    }

    inline SIMD ADSR::getOutput() {
        return output;
    }

    inline void ADSR::setAttackRate(Float rate) {
        attackRate = rate;
        attackCoef = calcCoef(rate, targetRatioA);
        attackBase = (1.0 + targetRatioA) * (1.0 - attackCoef);
    }

    inline void ADSR::setDecayRate(Float rate) {
        decayRate = rate;
        decayCoef = calcCoef(rate, targetRatioDR);
        decayBase = (sustainLevel - targetRatioDR) * (1.0 - decayCoef);
    }

    inline void ADSR::setReleaseRate(Float rate) {
        releaseRate = rate;
        releaseCoef = calcCoef(rate, targetRatioDR);
        releaseBase = -targetRatioDR * (1.0 - releaseCoef);
    }

    inline Float ADSR::calcCoef(Float rate, Float targetRatio) {
        return (rate <= 0) ? 0.0 : exp(-log((1.0 + targetRatio) / targetRatio) / rate);
    }

    inline void ADSR::setSustainLevel(Float level) {
        sustainLevel = level;
        decayBase = (sustainLevel - targetRatioDR) * (1.0 - decayCoef);
    }

    inline void ADSR::setTargetRatioA(Float targetRatio) {
        if (targetRatio < 0.000000001)
            targetRatio = 0.000000001;  // -180 dB
        targetRatioA = targetRatio;
        attackCoef = calcCoef(attackRate, targetRatioA);
        attackBase = (1.0 + targetRatioA) * (1.0 - attackCoef);
    }

    inline void ADSR::setTargetRatioDR(Float targetRatio) {
        if (targetRatio < 0.000000001)
            targetRatio = 0.000000001;  // -180 dB
        targetRatioDR = targetRatio;
        decayCoef = calcCoef(decayRate, targetRatioDR);
        releaseCoef = calcCoef(releaseRate, targetRatioDR);
        decayBase = (sustainLevel - targetRatioDR) * (1.0 - decayCoef);
        releaseBase = -targetRatioDR * (1.0 - releaseCoef);
    }
}
