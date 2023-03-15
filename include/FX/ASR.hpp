#pragma once

#include "SoundObject.hpp"


namespace Envelopes
{
    // todo: add ports
    class ASR : public GeneratorProcessor {
    public:
        
        DspFloatType min=0.0;
        DspFloatType max=1.0;
        
        ASR(DspFloatType sample_rate=44100.0f) : GeneratorProcessor()
        {        
            sr = sample_rate;
            reset();
            setAttackRate(0);            
            setReleaseRate(0);
            setSustainLevel(1.0);
            setTargetRatioA(0.3);
            setTargetRatioDR(0.0001);        
        }

        ~ASR(void) {

        }


        ASR(DspFloatType a, DspFloatType d, DspFloatType s, DspFloatType r, DspFloatType sample_rate=44100.0f, DspFloatType cutoff=10.0f) 
        : GeneratorProcessor()
        {     
            sr = sample_rate;
            reset();
            setAttackRate(a * sample_rate);
            setSustainLevel(s);
            setReleaseRate(r * sample_rate);                
            setTargetRatioA(0.3);
            setTargetRatioDR(0.0001);
        }
        
        DspFloatType process(void);
        DspFloatType getOutput(void);
        int getState(void);
        void gate(int on);
        
        void setAllTimes(DspFloatType a, DspFloatType d, DspFloatType s, DspFloatType r) {
            reset();
            setAttackTime(a);
            setDecayTime(d);
            setSustainLevel(s);
            setReleaseTime(r);
        }

        enum {
            PORT_ATKTIME,
            PORT_RELTIME,
            PORT_ATKRATE,
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
            printf("PORT_RELTIME=2\n");
            printf("PORT_ATKRATE=3\n");
            printf("PORT_RELRATE=5\n");
            printf("PORT_SUSTAIN=6\n");
            printf("PORT_RATIOA=7\n");
            printf("PORT_RADIODR=8\n");
            printf("PORT_RESET=9\n");
            printf("PORT_NOTEON=10\n");
            printf("PORT_NOTEOFF=11\n");
        }
        void setPort(int port, DspFloatType v) {
            switch(port)
            {
                case PORT_ATKTIME: setAttackTime(v); break;
                case PORT_RELTIME: setReleaseTime(v); break;
                case PORT_ATKRATE: setAttackRate(v); break;
                case PORT_RELRATE: setReleaseRate(v); break;
                case PORT_SUSTAIN: setSustainLevel(v); break;
                case PORT_RATIOA: setTargetRatioA(v); break;
                case PORT_RATIODR: setTargetRatioDR(v); break;
                case PORT_NOTEON: if(v != 0.0) noteOn(); break;
                case PORT_NOTEOFF: if(v != 0.0) noteOff(); break;
                case PORT_RESET: if(v != 0.0) reset(); break;                
            }
        }

        void setAttackTime(DspFloatType rate)  { setAttackRate(rate*sr);}
        void setReleaseTime(DspFloatType rate) { setReleaseRate(rate*sr);}

        void setAttackRate(DspFloatType rate);
        void setReleaseRate(DspFloatType rate) ;
        void setSustainLevel(DspFloatType level);

        void setTargetRatioA(DspFloatType targetRatio);
        void setTargetRatioDR(DspFloatType targetRatio);

        void noteOn() { gate(true); }
        void noteOff() { gate(false); }
        void reset(void);

        
        DspFloatType Tick(DspFloatType I=1, DspFloatType A=1,DspFloatType X=1, DspFloatType Y=1) {            
            DspFloatType r = process();            
            return A*r;
        }

        
        void ProcessSIMD(size_t n, DspFloatType * in, DspFloatType* out) {
            #pragma omp simd aligned(in,out)
            for(size_t i = 0; i < n; i++)
            {
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
                    case env_sustain:
                        break;
                    case env_release:
                        output = releaseBase + output * releaseCoef;
                        if (output <= 0.0) {
                            output = 0.0;
                            state = env_idle;
                        }
                }
                out[i] =  output*(max-min) + min;                    
                if(in) out[i] *= in[i];
            }
        }
        void ProcessBlock(size_t n, DspFloatType * in, DspFloatType* out) {
			ProcessSIMD(n,in,out);
		}
		void ProcessInplace(size_t n, DspFloatType* out) {
			ProcessSIMD(n,nullptr,out);
		}
        
        enum envState {
            env_idle = 0,
            env_attack,
            env_sustain,
            env_release
        };

    protected:    
        int state;
        DspFloatType output;
        DspFloatType attackRate;
        DspFloatType releaseRate;
        DspFloatType attackCoef;
        DspFloatType releaseCoef;
        DspFloatType sustainLevel;
        DspFloatType targetRatioA;
        DspFloatType targetRatioDR;
        DspFloatType attackBase;
        DspFloatType releaseBase;
        DspFloatType sr;

        DspFloatType calcCoef(DspFloatType rate, DspFloatType targetRatio);
    };

    inline DspFloatType ASR::process() {
        switch (state) {
            case env_idle:
                break;
            case env_attack:
                output = attackBase + output * attackCoef;            
                if (output >= 1.0) {
                    output = 1.0;
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
        return output*(max-min) + min;    
    }

    inline void ASR::gate(int gate) {
        if (gate)
            state = env_attack;
        else if (state != env_idle)
            state = env_release;
    }

    inline int ASR::getState() {
        return state;
    }

    inline void ASR::reset() {
        state = env_idle;
        output = 0.0;
    }

    inline DspFloatType ASR::getOutput() {
        return output;
    }

    inline void ASR::setAttackRate(DspFloatType rate) {
        attackRate = rate;
        attackCoef = calcCoef(rate, targetRatioA);
        attackBase = (1.0 + targetRatioA) * (1.0 - attackCoef);
    }

    inline void ASR::setReleaseRate(DspFloatType rate) {
        releaseRate = rate;
        releaseCoef = calcCoef(rate, targetRatioDR);
        releaseBase = -targetRatioDR * (1.0 - releaseCoef);
    }

    inline DspFloatType ASR::calcCoef(DspFloatType rate, DspFloatType targetRatio) {
        return (rate <= 0) ? 0.0 : exp(-log((1.0 + targetRatio) / targetRatio) / rate);
    }

    inline void ASR::setSustainLevel(DspFloatType level) {
        sustainLevel = level;
    }

    inline void ASR::setTargetRatioA(DspFloatType targetRatio) {
        if (targetRatio < 0.000000001)
            targetRatio = 0.000000001;  // -180 dB
        targetRatioA = targetRatio;
        attackCoef = calcCoef(attackRate, targetRatioA);
        attackBase = (1.0 + targetRatioA) * (1.0 - attackCoef);
    }

    inline void ASR::setTargetRatioDR(DspFloatType targetRatio) {
        if (targetRatio < 0.000000001)
            targetRatio = 0.000000001;  // -180 dB
        targetRatioDR = targetRatio;
        releaseCoef = calcCoef(releaseRate, targetRatioDR);
        releaseBase = -targetRatioDR * (1.0 - releaseCoef);
    }
}
