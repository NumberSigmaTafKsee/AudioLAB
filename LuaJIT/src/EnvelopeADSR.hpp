#pragma once



class ADSR {
public:
	
    ADSR(double sample_rate=44100.0f)
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


    ADSR(double a, double d, double s, double r, double sample_rate=44100.0f, double cutoff=10.0f) {     
        sr = sample_rate;
        reset();
        setAttackRate(a * sample_rate);
        setDecayRate(d * sample_rate);
        setSustainLevel(s);
        setReleaseRate(r * sample_rate);                
        setTargetRatioA(0.3);
        setTargetRatioDR(0.0001);
    }
	
	double process(void);
    double getOutput(void);
    int getState(void);
	void gate(int on);
    
    void setAllTimes(double a, double d, double s, double r) {
        reset();
        setAttackTime(a);
        setDecayTime(d);
        setSustainLevel(s);
        setReleaseTime(r);
    }

    void setAttackTime(double rate)  { setAttackRate(rate*sr);}
    void setDecayTime(double rate)   { setDecayRate(rate*sr); }
    void setReleaseTime(double rate) { setReleaseRate(rate*sr);}

    void setAttackRate(double rate);
    void setDecayRate(double rate);
    void setReleaseRate(double rate) ;
	void setSustainLevel(double level);

    void setTargetRatioA(double targetRatio);
    void setTargetRatioDR(double targetRatio);

    void noteOn() { gate(true); }
    void noteOff() { gate(false); }
    void reset(void);

    double Tick(double in) {
        return process() * in;
    }
    double Tick(double I=1, double X=1, double Y=1, double Z=1) {
        double at = attackCoef;
        double dt = decayCoef;
        double rt = releaseCoef;
        double st = sustainLevel;
        attackCoef  *= X;
        decayCoef   *= Y;
        releaseCoef *= Z;
        sustainLevel *= X*Y*Z;
        double r = process();
        attackCoef = at;
        decayCoef  = dt;
        releaseCoef= rt;
        sustainLevel = st;
        return I*r;
    }

    template<typename T>
    void Process(size_t n, T * input, T* output) {
        for(size_t i = 0; i < n; i++)
        {
            output[i] = process()*input[i];
        }
    }
    
    template<typename T>
    void InplaceProcess(size_t n, T * samples) {
        for(size_t i = 0; i < n ; i++)
            samples[i] = process();
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
	double output;
	double attackRate;
	double decayRate;
	double releaseRate;
	double attackCoef;
	double decayCoef;
	double releaseCoef;
	double sustainLevel;
    double targetRatioA;
    double targetRatioDR;
    double attackBase;
    double decayBase;
    double releaseBase;
    double sr;

    double calcCoef(double rate, double targetRatio);
};

inline double ADSR::process() {
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

inline double ADSR::getOutput() {
	return output;
}

inline void ADSR::setAttackRate(double rate) {
    attackRate = rate;
    attackCoef = calcCoef(rate, targetRatioA);
    attackBase = (1.0 + targetRatioA) * (1.0 - attackCoef);
}

inline void ADSR::setDecayRate(double rate) {
    decayRate = rate;
    decayCoef = calcCoef(rate, targetRatioDR);
    decayBase = (sustainLevel - targetRatioDR) * (1.0 - decayCoef);
}

inline void ADSR::setReleaseRate(double rate) {
    releaseRate = rate;
    releaseCoef = calcCoef(rate, targetRatioDR);
    releaseBase = -targetRatioDR * (1.0 - releaseCoef);
}

inline double ADSR::calcCoef(double rate, double targetRatio) {
    return (rate <= 0) ? 0.0 : exp(-log((1.0 + targetRatio) / targetRatio) / rate);
}

inline void ADSR::setSustainLevel(double level) {
    sustainLevel = level;
    decayBase = (sustainLevel - targetRatioDR) * (1.0 - decayCoef);
}

inline void ADSR::setTargetRatioA(double targetRatio) {
    if (targetRatio < 0.000000001)
        targetRatio = 0.000000001;  // -180 dB
    targetRatioA = targetRatio;
    attackCoef = calcCoef(attackRate, targetRatioA);
    attackBase = (1.0 + targetRatioA) * (1.0 - attackCoef);
}

inline void ADSR::setTargetRatioDR(double targetRatio) {
    if (targetRatio < 0.000000001)
        targetRatio = 0.000000001;  // -180 dB
    targetRatioDR = targetRatio;
    decayCoef = calcCoef(decayRate, targetRatioDR);
    releaseCoef = calcCoef(releaseRate, targetRatioDR);
    decayBase = (sustainLevel - targetRatioDR) * (1.0 - decayCoef);
    releaseBase = -targetRatioDR * (1.0 - releaseCoef);
}
