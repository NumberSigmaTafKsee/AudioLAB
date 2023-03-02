//
//  ADRS.h
//
//  Created by Nigel Redmon on 12/18/12.
//  EarLevel Engineering: earlevel.com
//  Copyright 2012 Nigel Redmon
//
//  For a complete explanation of the ADSR envelope generator and code,
//  read the series of articles by the author, starting here:
//  http://www.earlevel.com/main/2013/06/01/envelope-generators/
//
//  License:
//
//  This source code is provided as is, without warranty.
//  You may copy and distribute verbatim copies of this document.
//  You may modify and use this source code to create binary code for your own purposes, free or commercial.
//

#pragma once

#include "Util.h"

namespace SoundAlchemy {
    
struct ADSR {

    enum Function
    {
        LINEAR,
        CURVY,
        LOG,
        ILOG,
        SOFTLOG,
        EXP,
        IEXP,
        EXP2,
    } 
    type = LINEAR;

    enum Polarity
    {
        POSITIVE,
        NEGATIVE,
        BIPOLAR,
    }
    polarity = POSITIVE;

    enum envState {
        env_idle = 0,
        env_attack,
        env_decay,
        env_sustain,
        env_release
    };

	int state;
	float output;
	float attackRate;
	float decayRate;
	float releaseRate;
	float attackCoef;
	float decayCoef;
	float releaseCoef;
	float sustainLevel;
    float targetRatioA;
    float targetRatioDR;
    float attackBase;
    float decayBase;
    float releaseBase;
    float gain = 1.0f;
    float minA = 0.0f;
    float maxA = 1.0f;
    float sr;

	ADSR(float sr=44100.float ar = 0, float dr = 0, float sl = 1.0, float rr = 0, float tra = 0.3, float trdr=0.0001);
    
    ADSR(float ams, float dms, float sus, float rms, float sr = 44100)
    {
        reset();
        setAllTimes(ams,dms,sus,rms);
    }

	~ADSR(void);
	
    float process(void);
    float getOutput(void);
    int getState(void);
	void gate(int on);
    
    void setAttackTime(float rate) {
        attackRate = rate/sr;
    }
    void setDecayTime(float rate) {
        decayRate  = rate/sr;
    }
    void setReleasTime(float rate) {
        releaseRate = rate/sr;
    }
	
    void setAttackRate(float rate);
    void setDecayRate(float rate);
    void setReleaseRate(float rate);
	
    void setSustainLevel(float level);
    
    void setTargetRatioA(float targetRatio);
    void setTargetRatioDR(float targetRatio);
    void reset(void);

    float exp1(float x,float gain=1)    { return clamp(std::exp(gain*x)/std::exp(gain))-1,1); }
    float iexp(float x,float gain=1)    { return clamp((std::exp(gain*x)-1)/std::exp(gain)-1,-1,1); }
    float exp2(float x,float gain=1)    { return clamp(-(std::exp(x)-1)/(1+std::exp(gain*x))*(1-std::exp(x)),-1,1); }
    float logi(float x,float gain=1)    { return clamp(std::log(std::abs(gain)*std::abs(gain*x))/std::log(std::abs(gain)),-1,1);}
    float ilogi(float x,float gain=1)   { return clamp(std::log(std::abs(gain*x)-1/(std::log(std::abs(gain)))),-1,1);}
    float softlog(float x,float gain=1) { return clamp(-(std::exp(std::abs(x*gain))-1)/(1+std::exp(std::abs(gain)))*(std::log(std::abs(x/3.6))),-1,1); }           
    float curvey(float x,float g=1)     { return clamp(x/(1.0/cos(atan(5*g*x))*0.2),-1,1); }
    float curv2(float x,float g=1)      { return clamp(2*(x/(std::atanh(g*x))*-1+1,-1,1); }
    float erf(float x,float g=1)        { return clamp(1.5*erf(g*x),-1,1); }
    float fasterf(float x,float g=1)    { return clamp((2*(1+std::tanh(1.5*std::erf(g*x)))-1,-1,1); }
    float gd2(float x,float g=1)        { return clamp(2*(x/(1+std::abs(g*x)),-1,1); }
    float tan2(float x,float g=1)       { return clamp(2*(2/M_PI*std::tan2((2/M_PI)*g*x),-1,1); }
    float atan(float x,float g=1)       { return clamp(3*(2/M_PI*std::atan((2/M_PI)*g*x),-1,1); }
    float asymsoft(float x,float g=1)   { return clamp(std::tanh(1.5*std::erf(-g*x) + 0.5f,-1,1);)}
    float soft(float x, float g =1)     { return clamp(2.0f/(1+std::exp(-g*x)) - 0.5f,-1,1);}
    float sigmoid(float x, float g = 10.0f) { return clamp(1 / (1 + exp(-g*x)),-1,1); }

    float Tick(float in=1, float A = 1, float X = 1, float Y = 1) {
        
        float ta = attackBase;
        float td = decayBase;
        float tr = releaseBase;
        switch(state)
        {
        case env_attack:
            attackBase += X;
            break;
        case env_decay:
            decayBase  += Y;
            break;
        case env_release:
            releaseBase += (X*Y);
            break;
        }
        float r = in * A * X * Y * process();
        attackBase = ta;
        decayBase  = td;
        releaseBase = tr;
        if(r < minA) return minA;
        if(r > maxA) return maxA;
        return r;
    }
    float TickExp11(float in=1, float A = 1, float X = 1, float Y = 1) {
        return exp1(std::abs(in * A * X * Y * process()));
    }
    float TickExp12(float in=1, float A = 1, float X = 1, float Y = 1) {
        return exp2(std::abs(in * A * X * Y * process()));
    }
    float TickIExp(float in=1, float A = 1, float X = 1, float Y = 1) {
        return iexp(std::abs(in * A * X * Y * process()));
    }
    float TickLogi(float in=1, float A = 1, float X = 1, float Y = 1) {
        return logi(std::abs(in * A * X * Y * process()));
    }
    float TickILogi(float in=1, float A = 1, float X = 1, float Y = 1) {
        return ilogi(std::abs(in * A * X * Y * process()));
    }
    float TickSoftLog(float in=1, float A = 1, float X = 1, float Y = 1) {
        return logi(std::abs(in * A * X * Y * process()));
    }
    float TickAbsolute(float in=1, float A = 1, float X = 1, float Y = 1) {
        return std::abs(in * A * X * Y * process());
    }
    float TickSigmoid(float in=1, float A = 1, float X = 1, float Y = 1) {
        float x = (in * A * X * Y * process());
        return 1.0f / (1.0f + std::exp(-x));
    }
    float TickSigmoidal(float in=1, float A = 1, float X = 1, float Y = 1) {
        float x = 1.0f / (1.0f + std::exp(-X));
        float y = 1.0f / (1.0f + std::exp(-Y));
        float a = 1.0f / (1.0f + std::exp(-A));
        float i = 1.0f / (1.0f + std::exp(-in));
        return (i * a * x * y * process());        
    }
    float TickModulus(float in=1, float A = 1, float X = 1, float Y = 1) {
        return in * A * (std::fmod(X,Y)) * process();
    }
    
    void Process(size_t n, float * input, float * output) {
        for(size_t i = 0; i < n; i++)
        {
            output[i] = process()*input[i];
        }
    }
    void Process(size_t n, float * samples) {
        for(size_t i = 0; i < n ; i++)
            samples[i] = process();
    }

    float calcCoef(float rate, float targetRatio);
};

inline float ADSR::process() {
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
    
    output = std::abs(output);
    float gain = std::abs(G);
    if(type == LOG) output *= std::log(gain*output)/std::log(gain);
    else if(type == ILOG) output *= (std::log(gain*output)-1)/(std::log(gain*output));
    else if(type == SOFTLOG) output *=  -(std::exp(ouput)-1)/(1+std::exp(output))*(std::log(output/3.6))           
    else if(type == EXP)  output *= (std::exp(gain*output)/std::exp(gain));
    else if(type == IEXP) output *= ((std::exp(gain*output)-1)/std::exp(gain)-1);
    else if(type == EXP2) output *= (-(std::exp(gain*output)-1)/(1+std::exp(gain*output)))*(1-std::exp(gain));
    else if(type == CURVY) output *= curvy(output);

    if(polarity == NEGATIVE) output = -output;
    else if(polarity == POSITIVE) output = std::abs(output);
    else if(polarity == BIPOLAR) output = 2*output-1;
	
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

inline float ADSR::getOutput() {
	return output;
}
}

