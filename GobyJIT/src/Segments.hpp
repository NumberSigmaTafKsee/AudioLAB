#pragma once

struct Segment
{

};

struct Ramp : public Segment
{
    double rate,sr,phase;

    Ramp(float F, float Fs = 44100.0f)
    {
        rate = F/Fs;
        phase = 0;
        sr = Fs;
    }    
    void setFrequency(float f) {
        rate = f / sr;
    }

    float Tick()
    {
        float r = phase;
        phase = fmod(phase + rate,1);
        return serpent_curve(r);
    }    
};

struct InverseRamp : public Segment
{
    double rate,sr,phase;

    InverseRamp(float F, float Fs = 44100.0f)
    {
        rate = F/Fs;
        phase = 0;
        sr = Fs;
    }    
    void setFrequency(float f) {
        rate = f / sr;
    }

    float Tick()
    {
        float r = 1.0 - phase;
        phase = fmod(phase + rate,1);
        return serpent_curve(r);
    }    
};


struct Attack : public Segment
{
    double rate,sr;
    double output;
    double target = 1.0;
    double ctr;

    Attack(double time, double Fs = 44100.0f) {
        sr = Fs;
        rate = time * sr;
        output = 0;
        ctr = 0;
    }
    void setRate(double time) {
        rate = time * sr;
    }
    double Tick()
    {
        output = ctr/rate;
        ctr++;
        if(output > target) output = target;
        return serpent_curve(output);
    }
    void Reset() {
        ctr = 0;
    }
    bool Finished()
    {
        return output >= target;
    }
};

struct Decay : public Segment
{
    double rate,sr;
    double output;    
    double target = 0.0;
    double ctr;

    Decay(double time, double level, double Fs = 44100.0f) {
        sr = Fs;
        rate = time * sr;
        output = level;
        ctr = 0;
    }
    void setRate(double time) {
        rate = time * sr;
    }
    double Tick()
    {
        output = (rate-ctr)/rate;
        ctr++;
        if(output < target) output = target;
        return serpent_curve(output);
    }
    void Reset() {
        ctr = 0;
    }
    bool Finished()
    {
        return output <= target;
    }
};

// this is a timed sustain level
struct Sustain : public Segment
{
    double rate,sr;
    double output;    
    double target = 0.0;
    double ctr;

    Sustain(double time, double level, double Fs = 44100.0f) {
        sr = Fs;
        rate = time * sr;
        output = level;
        ctr = 0;
    }
    void setRate(double time) {
        rate = time * sr;
    }
    double Tick()
    {        
        ctr++;        
        return serpent_curve(output);
    }
    void Reset() {
        ctr = 0;
    }
    bool Finished()
    {
        return ctr >= rate;
    }
};

struct ASR
{
    Attack  attack;    
    Decay   release;
    double  level;
    int     state = 0;

    ASR(float a, float s, float r, float fs = 44100.0f)
    : attack(a,fs), level(s), release(r,fs)
    {

    }
    void noteOn()
    {
        state = 0;
        attack.Reset();
        release.Reset();
    }
    void noteOff()
    {
        state = 2;
    }    
    double Tick()
    {
        if(state == -1) return 0;
        if(state == 0)
        {
            float r = attack.Tick();
            if(attack.Finished()) state = 1;
            return r;
        }
        if(state == 1)
        {
            return level;
        }
        if(state == 2)
        {
            float r = release.Tick();
            if(release.Finished()) state = -1;
            return r;
        }
        return 0;
    }
};


struct YangEnvelope
{
    std::vector<Segment*> segments;
    int current_segment;
};
