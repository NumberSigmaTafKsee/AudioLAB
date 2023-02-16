inline float ThirdInterp(const float x,const float L1,const float L0,const
float H0,const float H1)
{
    return
    L0 +
    .5f*
    x*(H0-L1 +
       x*(H0 + L0*(-2) + L1 +
          x*( (H0 - L0)*9 + (L1 - H1)*3 +
             x*((L0 - H0)*15 + (H1 -  L1)*5 +
                x*((H0 - L0)*6 + (L1 - H1)*2 )))));
}

// Quantum Super Sampling
// Pops too much right now
struct SuperADSR
{
    std::array<ADSR*,4> super_array;

    SuperADSR(float a, float d, float s, float r, float sr)
    {
        
        for(size_t i = 0; i < 4; i++)
        {
            float ar = noise.uniform_real_distribution(-1,1) * 0.01;
            float dr = noise.uniform_real_distribution(-1,1) * 0.01;
            float sr = noise.uniform_real_distribution(-1,1) * 0.01;
            float rr = noise.uniform_real_distribution(-1,1) * 0.01;

            super_array[i] = new ADSR(a +ar,d+dr,s+sr,r+rr,sr);
        }
    }
    void noteOn() {
        for(size_t i = 0; i < 4; i++) super_array[i]->noteOn();
    }
    void noteOff() {
        for(size_t i = 0; i < 4; i++) super_array[i]->noteOff();
    }
    float Tick()
    {        
        float out[4];
        for(size_t i = 0; i < 4; i++) out[i] = super_array[i]->Tick(); 
        return ThirdInterp(out[3] - out[0], out[0],out[1],out[2],out[3]);   
    }
};

// Quantum Super Sampling Oscillator
struct SuperOSC
{
    std::array<PolyBLEP*,4> super_array;
    size_t select = 0;
    float phase[4];
    float freq[4];

    SuperOSC(float sampleRate) 
    {        
         for(size_t i = 0; i < 4; i++) super_array[i] = new PolyBLEP(sampleRate,PolyBLEP::SAWTOOTH);                                   
    }
    void setFrequency(float f)
    {
        for(size_t i = 0; i < 4; i++)
        {
            // if this is too high it will sound like noise kaka
            float r = noise.uniform_real_distribution(-1,1) * 0.1;
            super_array[i]->setFrequency(f+r);
            phase[i] = 0.005*noise.rand();
            //super_array[i]->setPhase(phase[i]);            
            freq[i] = 0.01*noise.rand();            
        }        
    }    
    float Tick()
    {
        
        float out[4];
        float X[4];
        float Y[4];
        for(size_t i = 0; i < 4; i++) {
            X[i] = std::sin(2*M_PI*phase[i]);
            Y[i] = std::cos(2*M_PI*phase[i]);            
            phase[i] = fmod(phase[i] + freq[i]/sampleRate,1);
        }
        for(size_t i = 0; i < 4; i++) out[i] = Diode(super_array[i]->Tick(1,0.5+0.5*noise.rand(),X[i],Y[i]));        
        return ThirdInterp(out[3] - out[0], out[0],out[1],out[2],out[3]);
    }
};

