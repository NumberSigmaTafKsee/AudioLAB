#include "Std\StdRandom.h"

class Random : public Oszillator, public Modulator {
        
    public:
        
        Random(float sampleRate, int bufferSize);
        
        virtual float process() override;
        virtual float getOutput() override;
        virtual void setFrequency(double frequency) override;
        virtual void reset() override;
        virtual void setFine(float fine) override{};
        virtual float getFine() const override { return 0;};
    
    private:
        juce::Random random;
        float lastValue;
        int currentSample = 0;
        long currentTime = 0;
        int buffersize;
    };

Synthlab::Random::Random(float sampleRate, int bufferSize) : Oszillator(sampleRate) {
    this->buffersize = bufferSize;
    currentTime = Time::getMillisecondCounterHiRes();
}

float Synthlab::Random::process() {
    
    float period = (1.0f/frequency) * 1000;
    
    if (Time::getMillisecondCounterHiRes() - currentTime >= period) {
        currentTime = Time::getMillisecondCounterHiRes();
        
        lastValue = random.nextFloat();
    }
    
    return lastValue * volume;
}

float Synthlab::Random::getOutput() {
    return 0;
}

void Synthlab::Random::setFrequency(double frequency){
    this->frequency = frequency;
}

void Synthlab::Random::reset() {
    currentTime = 0;
    lastValue = 0;
}
