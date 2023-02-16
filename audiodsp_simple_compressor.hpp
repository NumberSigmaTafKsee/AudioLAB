// https://github.com/SneezeAttack/JuceCompressor
#include <cmath>
#include <algorithm>

class EnvelopeShaper{

    public:
    
	EnvelopeShaper();
    void processAudioSample(float& sample);
    void setAttack(float attack);
    void setRelease(float release);
	void setSampleRate(float sampleRate);
    
    private:
    
    float m_Envelope;
    float m_SampleRate;
    float m_AttackInMilliSeconds;
    float m_ReleaseInMilliSeconds;
    float m_Attack;
    float m_Release;
    
    void update();
    float calculate(float time);

};

class Compressor{
 
    
    private:
    
    float m_Threshold;
    float m_Ratio;
    
    EnvelopeShaper m_EnvelopeShaper;
    
    const float BOUND_LOG = -96.f;
    const float BOUND_LIN = decibelToAmplitude(BOUND_LOG);
    
    float AmplitudeToDecibel(float amplitude);
    float decibelToAmplitude(float db);
    
    
    public:
    
	Compressor();

    void setRatio(float ratio);
    void setTreshold(float thresh);
    void setAttack(float attack);
    void setRelease(float release);

	void setSampleRate(float sampleRate);
    
    void processAudioSample(float& sample);
};


EnvelopeShaper::EnvelopeShaper(){

    m_Envelope = 0.f;
	m_SampleRate = 44100.f;
	m_AttackInMilliSeconds = 10.f;
	m_ReleaseInMilliSeconds = 100.f;
    m_Attack = 0.f;
    m_Release = 0.f;
}

void EnvelopeShaper::processAudioSample(float& sample){

    if(sample > m_Envelope){
    
        m_Envelope += m_Attack * (sample - m_Envelope);
    }
    else if(sample < m_Envelope){
    
        m_Envelope += m_Release * (sample - m_Envelope);
    }
    sample = m_Envelope;

}

void EnvelopeShaper::setAttack(float attack){

    m_AttackInMilliSeconds = attack;
    update();

}


void EnvelopeShaper::setRelease(float release){

    m_ReleaseInMilliSeconds = release;
    update();

}


void EnvelopeShaper::update(){

    m_Attack = calculate(m_AttackInMilliSeconds);
    m_Release = calculate(m_ReleaseInMilliSeconds);
    
}

float EnvelopeShaper::calculate(float time){

    if(time <= 0.f||m_SampleRate <= 0.f){
        return 1;
    }
    return 1-exp(-1.f/(time*0.001f*m_SampleRate));


}

void EnvelopeShaper::setSampleRate(float sampleRate) {

	m_SampleRate = sampleRate;
}

Compressor::Compressor(){

    m_Threshold = 0.f;
    m_Ratio = 1.f;

}

void Compressor::processAudioSample(float& sample){

    float detectionSignal = sample;
    
    detectionSignal = fabs(detectionSignal);
    
    m_EnvelopeShaper.processAudioSample(detectionSignal);
    
    detectionSignal = AmplitudeToDecibel(detectionSignal);
    
    if(detectionSignal > m_Threshold){
     
        float scale = 1.f-(1.f/m_Ratio);
        float gain = scale * (m_Threshold - detectionSignal);
        
        gain = decibelToAmplitude(gain);
        
        sample *= gain;
    }

}

float Compressor::decibelToAmplitude(float db) {

    return pow(10.f, db / 20.f);

}


float Compressor::AmplitudeToDecibel(float amplitude) {

    amplitude = std::max(amplitude,BOUND_LIN);
    return 20.f*log10(amplitude);

}

void Compressor::setTreshold(float thresh) {

    m_Threshold = thresh;
}

void Compressor::setRatio(float ratio){

    m_Ratio = ratio;

}

void Compressor::setAttack(float attack){

    m_EnvelopeShaper.setAttack(attack);

}


void Compressor::setRelease(float release){

    m_EnvelopeShaper.setRelease(release);

}

void Compressor::setSampleRate(float sampleRate) {

	m_EnvelopeShaper.setSampleRate(sampleRate);
}