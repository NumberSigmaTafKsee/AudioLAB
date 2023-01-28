
float chorus_effect(float in, std::vector<float>& buffer, float Fs, int n,float depth, float rate,float predelay,float wet)
{    
    float t = (n-1)/Fs;
    float lfoMS = depth * std::sin(2*M_PI*rate*t) + predelay;
    float lfoSamples = (lfoMS/1000)*Fs;

    float mixPercent = wet;  
    float mix = mixPercent/100;

    float fracDelay = lfoSamples; 
    int   intDelay = floor(fracDelay); 
    float frac = fracDelay - intDelay; 
      
    float drySig = in; 
    float wetSig = (1.0f-frac)*buffer[intDelay] + (frac)*buffer[intDelay+1];
    
    float out = (1.0f-mix)*drySig + mix*wetSig;

    for(size_t i = 1; i < buffer.size()-1; i++)
        buffer[i] = buffer[i-1];
    buffer[0] = in;    
    return out;
}