/** ========================================================================= */
//  
class DelayLine {
/** ========================================================================= */
//  public methods
public:
    DelayLine(int maxBufferLength);
    forcedinline float process(float inputSample, float delayLen);
    void processBlock(float* buffer, float delayLen, int numSamples);
    void processBlock(float* buffer, const float* delayLenBuffer, int numSamples);
    void reset();

/** ========================================================================= */
//  public properties
public:
    // the max delay length in samples
    const int maxBufferLength;
    
    // linear feedback multiplier
    float feedback = 0.0f;

/** ========================================================================= */
//  private properties
private:
    // the delay buffer
    AudioBuffer<float> buffer;
    
    //  the write index
    //  i.e. the oldest sample in the buffer
    //       ---> becomes newest sample after writing
    int writeIdx = 0;

    float* data = nullptr;
};