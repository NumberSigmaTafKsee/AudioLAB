
// https://github.com/B4phomet/Juce-DOOM-Compressor
class CircularBuffer
{
public:
    CircularBuffer();
    CircularBuffer(int bufferSize, int delayLength);
    float getData();
    void setData(float data);
    void nextSample();
private:
    juce::AudioSampleBuffer buffer;
    int writeIndex;
    int readIndex;
    int delayLength;
};

class Compressor
{
public:
    Compressor();
    float compressSample(float data, float attack, float release);
private:
    CircularBuffer buffer;
    float tav, rms, gain;
};

CircularBuffer::CircularBuffer()
{
    buffer = juce::AudioSampleBuffer();
    writeIndex = readIndex = delayLength = 0;
}

CircularBuffer::CircularBuffer(int bufferSize, int delayLength)
{
    buffer = juce::AudioSampleBuffer(1, bufferSize);
    buffer.clear();
    writeIndex = delayLength;
    readIndex = 0;
    this->delayLength = delayLength;
}

float CircularBuffer::getData()
{
    return buffer.getSample(0, readIndex);
}

void CircularBuffer::setData(float data)
{
    buffer.setSample(0, writeIndex, data);
}

void CircularBuffer::nextSample()
{
    int bufferLength = buffer.getNumSamples();
    readIndex = ((bufferLength + writeIndex) - delayLength) % bufferLength;
    writeIndex = (writeIndex + 1) % bufferLength;
}


Compressor::Compressor()
{
    buffer = CircularBuffer(150, 20);
    tav = 0.01;
    rms = 0;
    gain = 1;
}

float Compressor::compressSample(float data, float attack, float release)
{
    //calculate RMS values of input
    rms = (1 - tav) * rms + tav * std::pow(data, 2.0f);
    float dbRMS = 10 * std::log10(rms);

    //calculate the compressor gain
    float slope = 1 - (1 / INFINITY);
    float dbGain = std::min(0.0f, (slope * (0.0f - dbRMS)));
    float newGain = std::pow(10, dbGain / 20);
    
    //apply attack/release curve to comp gain
    float coeff;
    if (newGain < gain) coeff = attack;
    else coeff = release;
    gain = (1 - coeff) * gain + coeff * newGain;

    //compress the sample
    float compressedSample = gain * buffer.getData();
    buffer.setData(data);
    buffer.nextSample();
    return compressedSample;
}
