
enum Waveform
{
    kWaveformSine = 0,
    kWaveformTriangle,
    kWaveformSquare,
    kWaveformSquareSlopedEdges,
    kNumWaveforms
} waveform = kWaveformSine;

float  lfoPhase;             
float  lfoFreqHz=0.1f;
double sampleRate        = 44100.0f;
double inverseSampleRate = 1/44100.0f;
int TotalNumInputChannels = 2;
int TotalNumOutputChannels = 2;


// Function for calculating "biased" LFO waveforms with output range [0, 1].
// Phase range [0, 1], output also [0, 1] (not [-1, +1] as for the ordinary Sine function).
float LFO_GetSample(float phase, Waveform waveform)
{
    switch (waveform)
    {
    case kWaveformTriangle:
        if (phase < 0.25f)
            return 0.5f + 2.0f*phase;
        else if (phase < 0.75f)
            return 1.0f - 2.0f*(phase - 0.25f);
        else
            return 2.0f*(phase - 0.75f);
    case kWaveformSquare:
        if (phase < 0.5f)
            return 1.0f;
        else
            return 0.0f;
    case kWaveformSquareSlopedEdges:
        if (phase < 0.48f)
            return 1.0f;
        else if (phase < 0.5f)
            return 1.0f - 50.0f*(phase - 0.48f);
        else if (phase < 0.98f)
            return 0.0f;
        else
            return 50.0f*(phase - 0.98f);
    case kWaveformSine:
    default:
        return 0.5f + 0.5f*sinf(TWOPI_F * phase);
    }
}


// Process one buffer ("block") of data
void Tremolo_ProcessBlock (size_t n, float * inputs, float * outputs)
{
    Undenormal denormal;

    // LFO phase starts at current LFO phase for each channel, advances by deltaPhi for each sample
    float phi = lfoPhase;
    float deltaPhi = float(lfoFreqHz * inverseSampleRate);

    // apply the same modulation to all input channels for which there is an output channel
    for (int channelIndex = 0; channelIndex < TotalNumInputChannels; channelIndex++)
    {
        // restart the phase sequence
        phi = lfoPhase;

        float* pIn  = inputs  + channelIndex;
        float* pOut = outputs + channelIndex;

        for (int i = 0; i < n; i++)
        {
            float modAmount = LFO_GetSample(phi, parameters.lfoWaveform);
            float x = (1.0f - parameters.modDepth * modAmount);
            *pIn *= x;
            *pOut *= x;
            pIn  += TotalNumInputChannels;
            pOut += TotalNumOutputChannels;

            // Update LFO phase, keeping in range [0, 1]
            phi += deltaPhi;
            while (phi >= 1.0) phi -= 1.0;
        }
    }
    // update the main LFO phase state variable, ready for the next processBlock() call
    lfoPhase = phi;    
}
