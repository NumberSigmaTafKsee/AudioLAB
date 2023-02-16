#include <iostream>
#include <vector>
#include <cmath>
#include <cstring>
#include "AudioDevice.h"
#include "PolyBLEP.h"
#include "ResonantLowpassFilter.h"
#include "Undenormal.hpp"
//#include "BasicDelayLine.h"

using std::cout;
using std::cin;
using std::endl;

PolyBLEP * X;
PolyBLEP * Y;
ResonantLowpassFilter2 filter(44100,440,5,0.5);
//BasicDelayLine delay(44100,200,90,90.0f);

enum Waveform
{
    kWaveformSine = 0,
    kWaveformTriangle,
    kWaveformSquare,
    kWaveformSquareSlopedEdges,
    kNumWaveforms
};

enum Interpolation
{
    kInterpolationNearestNeighbour = 1,
    kInterpolationLinear,
    kInterpolationCubic,
    kNumInterpolations
};

Waveform lfoWaveform = kWaveformSine;
float   lfoPhase;             
float   lfoFreqHz=1.f;
double  sampleRate        = 44100.0f;
double  inverseSampleRate = 1/44100.0f;
int     TotalNumInputChannels = 1;
int     TotalNumOutputChannels = 1;
float   modDepth = 0.5f;

float amp   = -6.0f;
float gain = pow(10,amp/20.0f);

///////////////////////////////////////////////////////////////////////////////
// LFO
// Function for calculating "biased" LFO waveforms with output range [0, 1].
// Phase range [0, 1], output also [0, 1] (not [-1, +1] as for the ordinary Sine function).
///////////////////////////////////////////////////////////////////////////////
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
        return 0.5f + 0.5f*sinf(2*M_PI*phase);
    }
}

void Interleave(size_t n, size_t channel, size_t stride, float * in, float * output) {
    size_t x = 0;
    for(size_t i = channel; i < n*stride; i += stride) {
        output[x++] = in[i];
    }
}
void Deinterleave(size_t n, size_t channel, size_t stride, float * in, float * out) {
    size_t x = channel;
    for(size_t i = 0; i < n; i++) {
        out[i] = in[x += stride];
    }
}


///////////////////////////////////////////////////////////////////////////////
// simple amp
///////////////////////////////////////////////////////////////////////////////
void Amplifier_ProcessBlock(size_t n, float * inputs, float * outputs)
{
    for(size_t i = 0; i < n; i++) {
        outputs[i] = gain * inputs[i];
    }
}

///////////////////////////////////////////////////////////////////////////////
// Chorus
///////////////////////////////////////////////////////////////////////////////
/* Stereo
void ChorusAudioProcessor::processBlock (AudioSampleBuffer& buffer, MidiBuffer& midiMessages)
{
    // Helpful information about this block of samples:
    const int numInputChannels = getNumInputChannels();     // How many input channels for our effect?
    const int numOutputChannels = getNumOutputChannels();   // How many output channels for our effect?
    const int numSamples = buffer.getNumSamples();          // How many samples in the buffer for this block?
    
    int channel, dpw; // dpr = delay read pointer; dpw = delay write pointer
    float dpr, currentDelay, ph;
    
    // Go through each channel of audio that's passed in. In this example we apply identical
    // effects to each channel, regardless of how many input channels there are. For some effects, like
    // a stereo chorus or panner, you might do something different for each channel.
    
    for (channel = 0; channel < numInputChannels; ++channel)
    {
        // channelData is an array of length numSamples which contains the audio for one channel
        float* channelData = buffer.getWritePointer(channel);
        
        // delayData is the circular buffer for implementing delay on this channel
        float* delayData = delayBuffer_.getWritePointer (jmin (channel, delayBuffer_.getNumChannels() - 1));
        
        // Make a temporary copy of any state variables declared in PluginProcessor.h which need to be
        // maintained between calls to processBlock(). Each channel needs to be processed identically
        // which means that the activity of processing one channel can't affect the state variable for
        // the next channel.
        
        dpw = delayWritePosition_;
        ph = lfoPhase_;
        
        for (int i = 0; i < numSamples; ++i)
        {
            const float in = channelData[i];
            float interpolatedSample = 0.0;
            float phaseOffset = 0.0;
            float weight;
            
            // Chorus can have more than 2 voices (where the original, undelayed signal counts as a voice).
            // In this implementation, all voices use the same LFO, but with different phase offsets. It
            // is also possible to use different waveforms and different frequencies for each voice.
            
            for(int j = 0; j < numVoices_ - 1; ++j)
            {
                if(stereo_ != 0 && numVoices_ > 2)
                {
                    // A stereo chorus pans each voice to a different location in the stereo field.
                    // How this is done depends on the number of voices:
                    // -- 2 voices: N/A (need at least 2 delayed voices for stereo chorus)
                    // -- 3 voices: 1 voice left, 1 voice right (0, 1)
                    // -- 4 voices: 1 voice left, 1 voice centre, 1 voice right (0, 0.5, 1)
                    // -- 5 voices: 1 voice left, 1 voice left-centre,
                    //              1 voice right-centre, 1 voice right (0, 0.33, 0.66, 1)
                    
                    weight = (float)j/(float)(numVoices_ - 2);
                    
                    // Left and right channels are mirrors of each other in weight
                    if(channel != 0)
                        weight = 1.0 - weight;
                }
                else
                    weight = 1.0;

                // Add the voice to the mix if it has nonzero weight
                if(weight != 0.0)
                {
                    // Recalculate the read pointer position with respect to the write pointer. A more efficient
                    // implementation might increment the read pointer based on the derivative of the LFO without
                    // running the whole equation again, but this format makes the operation clearer.
                    
                    currentDelay = delay_ + sweepWidth_*lfo(fmodf(ph + phaseOffset, 1.0f), waveform_);
                    dpr = fmodf((float)dpw - (float)(currentDelay * getSampleRate()) + (float)delayBufferLength_,
                                (float)delayBufferLength_);
                    
                    // In this example, the output is the input plus the contents of the delay buffer (weighted by delayMix)
                    // The last term implements a tremolo (variable amplitude) on the whole thing.
          
                    if(interpolation_ == kInterpolationLinear)
                    {
                        // Find the fraction by which the read pointer sits between two
                        // samples and use this to adjust weights of the samples
                        float fraction = dpr - floorf(dpr);
                        int previousSample = (int)floorf(dpr);
                        int nextSample = (previousSample + 1) % delayBufferLength_;
                        interpolatedSample = fraction*delayData[nextSample]
                            + (1.0f-fraction)*delayData[previousSample];
                    }
                    else if(interpolation_ == kInterpolationCubic)
                    {
                        // Cubic interpolation will produce cleaner results at the expense
                        // of more computation. This code uses the Catmull-Rom variant of
                        // cubic interpolation. To reduce the load, calculate a few quantities
                        // in advance that will be used several times in the equation:
                        
                        int sample1 = (int)floorf(dpr);
                        int sample2 = (sample1 + 1) % delayBufferLength_;
                        int sample3 = (sample2 + 1) % delayBufferLength_;
                        int sample0 = (sample1 - 1 + delayBufferLength_) % delayBufferLength_;
                        
                        float fraction = dpr - floorf(dpr);
                        float frsq = fraction*fraction;
                        
                        float a0 = -0.5f*delayData[sample0] + 1.5f*delayData[sample1]
                                    - 1.5f*delayData[sample2] + 0.5f*delayData[sample3];
                        float a1 = delayData[sample0] - 2.5f*delayData[sample1]
                                    + 2.0f*delayData[sample2] - 0.5f*delayData[sample3];
                        float a2 = -0.5f*delayData[sample0] + 0.5f*delayData[sample2];
                        float a3 = delayData[sample1];
                        
                        interpolatedSample = a0*fraction*frsq + a1*frsq + a2*fraction + a3;
                    }
                    else // Nearest neighbour interpolation
                    {
                        // Find the nearest input sample by rounding the fractional index to the
                        // nearest integer. It's possible this will round it to the end of the buffer,
                        // in which case we need to roll it back to the beginning.
                        int closestSample = (int)floorf(dpr + 0.5f);
                        if(closestSample == delayBufferLength_)
                            closestSample = 0;
                        interpolatedSample = delayData[closestSample];
                    }

                    // Store the output sample in the buffer, which starts by containing the input sample
                    channelData[i] += depth_ * weight * interpolatedSample;
                }
                
                // 3-voice chorus uses two voices in quadrature phase (90 degrees apart). Otherwise,
                // spread the voice phases evenly around the unit circle. (For 2-voice chorus, this
                // code doesn't matter since the loop only runs once.)
                if(numVoices_ < 3)
                    phaseOffset += 0.25f;
                else
                    phaseOffset += 1.0f / (float)(numVoices_ - 1);
            }
            
            // Store the current input in the delay buffer (no feedback in a chorus, unlike a flanger).
            delayData[dpw] = in;
            
            // Increment the write pointer at a constant rate. The read pointer will move at different
            // rates depending on the settings of the LFO, the delay and the sweep width.
            
            if (++dpw >= delayBufferLength_)
                dpw = 0;

            // Update the LFO phase, keeping it in the range 0-1
            ph += frequency_*inverseSampleRate_;
            if(ph >= 1.0)
                ph -= 1.0;
        }
    }
    
    // Having made a local copy of the state variables for each channel, now transfer the result
    // back to the main state variable so they will be preserved for the next call of processBlock()
    
    delayWritePosition_ = dpw;
    lfoPhase_ = ph;
    
    // In case we have more outputs than inputs, we'll clear any output
    // channels that didn't contain input data, (because these aren't
    // guaranteed to be empty - they may contain garbage).
    for (int i = numInputChannels; i < numOutputChannels; ++i)
    {
        buffer.clear (i, 0, buffer.getNumSamples());
    }
}
*/



std::vector<float> chorus_delay_buffer;

int   chorus_num_voices=4;
float chorus_delay;
float chorus_sweep_width;
float chorus_depth;
float chorus_frequency;  // LFO frequency (Hz)

Waveform        chorus_waveform;   // What shape should be used for the LFO
Interpolation   chorus_interpolation; // What type of interpolation to use

const float kMaximumDelay = 0.05;
const float kMaximumSweepWidth = 0.05;

int chorus_delayBufferLength_;
int chorus_delayWritePosition_;
float  chorus_lfoPhase_;   // Phase of the low-frequency oscillator
        


void Chorus_Init(float delay=30, float sweep=0.02,float depth=1.0, float freq=0.2f,Waveform waveform=kWaveformSine,
                    Interpolation interp = kInterpolationLinear,int voices=2)
{
    chorus_delay = delay/1000.0f;
    chorus_sweep_width = sweep;
    chorus_depth = depth;
    chorus_waveform = waveform;
    chorus_interpolation = interp;
    chorus_num_voices = voices;
    chorus_frequency = freq;

    
    chorus_lfoPhase_ = 0;
    chorus_delayWritePosition_ = 0;

    int len = (int)((kMaximumDelay + kMaximumSweepWidth)*sampleRate)+3;    
    chorus_delay_buffer.resize(len);
    chorus_delayBufferLength_ = len;
    memset(chorus_delay_buffer.data(),0,len*sizeof(float));    
}
void Chorus_ProcessBlock (size_t n, float * inputs, float * outputs)
{
    
    int channel, dpw; // dpr = delay read pointer; dpw = delay write pointer
    float dpr, currentDelay, ph;
    
    // delayData is the circular buffer for implementing delay on this channel
    float* delayData = &chorus_delay_buffer[0];
             
    dpw = chorus_delayWritePosition_;
    ph  = chorus_lfoPhase_;
        
    for (int i = 0; i < n; ++i)
    {
        float in = inputs[i];
        float interpolatedSample = 0.0;
        float phaseOffset = 0.0;
                    
        outputs[i] = in;
        for(int j = 0; j < chorus_num_voices-1; ++j)
        {                
            currentDelay = chorus_delay + chorus_sweep_width* LFO_GetSample(fmodf(ph + phaseOffset, 1.0f), chorus_waveform);
            dpr = fmodf((float)dpw - (float)(currentDelay * sampleRate) + (float)chorus_delayBufferLength_, (float)chorus_delayBufferLength_);
             
            if(chorus_interpolation == kInterpolationLinear)
            {
                // Find the fraction by which the read pointer sits between two
                // samples and use this to adjust weights of the samples
                float fraction = dpr - floorf(dpr);
                int previousSample = (int)floorf(dpr);
                int nextSample = (previousSample + 1) % chorus_delayBufferLength_;
                interpolatedSample = fraction*delayData[nextSample]
                    + (1.0f-fraction)*delayData[previousSample];
            }
            else if(chorus_interpolation == kInterpolationCubic)
            {
                // Cubic interpolation will produce cleaner results at the expense
                // of more computation. This code uses the Catmull-Rom variant of
                // cubic interpolation. To reduce the load, calculate a few quantities
                // in advance that will be used several times in the equation:
                
                int sample1 = (int)floorf(dpr);
                int sample2 = (sample1 + 1) % chorus_delayBufferLength_;
                int sample3 = (sample2 + 1) % chorus_delayBufferLength_;
                int sample0 = (sample1 - 1 + chorus_delayBufferLength_) % chorus_delayBufferLength_;
                
                float fraction = dpr - floorf(dpr);
                float frsq = fraction*fraction;
                
                float a0 = -0.5f*delayData[sample0] + 1.5f*delayData[sample1]
                            - 1.5f*delayData[sample2] + 0.5f*delayData[sample3];
                float a1 = delayData[sample0] - 2.5f*delayData[sample1]
                            + 2.0f*delayData[sample2] - 0.5f*delayData[sample3];
                float a2 = -0.5f*delayData[sample0] + 0.5f*delayData[sample2];
                float a3 = delayData[sample1];
                
                interpolatedSample = a0*fraction*frsq + a1*frsq + a2*fraction + a3;
            }
            else // Nearest neighbour interpolation
            {
                // Find the nearest input sample by rounding the fractional index to the
                // nearest integer. It's possible this will round it to the end of the buffer,
                // in which case we need to roll it back to the beginning.
                int closestSample = (int)floorf(dpr + 0.5f);
                if(closestSample == chorus_delayBufferLength_)
                    closestSample = 0;
                interpolatedSample = delayData[closestSample];
            }
            
            // Store the output sample in the buffer, which starts by containing the input sample
            outputs[i] += chorus_depth * interpolatedSample;
        
            // 3-voice chorus uses two voices in quadrature phase (90 degrees apart). Otherwise,
            // spread the voice phases evenly around the unit circle. (For 2-voice chorus, this
            // code doesn't matter since the loop only runs once.)
            if(chorus_num_voices < 3)
                phaseOffset += 0.25f;
            else
                phaseOffset += 1.0f / (float)(chorus_num_voices - 1);        
        }
        
        
        // Store the current input in the delay buffer (no feedback in a chorus, unlike a flanger).
        delayData[dpw] = in;
        
        
        // Increment the write pointer at a constant rate. The read pointer will move at different
        // rates depending on the settings of the LFO, the delay and the sweep width.
        
        if (++dpw >= chorus_delayBufferLength_)
            dpw = 0;

        // Update the LFO phase, keeping it in the range 0-1
        ph += chorus_frequency*inverseSampleRate;
        if(ph >= 1.0)
            ph -= 1.0;
    }    
    // Having made a local copy of the state variables for each channel, now transfer the result
    // back to the main state variable so they will be preserved for the next call of processBlock()    
    chorus_delayWritePosition_ = dpw;
    chorus_lfoPhase_ = ph;
    
}

///////////////////////////////////////////////////////////////////////////////
/* Ring Modulator
///////////////////////////////////////////////////////////////////////////////
// Process one buffer ("block") of data
void RingModProcessor::processBlock (AudioBuffer<float>& buffer, MidiBuffer&)
{
    ScopedNoDenormals noDenormals;

    // local copies of state variables carrierPhase and lfoPhase
    float cphi = carrierPhase;
    float lphi = lfoPhase;
    float dlphi = float(parameters.lfoFreqHz * inverseSampleRate);

    // apply the same modulation to all input channels for which there is an output channel
    int channelIndex = 0;
    for (; channelIndex < getTotalNumInputChannels(); channelIndex++)
    {
        // restart the phase sequence
        cphi = carrierPhase;
        lphi = lfoPhase;

        const float* pIn = buffer.getReadPointer(channelIndex);
        float* pOut = buffer.getWritePointer(channelIndex);

        for (int i = 0; i < buffer.getNumSamples(); i++)
        {
            // Carrier oscillator is a simple sine wave
            const float twoPi = 6.283185f;
            float carrier = sin(twoPi * cphi);

            // Ring modulation is just simple multiplication
            *pOut++ = *pIn++ * carrier;

            // Update carrier phase with FM, keeping in range [0, 1]
            float lfo = RingModLFO::getSample(lphi, parameters.lfoWaveform);
            float deltaCarrierHz = parameters.lfoWidthHz * lfo;
            float dcphi = float((parameters.carrierFreqHz + deltaCarrierHz) * inverseSampleRate);
            cphi += dcphi;
            while (cphi >= 1.0) cphi -= 1.0;

            // Update LFO phase, keeping in range [0, 1]
            lphi += dlphi;
            while (lphi >= 1.0) lphi -= 1.0;
        }
    }

    // update the main phase state variables, ready for the next processBlock() call
    carrierPhase = cphi;
    lfoPhase = lphi;

    // clear any remaining/excess output channels to zero
    for (; channelIndex < getTotalNumOutputChannels(); channelIndex++)
    {
        buffer.clear(channelIndex, 0, buffer.getNumSamples());
    }
}
*/

// Process one buffer ("block") of data
void Tremolo_ProcessBlock (size_t n, float * inputs, float * outputs)
{
    Undenormal denormal;

    // LFO phase starts at current LFO phase for each channel, advances by deltaPhi for each sample
    float phi = lfoPhase;
    float deltaPhi = float(lfoFreqHz * inverseSampleRate);

    // restart the phase sequence
    phi = lfoPhase;

    float* pIn  = inputs;
    float* pOut = outputs;

    for (int i = 0; i < n; i++)
    {
        float modAmount = LFO_GetSample(phi, lfoWaveform);
        float x = (1.0f - modDepth * modAmount);
        *pIn++ *= x;
        *pOut++ *= x;
        // Update LFO phase, keeping in range [0, 1]
        phi += deltaPhi;
        while (phi >= 1.0) phi -= 1.0;
    }

    // update the main LFO phase state variable, ready for the next processBlock() call
    lfoPhase = phi;    
}

// Process one buffer ("block") of data
void Tremolo_StereoProcessBlock (size_t n, float * inputs, float * outputs)
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
            float modAmount = LFO_GetSample(phi, lfoWaveform);
            float x = (1.0f - modDepth * modAmount);
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


///////////////////////////////////////////////////////////////////////////////
// Compressor
///////////////////////////////////////////////////////////////////////////////
bool autoTime=true;
float x_g[1024], x_l[1024], y_g[1024], y_l[1024], c[1024]; // input, output, control
float yL_prev=0;
int samplerate;

float compressionRatio = 10.0f;
float attackTimeMs = 15.0f, releaseTimeMs=100.0f;
float threshold=-24.0f; //dB
float makeUpGain=0.0f;
float currentGain = 0;

float ratio=10.0f,tauAttack=15.0f,tauRelease=100.0f,alphaAttack=15,alphaRelease=100;



void computeCompressionGain(size_t bufferSize, float * buffer)
{
    float alphaAttack = exp(-1.0f / (0.001f * sampleRate * attackTimeMs));
    float alphaRelease = exp(-1.0f / (0.001f * sampleRate * releaseTimeMs));
    float yl_avg = 0.0f;

    for (int i = 0; i < bufferSize; ++i)
    {
        // Level detection- estimate level using peak detector
        if (fabs(buffer[i]) < 0.000001f) x_g[i] = -120;
        else x_g[i] = 20 * log10(fabs(buffer[i]));

        // Gain computer- static apply input/output curve
        if (x_g[i] >= threshold)
            y_g[i] = threshold + (x_g[i] - threshold) / compressionRatio;
        else
            y_g[i] = x_g[i];

        x_l[i] = x_g[i] - y_g[i];

        // Ballistics- smoothing of the gain 
        if (x_l[i] > yL_prev)
            y_l[i] = alphaAttack * yL_prev + (1 - alphaAttack) * x_l[i];
        else
            y_l[i] = alphaRelease * yL_prev + (1 - alphaRelease) * x_l[i];

        // Accumulate averaged gain for the whole buffer (for GUI display)
        yl_avg += y_l[i];

        // find control
        c[i] = pow(10.0f, (makeUpGain - y_l[i]) / 20.0f);
        yL_prev = y_l[i];
    }

    yl_avg /= bufferSize;
    currentGain = pow(10.0f, -yl_avg / 20.0f);    
}

/* Stereo
void CompressorAudioProcessor::processBlock (AudioSampleBuffer& buffer, MidiBuffer& midiMessages)
{
	if (compressorONOFF)
	{
		inputBuffer.setSize(M,bufferSize);
		inputBuffer.clear();
		for (int m = 0 ; m < M ; ++m)
		{
			if ( (threshold< 0) )
			{
				inputBuffer.clear(m,0,bufferSize);
				// Mix down left-right to analyse the input		
				inputBuffer.addFrom(m,0,buffer,m*2,0,bufferSize,0.5);
				inputBuffer.addFrom(m,0,buffer,m*2+1,0,bufferSize,0.5);
				// compression : calculates the control voltage
				compressor(inputBuffer,m);
				// apply control voltage to the audio signal
				for (int i = 0 ; i < bufferSize ; ++i)
				{
					buffer.getWritePointer(2*m+0)[i] *= c[i];
					buffer.getWritePointer(2*m+1)[i] *= c[i];
				}
				inputBuffer.clear(m,0,bufferSize);
				// Mix down left-right to analyse the output
				inputBuffer.addFrom(m,0,buffer,m*2,0,bufferSize,0.5);
				inputBuffer.addFrom(m,0,buffer,m*2+1,0,bufferSize,0.5);
			}
		}
	}
}
// compressor functions
void CompressorAudioProcessor::compressor(AudioSampleBuffer &buffer, int m)
{
	alphaAttack = exp(-1/(0.001 * samplerate * tauAttack));
	alphaRelease= exp(-1/(0.001 * samplerate * tauRelease));
	for (int i = 0 ; i < bufferSize ; ++i)
	{
		//Level detection- estimate level using peak detector
		if (fabs(buffer.getWritePointer(m)[i]) < 0.000001) x_g[i] =-120;
		else x_g[i] =20*log10(fabs(buffer.getWritePointer(m)[i]));
		//Gain computer- static apply input/output curve
		if (x_g[i] >= threshold) y_g[i] = threshold+ (x_g[i] - threshold) / ratio;
		else y_g[i] = x_g[i];
		x_l[i] = x_g[i] - y_g[i];
		//Ballistics- smoothing of the gain 
		if (x_l[i]>yL_prev)  y_l[i]=alphaAttack * yL_prev+(1 - alphaAttack ) * x_l[i] ; 
		else				 y_l[i]=alphaRelease* yL_prev+(1 - alphaRelease) * x_l[i] ;
		//find control
		c[i] = pow(10,(makeUpGain - y_l[i])/20);
		yL_prev=y_l[i];
	}
}
*/

/* Compressor 2 Stereo
// Process one buffer ("block") of data
void CompressorProcessor::processBlock (AudioBuffer<float>& buffer, MidiBuffer&)
{
    ScopedNoDenormals noDenormals;

    int M = getTotalNumInputChannels() / 2;
    inputBuffer.setSize(M, bufferSize);

    for (int m = 0; m < M; ++m)
    {
        if ((parameters.threshold < 0.0f))
        {
            // Mix down left-right to analyse the input (feedforward)
            inputBuffer.clear(m, 0, bufferSize);
            inputBuffer.addFrom(m, 0, buffer, m * 2, 0, bufferSize, 0.5f);
            inputBuffer.addFrom(m, 0, buffer, m * 2 + 1, 0, bufferSize, 0.5f);

            // compression : calculates the control voltage
            computeCompressionGain(inputBuffer, m);

            // apply control voltage to the audio signal
            const float *pIL = buffer.getReadPointer(2 * m + 0);
            const float *pIR = buffer.getReadPointer(2 * m + 1);
            float* pOL = buffer.getWritePointer(2 * m + 0);
            float* pOR = buffer.getWritePointer(2 * m + 1);
            for (int i = 0; i < bufferSize; ++i)
            {
                float cv = c[i];
                *pOL++ = cv * *pIL++;
                *pOR++ = cv * *pIR++;
            }

            // Mix down left-right to analyse the output (feedback)
            //inputBuffer.clear(m, 0, bufferSize);
            //inputBuffer.addFrom(m, 0, buffer, m * 2, 0, bufferSize, 0.5f);
            //inputBuffer.addFrom(m, 0, buffer, m * 2 + 1, 0, bufferSize, 0.5f);
        }
    }
}

void CompressorProcessor::computeCompressionGain(AudioSampleBuffer &buffer, int m)
{
    float alphaAttack = exp(-1.0f / (0.001f * samplerate * parameters.attackTimeMs));
    float alphaRelease = exp(-1.0f / (0.001f * samplerate * parameters.releaseTimeMs));

    float yl_avg = 0.0f;

    for (int i = 0; i < bufferSize; ++i)
    {
        // Level detection- estimate level using peak detector
        if (fabs(buffer.getWritePointer(m)[i]) < 0.000001f) x_g[i] = -120;
        else x_g[i] = 20 * log10(fabs(buffer.getWritePointer(m)[i]));

        // Gain computer- static apply input/output curve
        if (x_g[i] >= parameters.threshold)
            y_g[i] = parameters.threshold + (x_g[i] - parameters.threshold) / parameters.compressionRatio;
        else
            y_g[i] = x_g[i];
        x_l[i] = x_g[i] - y_g[i];

        // Ballistics- smoothing of the gain 
        if (x_l[i] > yL_prev)
            y_l[i] = alphaAttack * yL_prev + (1 - alphaAttack) * x_l[i];
        else
            y_l[i] = alphaRelease * yL_prev + (1 - alphaRelease) * x_l[i];

        // Accumulate averaged gain for the whole buffer (for GUI display)
        yl_avg += y_l[i];

        // find control
        c[i] = pow(10.0f, (parameters.makeUpGain - y_l[i]) / 20.0f);
        yL_prev = y_l[i];
    }

    yl_avg /= bufferSize;
    currentGain = pow(10.0f, -yl_avg / 20.0f);
    sendChangeMessage();
}
*/

void compressor(size_t n, float * input)
{
	alphaAttack = exp(-1/(0.001 * samplerate * tauAttack));
	alphaRelease= exp(-1/(0.001 * samplerate * tauRelease));
	for (int i = 0 ; i <  n; ++i)
	{
		//Level detection- estimate level using peak detector
		if (fabs(input[i]) < 0.000001) x_g[i] =-120;
		else x_g[i] =20*log10(fabs(input[i]));
		//Gain computer- static apply input/output curve
		if (x_g[i] >= threshold) y_g[i] = threshold+ (x_g[i] - threshold) / ratio;
		else y_g[i] = x_g[i];
		x_l[i] = x_g[i] - y_g[i];
		//Ballistics- smoothing of the gain 
		if (x_l[i]>yL_prev)  y_l[i]=alphaAttack * yL_prev+(1 - alphaAttack ) * x_l[i] ; 
		else				 y_l[i]=alphaRelease* yL_prev+(1 - alphaRelease) * x_l[i] ;
		//find control
		c[i] = pow(10,(makeUpGain - y_l[i])/20);
		yL_prev=y_l[i];
	}
}
void Compressor_ProcessBlock (size_t n, float * input, float * output)
{
    Undenormal denormal;
    if (threshold < 0.0f)
    {        
        // compression : calculates the control voltage
        //computeCompressionGain(n,input);
        compressor(n,input);
        // apply control voltage to the audio signal        
        for (int i = 0; i < n; ++i)
        {
            float cv = c[i];
            output[i] = cv * input[i];            
        }        
    }
}


struct StreamCallbackData {
    const float * input;
    float * output;
    unsigned long frameCount;
    const PaStreamCallbackTimeInfo *timeinfo;
    StreamCallbackFlags flags;
    void * userData;
};

int ProcessFunction(StreamCallbackData * ptr) 
{
    float max=-999;
    for(size_t i = 0; i < ptr->frameCount; i++)    
    {
        ptr->output[i] = filter.Tick(X->Tick());
        if(ptr->output[i] > max) max=ptr->output[i];
    }        
    Chorus_ProcessBlock(ptr->frameCount,ptr->output,ptr->output);   
    //Tremolo_ProcessBlock(ptr->frameCount,ptr->output,ptr->output);   
    //Amplifier_ProcessBlock(ptr->frameCount,ptr->output,ptr->output);   
    //Compressor_ProcessBlock(ptr->frameCount,ptr->output,ptr->output);      
    return 0;
}

int StreamCallback(const void * input, void * output, unsigned long frameCount, const PaStreamCallbackTimeInfo * timeinfo, PaStreamCallbackFlags statusFlags, void * userData) {
    StreamCallbackData data;
    data.input = (float*)input;
    data.output = (float*)output;
    data.frameCount = frameCount;
    data.timeinfo = timeinfo;
    data.flags = (StreamCallbackFlags)statusFlags;
    data.userData = userData;
    return ProcessFunction(&data);
}

int main()
{
    AudioDevice device;
    std::cout << device.GetDeviceCount() << std::endl;
    for(size_t i = 0; i < device.GetHostApiCount(); i++) {
        AudioHostApiInfo api(i);
        cout << "API: " << api.name() << endl;
    }
    int device_id = 14;
    for(size_t i = 0; i < device.GetDeviceCount(); i++)
    {
        AudioDeviceInfo info(i);
        cout << "Device # " << i << ": " << info.name() << endl;
        if(!strcmp(info.name(),"pulse")) device_id = i;
    }

    Chorus_Init();   
    X = new PolyBLEP(44100,PolyBLEP::SAWTOOTH);
    Y = new PolyBLEP(44100,PolyBLEP::SAWTOOTH);
    AudioStreamParameters input(device_id,1,SAMPLE_FLOAT32,-1.0);
    AudioStreamParameters output(device_id,1,SAMPLE_FLOAT32,-1.0);
    AudioStream stream(&input,&output,44100,256,STREAM_NOFLAG,StreamCallback,NULL);
    stream.StartStream();
    Sleep(5000);
    stream.StopStream();
    
}