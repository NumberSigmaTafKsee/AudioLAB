#pragma once

#include "CEnvelopeDetector.hpp"
#include "Undenormal.hpp"
#include <algorithm>
#include <cmath>

extern double sampleRate;

struct Compreezor
{

    Compreezor(double sampleRate, int samplesPerBlock);
        
	float DetGain=1;  //Input Gain in dB    
	float Threshold=0; //Compressor Threshold
	float AttackTime=10; //Attack Time in Milliseconds
	float ReleaseTime=200; //Release Time in Milliseconds
	float Ratio=4; //Compression Ratio
	float OutputGain=1; //Output Gain in dB
	float KneeWidth=0; //Compressor Knee Width
	bool  DigitalAnalogue = false;; //Digital/Analogue style compression


	int DETECT_MODE_PEAK = 0;
	int DETECT_MODE_MS = 1;
	int DETECT_MODE_RMS = 2;

	CEnvelopeDetector  m_LeftDetector;
	CEnvelopeDetector m_RightDetector;

    float calcCompressorGain(float fDetectorValue, float fThreshold,float fRatio, float fKneeWidth, bool bLimit);
    void ProcessBlock(size_t numSamples, float ** inputs, float ** outputs);
};

void Compreezor::ProcessBlock(size_t numSamples, float ** inputs, float ** outputs)
{
	Undenormal noDenormals;
	const int totalNumInputChannels = 2;
	const int totalNumOutputChannels = 2;
	float fInputGain = std::pow(10.0, DetGain / 20.0);
	float fOutputGain = std::pow(10.0, OutputGain / 20.0);

	
    float* channelDataL = inputs[0];
    float* channelDataR = inputs[1];
    float* channelOutL  = outputs[0];
    float* channelOutR  = outputs[1];

    for (int sample = 0; sample < numSamples; ++sample)
    {       // ..do something to the data...			    

        channelOutL[sample] = channelDataL[sample] * DetGain;
        channelOutR[sample] = channelDataR[sample] * DetGain;
        
        // detect left channel
        float fLeftDetector = m_LeftDetector.detect(channelOutL[sample]);
        // gain calc
        float fGn;
        // branch

        //if (m_uProcessorType == COMP) //always true for this project
        fGn = calcCompressorGain(fLeftDetector, Threshold, Ratio, KneeWidth, false);
        // form output and apply make up gain
        channelOutL[sample] = fGn*channelOutL[sample] * OutputGain;

        float fRightDetector = m_RightDetector.detect(channelOutR[sample]);

        fGn = calcCompressorGain(fRightDetector, Threshold, Ratio, KneeWidth, false);
        // form output and apply make up gain
        channelOutR[sample] = fGn*channelOutR[sample] * OutputGain;
    }

}

float Compreezor::calcCompressorGain(float fDetectorValue, float fThreshold,
	float fRatio, float fKneeWidth, bool bLimit)
{
	// slope variable
	float CS = 1.0 - 1.0 / fRatio; // [ Eq. 13.1 ]
								   // limiting is infinite ratio thus CS->1.0
								   //if (bLimit)
								   //CS = 1;

								   // soft-knee with detection value in range?
	if (fKneeWidth > 0 && fDetectorValue > (fThreshold - fKneeWidth / 2.0) &&
		fDetectorValue < fThreshold + fKneeWidth / 2.0)
	{
		// setup for Lagrange
		double x[2];
		double y[2];
		x[0] = fThreshold - fKneeWidth / 2.0;
		x[1] = fThreshold + fKneeWidth / 2.0;
		x[1] = std::min(0.0, x[1]); // top limit is 0dBFS
		y[0] = 0; // CS = 0 for 1:1 ratio
		y[1] = CS; // current CS

				   // interpolate & overwrite CS
		CS = lagrpol(&x[0], &y[0], 2, fDetectorValue);
	}
	// compute gain; threshold and detection values are in dB
	float yG = CS * (fThreshold - fDetectorValue); // [ Eq. 13.1 ]
												   // clamp; this allows ratios of 1:1 to still operate
	yG = std::min(0.0f, yG);

	// convert back to linear
	return std::pow(10.0, yG / 20.0);
}

//==============================================================================
Compreezor::Compreezor(double sampleRate, int samplesPerBlock)
{
	// Use this method as the place to do any pre-playback
	// initialisation that you need..
	// init the envelope detectors
	// set all params at once with this function; see function definition

	

	if (DigitalAnalogue == true) //Digit
	{
		m_LeftDetector.init((float)sampleRate, AttackTime = 0.0, ReleaseTime,
			false, DETECT_MODE_RMS, true);
		m_RightDetector.init((float)sampleRate, AttackTime,
			ReleaseTime, false, DETECT_MODE_RMS, true);
	}
	else
	{
		m_LeftDetector.init((float)sampleRate, AttackTime = 0.0, ReleaseTime,
			true, DETECT_MODE_RMS, true);
		m_RightDetector.init((float)sampleRate, AttackTime,
			ReleaseTime, true, DETECT_MODE_RMS, true);
	}

}
