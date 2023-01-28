/*
  ==============================================================================

    PitchShifter.h
    Created: 24 May 2020 12:49:44pm
    Author:  mw

  ==============================================================================
*/
#pragma once
#define M_PI 3.1415926535

class PitchShifter
{
public:
    PitchShifter(unsigned int sample_rate = 44100) :
        sample_rate_(sample_rate)
    {
        setPitch(1.0);
    }

    float processSample(float sample) 
    {
		smbPitchShift(pitch_, 1, 1024, 4, sample_rate_, &sample, &sample);
        return sample;
    }

    /**
    The pitch must be between 0.5 (octave lower) and 2 (octave higher). 1 is no pitch change.
    */
    bool setPitch(float pitch) 
    {
        if (!(0.5f <= pitch <= 2.f)) 
        {
            return false;
        }

		this->pitch_ = pitch;
        return true;
    }

    bool setSampleRate(unsigned int sample_rate)
    {
        if (!(22050 <= sample_rate <= 192000))
        {
            return false;
        }

        sample_rate_ = sample_rate;
        return true;
    }

private:
    void smbFft(float* fftBuffer, long fftFrameSize, long sign);
    double smbAtan2(double x, double y);
    void smbPitchShift(float pitchShift, long numSampsToProcess, long fftFrameSize, long osamp, float sampleRate, float* indata, float* outdata);

    static constexpr size_t max_frame_length = 8192;
    float inFIFO_[max_frame_length];
    float outFIFO_[max_frame_length];
    float FFTworkspace_[2 * max_frame_length];
    float lastPhase_[max_frame_length / 2 + 1];
    float phaseSum_[max_frame_length / 2 + 1];
    float outputAccumulator_[2 * max_frame_length];
    float analyzedFrequency_[max_frame_length];
    float analyzedMagnitude_[max_frame_length];
    float synthesizedFrequency_[max_frame_length];
    float synthesizedMagnitude_[max_frame_length];

    long range_over_ = 0;
    bool initialized_ = false;
	float pitch_ = 1.f;
    unsigned int sample_rate_ = 44100;
};