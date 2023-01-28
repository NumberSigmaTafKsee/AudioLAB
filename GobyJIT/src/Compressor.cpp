
#include "StereoVector.hpp"
#include "Undenormal.hpp"

struct Compressor
{
    StereoVector mixedDownInput;
    float xl;
    float yl;
    float xg;
    float yg;
    float control;

    float inputLevel;
    float ylPrev;

    float inverseSampleRate;
    float inverseE;
    float calculateAttackOrRelease (float value);

    bool bypassSwitch;
    float paramThreshold;
    float paramRatio;
    float paramAttack;
    float paramRelease;
    float paramMakeupGain;

    int numInputChannels;
    int numOutputChannels;

    enum CompressorType {
        COMPRESSOR_LIMITER=0,
        EXPANDER_NOISEGATE=1
    } 
    type;

    // mixer -> wet/dry -> processors -> output
    // sidechain -> processors -> compressor input/mix

    Compressor(float sr, int channels=2, CompressorType Mode=COMPRESSOR_LIMITER, float ThresholdDB=-24.0f, float Ratio = 50.0f,float AttackMS=2.0f, float ReleaseMS=300.0f, float MakeupGain=0.0f)
    {
        const double smoothTime = 1e-3;
        paramThreshold = ThresholdDB;
        paramRatio = Ratio;
        paramAttack = AttackMS;
        paramRelease = ReleaseMS;
        paramMakeupGain = MakeupGain;
        paramBypass = false;

        mixedDownInput.resize(samplesPerBlock);

        inputLevel = 0.0f;
        ylPrev = 0.0f;
        inverseSampleRate = 1.0f / (float)getSampleRate();
        inverseE = 1.0f / M_E;}

        numInputChannels = channels;
        numOutputChannels = channels;
    }


    void Process(StereoVector& input, StereoVector& output)
    {
        Undenormal noDenormals;

        const int numInputChannels  = numChannels;
        const int numOutputChannels = numChannels;
        const int numSamples = input.size();

        //======================================

        if ((bool)paramBypass.getTargetValue())
            return;

        //======================================

        mixedDownInput.clear();
        for (int channel = 0; channel < numInputChannels; ++channel)
        {
                for(size_t i = 0; i < input[channel].size(); i++)
                    mixedDownInput[channel] += input[channel][i] / numChannels;
        }
            
        for (int sample = 0; sample < numSamples; ++sample) {
            bool expander = (bool)paramMode.getTargetValue();
            float T = paramThreshold.getNextValue();
            float R = paramRatio.getNextValue();
            float alphaA = calculateAttackOrRelease (paramAttack.getNextValue());
            float alphaR = calculateAttackOrRelease (paramRelease.getNextValue());
            float makeupGain = paramMakeupGain.getNextValue();

            float inputSquared = powf (mixedDownInput.getSample (0, sample), 2.0f);
            if (expander) {
                const float averageFactor = 0.9999f;
                inputLevel = averageFactor * inputLevel + (1.0f - averageFactor) * inputSquared;
            } else {
                inputLevel = inputSquared;
            }
            xg = (inputLevel <= 1e-6f) ? -60.0f : 10.0f * log10f (inputLevel);

            // Expander
            if (expander) {
                if (xg > T)
                    yg = xg;
                else
                    yg = T + (xg - T) * R;

                xl = xg - yg;

                if (xl < ylPrev)
                    yl = alphaA * ylPrev + (1.0f - alphaA) * xl;
                else
                    yl = alphaR * ylPrev + (1.0f - alphaR) * xl;

            // Compressor
            } else {
                if (xg < T)
                    yg = xg;
                else
                    yg = T + (xg - T) / R;

                xl = xg - yg;

                if (xl > ylPrev)
                    yl = alphaA * ylPrev + (1.0f - alphaA) * xl;
                else
                    yl = alphaR * ylPrev + (1.0f - alphaR) * xl;
            }

            control = powf (10.0f, (makeupGain - yl) * 0.05f);
            ylPrev = yl;

            for (int channel = 0; channel < numInputChannels; ++channel) {
                float newValue = buffer.getSample (channel, sample) * control;
                buffer.setSample (channel, sample, newValue);
            }
        }
    }

    float calculateAttackOrRelease (float value)
    {
        if (value == 0.0f)
            return 0.0f;
        else
            return pow (inverseE, inverseSampleRate / value);
    }

}