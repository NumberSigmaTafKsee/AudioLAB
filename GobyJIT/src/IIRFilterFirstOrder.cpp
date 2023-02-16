
//#include "../JuceLibraryCode/JuceHeader.h"

#include <math.h>
#include "IIRFilterFirstOrder.h"

 #define SNAP_TO_ZERO(n)    if (! (n < -1.0e-8f || n > 1.0e-8f)) n = 0;

IIRCoefficientsSinglePole::IIRCoefficientsSinglePole() noexcept
{
    //zeromem (coefficients, sizeof (coefficients));
    
    
    //I don't know why memset won't work
    coefficients[0] = 0.f;
    coefficients[1] = 0.f;
    coefficients[2] = 0.f;
}

IIRCoefficientsSinglePole::~IIRCoefficientsSinglePole() noexcept {}

IIRCoefficientsSinglePole::IIRCoefficientsSinglePole (double c0, double c1, double c2, double c3) noexcept
{
    auto a = 1.0 / c2;
    
    coefficients[0] = (float) (c0 * a);
    coefficients[1] = (float) (c1 * a);
    coefficients[2] = (float) (c3 * a);
}

IIRCoefficientsSinglePole IIRCoefficientsSinglePole::makeLowPass (double sampleRate,
                                              double frequency) noexcept
{
    //return makeLowPass (sampleRate, frequency, 1.0 / MathConstants<double>::sqrt2);
    
    return makeLowPass (sampleRate, frequency, 1.0 / sqrt(2.0));

}

IIRCoefficientsSinglePole IIRCoefficientsSinglePole::makeLowPass (double sampleRate,
                                              double frequency,
                                              double Q) noexcept
{
    
    //jassert (sampleRate > 0.0);
    //jassert (frequency > 0 && frequency <= static_cast<float> (sampleRate * 0.5));
    
    //auto n = std::tan (MathConstants<double>::pi * frequency / sampleRate);
    
    //ensure valid sample rate and freq and we're below nyquist
    if (sampleRate <= 0.0 || frequency <= 0.0 || frequency > static_cast<float>(sampleRate * 0.5))
    {
        frequency = 140.0;
        sampleRate = 48000;
    }
    
    
    auto n = tan(M_PI * frequency / sampleRate);
    
    
    return IIRCoefficientsSinglePole(n, n, n + 1, n - 1);
}

IIRCoefficientsSinglePole IIRCoefficientsSinglePole::makeHighPass (double sampleRate,
                                               double frequency) noexcept
{
    return makeHighPass (sampleRate, frequency, 1.0 / sqrt(2.0));
}

IIRCoefficientsSinglePole IIRCoefficientsSinglePole::makeHighPass (double sampleRate,
                                               double frequency,
                                               double Q) noexcept
{
    
    //jassert (sampleRate > 0.0);
    //jassert (frequency > 0 && frequency <= static_cast<float> (sampleRate * 0.5));
    
    //auto n = std::tan (MathConstants<double>::pi * frequency / sampleRate);
    
    //ensure valid sample rate and freq and we're below nyquist
    if (sampleRate <= 0.0 || frequency <= 0.0 || frequency > static_cast<float>(sampleRate * 0.5))
    {
        frequency = 140.0;
        sampleRate = 48000;
    }
    
    auto n = tan(M_PI * frequency / sampleRate);

    
    return IIRCoefficientsSinglePole(1, -1, n + 1, n - 1);
}

//==============================================================================
IIRFilterFirstOrder::IIRFilterFirstOrder() noexcept {}

IIRFilterFirstOrder::~IIRFilterFirstOrder() noexcept {}

void IIRFilterFirstOrder::setCoefficients (const IIRCoefficientsSinglePole& newCoefficients) noexcept
{
    coefficients = newCoefficients;
}

//==============================================================================
void IIRFilterFirstOrder::reset() noexcept
{
    v1 = 0.f;
}

float IIRFilterFirstOrder::processSingleSampleRaw (float in) noexcept
{
    auto out = coefficients.coefficients[0] * in + v1;
    
    //JUCE_SNAP_TO_ZERO (out);
    
    SNAP_TO_ZERO(out);
    
    v1 = coefficients.coefficients[1] * in - coefficients.coefficients[2] * out;
    
    return out;
}
