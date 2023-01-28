#include <cmath>

//#define PI 3.14159265358979323846
#define TAU 6.28318530717958647692

namespace FX::Distortion
{
    class Distortion
    {
    public:
        struct Controls {
            // Distortion mode, 0 = bypass, 1 = soft-clip, 2 = hard-clip
            int mode;
            // Drive, [1., ?), the amount of gain prior to the non-linearity
            float drive;
            // Threshold, [0., 1], the amplitude at which to start fold-back distortion
            float threshold;
            // Mix, [0., 1.] ratio between a dry and wet signal
            float mix;
        } controls;
        
        Distortion();
        ~Distortion();
        float processSample(float sample);
        
    private:
        // Intermediate values
        float input, output = 0.f;
        float softClipThreshold = 2.f / 3.f;
        
        // Nonlinearities
        float softClip(float sample);
        float arctangent(float sample, float alpha);
        float hardClip(float sample);
        float squareLaw(float sample, float alpha);
        float cubicWaveShaper(float sample);
        
        float foldback(float sample);
        float waveShaper1(float sample, float alpha);
        float waveShaper2(float sample, float alpha);
        float waveShaper3(float sample, float alpha);
        
        float gloubiBoulga(float sample);
        float gloubiApprox(float sample);
    };

    Distortion::Distortion() {
        controls.mode = 0;
        controls.drive = 1.f;
        controls.mix = 0.f;
    }

    Distortion::~Distortion() {}

    float Distortion::processSample(float sample)
    {
        input = sample;
        output = input * controls.drive;
        
        switch (controls.mode) {
            case 1:
                output = softClip(output);
                break;
            case 2:
                output = arctangent(input, controls.drive);
                break;
            case 3:
                output = hardClip(output);
                break;
            case 4:
                output = squareLaw(input, controls.drive);
                break;
            case 5:
                output = cubicWaveShaper(output);
                break;
            case 6:
                output = foldback(output);
                break;
            case 7:
                output = gloubiApprox(output);
                break;
            case 8:
                output = gloubiBoulga(output);
                break;
            default:
                output = input;
                break;
        }
        
        return (1.f - controls.mix) * input + controls.mix * output;
    }

    /** Cubic soft-clipping nonlinearity
     
        Use 3x oversampling to eliminate aliasing
    */
    float Distortion::softClip(float sample)
    {
        if (sample < -1.f) {
            return -softClipThreshold;
        }
        else if (sample > 1.f) {
            return softClipThreshold;
        }
        else {
            return sample - ((sample * sample * sample) / 3.f);
        }
    }

    // Arctangent nonlinearity
    float Distortion::arctangent(float sample, float alpha)
    {
        // f(x) = (2 / PI) * arctan(alpha * x[n]), where alpha >> 1 (drive param)
        return (2.f / M_PI)* std::atan(alpha * sample);
    }

    // Hard-clipping nonlinearity
    float Distortion::hardClip(float sample)
    {
        if (sample < -1.f) {
            return -1.f;
        }
        else if (sample > 1.f) {
            return 1.f;
        }
        else {
            return sample;
        }
    }

    // Square law series expansion
    float Distortion::squareLaw(float sample, float alpha)
    {
        return sample + alpha * sample * sample;
    }

    /** A cubic nonlinearity, input range: [-1, 1]?
     
        Use 3x oversampling to eliminate aliasing
    */
    float Distortion::cubicWaveShaper(float sample)
    {
        return 1.5f * sample - 0.5f * sample * sample * sample;
    }

    // Foldback nonlinearity, input range: (-inf, inf)
    float Distortion::foldback(float sample)
    {
        // Threshold should be > 0.f
        if (sample > controls.threshold || sample < -controls.threshold) {
            sample = fabs(fabs(fmod(sample - controls.threshold,
                                    controls.threshold * 4))
                        - controls.threshold * 2) - controls.threshold;
        }
        return sample;
    }

    // A nonlinearity by Partice Tarrabia and Bram de Jong
    float Distortion::waveShaper1(float sample, float alpha)
    {
        const float k = 2.f * alpha / (1.f - alpha);
        return (1.f + k) * sample / (1.f + k * std::fabs(sample));
    }

    // A nonlinearity by Jon Watte
    float Distortion::waveShaper2(float sample, float alpha)
    {
        const float z = M_PI * alpha;
        const float s = 1.f / std::sin(z);
        const float b = 1.f / alpha;
        
        if (sample > b) {
            return 1.f;
        }
        else {
            return std::sin(z * sample) * s;
        }
    }

    // A nonlinearity by Bram de Jong, input range: [-1, 1]
    float Distortion::waveShaper3(float sample, float alpha)
    {
        // original design requires sample be positive
        // alpha: [0, 1]
        bool isNegative = false;
        float output = sample;
        if (sample < 0.f) {
            isNegative = true;
            output = -output;
        }
        
        if (output > alpha) {
            output = alpha + (output - alpha)
                / (1.f + std::pow(((output - alpha) / (1.f - alpha)), 2.f));
        }
        if (output > 1.f) {
            output = (alpha + 1.f) / 2.f;
        }
        
        if (isNegative) {
            output = -output;
        }
        
        return output;
    }


    /** A nonlinearity by Laurent de Soras (allegedily)
        This is very expensive, and someone recommended using
        f(x) = x - 0.15 * x^2 - 0.15 * x^3 for a fast approximation.
    */
    float Distortion::gloubiBoulga(float sample)
    {
        const double x = sample * 0.686306;
        const double a = 1 + std::exp(std::sqrt(fabs(x)) * -0.75);
        return (std::exp(x) - std::exp(-x * a)) / (std::exp(x) + std::exp(-x));
    }

    // Approximation based on description in gloubiBoulga
    float Distortion::gloubiApprox(float sample)
    {
        return sample - (0.15f * sample * sample) - (0.15f * sample * sample * sample);
    }
}