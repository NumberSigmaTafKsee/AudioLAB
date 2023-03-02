
#include "AudioProcessor.hpp"

// AMPLIFIER        - AudioTK, Tubes, Diodes, Distortion, Non-Linear, Wave Digital Filters
// WDF-9000         - chow wdf, rt-wdf, circuit loader, chow circuit creator, scatter matrix
// FUZZ-9000        - Fuzzy Logic Signal Processor, Fuzzy Numbers, Exp Fuzzy Numbers
// MDFA-9000        - Mathemtaical Distortion Morphing Amplifier, Piecewise, BIP/QUAD,, Time Varying, Curves/Beziers/Splines

#include "Analog/LiquidMoog.hpp"
#include "Distortion/SstWaveshaper.hpp"
#include "Distortion/Waveshaping.hpp"
#include "Distortion/ChebyDistortion.hpp"
#include "Analog/MDFM-1000.hpp"

struct DistortionProcessor : public MonoAudioProcessor
{
    Distortion::AmplifierN<5> amplifier;
    Distortion::QuadAmplifier quadamp;
    Distortion::BipolarAmplifier biamp;
    Distortion::BipolarAmplifier2 biamp2;
    Distortion::TwinAmplifier tamp;
    Distortion::RangeAmplifier ramp;
    Distortion::Range4Amplifier ramp2;
    Distortion::ChebyDistortion<4> chebyd;

    // FM5150 FM/Phase Audio Distortion
    double fmdistortion(double x, double A=1, double B=1, double X=0, double Y=0)
    {
        return sin(2*M_PI*(A*sin(2*M_PI*(B*x+Y))+X));
    }
    double pmdistortion(double x, double X=0)
    {
        return sin(2*M_PI*(x+X));
    }


};