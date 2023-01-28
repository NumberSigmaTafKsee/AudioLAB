#include <cmath>
#include <cfloat>
#include <cstdio>

#include <dsp/DSPVector.hh>

#define VECLEN      1024    // vector length
#define FREQ        0.1     // normalized frequency


void sig_gen (clDSPVector<float> &dspvSamples)
{
    int iCycles;
    int iSampleCntr;
    float fPI;
    float fDeltaF;
    
    iCycles = (int) (1.0 / FREQ + 0.5);
    fPI = acosf(-1);
    fDeltaF = fPI * 2.0f / iCycles;
    for (iSampleCntr = 0; iSampleCntr < dspvSamples.Size(); iSampleCntr++)
        dspvSamples[iSampleCntr] = sinf(iSampleCntr % iCycles * fDeltaF - fPI);
}


int main (int argc, char *argv[])
{
    int iCntr;
    FILE *fileSamples;
    FILE *fileSpectrum;
    clDSPVector<float> dspvSamples(VECLEN);
    clDSPVector<stSCplx> dspvSpectrum;
    clDSPVector<float> dspvMagnitudeSpectrum;

    sig_gen(dspvSamples);
    
    fileSamples = fopen("example1.samples.txt", "wt");
    for (iCntr = 0; iCntr < dspvSamples.Size(); iCntr++)
        fprintf(fileSamples, "%g\n", dspvSamples[iCntr]);
    fclose(fileSamples);

    dspvSpectrum.FFTInitialize(dspvSamples.Size(), true);
    dspvSpectrum.FFT(dspvSamples);
    dspvMagnitudeSpectrum.Magnitude(dspvSpectrum);
    
    fileSpectrum = fopen("example1.spectrum.txt", "wt");
    for (iCntr = 0; iCntr < dspvMagnitudeSpectrum.Size(); iCntr++)
        fprintf(fileSpectrum, "%g\n", dspvMagnitudeSpectrum[iCntr]);
    fclose(fileSpectrum);
    
    return 0;
}
