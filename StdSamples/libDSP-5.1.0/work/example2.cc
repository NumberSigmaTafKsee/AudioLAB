#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cfloat>
#include <ctime>

#include <dsp/DSPOp.hh>
#include <dsp/Filter2.hh>

#define VECLEN      1024    // vector length


int main (int argc, char *argv[])
{
    int iRandomValue;
    long lSampleCntr;
    float fScale;
    float fpInData[VECLEN];
    float fpOutData[VECLEN];
    stSCplx sValue;
    stSCplx spInData[VECLEN / 2 + 1];
    FILE *fileSamples;
    clDSPOp DSP;
    clFilter2 Filter;
    
    srand(time(NULL));
    fScale = 1.0f / RAND_MAX;
    DSP.FFTInitialize(VECLEN, true);
    sValue.R = 2.0f / (VECLEN / 2 + 1);
    sValue.I = 0;

    //Filter.InitializeLP(0.45f, 0.55f, 72.0f);
    //Filter.InitializeLP(0.45f, 0.55f, 144.49f, 0.75f);
    Filter.InitializeLP(0.45f, 0.55f, 192.66f, 0.75f);
    printf("Delay for filter 1 is %li\n", Filter.GetDelay());
    while (!Filter.Get(fpOutData, VECLEN))
    {
        /*for (lSampleCntr = 0; lSampleCntr < VECLEN; lSampleCntr++)
        {
            iRandomValue = rand();
            fpInData[lSampleCntr] = (iRandomValue & 1) ?
                (-iRandomValue * fScale) : (iRandomValue * fScale);
        }*/

        DSP.Set(spInData, sValue, VECLEN / 2 + 1);
        DSP.IFFTo(fpInData, spInData);

        Filter.Put(fpInData, VECLEN);
    }
    Filter.Uninitialize();

    fileSamples = fopen("example2.samples1.txt", "wt");
    for (lSampleCntr = 0; lSampleCntr < VECLEN; lSampleCntr++)
        fprintf(fileSamples, "%g\n", fpOutData[lSampleCntr]);
    fclose(fileSamples);

    //Filter.InitializeHP(0.55f, 0.45f, 72.0f);
    //Filter.InitializeHP(0.11f, 0.10f, 72.0f);
    //Filter.InitializeHP(0.15f, 0.05f, 72.0f);
    Filter.Initialize(1024, (float *) NULL, 0.75f, DSP_FILT_DEF_BETAF);
    float fCorner = 0.105f;
    Filter.DesignHP(&fCorner);
    printf("Delay for filter 2 is %li\n", Filter.GetDelay());
    while (!Filter.Get(fpOutData, VECLEN))
    {
        /*for (lSampleCntr = 0; lSampleCntr < VECLEN; lSampleCntr++)
        {
            iRandomValue = rand();
            fpInData[lSampleCntr] = (iRandomValue & 1) ?
                (-iRandomValue * fScale) : (iRandomValue * fScale);
        }*/

        DSP.Set(spInData, sValue, VECLEN / 2 + 1);
        DSP.IFFTo(fpInData, spInData);

        Filter.Put(fpInData, VECLEN);
    }
    Filter.Uninitialize();

    fileSamples = fopen("example2.samples2.txt", "wt");
    for (lSampleCntr = 0; lSampleCntr < VECLEN; lSampleCntr++)
        fprintf(fileSamples, "%g\n", fpOutData[lSampleCntr]);
    fclose(fileSamples);

    return 0;
}
