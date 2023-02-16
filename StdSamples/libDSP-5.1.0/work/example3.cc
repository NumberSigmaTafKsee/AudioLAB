#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cfloat>

#include <dsp/DSPOp.hh>

#define VECLEN      1024    // vector length
//#define VECPOS      (VECLEN / 128)
#define VECPOS      13


int main (int argc, char *argv[])
{
    long lSampleCntr;
    float fpOutData[VECLEN];
    stSCplx spInData[VECLEN / 2 + 1];
    clDSPOp DSP;

    for (lSampleCntr = 0; lSampleCntr < (VECLEN / 2 + 1); lSampleCntr++)
    {
        if (lSampleCntr % VECPOS == 0)
        {
            spInData[lSampleCntr].R = 1.0;
            spInData[lSampleCntr].I = 0.0;
        }
        else
        {
            spInData[lSampleCntr].R = 0.0;
            spInData[lSampleCntr].I = 0.0;
        }
        /*spInData[lSampleCntr].R = 1.0;
        spInData[lSampleCntr].I = 0.0;*/
    }
    DSP.FFTInitialize(VECLEN, true);
    DSP.IFFTo(fpOutData, spInData);
    DSP.FFTUninitialize();
    
    for (lSampleCntr = 0; lSampleCntr < VECLEN; lSampleCntr++)
    {
        printf("%g\n", fpOutData[lSampleCntr]);
    }
    
    return 0;
}
