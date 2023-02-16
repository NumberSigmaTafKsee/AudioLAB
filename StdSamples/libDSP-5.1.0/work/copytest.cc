#include <cstdio>

#include <MersenneTwister.h>

#include <dsp/DSPOp.hh>


int main (int argc, char *argv[])
{
    int iCntr;
    float *fpSrcData;
    float *fpDestData;
    double *dpSrcData;
    double *dpDestData;
    MTRand MTR;
    clDSPOp DSP;
    
    fpSrcData = new float[1048576 - 1];
    fpDestData = new float[1048576 - 1];
    dpSrcData = new double[1048576 - 1];
    dpDestData = new double[1048576 - 1];
    
    for (iCntr = 0; iCntr < 1048576 - 1; iCntr++)
    {
        fpSrcData[iCntr] = MTR.rand() * 2.0 + -1.0;
        dpSrcData[iCntr] = MTR.rand() * 2.0 + -1.0;
    }
    
    DSP.Copy(fpDestData, fpSrcData, 1048576 - 1);
    DSP.Copy(dpDestData, dpSrcData, 1048576 - 1);
    
    for (iCntr = 0; iCntr < 1048576 - 1; iCntr++)
    {
        if ((fpDestData[iCntr] != fpSrcData[iCntr]) ||
            (dpDestData[iCntr] != dpSrcData[iCntr]))
        {
            puts("Error!");
            break;
        }
    }

    delete fpSrcData;
    delete fpDestData;
    delete dpSrcData;
    delete dpDestData;
    
    return 0;
}
