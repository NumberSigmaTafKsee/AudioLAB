#include <cstdio>

#include <dsp/DSPVector.hh>


int main (int argc, char *argv[])
{
    float fpTestData[] = { 1.0, 2.0, 3.0, 4.0, 5.0 };
    clDSPVector<float> dspvTest;
    clDSPVector<float> dspvTest2;

    dspvTest.Put(fpTestData, 5);
    dspvTest2 = dspvTest;
    dspvTest2 += dspvTest;
    dspvTest2 += 1.0;
    
    printf("%li\n", dspvTest2.Size());
    printf("%f %f %f %f %f\n", 
        dspvTest2[0],
        dspvTest2[1],
        dspvTest2[2],
        dspvTest2[3],
        dspvTest2[4]);

    float fMin, fMax;
    float fStdDev, fMean;
    dspvTest2.MinMax(fMin, fMax);
    printf("%f %f\n", fMin, fMax);
    dspvTest2.StdDev(fStdDev, fMean);
    printf("%f %f\n", fStdDev, fMean);

    stSCplx spCplxTestData[] = { 
        { 1.0, -1.0 }, 
        { 0.5, -0.5 },
        { 1.0, 0.0 },
        { 0.0, -1.0 },
        { 1.0, 1.0 } };
    clDSPVector<stSCplx> dspvCTest;
    clDSPVector<stSCplx> dspvCTest2;

    dspvCTest.Put(spCplxTestData, 5);
    dspvCTest2 = dspvCTest;

    dspvCTest2 *= dspvCTest;

    printf("%li\n", dspvCTest2.Size());
    printf("(%f, %f) (%f, %f) (%f, %f) (%f, %f) (%f, %f)\n",
        dspvCTest2[0].R, dspvCTest2[0].I,
        dspvCTest2[1].R, dspvCTest2[1].I,
        dspvCTest2[2].R, dspvCTest2[2].I,
        dspvCTest2[3].R, dspvCTest2[3].I,
        dspvCTest2[4].R, dspvCTest2[4].I);

    clDSPVector<float> dspvTest3;

    dspvTest3.WinBartlett(5);
    //dspvTest3.WinDolphChebyshev(0.005623413f, 5);

    printf("%li\n", dspvTest3.Size());
    printf("%f %f %f %f %f\n",
        dspvTest3[0],
        dspvTest3[1],
        dspvTest3[2],
        dspvTest3[3],
        dspvTest3[4]);

    dspvTest3.WinDolphChebyshev(0.005623413f, 5);

    printf("%li\n", dspvTest3.Size());
    printf("%f %f %f %f %f\n",
        dspvTest3[0],
        dspvTest3[1],
        dspvTest3[2],
        dspvTest3[3],
        dspvTest3[4]);

    return 0;
}
