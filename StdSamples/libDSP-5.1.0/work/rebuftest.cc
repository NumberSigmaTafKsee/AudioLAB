#include <cstdio>
#include "dsp/ReBuffer.hh"

int main (int argc, char *argv[])
{
    float fpIn[] = { 0, 1, 2, 3, 4 };
    float fpOut[10];
    clReBuffer RB;
    
    RB.Put(fpIn, 5);
    printf("%i ", (int) RB.Get(fpOut, 2));
    printf("%g %g\n", fpOut[0], fpOut[1]);
    printf("%i ", (int) RB.Get(fpOut, 2));
    printf("%g %g\n", fpOut[0], fpOut[1]);
    printf("%i\n", (int) RB.Get(fpOut, 5));

    RB.Put(fpIn, 5);
    printf("%i ", (int) RB.Get(fpOut, 5));
    printf("%g %g %g %g %g\n", fpOut[0], fpOut[1], fpOut[2], fpOut[3], fpOut[4]);

    printf("size: %i\n", RB.GetCount());

    RB.Put(fpIn, 5);
    RB.Put(fpIn, 5);
    printf("%i ", (int) RB.Get(fpOut, 10));
    printf("%g %g %g %g %g ", fpOut[0], fpOut[1], fpOut[2], fpOut[3], fpOut[4]);
    printf("%g %g %g %g %g\n", fpOut[5], fpOut[6], fpOut[7], fpOut[8], fpOut[9]);
    
    RB.Put(fpIn, 5);
    printf("%i ", (int) RB.Get(fpOut, 5));
    printf("%g %g %g %g %g\n", fpOut[0], fpOut[1], fpOut[2], fpOut[3], fpOut[4]);

    RB.Put(fpIn, 5);
    printf("%i ", (int) RB.Get(fpOut, 5));
    printf("%g %g %g %g %g\n", fpOut[0], fpOut[1], fpOut[2], fpOut[3], fpOut[4]);

        
    clReBuffer RB2;
    
    puts("---");

    RB2.Put(fpIn, 4);
    RB2.Put(fpIn, 4);
    printf("%i ", (int) RB2.Get(fpOut, 8));
    printf("%g %g %g %g ", fpOut[0], fpOut[1], fpOut[2], fpOut[3]);
    printf("%g %g %g %g\n", fpOut[4], fpOut[5], fpOut[6], fpOut[7]);
    RB2.Put(fpIn, 4);
    
    
    int i;
    int iValue = 0;
    int iValue2 = 0;
    float fpIn2[512];
    float fpOut2[819];
    clReBuffer RB3;
    
    puts("---");
    while (iValue < 4096)
    {
        for (i = 0; i < 512; i++)
            fpIn2[i] = (float) iValue++;
        RB3.Put(fpIn2, 512);
        while (RB3.Get(fpOut2, 819))
        {
            for (i = 0; i < 819; i++)
            {
                if (((int) fpOut2[i]) != iValue2)
                {
                    printf("[%i] == %i != %i\n", i, (int) fpOut2[i], iValue2);
                    return 1;
                }
                iValue2++;
            }
        }
    }

    return 0;
}
