#include <stdio.h>

#include "dsp/dspop.h"


int main (int argc, char *argv[])
{
    stSCplx sV1, sV2, sV3;

    dsp_init();

    sV1.R = 1.0f;
    sV1.I = 0.0f;
    sV2.R = 0.5f;
    sV2.I = 0.2f;

    dsp_caddf(&sV1, sV2, 1);

    printf("%g + %gi\n", sV1.R, sV1.I);

    sV1.R = 1.0f;
    sV1.I = 0.0f;
    sV2.R = 0.5f;
    sV2.I = 0.2f;
    sV3.R = 0.0f;
    sV3.I = 0.0f;

    dsp_cmuladdf(&sV3, &sV1, &sV2, 1);

    printf("%g + %gi\n", sV3.R, sV3.I);

    return 0;
}
