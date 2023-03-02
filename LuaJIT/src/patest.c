#include <portaudio.h>
#include <stdio.h>

void main()
{
    PaStreamParameters outputParameters;
    PaStream *stream;
    PaError err;
    int i;

    err = Pa_Initialize();
    printf("%d\n",err);

}