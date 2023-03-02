#include <fftw3.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

typedef struct
{
    int inFrameSize;
    int inWindowSize;
    int inSampleRate;
    float *inWindowing;
    fftwf_plan inPlan;
    int outFrameSize;
    int outWindowSize;
    int outSampleRate;
    float *outWindowing;
    fftwf_plan outPlan;
    float *inFifo;
    float *synthesisMem;
    fftwf_complex *samples;
    int pos;
} FFT_Resampler_Handle;

#define max(a,b) (a > b? a:b)

void FFT_Resampler_Free(FFT_Resampler_Handle *handle)
{
    if (handle) {
        if (handle->inFifo) {
            free(handle->inFifo);
            handle->inFifo = NULL;
        }

        if (handle->inPlan) {
            fftwf_destroy_plan(handle->inPlan);
            handle->inPlan = NULL;
        }

        if (handle->outPlan) {
            fftwf_destroy_plan(handle->outPlan);
            handle->outPlan = NULL;
        }

        if (handle->samples) {
            fftw_free(handle->samples);
            handle->samples = NULL;
        }

        if (handle->synthesisMem) {
            fftw_free(handle->synthesisMem);
            handle->synthesisMem = NULL;
        }

        if (handle->inWindowing) {
            fftw_free(handle->inWindowing);
            handle->inWindowing = NULL;
        }

        if (handle->outWindowing) {
            fftw_free(handle->outWindowing);
            handle->outWindowing = NULL;
        }
    }
}

void FFT_Resampler_reset(FFT_Resampler_Handle *handle)
{
    if (handle)
        handle->pos = 0;
}

int FFT_Resampler_Init(FFT_Resampler_Handle *handle, size_t inSampleRate, size_t outSampleRate, size_t nFFT)
{
    if (handle) {
        handle->pos = 0;
        if (outSampleRate < inSampleRate) {
            nFFT = nFFT * inSampleRate * 128 / outSampleRate;
        }
        else {
            nFFT = nFFT * outSampleRate * 128 / inSampleRate;
        }
        nFFT += (nFFT % 2);

        handle->inFrameSize = nFFT;
        handle->inWindowSize = nFFT * 2;
        handle->inSampleRate = inSampleRate;
        handle->outSampleRate = outSampleRate;
        handle->outFrameSize = handle->inFrameSize * outSampleRate / inSampleRate;
        handle->outWindowSize = handle->outFrameSize * 2;        

        handle->outWindowing = (float *) fftwf_alloc_real(handle->outFrameSize);
        handle->inFifo = (float *) fftwf_alloc_real(max(handle->inWindowSize, handle->outWindowSize));
        handle->samples = (fftwf_complex *) fftwf_alloc_complex(max(handle->inWindowSize, handle->outWindowSize));
        handle->inWindowing = (float *) fftwf_alloc_real(handle->inFrameSize);
        handle->synthesisMem = (float *) fftwf_alloc_real(handle->outFrameSize);
                
        handle->inPlan = fftwf_plan_dft_r2c_1d(handle->inWindowSize,handle->inFifo,handle->samples,FFTW_ESTIMATE);        
        handle->outPlan = fftwf_plan_dft_c2r_1d(handle->outWindowSize,handle->samples,handle->synthesisMem,FFTW_ESTIMATE);
        
        
        if ((handle->inFifo == NULL) || (handle->inPlan == NULL) || (handle->outPlan == NULL)
            || (handle->samples == NULL)
            || (handle->synthesisMem == NULL) || (handle->inWindowing == NULL) || (handle->outWindowing == NULL)
            ) {
            FFT_Resampler_Free(handle);
            return 0;
        }
        float norm = 1.0f / handle->inWindowSize;
        for (size_t i = 0; i < handle->inFrameSize; i++) {
            double t = sin(.5 * M_PI * (i + .5) / handle->inFrameSize);
            handle->inWindowing[i] = (float) sin(.5 * M_PI * t * t) * norm;
        }
        for (size_t i = 0; i < handle->outFrameSize; i++) {
            double t = sin(.5 * M_PI * (i + .5) / handle->outFrameSize);
            handle->outWindowing[i] = (float) sin(.5 * M_PI * t * t);
        }
        return 1;
    }
    return 0;
}

int FFT_Resampler_Proc(FFT_Resampler_Handle *handle, const float *input, float *output)
{
    if ((input == NULL) || (handle == NULL) || (output == NULL)) {
        return -1;
    }
    float *inFifo = handle->inFifo;
    float *synthesis_mem = handle->synthesisMem;
    for (size_t i = 0; i < handle->inFrameSize; i++) {
        inFifo[i] *= handle->inWindowing[i];
        inFifo[handle->inWindowSize - 1 - i] = input[handle->inFrameSize - 1 - i] * handle->inWindowing[i];
    }
    fftwf_execute(handle->inPlan);
    if (handle->outWindowSize < handle->inWindowSize) {
        int half_output = (handle->outWindowSize / 2);
        int diff_size = handle->inWindowSize - handle->outWindowSize;
        memset(handle->samples + half_output, 0, diff_size * sizeof(fftw_complex));
    }
    else if (handle->outWindowSize > handle->inWindowSize) {
        int half_input = handle->inWindowSize / 2;
        int diff_size = handle->outWindowSize - handle->inWindowSize;
        memmove(handle->samples + half_input + diff_size, handle->samples + half_input,
                half_input * sizeof(fftw_complex));
        memset(handle->samples + half_input, 0, diff_size * sizeof(fftw_complex));
    }
    fftwf_execute(handle->outPlan);
    for (size_t i = 0; i < handle->outFrameSize; i++) {
        output[i] = inFifo[i] * handle->outWindowing[i] + synthesis_mem[i];
        inFifo[handle->outWindowSize - 1 - i] *= handle->outWindowing[i];
    }
    memcpy(synthesis_mem, inFifo + handle->outFrameSize, handle->outFrameSize * sizeof(float));
    memcpy(inFifo, input, handle->inFrameSize * sizeof(float));
    if (handle->pos == 0) {
        handle->pos++;
        return 0;
    }
    handle->pos++;
    return 1;
}