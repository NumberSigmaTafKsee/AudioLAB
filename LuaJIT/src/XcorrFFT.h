#pragma once

#include <complex.h>

#include "W:\abid\USCT\SignalDetection\SignalDetection\fftw\include\fftw3.h"
#include "stdafx.h"


#define REAL 0
#define IMAG 1

//void fft(fftw_complex *input, fftw_complex *output, int N);
//void ifft(fftw_complex *input, fftw_complex* output, int N);
//void displayComplex(fftw_complex *y, int N);
//void displayReal(fftw_complex *y, int N);

//void MyXcorr_fft_1D(double *tab1_r, double *tab2_r, int n1, double *Result_r, double *lags);
void MyXcorr_fft(float *tab1_r, float *tab2_r, int n1, int m, float *Result_r, int nthreads);
void Myinterp_fft_Debug(float *tab1_r, float *tab2_r, int n1, int m, float *Result_r, int nthreads, int RESAMPLE_FACTOR);
void Myinterp_fft_SimpleVersion(float *tab1_r, float *tab2_r, int n1, int m, float *Result_r, int nthreads, int RESAMPLE_FACTOR);
//search for the maximum
void max_index_2Darray_according_1D(float *a, int n, int m, int * start, int *end, int * max_i);

void Myinterp_fft_DebugTHREADS(float *tab1_r, float *tab2_r, int n1, int m, float *Result_r, int nthreads, int RESAMPLE_FACTOR);
