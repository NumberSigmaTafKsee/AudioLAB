/*
MDCT/IMDCT is the most widely used filterbank in digital audio coding, e.g. MP3, AAC, WMA,
OGG Vorbis, ATRAC.

suppose input x and N=size(x,1)/2. the MDCT transform matrix is
    C=cos(pi/N*([0:2*N-1]'+.5+.5*N)*([0:N-1]+.5));
then MDCT spectrum for input x is
    y=C'*x;
A well known fast algorithm is based on FFT :
 (1) fold column-wisely the 2*N rows into N rows
 (2) complex arrange the N rows into N/2 rows
 (3) pre-twiddle, N/2-point complex fft, post-twiddle
 (4) reorder to form the MDCT spectrum
in fact, (2)-(4) is a fast DCT-IV algorithm.

Implementation of the algorithm can be found in faac, but a little bit mess to extract for
standalone use, and I ran into that problem. So I wrote some c codes to implement
MDCT/IMDCT for any length that is a multiple of 4. Hopefully they will be useful to people
here.

I benchmarked the codes using 3 FFT routines, FFT in faac, kiss_fft, and the awful FFTW.
MDCT based on FFTW is the fastest, 2048-point MDCT single precision 10^5 times in 1.54s,
about 50% of FFT in faac on my Petium IV 3G Hz.

An author of the FFTW, Steven G. Johnson, has a hard-coded fixed size MDCT of 256 input
samples(http://jdj.mit.edu/~stevenj/mdct_128nr.c). My code is 13% slower than his.


Using the codes is very simple:
 (1) init (declare first "extern void* mdctf_init(int)")
       void* m_plan = mdctf_init(N);
 (2) run mdct/imdct as many times as you wish
       mdctf(freq, time, m_plan);
 (3) free
       mdctf_free(m_plan);

Of course you need the the fftw library. On Linux, gcc options are "-O2 -lfftw3f -lm".
This is single precision.


Enjoy :)
*/

/*********************************************************
  MDCT/IMDCT of 4x length, Single Precision, based on FFTW
  shuhua dot zhang at gmail dot com
  Dept. of E.E., Tsinghua University
 *********************************************************/

  #include <stdio.h>
  #include <stdlib.h>
  #include <math.h>
  #include <fftw3.h>
  
  
  typedef struct {
      int              N;              // Number of time data points
      float*           twiddle;        // Twiddle factor
      fftwf_complex*   fft_in;         // fft workspace, input
      fftwf_complex*   fft_out;        // fft workspace, output
      fftwf_plan       fft_plan;       // fft configuration
  } mdctf_plan;
  
  
  mdctf_plan* mdctf_init(int N);
  void mdctf_free(mdctf_plan* m_plan);
  void mdctf(float* mdct_line, float* time_signal, mdctf_plan* m_plan);
  void imdctf(float* time_signal, float* mdct_line, mdctf_plan* m_plan);
  
  
  mdctf_plan* mdctf_init(int N)
  {
      mdctf_plan* m_plan;
      double alpha, omiga, scale;
      int    n;
  
      if( 0x00 != (N & 0x03))
      {
          fprintf(stderr, " Expecting N a multiple of 4\n");
          return NULL;
      }
  
      m_plan = (mdctf_plan*) malloc(sizeof(mdctf_plan));
  
      m_plan->N = N;
  
      m_plan->twiddle = (float*) malloc(sizeof(float) * N >> 1);
      alpha = 2.f * M_PI / (8.f * N);
      omiga = 2.f * M_PI / N;
      scale = sqrt(sqrt(2.f / N));
      for(n = 0; n < (N >> 2); n++)
      {
          m_plan->twiddle[2*n+0] = (float) (scale * cos(omiga * n + alpha));
          m_plan->twiddle[2*n+1] = (float) (scale * sin(omiga * n + alpha));
      }
  
      m_plan->fft_in   = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * N >> 2);
      m_plan->fft_out  = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * N >> 2);
      m_plan->fft_plan = fftwf_plan_dft_1d(N >> 2,
                                           m_plan->fft_in,
                                           m_plan->fft_out,
                                           FFTW_FORWARD,
                                           FFTW_MEASURE);
  
      return m_plan;
  
  }
  
  
  void mdctf_free(mdctf_plan* m_plan)
  {
      fftwf_destroy_plan(m_plan->fft_plan);
      fftwf_free(m_plan->fft_in);
      fftwf_free(m_plan->fft_out);
      free(m_plan->twiddle);
      free(m_plan);
  }
  
  
  void mdctf(float* mdct_line, float* time_signal, mdctf_plan* m_plan)
  {
      float  *xr, *xi, r0, i0;
      float  *cos_tw, *sin_tw, c, s;
      int     N4, N2, N34, N54, n;
  
      N4  = (m_plan->N) >> 2;
      N2  = 2 * N4;
      N34 = 3 * N4;
      N54 = 5 * N4;
  
      cos_tw = m_plan->twiddle;
      sin_tw = cos_tw + 1;
  
      /* odd/even folding and pre-twiddle */
      xr = (float*) m_plan->fft_in;
      xi = xr + 1;
      for(n = 0; n < N4; n += 2)
      {
          r0 = time_signal[N34-1-n] + time_signal[N34+n];
          i0 = time_signal[N4+n]    - time_signal[N4-1-n];
  
          c = cos_tw[n];
          s = sin_tw[n];
  
          xr[n] = r0 * c + i0 * s;
          xi[n] = i0 * c - r0 * s;
      }
  
      for(; n < N2; n += 2)
      {
          r0 = time_signal[N34-1-n] - time_signal[-N4+n];
          i0 = time_signal[N4+n]    + time_signal[N54-1-n];
  
          c = cos_tw[n];
          s = sin_tw[n];
  
          xr[n] = r0 * c + i0 * s;
          xi[n] = i0 * c - r0 * s;
      }
  
      /* complex FFT of N/4 long */
      fftwf_execute(m_plan->fft_plan);
  
      /* post-twiddle */
      xr = (float*) m_plan->fft_out;
      xi = xr + 1;
      for(n = 0; n < N2; n += 2)
      {
          r0 = xr[n];
          i0 = xi[n];
  
          c = cos_tw[n];
          s = sin_tw[n];
  
          mdct_line[n]      = - r0 * c - i0 * s;
          mdct_line[N2-1-n] = - r0 * s + i0 * c;
      }
  }
  
  
  void imdctf(float* time_signal, float* mdct_line, mdctf_plan* m_plan)
  {
      float  *xr, *xi, r0, i0, r1, i1;
      float  *cos_tw, *sin_tw, c, s;
      int     N4, N2, N34, N54, n;
  
      N4  = (m_plan->N) >> 2;
      N2  = 2 * N4;
      N34 = 3 * N4;
      N54 = 5 * N4;
  
      cos_tw = m_plan->twiddle;
      sin_tw = cos_tw + 1;
  
         /* pre-twiddle */
      xr = (float*) m_plan->fft_in;
      xi = xr + 1;
      for(n = 0; n < N2; n += 2)
      {
          r0 =  mdct_line[n];
          i0 =  mdct_line[N2-1-n];
  
          c = cos_tw[n];
          s = sin_tw[n];
  
              xr[n] = -2.f * (i0 * s + r0 * c);
          xi[n] = -2.f * (i0 * c - r0 * s);
      }
  
      /* complex FFT of N/4 long */
      fftwf_execute(m_plan->fft_plan);
  
      /* odd/even expanding and post-twiddle */
      xr = (float*) m_plan->fft_out;
      xi = xr + 1;
      for(n = 0; n < N4; n += 2)
      {
          r0 = xr[n];
          i0 = xi[n];
  
          c = cos_tw[n];
          s = sin_tw[n];
  
          r1 = r0 * c + i0 * s;
          i1 = r0 * s - i0 * c;
  
          time_signal[N34-1-n] =  r1;
          time_signal[N34+n]   =  r1;
          time_signal[N4+n]    =  i1;
          time_signal[N4-1-n]  = -i1;
      }
  
      for(; n < N2; n += 2)
      {
          r0 = xr[n];
          i0 = xi[n];
  
          c = cos_tw[n];
          s = sin_tw[n];
  
          r1 = r0 * c + i0 * s;
          i1 = r0 * s - i0 * c;
  
          time_signal[N34-1-n] =  r1;
          time_signal[-N4+n]   = -r1;
          time_signal[N4+n]    =  i1;
          time_signal[N54-1-n] =  i1;
      }
  }