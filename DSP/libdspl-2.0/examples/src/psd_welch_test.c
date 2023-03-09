#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"


#define N         8192
#define FS        1.0



void psd_plot(int argc, char* argv[], char* fn_png, char* fn_data)
{
    char cmd[512] = {0};
    void* hplot;  /* GNUPLOT handle               */
    /* plotting by GNUPLOT */
    gnuplot_create(argc, argv, 680, 480, fn_png, &hplot);
    gnuplot_cmd(hplot, "unset key");
    gnuplot_cmd(hplot, "set grid");
    gnuplot_cmd(hplot, "set xlabel 'frequency, Hz'");
    gnuplot_cmd(hplot, "set ylabel 'Power Spectral Density, [dB/Hz]'");
    gnuplot_cmd(hplot, "set yrange [-60: 40]");
    
    sprintf(cmd, "plot '%s' w l lt -1", fn_data);
    
    gnuplot_cmd(hplot, cmd);
    gnuplot_close(hplot);
}


int main(int argc, char* argv[])
{
    void* hdspl;  /* DSPL handle                  */

    hdspl = dspl_load();    /* Load DSPL function */
    random_t rnd = {0};     /* random structure   */
    
    double *x   = NULL;
    double *psd = NULL;
    double *frq = NULL;

    int k, err;

    x   = (double*)malloc(N*sizeof(double));
    psd = (double*)malloc(N*sizeof(double));
    frq = (double*)malloc(N*sizeof(double));
    
    /* input signal fill as noise -20 dB/Hz power spectrum density */
    random_init(&rnd, RAND_TYPE_MRG32K3A, NULL);
    randn(x, N, 0.0, 0.1, &rnd);
    
    /* x[k] = 0.1 * cos(2 * pi * 0.2 * k) + cos(2 * pi * 0.26 * k) + noise */
    for(k = 0; k < N; k++)
        x[k] += cos(M_2PI * 0.26 * (double)k) + 0.1*cos(M_2PI*0.2* (double)k); 
               
    /* Twosided PSD logscale magnitude n = 8192, nfft = 8192. 
       This case is periodogram */
       
    err = psd_welch(x, N, DSPL_WIN_BLACKMAN, 0, 8192, 4096,  NULL, FS, 
                       DSPL_FLAG_LOGMAG | DSPL_FLAG_PSD_TWOSIDED, psd, frq);
    writetxt(frq, psd, 8192, "dat/psd_welch_8192.txt");
    psd_plot(argc, argv,     "img/psd_welch_8192.png", 
                             "dat/psd_welch_8192.txt");


    err = psd_welch(x, N, DSPL_WIN_BLACKMAN, 0, 1024, 512,  NULL, FS, 
                       DSPL_FLAG_LOGMAG | DSPL_FLAG_PSD_TWOSIDED, psd, frq);
    writetxt(frq, psd, 1024, "dat/psd_welch_1024.txt");
    psd_plot(argc, argv,     "img/psd_welch_1024.png", 
                             "dat/psd_welch_1024.txt");

    err = psd_welch(x, N, DSPL_WIN_BLACKMAN, 0, 256, 128,  NULL, FS, 
                       DSPL_FLAG_LOGMAG | DSPL_FLAG_PSD_TWOSIDED, psd, frq);
    writetxt(frq, psd, 256, "dat/psd_welch_256.txt");
    psd_plot(argc, argv,    "img/psd_welch_256.png", 
                            "dat/psd_welch_256.txt");

    dspl_free(hdspl);      /* free dspl handle */
    if(x)
      free(x);
    if(psd)
      free(psd);
    if(frq)
      free(frq);
  return 0;
}


