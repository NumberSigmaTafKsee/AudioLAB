#include <cstdio>
#include <cmath>

#include <dsp/DSPVector.hh>


//#define WINSIZE 1024
//#define WINSIZE 8
#define WINSIZE 32


int main (int argc, char *argv[])
{
    clDSPVector<double> dspvWin;

    //dspvWin.WinBartlett(WINSIZE);
    //dspvWin.WinHanning(WINSIZE);
    //dspvWin.WinHamming(WINSIZE);
    //dspvWin.WinBlackman(WINSIZE, 0.0);
    //dspvWin.WinExactBlackman(WINSIZE);
    //dspvWin.WinBlackmanHarris(WINSIZE);
    //dspvWin.WinCosTapered(WINSIZE);
    //dspvWin.WinFlatTop(WINSIZE);
    //dspvWin.WinTukey(WINSIZE);
    //dspvWin.WinKaiser(7.04f, WINSIZE);
    //dspvWin.WinKaiser(14.96454265f, WINSIZE);
    //dspvWin.WinKaiser(sqrtf(2.0f), WINSIZE);
    //dspvWin.WinKaiserBessel(4.0f, WINSIZE);
    //dspvWin.WinDolphChebyshev(1.0 / pow(2.0, 24.0), WINSIZE); //0.005623413
    dspvWin.WinDolphChebyshev(1.0 / pow(10.0, 3.0 * WINSIZE / 20.0), WINSIZE);

    for (int i = 0; i < WINSIZE; i++)
    {
        fprintf(stdout, " %g", dspvWin[i]);
    }
    fprintf(stdout, "\n");

    return 0;
}
