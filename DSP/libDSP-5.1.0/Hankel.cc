/*

    Hankel transform (Abel transform followed by FFT)
    Copyright (C) 1995 Colorado School of Mines. All rights reserved.
    Copyright (C) 2001-2003 Jussi Laako

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

    Source converted from cwplib (by Dave Hale and Lydia Deng / Colorado 
    School of Mines).

    References:
    Hansen, E. W., 1985, Fast Hankel transform algorithm:  IEEE Trans. on
    Acoustics, Speech and Signal Processing, v. ASSP-33, n. 3, p. 666-671.
    (Beware of several errors in the equations in this paper!)

*/


#ifdef USE_INTEL_MATH
    #include <mathimf.h>
#else
    #include <math.h>
#endif
#include <float.h>

#include "dsp/Hankel.hh"


#define ABEL_NSE 9
const float fpH[ABEL_NSE] = {
    1.000000000000000000f,
    0.610926299405048390f,
    0.895089852938535935f,
    1.34082948787002865f,
    2.02532848558443890f,
    3.18110895533701843f,
    5.90898360396353794f,
    77.6000213494180286f,
    528.221800846070892f
};
const double dpH[ABEL_NSE] = {
    1.000000000000000000,
    0.610926299405048390,
    0.895089852938535935,
    1.34082948787002865,
    2.02532848558443890,
    3.18110895533701843,
    5.90898360396353794,
    77.6000213494180286,
    528.221800846070892
};
const float fpLambda[ABEL_NSE] = {
    0.000000000000000000f,
    -2.08424632126539366f,
    -5.78928630565552371f,
    -14.6268676854951032f,
    -35.0617158334443104f,
    -83.3258406398958158f,
    -210.358805421311445f,
    -6673.64911325382036f,
    -34897.7050244132261f
};
const double dpLambda[ABEL_NSE] = {
    0.000000000000000000,
    -2.08424632126539366,
    -5.78928630565552371,
    -14.6268676854951032,
    -35.0617158334443104,
    -83.3258406398958158,
    -210.358805421311445,
    -6673.64911325382036,
    -34897.7050244132261
};


void clHankel::InitAbel (const float *fpNullPtr)
{
    long i, j;
    float *a, *b0, *b1, fi, hj, lambdaj, scale, temp;

    a = (float *) A.Size(ABEL_NSE * lSize * sizeof(float));
    b0 = (float *) B0.Size(ABEL_NSE * lSize * sizeof(float));
    b1 = (float *) B1.Size(ABEL_NSE * lSize * sizeof(float));
    for (i = 1; i < lSize; ++i) 
    {
        fi = i + 1.0f;
        for (j = 0; j < ABEL_NSE; ++j) 
        {
            hj = fpH[j];
            lambdaj = fpLambda[j];
            a[i * ABEL_NSE + j] = temp = 
                powf(fi / (fi - 1.0f), lambdaj);
            temp *= fi / (fi - 1.0f);
            scale = 2.0f * hj * (fi - 1.0f) /
                ((lambdaj + 1.0f) * (lambdaj + 2.0f));
            b0[i * ABEL_NSE + j] = 
                scale * (fi - 1.0f + (lambdaj + 2.0f - fi) * temp);
            b1[i * ABEL_NSE + j] = 
                -scale * (lambdaj + 1.0f + fi - fi * temp);
        }
    }
}


void clHankel::InitAbel (const double *dpNullPtr)
{
    long i, j;
    double *a, *b0, *b1, fi, hj, lambdaj, scale, temp;

    a = (double *) A.Size(ABEL_NSE * lSize * sizeof(double));
    b0 = (double *) B0.Size(ABEL_NSE * lSize * sizeof(double));
    b1 = (double *) B1.Size(ABEL_NSE * lSize * sizeof(double));
    for (i = 1; i < lSize; ++i)
    {
        fi = i + 1.0;
        for (j = 0; j < ABEL_NSE; ++j)
        {
            hj = fpH[j];
            lambdaj = fpLambda[j];
            a[i * ABEL_NSE + j] = temp = 
                pow(fi / (fi - 1.0), lambdaj);
            temp *= fi / (fi - 1.0);
            scale = 2.0 * hj * (fi - 1.0) /
                ((lambdaj + 1.0) * (lambdaj + 2.0));
            b0[i * ABEL_NSE + j] = 
                scale * (fi - 1.0 + (lambdaj + 2.0 - fi) * temp);
            b1[i * ABEL_NSE + j] = 
                -scale * (lambdaj + 1.0 + fi - fi * temp);
        }
    }
}


void clHankel::DoAbel (float *fpRes, const float *fpSrc)
{
    long i, j;
    float sum, fi, fip1;
    float *a, *b0, *b1, xi[ABEL_NSE];

    a = (float *) A.GetPtr();
    b0 = (float *) B0.GetPtr();
    b1 = (float *) B1.GetPtr();
    fi = fpSrc[lSize - 1];
    fpRes[0] = 0.5f * fpSrc[0] + fi;
    sum = 0.0f;
    for (j = 0; j < ABEL_NSE; ++j) 
    {
        xi[j] = b1[(lSize - 1) * ABEL_NSE + j] * fi;
        sum += xi[j];
    }
    fpRes[lSize - 1] = sum;
    for (i = lSize - 2; i > 0; --i) 
    {
        fip1 = fi;
        fi = fpSrc[i];
        fpRes[0] += fi;
        sum = 0.0f;
        for (j = 0; j < ABEL_NSE; ++j) 
        {
            xi[j] = a[i * ABEL_NSE + j] * xi[j] + 
                b0[i * ABEL_NSE + j] * fip1 + 
                b1[i * ABEL_NSE + j] * fi;
            sum += xi[j];
        }
        fpRes[i] = sum;
    }
    fpRes[0] *= 2.0f;
}


void clHankel::DoAbel (double *dpRes, const double *dpSrc)
{
    long i, j;
    double sum, fi, fip1;
    double *a, *b0, *b1, xi[ABEL_NSE];

    a = (double *) A.GetPtr();
    b0 = (double *) B0.GetPtr();
    b1 = (double *) B1.GetPtr();
    fi = dpSrc[lSize - 1];
    dpRes[0] = 0.5 * dpSrc[0] + fi;
    sum = 0.0;
    for (j = 0; j < ABEL_NSE; ++j)
    {
        xi[j] = b1[(lSize - 1) * ABEL_NSE + j] * fi;
        sum += xi[j];
    }
    dpRes[lSize - 1] = sum;
    for (i = lSize - 2; i > 0; --i)
    {
        fip1 = fi;
        fi = dpSrc[i];
        dpRes[0] += fi;
        sum = 0.0;
        for (j = 0; j < ABEL_NSE; ++j)
        {
            xi[j] = a[i * ABEL_NSE + j] * xi[j] + 
                b0[i * ABEL_NSE + j] * fip1 + 
                b1[i * ABEL_NSE + j] * fi;
            sum += xi[j];
        }
        dpRes[i] = sum;
    }
    dpRes[0] *= 2.0;
}


void clHankel::UninitAbel ()
{
    A.Free();
    B0.Free();
    B1.Free();
}


clHankel::clHankel ()
{
    lSize = 0;
}


clHankel::~clHankel ()
{
    if (lSize != 0) Uninitialize();
}


void clHankel::Initialize (long lNewSize, const float *fpNullPtr)
{
    lSize = lNewSize / 2 + 1;
    lFFTSize = lNewSize;
    InitAbel(fpNullPtr);
    GX.Size(lFFTSize * sizeof(float));
    //GK.Size(lSize * sizeof(float) * 2);  // for PFAFFT
    GK.Size(lFFTSize * sizeof(stSCplx));  // for "normal" FFT
    fOutScale0 = (float) (1.0 / (2.0 * acos(-1.0)));
    fOutScale1 = (float) (-1.0 / (2.0 * acos(-1.0)));
    DSP.FFTInitialize(lFFTSize, false);
}


void clHankel::Initialize (long lNewSize, const double *dpNullPtr)
{
    lSize = lNewSize / 2 + 1;
    lFFTSize = lNewSize;
    InitAbel(dpNullPtr);
    GX.Size(lFFTSize * sizeof(double));
    //GK.Size(lSize * sizeof(double) * 2);  // for PFAFFT
    GK.Size(lFFTSize * sizeof(stDCplx));  // for "normal" FFT
    dOutScale0 = 1.0 / (2.0 * acos(-1.0));
    dOutScale1 = -1.0 / (2.0 * acos(-1.0));
    DSP.FFTInitialize(lFFTSize, false);
}


void clHankel::Uninitialize ()
{
    UninitAbel();
    GX.Free();
    GK.Free();
    lSize = 0;
}


void clHankel::Process0 (float *fpRes, const float *fpSrc)
{
    long i;
    float fPreScale;
    float *gx;
    //float *gk;
    stpSCplx gk;
#   ifdef __GNUG__
    stSCplx spTemp[lFFTSize];
#   else
    clDSPAlloc Temp(lFFTSize * sizeof(stSCplx));
    stpSCplx spTemp = Temp;
#   endif

    gx = GX;
    gk = GK;
    DoAbel(gx, fpSrc);
    /*for (i = lSize; i < lFFTSize; ++i)
        gx[i] = gx[lFFTSize - i];
    FFT.pfarc(-1, lFFTSize, gx, gk);
    for (i = 0; i < lSize; ++i)
        fpRes[i] = fOutScale0 * gk[i << 1];*/
    fPreScale = 2.0f / lSize;
    for (i = 0; i < (lFFTSize >> 1); i++)
    {
        spTemp[i].R = gx[i] * fPreScale;
        spTemp[i].I = 0;
    }
    for (i = (lFFTSize >> 1); i < lFFTSize; i++)
    {
        spTemp[i].R = gx[lFFTSize - i] * fPreScale;
        spTemp[i].I = 0;
    }
    DSP.IFFTo(gk, spTemp);
    for (i = 0; i < lSize; i++)
    {
        fpRes[i] = sqrt(gk[i].R * gk[i].R + gk[i].I * gk[i].I);
    }
}


void clHankel::Process0 (double *dpRes, const double *dpSrc)
{
    long i;
    double dPreScale;
    double *gx;
    //double *gk;
    stpDCplx gk;
#   ifdef __GNUG__
    stDCplx spTemp[lFFTSize];
#   else
    clDSPAlloc Temp(lFFTSize * sizeof(stDCplx));
    stpDCplx spTemp = Temp;
#   endif

    gx = GX;
    gk = GK;
    DoAbel(gx, dpSrc);
    /*for (i = lSize; i < lFFTSize; ++i)
        gx[i] = gx[lFFTSize - i];
    FFT.pfarc(-1, lFFTSize, gx, gk);
    for (i = 0; i < lSize; ++i)
        dpRes[i] = dOutScale0 * gk[i << 1];*/
    dPreScale = 2.0 / lSize;
    for (i = 0; i < (lFFTSize >> 1); i++)
    {
        spTemp[i].R = gx[i] * dPreScale;
        spTemp[i].I = 0;
    }
    for (i = (lFFTSize >> 1); i < lFFTSize; i++)
    {
        spTemp[i].R = gx[lFFTSize - i] * dPreScale;
        spTemp[i].I = 0;
    }
    DSP.IFFTo(gk, spTemp);
    for (i = 0; i < lSize; i++)
    {
        dpRes[i] = sqrt(gk[i].R * gk[i].R + gk[i].I * gk[i].I);
    }
}


void clHankel::Process1 (float *fpRes, const float *fpSrc)
{
    long i;
    float fPreScale;
    float *gx;
    //float *gk;
    stpSCplx gk;
#   ifdef __GNUG__
    stSCplx spTemp[lFFTSize];
#   else
    clDSPAlloc Temp(lFFTSize * sizeof(stSCplx));
    stpSCplx spTemp = Temp;
#   endif

    gx = GX;
    gk = GK;
    for (i = 1; i < lSize; ++i)
        gx[i] = fpSrc[i] / i;
    DoAbel(gx, gx);
    for (i = 0; i < lSize; ++i)
        gx[i] *= i;
    for (i = lSize; i < lFFTSize; ++i)
        gx[i] = -gx[lFFTSize - i];
    gx[lFFTSize / 2] = 0.0f;
    /*FFT.pfarc(-1, lFFTSize, gx, gk);
    for (i = 0; i < lSize; ++i)
        fpRes[i] = fOutScale1 * gk[(i << 1) + 1];*/
    fPreScale = 2.0f / lSize;
    for (i = 0; i < lFFTSize; i++)
    {
        spTemp[i].R = gx[i] * fPreScale;
        spTemp[i].I = 0;
    }
    DSP.IFFTo(gk, spTemp);
    for (i = 0; i < lSize; i++)
    {
        fpRes[i] = sqrt(gk[i].R * gk[i].R + gk[i].I * gk[i].I);
    }
}


void clHankel::Process1 (double *dpRes, const double *dpSrc)
{
    long i;
    double dPreScale;
    double *gx;
    //double *gk;
    stpDCplx gk;
#   ifdef __GNUG__
    stDCplx spTemp[lFFTSize];
#   else
    clDSPAlloc Temp(lFFTSize * sizeof(stDCplx));
    stpDCplx spTemp = Temp;
#   endif

    gx = GX;
    gk = GK;
    for (i = 1; i < lSize; ++i)
        gx[i] = dpSrc[i] / i;
    DoAbel(gx, gx);
    for (i = 0; i < lSize; ++i)
        gx[i] *= i;
    for (i = lSize; i < lFFTSize; ++i)
        gx[i] = -gx[lFFTSize - i];
    gx[lFFTSize / 2] = 0.0;
    /*FFT.pfarc(-1, lFFTSize, gx, gk);
    for (i = 0; i < lSize; ++i)
        dpRes[i] = dOutScale1 * gk[(i << 1) + 1];*/
    dPreScale = 2.0 / lSize;
    for (i = 0; i < lFFTSize; i++)
    {
        spTemp[i].R = gx[i] * dPreScale;
        spTemp[i].I = 0;
    }
    DSP.IFFTo(gk, spTemp);
    for (i = 0; i < lSize; i++)
    {
        dpRes[i] = sqrt(gk[i].R * gk[i].R + gk[i].I * gk[i].I);
    }
}
