/*

    FFT filter
    Copyright (C) 1999-2003 Jussi Laako

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

*/


#include <cstdio>
#include <cstdlib>
#ifdef USE_INTEL_MATH
    #include <mathimf.h>
#else
    #include <math.h>
#endif
#include <float.h>

#include "dsp/Filter.hh"


void clFilter::InitCoeffsS ()
{
    long lCoeffCntr;
    stpSCplx spSCCoeffs = CCoeffs;

    for (lCoeffCntr = 0; lCoeffCntr < lSpectPoints; lCoeffCntr++)
    {
        spSCCoeffs[lCoeffCntr].R = 1.0f;
        spSCCoeffs[lCoeffCntr].I = 0.0f;
    }
    ReadyFilterS();
}


void clFilter::InitCoeffsD ()
{
    long lCoeffCntr;
    stpDCplx spDCCoeffs = CCoeffs;

    for (lCoeffCntr = 0; lCoeffCntr < lSpectPoints; lCoeffCntr++)
    {
        spDCCoeffs[lCoeffCntr].R = 1.0;
        spDCCoeffs[lCoeffCntr].I = 0.0;
    }
    ReadyFilterD();
}


void clFilter::ReadyFilterS ()
{
    float *fpCoeffWin = CoeffWin;
    float *fpProc = Proc;
    stpSCplx spSCCoeffs = CCoeffs;
    
    fpCoeffWin = CoeffWin;
    IFFTo(fpProc, spSCCoeffs);
    Mul(fpProc, fpCoeffWin, lFFTSize);
    FFTi(spSCCoeffs, fpProc);
}


void clFilter::ReadyFilterD ()
{
    double *dpCoeffWin = CoeffWin;
    double *dpProc = Proc;
    stpDCplx spDCCoeffs = CCoeffs;

    dpCoeffWin = CoeffWin;
    IFFTo(dpProc, spDCCoeffs);
    Mul(dpProc, dpCoeffWin, lFFTSize);
    FFTi(spDCCoeffs, dpProc);
}


float clFilter::GetKaiserBeta (float fAlpha)
{
    if (fAlpha < 21.0f)
    {
        return 0.0f;
    }
    else if (fAlpha >= 21.0f && fAlpha <= 50.0f)
    {
        return (0.5842f * powf(fAlpha - 21.0f, 0.4f) + 
            0.07886f * (fAlpha - 21.0f));
    }
    else
    {
        return (0.1102f * (fAlpha - 8.7f));
    }
}


double clFilter::GetKaiserBeta (double dAlpha)
{
    if (dAlpha < 21.0)
    {
        return 0.0;
    }
    else if (dAlpha >= 21.0 && dAlpha <= 50.0)
    {
        return (0.5842 * pow(dAlpha - 21.0, 0.4) + 
            0.07886 * (dAlpha - 21.0));
    }
    else
    {
        return (0.1102 * (dAlpha - 8.7));
    }
}


clFilter::clFilter ()
{
    bInitialized = false;
}


clFilter::~clFilter ()
{
    if (bInitialized) Uninitialize();
}


bool clFilter::Initialize (long lWindowSize, const float *fpFiltCoeffs,
    float fOverlap, float fBeta, int iSmoothWindow)
{
    float *fpCoeffWin;
    float *fpTemp;
    double *dpTemp2;
    clDSPAlloc Temp;
    clDSPAlloc Temp2;

    if (bInitialized) Uninitialize();
    bInitialized = true;
    lFFTSize = lWindowSize;
    lHalfSize = lFFTSize / 2;
    lOldSize = (long) ((float) lWindowSize * fOverlap + 0.5f);
    lNewSize = lWindowSize - lOldSize;
    lSpectPoints = lFFTSize / 2 + 1;
    Prev.Size(lOldSize * sizeof(float));
    Proc.Size(lFFTSize * sizeof(float));
    CCoeffs.Size(lSpectPoints * sizeof(stSCplx));
    CProc.Size(lSpectPoints * sizeof(stSCplx));
    CoeffWin.Size(lFFTSize * sizeof(float));

    fpTemp = (float *) Temp.Size(lFFTSize * sizeof(float));
    dpTemp2 = (double *) Temp2.Size(lFFTSize * sizeof(double));
    fpCoeffWin = CoeffWin;
    switch (iSmoothWindow)
    {
        case FILTER_SMOOTH_NONE:
            Set(dpTemp2, 1.0, lFFTSize);
            break;
        case FILTER_SMOOTH_KAISER:
            WinKaiser(dpTemp2, (double) fBeta, lFFTSize);
            break;
        case FILTER_SMOOTH_KAISER_BESSEL:
            WinKaiserBessel(dpTemp2, 4.0, lFFTSize);
            break;
        case FILTER_SMOOTH_DOLPH_CHEBYSHEV:
            WinDolphChebyshev(dpTemp2, 1.0 / pow(2.0, 32.0), lFFTSize);
            break;
    }
    Convert(fpTemp, dpTemp2, lFFTSize);
    Copy(fpCoeffWin, &fpTemp[lHalfSize], lHalfSize);
    Copy(&fpCoeffWin[lHalfSize], fpTemp, lHalfSize);

    FFTInitialize(lFFTSize, true);
    if (fpFiltCoeffs != NULL)
    {
        SetCoeffs(fpFiltCoeffs);
    }
    else
    {
        InitCoeffsS();
    }
    Zero((float *) Prev, lOldSize);

    return true;
}


bool clFilter::Initialize (long lWindowSize, const double *dpFiltCoeffs,
    double dOverlap, double dBeta, int iSmoothWindow)
{
    double *dpTemp;
    double *dpCoeffWin;
    clDSPAlloc Temp;

    if (bInitialized) Uninitialize();
    bInitialized = true;
    lFFTSize = lWindowSize;
    lHalfSize = lFFTSize / 2;
    lOldSize = (long) ((double) lWindowSize * dOverlap + 0.5);
    lNewSize = lWindowSize - lOldSize;
    lSpectPoints = lFFTSize / 2 + 1;
    Prev.Size(lOldSize * sizeof(double));
    Proc.Size(lFFTSize * sizeof(double));
    CCoeffs.Size(lSpectPoints * sizeof(stDCplx));
    CProc.Size(lSpectPoints * sizeof(stDCplx));
    CoeffWin.Size(lFFTSize * sizeof(double));

    dpTemp = (double *) Temp.Size(lFFTSize * sizeof(double));
    dpCoeffWin = CoeffWin;
    switch (iSmoothWindow)
    {
        case FILTER_SMOOTH_NONE:
            Set(dpTemp, 1.0, lFFTSize);
            break;
        case FILTER_SMOOTH_KAISER:
            WinKaiser(dpTemp, dBeta, lFFTSize);
            break;
        case FILTER_SMOOTH_KAISER_BESSEL:
            WinKaiserBessel(dpTemp, 4.0, lFFTSize);
            break;
        case FILTER_SMOOTH_DOLPH_CHEBYSHEV:
            WinDolphChebyshev(dpTemp, 1.0 / pow(2.0, 32.0), lFFTSize);
            break;
    }
    Copy(dpCoeffWin, &dpTemp[lHalfSize], lHalfSize);
    Copy(&dpCoeffWin[lHalfSize], dpTemp, lHalfSize);

    FFTInitialize(lFFTSize, true);
    if (dpFiltCoeffs != NULL)
    {
        SetCoeffs(dpFiltCoeffs);
    }
    else
    {
        InitCoeffsD();
    }
    Zero((double *) Prev, lOldSize);

    return true;
}


bool clFilter::InitFloat (long lWindowSize)
{
    float *fpNoCoeffs = NULL;

    return Initialize(lWindowSize, fpNoCoeffs);
}


bool clFilter::InitDouble (long lWindowSize)
{
    double *dpNoCoeffs = NULL;

    return Initialize(lWindowSize, dpNoCoeffs);
}


bool clFilter::InitializeLP (float fPassBand, float fStopBand,
    float fRippleRatio, float fOverlap)
{
    long lWindowSize;
    long lTwosExp;
    float fDeltaOmega;
    float fBeta;
    float fCorner;
    float *fpNullPtr = NULL;
    
    if (fPassBand >= fStopBand) return false;

    fCorner = (fPassBand + fStopBand) * 0.5f;
    fDeltaOmega = fStopBand - fPassBand;
    lTwosExp = (long)
        (logf((fRippleRatio - 8.0f) / (2.285f * fDeltaOmega)) / 
        logf(2.0f) + 0.5f);
    lWindowSize = (long) powf(2.0f, lTwosExp);
    fBeta = GetKaiserBeta(fRippleRatio);
    
    if (!Initialize(lWindowSize, fpNullPtr, fOverlap, fBeta, 
        FILTER_SMOOTH_KAISER)) return false;
    DesignLP(&fCorner);

    return true;
}


bool clFilter::InitializeLP (double dPassBand, double dStopBand,
    double dRippleRatio, double dOverlap)
{
    long lWindowSize;
    long lTwosExp;
    double dDeltaOmega;
    double dBeta;
    double dCorner;
    double *dpNullPtr = NULL;
    
    if (dPassBand >= dStopBand) return false;

    dCorner = (dPassBand + dStopBand) * 0.5;
    dDeltaOmega = dStopBand - dPassBand;
    lTwosExp = (long)
        (log((dRippleRatio - 8.0) / (2.285 * dDeltaOmega)) / 
        log(2.0) + 0.5);
    lWindowSize = (long) pow(2.0, lTwosExp);
    dBeta = GetKaiserBeta(dRippleRatio);
    
    if (!Initialize(lWindowSize, dpNullPtr, dOverlap, dBeta, 
        FILTER_SMOOTH_KAISER)) return false;
    DesignLP(&dCorner);

    return true;
}


bool clFilter::InitializeHP (float fPassBand, float fStopBand,
    float fRippleRatio, float fOverlap)
{
    long lWindowSize;
    long lTwosExp;
    float fDeltaOmega;
    float fBeta;
    float fCorner;
    float *fpNullPtr = NULL;
    
    if (fPassBand <= fStopBand) return false;

    fCorner = (fPassBand + fStopBand) * 0.5f;
    fDeltaOmega = fPassBand - fStopBand;
    lTwosExp = (long)
        (logf((fRippleRatio - 8.0f) / (2.285f * fDeltaOmega)) / 
        logf(2.0f) + 0.5f);
    lWindowSize = (long) powf(2.0f, lTwosExp);
    fBeta = GetKaiserBeta(fRippleRatio);
    
    if (!Initialize(lWindowSize, fpNullPtr, fOverlap, fBeta, 
        FILTER_SMOOTH_KAISER)) return false;
    DesignHP(&fCorner);

    return true;
}


bool clFilter::InitializeHP (double dPassBand, double dStopBand,
    double dRippleRatio, double dOverlap)
{
    long lWindowSize;
    long lTwosExp;
    double dDeltaOmega;
    double dBeta;
    double dCorner;
    double *dpNullPtr = NULL;
    
    if (dPassBand <= dStopBand) return false;

    dCorner = (dPassBand + dStopBand) * 0.5;
    dDeltaOmega = dPassBand - dStopBand;
    lTwosExp = (long)
        (log((dRippleRatio - 8.0) / (2.285 * dDeltaOmega)) / 
        log(2.0) + 0.5);
    lWindowSize = (long) pow(2.0, lTwosExp);
    dBeta = GetKaiserBeta(dRippleRatio);
    
    if (!Initialize(lWindowSize, dpNullPtr, dOverlap, dBeta, 
        FILTER_SMOOTH_KAISER)) return false;
    DesignHP(&dCorner);

    return true;
}


void clFilter::Uninitialize ()
{
    FFTUninitialize();
    CoeffWin.Free();
    Prev.Free();
    Proc.Free();
    CCoeffs.Free();
    CProc.Free();
    InBuf.Clear();
    OutBuf.Clear();
}


void clFilter::SetCoeffs (const float *fpFiltCoeffs)
{
    long lCoeffCntr;
    stpSCplx spSCCoeffs = CCoeffs;

    for (lCoeffCntr = 0; lCoeffCntr < lSpectPoints; lCoeffCntr++)
    {
        spSCCoeffs[lCoeffCntr].R = fpFiltCoeffs[lCoeffCntr];
        spSCCoeffs[lCoeffCntr].I = 0.0f;
    }
    ReadyFilterS();
}


void clFilter::SetCoeffs (const double *dpFiltCoeffs)
{
    long lCoeffCntr;
    stpDCplx spDCCoeffs = CCoeffs;

    for (lCoeffCntr = 0; lCoeffCntr < lSpectPoints; lCoeffCntr++)
    {
        spDCCoeffs[lCoeffCntr].R = dpFiltCoeffs[lCoeffCntr];
        spDCCoeffs[lCoeffCntr].I = 0.0;
    }
    ReadyFilterD();
}


void clFilter::SetCoeffs (const stpSCplx spSFiltCoeffs, bool bSmooth)
{
    stpSCplx spSCCoeffs = CCoeffs;

    Copy(spSCCoeffs, spSFiltCoeffs, lSpectPoints);
    if (bSmooth) ReadyFilterS();
}


void clFilter::SetCoeffs (const stpDCplx spDFiltCoeffs, bool bSmooth)
{
    stpDCplx spDCCoeffs = CCoeffs;

    Copy(spDCCoeffs, spDFiltCoeffs, lSpectPoints);
    if (bSmooth) ReadyFilterD();
}


void clFilter::GetCoeffs (float *fpFiltCoeffs)
{
    long lCoeffCntr;
    stpSCplx spSCCoeffs = CCoeffs;

    for (lCoeffCntr = 0; lCoeffCntr < lSpectPoints; lCoeffCntr++)
    {
        fpFiltCoeffs[lCoeffCntr] = spSCCoeffs[lCoeffCntr].R;
    }
}


void clFilter::GetCoeffs (double *dpFiltCoeffs)
{
    long lCoeffCntr;
    stpDCplx spDCCoeffs = CCoeffs;

    for (lCoeffCntr = 0; lCoeffCntr < lSpectPoints; lCoeffCntr++)
    {
        dpFiltCoeffs[lCoeffCntr] = spDCCoeffs[lCoeffCntr].R;
    }
}


void clFilter::GetCoeffs (stpSCplx spSFiltCoeffs)
{
    stpSCplx spSCCoeffs = CCoeffs;

    Copy(spSFiltCoeffs, spSCCoeffs, lSpectPoints);
}


void clFilter::GetCoeffs (stpDCplx spDFiltCoeffs)
{
    stpDCplx spDCCoeffs = CCoeffs;

    Copy(spDFiltCoeffs, spDCCoeffs, lSpectPoints);
}


long clFilter::GetDelay ()
{
    return (lOldSize / 2);
}


void clFilter::Put (const float *fpSrcData, long lSrcCount)
{
    float *fpPrev = Prev;
    float *fpProc = Proc;
    stpSCplx spSCCoeffs = CCoeffs;
    stpSCplx spSCProc = CProc;

    InBuf.Put(fpSrcData, lSrcCount);
    while (InBuf.Get(&fpProc[lOldSize], lNewSize))
    {
        Copy(fpProc, fpPrev, lOldSize);
        Copy(fpPrev, &fpProc[lNewSize], lOldSize);

        FFTi(spSCProc, fpProc);
        Mul(spSCProc, spSCCoeffs, lSpectPoints);
        IFFTo(fpProc, spSCProc);

        OutBuf.Put(&fpProc[lOldSize >> 1], lNewSize);
    }
}


void clFilter::Put (const double *dpSrcData, long lSrcCount)
{
    double *dpPrev = Prev;
    double *dpProc = Proc;
    stpDCplx spDCCoeffs = CCoeffs;
    stpDCplx spDCProc = CProc;

    InBuf.Put(dpSrcData, lSrcCount);
    while (InBuf.Get(&dpProc[lOldSize], lNewSize))
    {
        Copy(dpProc, dpPrev, lOldSize);
        Copy(dpPrev, &dpProc[lNewSize], lOldSize);

        FFTi(spDCProc, dpProc);
        Mul(spDCProc, spDCCoeffs, lSpectPoints);
        IFFTo(dpProc, spDCProc);

        OutBuf.Put(&dpProc[lOldSize >> 1], lNewSize);
    }
}


void clFilter::Put (const float *fpSrcData, long lSrcCount, 
    const stpSCplx spSCCoeffs)
{
    float *fpPrev = Prev;
    float *fpProc = Proc;
    stpSCplx spSCProc = CProc;

    InBuf.Put(fpSrcData, lSrcCount);
    while (InBuf.Get(&fpProc[lOldSize], lNewSize))
    {
        Copy(fpProc, fpPrev, lOldSize);
        Copy(fpPrev, &fpProc[lNewSize], lOldSize);

        FFTi(spSCProc, fpProc);
        Mul(spSCProc, spSCCoeffs, lSpectPoints);
        IFFTo(fpProc, spSCProc);

        OutBuf.Put(&fpProc[lOldSize >> 1], lNewSize);
    }
}


void clFilter::Put (const double *dpSrcData, long lSrcCount, 
    const stpDCplx spDCCoeffs)
{
    double *dpPrev = Prev;
    double *dpProc = Proc;
    stpDCplx spDCProc = CProc;

    InBuf.Put(dpSrcData, lSrcCount);
    while (InBuf.Get(&dpProc[lOldSize], lNewSize))
    {
        Copy(dpProc, dpPrev, lOldSize);
        Copy(dpPrev, &dpProc[lNewSize], lOldSize);

        FFTi(spDCProc, dpProc);
        Mul(spDCProc, spDCCoeffs, lSpectPoints);
        IFFTo(dpProc, spDCProc);

        OutBuf.Put(&dpProc[lOldSize >> 1], lNewSize);
    }
}


bool clFilter::Get (float *fpDestData, long lDestCount)
{
    return OutBuf.Get(fpDestData, lDestCount);
}


bool clFilter::Get (double *dpDestData, long lDestCount)
{
    return OutBuf.Get(dpDestData, lDestCount);
}


void clFilter::DesignLP (float *fpCorner, bool bDCBlock)
{
    long lPassCount;
    float *fpCoeffs;
    clDSPAlloc Coeffs(lSpectPoints * sizeof(float));

    fpCoeffs = Coeffs;
    lPassCount = (long) floor((float) lSpectPoints * (*fpCorner));
    *fpCorner = (float) lPassCount / (float) lSpectPoints;
    Zero(fpCoeffs, lSpectPoints);
    Set(fpCoeffs, 1.0f, lPassCount);
    if (bDCBlock) fpCoeffs[0] = 0.0f;
    SetCoeffs(fpCoeffs);
}


void clFilter::DesignLP (double *dpCorner, bool bDCBlock)
{
    long lPassCount;
    double *dpCoeffs;
    clDSPAlloc Coeffs(lSpectPoints * sizeof(double));

    dpCoeffs = Coeffs;
    lPassCount = (long) floor((double) lSpectPoints * (*dpCorner));
    *dpCorner = (double) lPassCount / (double) lSpectPoints;
    Zero(dpCoeffs, lSpectPoints);
    Set(dpCoeffs, 1.0, lPassCount);
    if (bDCBlock) dpCoeffs[0] = 0.0;
    SetCoeffs(dpCoeffs);
}


void clFilter::DesignLP (float *fpCorner, float fSampleRate, bool bDCBlock)
{
    float fNyquist;
    float fNormCorner;

    fNyquist = fSampleRate / 2.0f;
    fNormCorner = (*fpCorner) / fNyquist;
    DesignLP(&fNormCorner, bDCBlock);
    *fpCorner = fNormCorner * fNyquist;
}


void clFilter::DesignLP (double *dpCorner, double dSampleRate, bool bDCBlock)
{
    double dNyquist;
    double dNormCorner;

    dNyquist = dSampleRate / 2.0;
    dNormCorner = (*dpCorner) / dNyquist;
    DesignLP(&dNormCorner, bDCBlock);
    *dpCorner = dNormCorner * dNyquist;
}


void clFilter::DesignHP (float *fpCorner)
{
    long lStopCount;
    float *fpCoeffs;
    clDSPAlloc Coeffs(lSpectPoints * sizeof(float));

    fpCoeffs = Coeffs;
    lStopCount = (long) ceil((float) lSpectPoints * (*fpCorner));
    *fpCorner = (float) lStopCount / (float) lSpectPoints;
    Set(fpCoeffs, 1.0f, lSpectPoints);
    Zero(fpCoeffs, lStopCount);
    SetCoeffs(fpCoeffs);
}


void clFilter::DesignHP (double *dpCorner)
{
    long lStopCount;
    double *dpCoeffs;
    clDSPAlloc Coeffs(lSpectPoints * sizeof(double));

    dpCoeffs = Coeffs;
    lStopCount = (long) ceil((double) lSpectPoints * (*dpCorner));
    *dpCorner = (double) lStopCount / (double) lSpectPoints;
    Set(dpCoeffs, 1.0, lSpectPoints);
    Zero(dpCoeffs, lStopCount);
    SetCoeffs(dpCoeffs);
}


void clFilter::DesignHP (float *fpCorner, float fSampleRate)
{
    float fNyquist;
    float fNormCorner;

    fNyquist = fSampleRate / 2.0f;
    fNormCorner = (*fpCorner) / fNyquist;
    DesignHP(&fNormCorner);
    *fpCorner = fNormCorner * fNyquist;
}


void clFilter::DesignHP (double *dpCorner, double dSampleRate)
{
    double dNyquist;
    double dNormCorner;

    dNyquist = dSampleRate / 2.0;
    dNormCorner = (*dpCorner) / dNyquist;
    DesignHP(&dNormCorner);
    *dpCorner = dNormCorner * dNyquist;
}


void clFilter::DesignBP (float *fpLowCorner, float *fpHighCorner)
{
    long lLowIdx;
    long lHighIdx;
    float *fpCoeffs;
    clDSPAlloc Coeffs(lSpectPoints * sizeof(float));

    fpCoeffs = Coeffs;
    lLowIdx = (long) ceil((float) lSpectPoints * (*fpLowCorner));
    lHighIdx = (long) floor((float) lSpectPoints * (*fpHighCorner));
    *fpLowCorner = (float) lLowIdx / (float) lSpectPoints;
    *fpHighCorner = (float) lHighIdx / (float) lSpectPoints;
    Zero(fpCoeffs, lSpectPoints);
    Set(&fpCoeffs[lLowIdx], 1.0f, lHighIdx - lLowIdx);
    SetCoeffs(fpCoeffs);
}


void clFilter::DesignBP (double *dpLowCorner, double *dpHighCorner)
{
    long lLowIdx;
    long lHighIdx;
    double *dpCoeffs;
    clDSPAlloc Coeffs(lSpectPoints * sizeof(double));

    dpCoeffs = Coeffs;
    lLowIdx = (long) ceil((double) lSpectPoints * (*dpLowCorner));
    lHighIdx = (long) floor((double) lSpectPoints * (*dpHighCorner));
    *dpLowCorner = (double) lLowIdx / (double) lSpectPoints;
    *dpHighCorner = (double) lHighIdx / (double) lSpectPoints;
    Zero(dpCoeffs, lSpectPoints);
    Set(&dpCoeffs[lLowIdx], 1.0, lHighIdx - lLowIdx);
    SetCoeffs(dpCoeffs);
}


void clFilter::DesignBP (float *fpLowCorner, float *fpHighCorner, 
    float fSampleRate)
{
    float fNyquist;
    float fNormLowCorner;
    float fNormHighCorner;

    fNyquist = fSampleRate / 2.0f;
    fNormLowCorner = (*fpLowCorner) / fNyquist;
    fNormHighCorner = (*fpHighCorner) / fNyquist;
    DesignBP(&fNormLowCorner, &fNormHighCorner);
    *fpLowCorner = fNormLowCorner * fNyquist;
    *fpHighCorner = fNormHighCorner * fNyquist;
}


void clFilter::DesignBP (double *dpLowCorner, double *dpHighCorner,
    double dSampleRate)
{
    double dNyquist;
    double dNormLowCorner;
    double dNormHighCorner;

    dNyquist = dSampleRate / 2.0;
    dNormLowCorner = (*dpLowCorner) / dNyquist;
    dNormHighCorner = (*dpHighCorner) / dNyquist;
    DesignBP(&dNormLowCorner, &dNormHighCorner);
    *dpLowCorner = dNormLowCorner * dNyquist;
    *dpHighCorner = dNormHighCorner * dNyquist;
}


void clFilter::DesignBR (float *fpLowCorner, float *fpHighCorner)
{
    long lLowIdx;
    long lHighIdx;
    float *fpCoeffs;
    clDSPAlloc Coeffs(lSpectPoints * sizeof(float));

    fpCoeffs = Coeffs;
    lLowIdx = (long) ceil((float) lSpectPoints * (*fpLowCorner));
    lHighIdx = (long) floor((float) lSpectPoints * (*fpHighCorner));
    *fpLowCorner = (float) lLowIdx / (float) lSpectPoints;
    *fpHighCorner = (float) lHighIdx / (float) lSpectPoints;
    Set(fpCoeffs, 1.0f, lSpectPoints);
    Zero(&fpCoeffs[lLowIdx], lHighIdx - lLowIdx);
    SetCoeffs(fpCoeffs);
}


void clFilter::DesignBR (double *dpLowCorner, double *dpHighCorner)
{
    long lLowIdx;
    long lHighIdx;
    double *dpCoeffs;
    clDSPAlloc Coeffs(lSpectPoints * sizeof(double));

    dpCoeffs = Coeffs;
    lLowIdx = (long) ceil((double) lSpectPoints * (*dpLowCorner));
    lHighIdx = (long) floor((double) lSpectPoints * (*dpHighCorner));
    *dpLowCorner = (double) lLowIdx / (double) lSpectPoints;
    *dpHighCorner = (double) lHighIdx / (double) lSpectPoints;
    Set(dpCoeffs, 1.0, lSpectPoints);
    Zero(&dpCoeffs[lLowIdx], lHighIdx - lLowIdx);
    SetCoeffs(dpCoeffs);
}


void clFilter::DesignBR (float *fpLowCorner, float *fpHighCorner,
    float fSampleRate)
{
    float fNyquist;
    float fNormLowCorner;
    float fNormHighCorner;

    fNyquist = fSampleRate / 2.0f;
    fNormLowCorner = (*fpLowCorner) / fNyquist;
    fNormHighCorner = (*fpHighCorner) / fNyquist;
    DesignBR(&fNormLowCorner, &fNormHighCorner);
    *fpLowCorner = fNormLowCorner * fNyquist;
    *fpHighCorner = fNormHighCorner * fNyquist;
}


void clFilter::DesignBR (double *dpLowCorner, double *dpHighCorner,
    double dSampleRate)
{
    double dNyquist;
    double dNormLowCorner;
    double dNormHighCorner;

    dNyquist = dSampleRate / 2.0;
    dNormLowCorner = (*dpLowCorner) / dNyquist;
    dNormHighCorner = (*dpHighCorner) / dNyquist;
    DesignBR(&dNormLowCorner, &dNormHighCorner);
    *dpLowCorner = dNormLowCorner * dNyquist;
    *dpHighCorner = dNormHighCorner * dNyquist;
}
