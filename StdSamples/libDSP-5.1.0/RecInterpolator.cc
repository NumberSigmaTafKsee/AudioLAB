/*

    Recursive interpolation class
    Copyright (C) 2000-2003 Jussi Laako

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
#ifdef USE_INTEL_MATH
    #include <mathimf.h>
#else
    #include <math.h>
#endif
#include <float.h>

#include "dsp/RecInterpolator.hh"


void clRecInterpolator::InitHalves (double dCenter)
{
    bool bInverted = false;
    double dBandCenter = 0.5;
    double dBandWidth = 0.5;
    long lInitCntr;

    for (lInitCntr = 0; lInitCntr < lSubRounds; lInitCntr++)
    {
        dBandWidth *= 0.5;
        if (dCenter <= dBandCenter)
        {
            bpHalves[lInitCntr] = (!bInverted) ? false : true;
            dBandCenter = dBandCenter - dBandWidth;
        }
        else
        {
            bpHalves[lInitCntr] = (!bInverted) ? true : false;
            dBandCenter = dBandCenter + dBandWidth;
        }
        if (bpHalves[lInitCntr]) bInverted = !bInverted;
    }
}


clRecInterpolator::clRecInterpolator ()
{
    bInitialized = false;
    iType = FILTER_TYPE_FFT;
    lSubRounds = 0;
}


clRecInterpolator::~clRecInterpolator ()
{
    if (bInitialized) Uninitialize();
}


bool clRecInterpolator::Initialize (long lIntFact, long lFiltSize, 
    const float *fpNullPtr, float fBandCenter, int iFilterType)
{
    long lInitCntr;

    if (bInitialized) Uninitialize();
    iType = iFilterType;
    lFilterSize = abs((int) lFiltSize);
    lFactor = lIntFact;
    lIntSize = (long) (lFilterSize * (1.0f - DSP_FILT_DEF_OVERLAPF));
    lSubRounds = (long) (log(lFactor) / log(2) + 0.5);
    if (lSubRounds > RECINT_MAX_SUB_ROUNDS) return false;
    InitHalves(fBandCenter);
    for (lInitCntr = 0; lInitCntr < lSubRounds; lInitCntr++)
    {
        switch (iType)
        {
            case FILTER_TYPE_FFT:
                FFTIntBank[lInitCntr].Initialize(2, lFiltSize, fpNullPtr,
                    bpHalves[lInitCntr]);
                break;
            case FILTER_TYPE_FIR:
                FIRIntBank[lInitCntr].Initialize(2, fpNullPtr, 
                    bpHalves[lInitCntr]);
                break;
            case FILTER_TYPE_IIR:
                IIRIntBank[lInitCntr].Initialize(2, fpNullPtr,
                    bpHalves[lInitCntr]);
                break;
        }
    }
    IntBuf.Size(lIntSize * sizeof(float));
    bInitialized = true;
    return true;
}


bool clRecInterpolator::Initialize (long lIntFact, long lFiltSize,
    const double *dpNullPtr, double dBandCenter, int iFilterType)
{
    long lInitCntr;

    if (bInitialized) Uninitialize();
    iType = iFilterType;
    lFilterSize = lFiltSize;
    lFactor = lIntFact;
    lIntSize = (long) (lFilterSize * (1.0 - DSP_FILT_DEF_OVERLAP));
    lSubRounds = (long) (log(lFactor) / log(2) + 0.5);
    if (lSubRounds > RECINT_MAX_SUB_ROUNDS) return false;
    InitHalves(dBandCenter);
    for (lInitCntr = 0; lInitCntr < lSubRounds; lInitCntr++)
    {
        switch (iType)
        {
            case FILTER_TYPE_FFT:
                FFTIntBank[lInitCntr].Initialize(2, lFiltSize, dpNullPtr,
                    bpHalves[lInitCntr]);
                break;
            case FILTER_TYPE_FIR:
                FIRIntBank[lInitCntr].Initialize(2, dpNullPtr, 
                    bpHalves[lInitCntr]);
                break;
            case FILTER_TYPE_IIR:
                IIRIntBank[lInitCntr].Initialize(2, dpNullPtr,
                    bpHalves[lInitCntr]);
                break;
        }
    }
    IntBuf.Size(lIntSize * sizeof(double));
    bInitialized = true;
    return true;
}


bool clRecInterpolator::Initialize (long lIntFact, long lFiltSize, 
    const float *fpNullPtr, float fBandCenter, bool bUseFIR)
{
    if (!bUseFIR)
        return Initialize(lIntFact, lFiltSize, fpNullPtr, fBandCenter,
            FILTER_TYPE_FFT);
    else
        return Initialize(lIntFact, lFiltSize, fpNullPtr, fBandCenter,
            FILTER_TYPE_FIR);
}


bool clRecInterpolator::Initialize (long lIntFact, long lFiltSize, 
    const double *dpNullPtr, double dBandCenter, bool bUseFIR)
{
    if (!bUseFIR)
        return Initialize(lIntFact, lFiltSize, dpNullPtr, dBandCenter,
            FILTER_TYPE_FFT);
    else
        return Initialize(lIntFact, lFiltSize, dpNullPtr, dBandCenter,
            FILTER_TYPE_FIR);
}


void clRecInterpolator::Uninitialize ()
{
    long lUninitCntr;

    IntBuf.Free();
    for (lUninitCntr = 0; lUninitCntr < lSubRounds; lUninitCntr++)
    {
        switch (iType)
        {
            case FILTER_TYPE_FFT:
                FFTIntBank[lUninitCntr].Uninitialize();
                break;
            case FILTER_TYPE_FIR:
                FIRIntBank[lUninitCntr].Uninitialize();
                break;
            case FILTER_TYPE_IIR:
                IIRIntBank[lUninitCntr].Uninitialize();
                break;
        }
    }
    bInitialized = false;
}


void clRecInterpolator::Put (const float *fpSrcData, long lSrcCount)
{
    long lRecCntr;
    long lSubCount;
    float *fpIntBuf;

    fpIntBuf = IntBuf;
    switch (iType)
    {
        case FILTER_TYPE_FFT:
            FFTIntBank[0].Put(fpSrcData, lSrcCount);
            for (lRecCntr = 0; lRecCntr < (lSubRounds - 1); lRecCntr++)
            {
                while (FFTIntBank[lRecCntr].Get(fpIntBuf, lIntSize))
                    FFTIntBank[lRecCntr + 1].Put(fpIntBuf, lIntSize);
            }
            break;
        case FILTER_TYPE_FIR:
            FIRIntBank[0].Put(fpSrcData, lSrcCount);
            for (lRecCntr = 0; lRecCntr < (lSubRounds - 1); lRecCntr++)
            {
                lSubCount = (lIntSize << lRecCntr);
                while (FIRIntBank[lRecCntr].Get(fpIntBuf, lSubCount))
                    FIRIntBank[lRecCntr + 1].Put(fpIntBuf, lSubCount);
            }
            break;
        case FILTER_TYPE_IIR:
            IIRIntBank[0].Put(fpSrcData, lSrcCount);
            for (lRecCntr = 0; lRecCntr < (lSubRounds - 1); lRecCntr++)
            {
                lSubCount = (lIntSize << lRecCntr);
                while (IIRIntBank[lRecCntr].Get(fpIntBuf, lSubCount))
                    IIRIntBank[lRecCntr + 1].Put(fpIntBuf, lSubCount);
            }
            break;
    }
}


void clRecInterpolator::Put (const double *dpSrcData, long lSrcCount)
{
    long lRecCntr;
    long lSubCount;
    double *dpIntBuf;

    dpIntBuf = IntBuf;
    switch (iType)
    {
        case FILTER_TYPE_FFT:
            FFTIntBank[0].Put(dpSrcData, lSrcCount);
            for (lRecCntr = 0; lRecCntr < (lSubRounds - 1); lRecCntr++)
            {
                while (FFTIntBank[lRecCntr].Get(dpIntBuf, lIntSize))
                    FFTIntBank[lRecCntr + 1].Put(dpIntBuf, lIntSize);
            }
            break;
        case FILTER_TYPE_FIR:
            FIRIntBank[0].Put(dpSrcData, lSrcCount);
            for (lRecCntr = 0; lRecCntr < (lSubRounds - 1); lRecCntr++)
            {
                lSubCount = (lIntSize << lRecCntr);
                while (FIRIntBank[lRecCntr].Get(dpIntBuf, lSubCount))
                    FIRIntBank[lRecCntr + 1].Put(dpIntBuf, lSubCount);
            }
            break;
        case FILTER_TYPE_IIR:
            IIRIntBank[0].Put(dpSrcData, lSrcCount);
            for (lRecCntr = 0; lRecCntr < (lSubRounds - 1); lRecCntr++)
            {
                lSubCount = (lIntSize << lRecCntr);
                while (IIRIntBank[lRecCntr].Get(dpIntBuf, lSubCount))
                    IIRIntBank[lRecCntr + 1].Put(dpIntBuf, lSubCount);
            }
            break;
    }
}


bool clRecInterpolator::Get (float *fpDestData, long lDestCount)
{
    switch (iType)
    {
        case FILTER_TYPE_FFT:
            return FFTIntBank[lSubRounds - 1].Get(fpDestData, lDestCount);
        case FILTER_TYPE_FIR:
            return FIRIntBank[lSubRounds - 1].Get(fpDestData, lDestCount);
        case FILTER_TYPE_IIR:
            return IIRIntBank[lSubRounds - 1].Get(fpDestData, lDestCount);
    }
    return false;
}


bool clRecInterpolator::Get (double *dpDestData, long lDestCount)
{
    switch (iType)
    {
        case FILTER_TYPE_FFT:
            return FFTIntBank[lSubRounds - 1].Get(dpDestData, lDestCount);
        case FILTER_TYPE_FIR:
            return FIRIntBank[lSubRounds - 1].Get(dpDestData, lDestCount);
        case FILTER_TYPE_IIR:
            return IIRIntBank[lSubRounds - 1].Get(dpDestData, lDestCount);
    }
    return false;
}
