/*

    Recursive decimation class
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

#include "dsp/RecDecimator.hh"


void clRecDecimator::InitHalves (double dCenter)
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


clRecDecimator::clRecDecimator ()
{
    bInitialized = false;
    iType = FILTER_TYPE_FFT;
    lSubRounds = 0;
}


clRecDecimator::~clRecDecimator ()
{
    if (bInitialized) Uninitialize();
}


bool clRecDecimator::Initialize (long lDecFact, long lFiltSize, 
    const float *fpNullPtr, float fBandCenter, int iFilterType)
{
    long lInitCntr;

    if (bInitialized) Uninitialize();
    iType = iFilterType;
    lFilterSize = abs((int) lFiltSize);
    lFactor = lDecFact;
    lDecSize = (long) (lFilterSize * (1.0f - DSP_FILT_DEF_OVERLAPF));
    lSubRounds = (long) (log(lFactor) / log(2) + 0.5);
    if (lSubRounds > RECDEC_MAX_SUB_ROUNDS) return false;
    InitHalves(fBandCenter);
    for (lInitCntr = 0; lInitCntr < lSubRounds; lInitCntr++)
    {
        switch (iType)
        {
            case FILTER_TYPE_FFT:
                FFTDecBank[lInitCntr].Initialize(2, lFiltSize, fpNullPtr,
                    bpHalves[lInitCntr]);
                break;
            case FILTER_TYPE_FIR:
                FIRDecBank[lInitCntr].Initialize(2, fpNullPtr, 
                    bpHalves[lInitCntr]);
                break;
            case FILTER_TYPE_IIR:
                IIRDecBank[lInitCntr].Initialize(2, fpNullPtr,
                    bpHalves[lInitCntr]);
                break;
        }
    }
    DecBuf.Size(lDecSize * sizeof(float));
    bInitialized = true;
    return true;
}


bool clRecDecimator::Initialize (long lDecFact, long lFiltSize,
    const double *dpNullPtr, double dBandCenter, int iFilterType)
{
    long lInitCntr;

    if (bInitialized) Uninitialize();
    iType = iFilterType;
    lFilterSize = lFiltSize;
    lFactor = lDecFact;
    lDecSize = (long) (lFilterSize * (1.0 - DSP_FILT_DEF_OVERLAP));
    lSubRounds = (long) (log(lFactor) / log(2) + 0.5);
    if (lSubRounds > RECDEC_MAX_SUB_ROUNDS) return false;
    InitHalves(dBandCenter);
    for (lInitCntr = 0; lInitCntr < lSubRounds; lInitCntr++)
    {
        switch (iType)
        {
            case FILTER_TYPE_FFT:
                FFTDecBank[lInitCntr].Initialize(2, lFiltSize, dpNullPtr,
                    bpHalves[lInitCntr]);
                break;
            case FILTER_TYPE_FIR:
                FIRDecBank[lInitCntr].Initialize(2, dpNullPtr, 
                    bpHalves[lInitCntr]);
                break;
            case FILTER_TYPE_IIR:
                IIRDecBank[lInitCntr].Initialize(2, dpNullPtr,
                    bpHalves[lInitCntr]);
                break;
        }
    }
    DecBuf.Size(lDecSize * sizeof(double));
    bInitialized = true;
    return true;
}


bool clRecDecimator::Initialize (long lDecFact, long lFiltSize, 
    const float *fpNullPtr, float fBandCenter, bool bUseFIR)
{
    if (!bUseFIR)
        return Initialize(lDecFact, lFiltSize, fpNullPtr, fBandCenter, 
            FILTER_TYPE_FFT);
    else
        return Initialize(lDecFact, lFiltSize, fpNullPtr, fBandCenter, 
            FILTER_TYPE_FIR);
}


bool clRecDecimator::Initialize (long lDecFact, long lFiltSize,
    const double *dpNullPtr, double dBandCenter, bool bUseFIR)
{
    if (!bUseFIR)
        return Initialize(lDecFact, lFiltSize, dpNullPtr, dBandCenter, 
            FILTER_TYPE_FFT);
    else
        return Initialize(lDecFact, lFiltSize, dpNullPtr, dBandCenter, 
            FILTER_TYPE_FIR);
}


void clRecDecimator::Uninitialize ()
{
    long lUninitCntr;

    DecBuf.Free();
    for (lUninitCntr = 0; lUninitCntr < lSubRounds; lUninitCntr++)
    {
        switch (iType)
        {
            case FILTER_TYPE_FFT:
                FFTDecBank[lUninitCntr].Uninitialize();
                break;
            case FILTER_TYPE_FIR:
                FIRDecBank[lUninitCntr].Uninitialize();
                break;
            case FILTER_TYPE_IIR:
                IIRDecBank[lUninitCntr].Uninitialize();
                break;
        }
    }
    bInitialized = false;
}


void clRecDecimator::Put (const float *fpSrcData, long lSrcCount)
{
    long lRecCntr;
    long lSubCount;
    float *fpDecBuf;

    fpDecBuf = DecBuf;
    switch (iType)
    {
        case FILTER_TYPE_FFT:
            FFTDecBank[0].Put(fpSrcData, lSrcCount);
            for (lRecCntr = 0; lRecCntr < (lSubRounds - 1); lRecCntr++)
            {
                while (FFTDecBank[lRecCntr].Get(fpDecBuf, lDecSize))
                    FFTDecBank[lRecCntr + 1].Put(fpDecBuf, lDecSize);
            }
            break;
        case FILTER_TYPE_FIR:
            FIRDecBank[0].Put(fpSrcData, lSrcCount);
            for (lRecCntr = 0; lRecCntr < (lSubRounds - 1); lRecCntr++)
            {
                lSubCount = (lDecSize >> lRecCntr);
                while (FIRDecBank[lRecCntr].Get(fpDecBuf, lSubCount))
                    FIRDecBank[lRecCntr + 1].Put(fpDecBuf, lSubCount);
            }
            break;
        case FILTER_TYPE_IIR:
            IIRDecBank[0].Put(fpSrcData, lSrcCount);
            for (lRecCntr = 0; lRecCntr < (lSubRounds - 1); lRecCntr++)
            {
                lSubCount = (lDecSize >> lRecCntr);
                while (IIRDecBank[lRecCntr].Get(fpDecBuf, lSubCount))
                    IIRDecBank[lRecCntr + 1].Put(fpDecBuf, lSubCount);
            }
            break;
    }
}


void clRecDecimator::Put (const double *dpSrcData, long lSrcCount)
{
    long lRecCntr;
    long lSubCount;
    double *dpDecBuf;

    dpDecBuf = DecBuf;
    switch (iType)
    {
        case FILTER_TYPE_FFT:
            FFTDecBank[0].Put(dpSrcData, lSrcCount);
            for (lRecCntr = 0; lRecCntr < (lSubRounds - 1); lRecCntr++)
            {
                while (FFTDecBank[lRecCntr].Get(dpDecBuf, lDecSize))
                    FFTDecBank[lRecCntr + 1].Put(dpDecBuf, lDecSize);
            }
            break;
        case FILTER_TYPE_FIR:
            FIRDecBank[0].Put(dpSrcData, lSrcCount);
            for (lRecCntr = 0; lRecCntr < (lSubRounds - 1); lRecCntr++)
            {
                lSubCount = (lDecSize >> lRecCntr);
                while (FIRDecBank[lRecCntr].Get(dpDecBuf, lSubCount))
                    FIRDecBank[lRecCntr + 1].Put(dpDecBuf, lSubCount);
            }
            break;
        case FILTER_TYPE_IIR:
            IIRDecBank[0].Put(dpSrcData, lSrcCount);
            for (lRecCntr = 0; lRecCntr < (lSubRounds - 1); lRecCntr++)
            {
                lSubCount = (lDecSize >> lRecCntr);
                while (IIRDecBank[lRecCntr].Get(dpDecBuf, lSubCount))
                    IIRDecBank[lRecCntr + 1].Put(dpDecBuf, lSubCount);
            }
            break;
    }
}


bool clRecDecimator::Get (float *fpDestData, long lDestCount)
{
    switch (iType)
    {
        case FILTER_TYPE_FFT:
            return FFTDecBank[lSubRounds - 1].Get(fpDestData, lDestCount);
        case FILTER_TYPE_FIR:
            return FIRDecBank[lSubRounds - 1].Get(fpDestData, lDestCount);
        case FILTER_TYPE_IIR:
            return IIRDecBank[lSubRounds - 1].Get(fpDestData, lDestCount);
    }
    return false;
}


bool clRecDecimator::Get (double *dpDestData, long lDestCount)
{
    switch (iType)
    {
        case FILTER_TYPE_FFT:
            return FFTDecBank[lSubRounds - 1].Get(dpDestData, lDestCount);
        case FILTER_TYPE_FIR:
            return FIRDecBank[lSubRounds - 1].Get(dpDestData, lDestCount);
        case FILTER_TYPE_IIR:
            return IIRDecBank[lSubRounds - 1].Get(dpDestData, lDestCount);
    }
    return false;
}
