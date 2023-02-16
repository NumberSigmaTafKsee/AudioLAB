/*

    FIR Decimation class
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

#include "dsp/FIRDecimator.hh"


clFIRDecimator::clFIRDecimator ()
{
}


clFIRDecimator::~clFIRDecimator ()
{
    Uninitialize();
}


void clFIRDecimator::Uninitialize ()
{
    DecBuf.Free();
    clFIRMultiRate::Uninitialize();
}


void clFIRDecimator::Put (const float *fpSrcData, long lSrcCount)
{
    InBuf.Put(fpSrcData, lSrcCount);
}


void clFIRDecimator::Put (const double *dpSrcData, long lSrcCount)
{
    InBuf.Put(dpSrcData, lSrcCount);
}


bool clFIRDecimator::Get (float *fpDestData, long lDestCount)
{
    long lWorkCount;
    float *fpDecBuf;

    lWorkCount = lDestCount * lFactor;
    fpDecBuf = (float *) DecBuf.Size(lWorkCount * sizeof(float));
    if (!InBuf.Get(fpDecBuf, lWorkCount))
        return false;
    FIR.FIRFilter(fpDecBuf, lWorkCount);
    DSP.Decimate(fpDestData, fpDecBuf, lFactor, lWorkCount);
    DSP.Mul(fpDestData, fGain, lDestCount);
    return true;
}


bool clFIRDecimator::Get (double *dpDestData, long lDestCount)
{
    long lWorkCount;
    double *dpDecBuf;

    lWorkCount = lDestCount * lFactor;
    dpDecBuf = (double *) DecBuf.Size(lWorkCount * sizeof(double));
    if (!InBuf.Get(dpDecBuf, lWorkCount))
        return false;
    FIR.FIRFilter(dpDecBuf, lWorkCount);
    DSP.Decimate(dpDestData, dpDecBuf, lFactor, lWorkCount);
    DSP.Mul(dpDestData, dGain, lDestCount);
    return true;
}
