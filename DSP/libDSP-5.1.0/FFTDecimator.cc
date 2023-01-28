/*

    FFT Decimation class
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

#include "dsp/FFTDecimator.hh"


clFFTDecimator::clFFTDecimator ()
{
}


clFFTDecimator::~clFFTDecimator ()
{
    if (bInitialized) Uninitialize();
}


void clFFTDecimator::Uninitialize ()
{
    DecBuf.Free();
    clFFTMultiRate::Uninitialize();
}


void clFFTDecimator::Put (const float *fpSrcData, long lSrcCount)
{
    Filter.Put(fpSrcData, lSrcCount);
}


void clFFTDecimator::Put (const double *dpSrcData, long lSrcCount)
{
    Filter.Put(dpSrcData, lSrcCount);
}


bool clFFTDecimator::Get (float *fpDestData, long lDestCount)
{
    long lWorkCount;
    float *fpDecBuf;

    lWorkCount = lDestCount * lFactor;
    fpDecBuf = (float *) DecBuf.Size(lWorkCount * sizeof(float));
    if (!Filter.Get(fpDecBuf, lWorkCount)) return false;
    DSP.Decimate(fpDestData, fpDecBuf, lFactor, lWorkCount);
    return true;
}


bool clFFTDecimator::Get (double *dpDestData, long lDestCount)
{
    long lWorkCount;
    double *dpDecBuf;

    lWorkCount = lDestCount * lFactor;
    dpDecBuf = (double *) DecBuf.Size(lWorkCount * sizeof(float));
    if (!Filter.Get(dpDecBuf, lWorkCount)) return false;
    DSP.Decimate(dpDestData, dpDecBuf, lFactor, lWorkCount);
    return true;
}
