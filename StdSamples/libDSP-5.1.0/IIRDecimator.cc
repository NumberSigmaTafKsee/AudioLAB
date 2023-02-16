/*

    IIR Decimation class
    Copyright (C) 2002-2003 Jussi Laako

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

#include "dsp/IIRDecimator.hh"


clIIRDecimator::clIIRDecimator ()
{
}


clIIRDecimator::~clIIRDecimator ()
{
    Uninitialize();
}


void clIIRDecimator::Uninitialize ()
{
    DecBuf.Free();
    clIIRMultiRate::Uninitialize();
}


void clIIRDecimator::Put (const float *fpSrcData, long lSrcCount)
{
    InBuf.Put(fpSrcData, lSrcCount);
}


void clIIRDecimator::Put (const double *dpSrcData, long lSrcCount)
{
    InBuf.Put(dpSrcData, lSrcCount);
}


bool clIIRDecimator::Get (float *fpDestData, long lDestCount)
{
    long lWorkCount;
    float *fpDecBuf;

    lWorkCount = lDestCount * lFactor;
    fpDecBuf = (float *) DecBuf.Size(lWorkCount * sizeof(float));
    if (!InBuf.Get(fpDecBuf, lWorkCount))
        return false;
    Process(fpDecBuf, lWorkCount);
    DSP.Decimate(fpDestData, fpDecBuf, lFactor, lWorkCount);
    return true;
}


bool clIIRDecimator::Get (double *dpDestData, long lDestCount)
{
    long lWorkCount;
    double *dpDecBuf;

    lWorkCount = lDestCount * lFactor;
    dpDecBuf = (double *) DecBuf.Size(lWorkCount * sizeof(double));
    if (!InBuf.Get(dpDecBuf, lWorkCount))
        return false;
    Process(dpDecBuf, lWorkCount);
    DSP.Decimate(dpDestData, dpDecBuf, lFactor, lWorkCount);
    return true;
}
