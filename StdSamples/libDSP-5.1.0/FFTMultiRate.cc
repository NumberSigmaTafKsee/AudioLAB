/*

    Base class for FFT based multirate filters
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

#include "dsp/FFTMultiRate.hh"


clFFTMultiRate::clFFTMultiRate ()
{
    bInitialized = false;
}


clFFTMultiRate::~clFFTMultiRate ()
{
    if (bInitialized) Uninitialize();
}


bool clFFTMultiRate::Initialize (long lFactorP, long lFiltSize, 
    const float *fpNullPtr, bool bHighPass)
{
    float fCorner;
    float fPassBand;
    float fStopBand;

    if (bInitialized) Uninitialize();
    lFilterSize = abs((int) lFiltSize);
    lFactor = lFactorP;
    fCorner = 1.0f / (float) lFactor;
    if (lFiltSize <= 0)
    {
        if (!bHighPass)
        {
            fPassBand = fCorner * (1.0f - FFTMULTIRATE_DELTAOMEGA / 2.0f);
            fStopBand = fCorner * (1.0f + FFTMULTIRATE_DELTAOMEGA / 2.0f);
            Filter.InitializeLP(fPassBand, fStopBand, 
                FFTMULTIRATE_RIPPLERATIO,
                FFTMULTIRATE_OVERLAP);
        }
        else
        {
            fPassBand = fCorner * (1.0f + FFTMULTIRATE_DELTAOMEGA / 2.0f);
            fStopBand = fCorner * (1.0f - FFTMULTIRATE_DELTAOMEGA / 2.0f);
            Filter.InitializeHP(fPassBand, fStopBand, 
                FFTMULTIRATE_RIPPLERATIO,
                FFTMULTIRATE_OVERLAP);
        }
    }
    else
    {
        Filter.Initialize(lFilterSize, fpNullPtr);
        if (!bHighPass)
            Filter.DesignLP(&fCorner);
        else
            Filter.DesignHP(&fCorner);
    }
    return true;
}


bool clFFTMultiRate::Initialize (long lFactorP, long lFiltSize,
    const double *dpNullPtr, bool bHighPass)
{
    double dCorner;
    double dPassBand;
    double dStopBand;

    if (bInitialized) Uninitialize();
    lFilterSize = abs((int) lFiltSize);
    lFactor = lFactorP;
    dCorner = 1.0 / (double) lFactor;
    if (lFiltSize <= 0)
    {
        if (!bHighPass)
        {
            dPassBand = dCorner * (1.0 - FFTMULTIRATE_DELTAOMEGA / 2.0);
            dStopBand = dCorner * (1.0 + FFTMULTIRATE_DELTAOMEGA / 2.0);
            Filter.InitializeLP(dPassBand, dStopBand, 
                (double) FFTMULTIRATE_RIPPLERATIO,
                (double) FFTMULTIRATE_OVERLAP);
        }
        else
        {
            dPassBand = dCorner * (1.0 + FFTMULTIRATE_DELTAOMEGA / 2.0);
            dStopBand = dCorner * (1.0 - FFTMULTIRATE_DELTAOMEGA / 2.0);
            Filter.InitializeHP(dPassBand, dStopBand, 
                (double) FFTMULTIRATE_RIPPLERATIO,
                (double) FFTMULTIRATE_OVERLAP);
        }
    }
    else
    {
        Filter.Initialize(lFilterSize, dpNullPtr);
        if (!bHighPass)
            Filter.DesignLP(&dCorner);
        else
            Filter.DesignHP(&dCorner);
    }
    return true;
}


void clFFTMultiRate::Uninitialize ()
{
    Filter.Uninitialize();
    bInitialized = false;
}
