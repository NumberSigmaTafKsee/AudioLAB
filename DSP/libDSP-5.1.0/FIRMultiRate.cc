/*

    Base class for FIR based multirate filters
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

#include "dsp/FIRMultiRate.hh"

/*#include "Dec2Filter2.h"
#include "Dec3Filter2.h"
#include "Dec4Filter2.h"
#include "Dec8Filter2.h"
#include "Dec2hpFilter2.h"
#include "Dec3hpFilter2.h"
#include "Dec4hpFilter2.h"
#include "Dec8hpFilter2.h"*/
#include "Dec2Filter3.h"
#include "Dec3Filter3.h"
#include "Dec4Filter3.h"
#include "Dec8Filter3.h"
#include "Dec2hpFilter3.h"
#include "Dec3hpFilter3.h"
#include "Dec4hpFilter3.h"
#include "Dec8hpFilter3.h"


clFIRMultiRate::clFIRMultiRate ()
{
}


clFIRMultiRate::~clFIRMultiRate ()
{
    Uninitialize();
}


bool clFIRMultiRate::Initialize (long lFactorP, const float *fpNullPtr,
    bool bHighPass)
{
    lFactor = lFactorP;
    switch (lFactor)
    {
        case 2:
            if (!bHighPass)
            {
                fGain = fDec2FilterGain;
                FIR.FIRAllocate(fpDec2FilterCoeffs, lDec2FilterLen);
            }
            else
            {
                fGain = fDec2hpFilterGain;
                FIR.FIRAllocate(fpDec2hpFilterCoeffs, lDec2hpFilterLen);
            }
            break;
        case 3:
            if (!bHighPass)
            {
                fGain = fDec3FilterGain;
                FIR.FIRAllocate(fpDec3FilterCoeffs, lDec3FilterLen);
            }
            else
            {
                fGain = fDec3hpFilterGain;
                FIR.FIRAllocate(fpDec3hpFilterCoeffs, lDec3hpFilterLen);
            }
            break;
        case 4:
            if (!bHighPass)
            {
                fGain = fDec4FilterGain;
                FIR.FIRAllocate(fpDec4FilterCoeffs, lDec4FilterLen);
            }
            else
            {
                fGain = fDec4hpFilterGain;
                FIR.FIRAllocate(fpDec4hpFilterCoeffs, lDec4hpFilterLen);
            }
            break;
        case 8:
            if (!bHighPass)
            {
                fGain = fDec8FilterGain;
                FIR.FIRAllocate(fpDec8FilterCoeffs, lDec8FilterLen);
            }
            else
            {
                fGain = fDec8hpFilterGain;
                FIR.FIRAllocate(fpDec8hpFilterCoeffs, lDec8hpFilterLen);
            }
            break;
        default:
            return false;
    }
    return true;
}


bool clFIRMultiRate::Initialize (long lFactorP, const double *dpNullPtr,
    bool bHighPass)
{
    lFactor = lFactorP;
    switch (lFactor)
    {
        case 2:
            if (!bHighPass)
            {
                dGain = dDec2FilterGain;
                FIR.FIRAllocate(dpDec2FilterCoeffs, lDec2FilterLen);
            }
            else
            {
                dGain = dDec2hpFilterGain;
                FIR.FIRAllocate(dpDec2hpFilterCoeffs, lDec2hpFilterLen);
            }
            break;
        case 3:
            if (!bHighPass)
            {
                dGain = dDec3FilterGain;
                FIR.FIRAllocate(dpDec3FilterCoeffs, lDec3FilterLen);
            }
            else
            {
                dGain = dDec3hpFilterGain;
                FIR.FIRAllocate(dpDec3hpFilterCoeffs, lDec3hpFilterLen);
            }
            break;
        case 4:
            if (!bHighPass)
            {
                dGain = dDec4FilterGain;
                FIR.FIRAllocate(dpDec4FilterCoeffs, lDec4FilterLen);
            }
            else
            {
                dGain = dDec4hpFilterGain;
                FIR.FIRAllocate(dpDec4hpFilterCoeffs, lDec4hpFilterLen);
            }
            break;
        case 8:
            if (!bHighPass)
            {
                dGain = dDec8FilterGain;
                FIR.FIRAllocate(dpDec8FilterCoeffs, lDec8FilterLen);
            }
            else
            {
                dGain = dDec8hpFilterGain;
                FIR.FIRAllocate(dpDec8hpFilterCoeffs, lDec8hpFilterLen);
            }
            break;
        default:
            return false;
    }
    return true;
}


void clFIRMultiRate::Uninitialize ()
{
    FIR.FIRFree();
}
