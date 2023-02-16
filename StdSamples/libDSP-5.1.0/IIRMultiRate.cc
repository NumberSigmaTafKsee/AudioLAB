/*

    Base class for IIR based multirate filters
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

#include "dsp/IIRMultiRate.hh"

/*#include "Dec2IIR.h"
#include "Dec2hpIIR.h"
#include "Dec3IIR.h"
#include "Dec3hpIIR.h"*/
/*#include "Dec2IIR2.h"
#include "Dec2hpIIR2.h"
#include "Dec3IIR2.h"
#include "Dec3hpIIR2.h"*/
#include "Dec2IIR3.h"
#include "Dec2hpIIR3.h"
#include "Dec3IIR3.h"
#include "Dec3hpIIR3.h"


clIIRMultiRate::clIIRMultiRate ()
{
    lFactor = 1;
}


clIIRMultiRate::~clIIRMultiRate ()
{
    Uninitialize();
}


bool clIIRMultiRate::Initialize (long lFactorP, const float *fpNullPtr,
    bool bHighPass)
{
    lFactor = lFactorP;
    switch (lFactor)
    {
        case 2:
            if (!bHighPass)
                clIIRCascade::Initialize(fpDec2IIRCoeffs, lDec2IIRSize);
            else
                clIIRCascade::Initialize(fpDec2hpIIRCoeffs, lDec2hpIIRSize);
            break;
        case 3:
            if (!bHighPass)
                clIIRCascade::Initialize(fpDec3IIRCoeffs, lDec3IIRSize);
            else
                clIIRCascade::Initialize(fpDec3hpIIRCoeffs, lDec3hpIIRSize);
            break;
        default:
            return false;
    };
    return true;
}


bool clIIRMultiRate::Initialize (long lFactorP, const double *dpNullPtr,
    bool bHighPass)
{
    lFactor = lFactorP;
    switch (lFactor)
    {
        case 2:
            if (!bHighPass)
                clIIRCascade::Initialize(dpDec2IIRCoeffs, lDec2IIRSize);
            else
                clIIRCascade::Initialize(dpDec2hpIIRCoeffs, lDec2hpIIRSize);
            break;
        case 3:
            if (!bHighPass)
                clIIRCascade::Initialize(dpDec3IIRCoeffs, lDec3IIRSize);
            else
                clIIRCascade::Initialize(dpDec3hpIIRCoeffs, lDec3hpIIRSize);
            break;
        default:
            return false;
    };
    return true;
}


void clIIRMultiRate::Uninitialize ()
{
    clIIRCascade::Uninitialize();
}
