/*

    DSP operations
    Copyright (C) 1998-2005 Jussi Laako

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
#include <climits>
#ifdef USE_INTEL_MATH
    #include <mathimf.h>
#else
    #include <math.h>
#endif
#include <float.h>
#ifdef DSP_IPP
    #include <ipp.h>
#endif
#ifdef __INTEL_COMPILER
    #include <xmmintrin.h>
    #include <pmmintrin.h>
#endif
#if (defined(DSP_X86) && (!defined(__INTEL_COMPILER)))
    #include <xmmintrin.h>
#endif

#include <commoncpp/Compilers.hh>

#include "dsp/DSPOp.hh"
#ifdef DSP_X86
#include "dsp/X86.h"
#endif
#if (defined(DSP_X86_64) && (!defined(DSP_X86)))
#include "dsp/X86-64.h"
#endif


#ifdef DSP_X86
bool bHave3DNow = false;
int iHaveSSE = 0;
unsigned int uiOldCSR = 0;
#endif


extern "C"
{

INLINE int FloatCompare (const void *vpValue1, const void *vpValue2)
{
    if (*(float *)vpValue1 == *(float *)vpValue2)
    {
        return 0;
    }
    if (*(float *)vpValue1 < *(float *)vpValue2)
    {
        return -1;
    }
    else
    {
        return 1;
    }
}


INLINE int DoubleCompare (const void *vpValue1, const void *vpValue2)
{
    if (*(double *)vpValue1 == *(double *)vpValue2)
    {
        return 0;
    }
    if (*(double *)vpValue1 < *(double *)vpValue2)
    {
        return -1;
    }
    else
    {
        return 1;
    }
}


INLINE int LongCompare (const void *vpValue1, const void *vpValue2)
{
    if (*(long *)vpValue1 == *(long *)vpValue2)
    {
        return 0;
    }
    if (*(long *)vpValue1 < *(long *)vpValue2)
    {
        return -1;
    }
    else
    {
        return 1;
    }
}

}


signed long clDSPOp::Round (float fValue)
{
    signed long slTemp;

    // nearbyintf() seems to be buggy in some glibc versions
    #if defined(USE_INTEL_MATH)
    slTemp = lrintf(fValue);
    #elif defined(_ISOC9X_SOURCE)
    slTemp = (signed long) nearbyintf(fValue);
    #else
    slTemp = (fValue >= 0.0F) ?
        ((signed long) (fValue + 0.5F)) : ((signed long) (fValue - 0.5F));
    #endif
    return slTemp;
}


signed long clDSPOp::Round (double dValue)
{
    signed long slTemp;

    // nearbyint() seems to be buggy in some glibc versions
    #if defined(USE_INTEL_MATH)
    slTemp = lrint(dValue);
    #elif defined(_ISOC9X_SOURCE)
    slTemp = (signed long) nearbyint(dValue);
    #else
    slTemp = (dValue >= 0.0) ?
        ((signed long) (dValue + 0.5)) : ((signed long) (dValue - 0.5));
    #endif
    return slTemp;
}


INLINE void clDSPOp::Cart2Polar (float *fpMagn, float *fpPhase,
    float fReal, float fImag)
{
    #ifndef _ISOC9X_SOURCE
    *fpMagn = sqrtf(fReal * fReal + fImag * fImag);
    #else
    *fpMagn = hypotf(fReal, fImag);
    #endif
    *fpPhase = atan2f(fImag, fReal);
}


INLINE void clDSPOp::Cart2Polar (double *dpMagn, double *dpPhase,
    double dReal, double dImag)
{
    #ifndef _ISOC9X_SOURCE
    *dpMagn = sqrt(dReal * dReal + dImag * dImag);
    #else
    *dpMagn = hypot(dReal, dImag);
    #endif
    *dpPhase = atan2(dImag, dReal);
}


INLINE void clDSPOp::Cart2Polar (float *fpMagn, float *fpPhase,
    const stpSCplx spCplx)
{
    #ifndef _ISOC9X_SOURCE
    *fpMagn = sqrtf(spCplx->R * spCplx->R + spCplx->I * spCplx->I);
    #else
    *fpMagn = hypotf(spCplx->R, spCplx->I);
    #endif
    *fpPhase = atan2f(spCplx->I, spCplx->R);
}


INLINE void clDSPOp::Cart2Polar (double *dpMagn, double *dpPhase,
    const stpDCplx spCplx)
{
    #ifndef _ISOC9X_SOURCE
    *dpMagn = sqrt(spCplx->R * spCplx->R + spCplx->I * spCplx->I);
    #else
    *dpMagn = hypot(spCplx->R, spCplx->I);
    #endif
    *dpPhase = atan2(spCplx->I, spCplx->R);
}


INLINE void clDSPOp::Cart2Polar (stpSPolar spPolar, const stpSCplx spCplx)
{
    #ifndef _ISOC9X_SOURCE
    spPolar->M = sqrtf(spCplx->R * spCplx->R + spCplx->I * spCplx->I);
    #else
    spPolar->M = hypotf(spCplx->R, spCplx->I);
    #endif
    spPolar->P = atan2f(spCplx->I, spCplx->R);
}


INLINE void clDSPOp::Cart2Polar (stpDPolar spPolar, const stpDCplx spCplx)
{
    #ifndef _ISOC9X_SOURCE
    spPolar->M = sqrt(spCplx->R * spCplx->R + spCplx->I * spCplx->I);
    #else
    spPolar->M = hypot(spCplx->R, spCplx->I);
    #endif
    spPolar->P = atan2(spCplx->I, spCplx->R);
}


INLINE void clDSPOp::Cart2Polar (utpSCoord upCoord)
{
    #ifndef _ISOC9X_SOURCE
    upCoord->P.M = sqrtf(upCoord->C.R * upCoord->C.R +
        upCoord->C.I * upCoord->C.I);
    #else
    upCoord->P.M = hypotf(upCoord->C.R, upCoord->C.I);
    #endif
    upCoord->P.P = atan2f(upCoord->C.I, upCoord->C.R);
}


INLINE void clDSPOp::Cart2Polar (utpDCoord upCoord)
{
    #ifndef _ISOC9X_SOURCE
    upCoord->P.M = sqrt(upCoord->C.R * upCoord->C.R +
        upCoord->C.I * upCoord->C.I);
    #else
    upCoord->P.M = hypot(upCoord->C.R, upCoord->C.I);
    #endif
    upCoord->P.P = atan2(upCoord->C.I, upCoord->C.R);
}


INLINE void clDSPOp::Polar2Cart (float *fpReal, float *fpImag,
    float fMagn, float fPhase)
{
    #ifndef _GNU_SOURCE
    *fpReal = fMagn * cosf(fPhase);
    *fpImag = fMagn * sinf(fPhase);
    #else
    sincosf(fPhase, fpImag, fpReal);
    *fpReal *= fMagn;
    *fpImag *= fMagn;
    #endif
}


INLINE void clDSPOp::Polar2Cart (double *dpReal, double *dpImag,
    double dMagn, double dPhase)
{
    #ifndef _GNU_SOURCE
    *dpReal = dMagn * cos(dPhase);
    *dpImag = dMagn * sin(dPhase);
    #else
    sincos(dPhase, dpImag, dpReal);
    *dpReal *= dMagn;
    *dpImag *= dMagn;
    #endif
}


INLINE void clDSPOp::Polar2Cart (stpSCplx spCplx, float fMagn, float fPhase)
{
    #ifndef _GNU_SOURCE
    spCplx->R = fMagn * cosf(fPhase);
    spCplx->I = fMagn * sinf(fPhase);
    #else
    sincosf(fPhase, &spCplx->I, &spCplx->R);
    spCplx->R *= fMagn;
    spCplx->I *= fMagn;
    #endif
}


INLINE void clDSPOp::Polar2Cart (stpDCplx spCplx, double dMagn, double dPhase)
{
    #ifndef _GNU_SOURCE
    spCplx->R = dMagn * cos(dPhase);
    spCplx->I = dMagn * sin(dPhase);
    #else
    sincos(dPhase, &spCplx->I, &spCplx->R);
    spCplx->R *= dMagn;
    spCplx->I *= dMagn;
    #endif
}


INLINE void clDSPOp::Polar2Cart (stpSCplx spCplx, const stpSPolar spPolar)
{
    #ifndef _GNU_SOURCE
    spCplx->R = spPolar->M * cosf(spPolar->P);
    spCplx->I = spPolar->M * sinf(spPolar->P);
    #else
    sincosf(spPolar->P, &spCplx->I, &spCplx->R);
    spCplx->R *= spPolar->M;
    spCplx->I *= spPolar->M;
    #endif
}


INLINE void clDSPOp::Polar2Cart (stpDCplx spCplx, const stpDPolar spPolar)
{
    #ifndef _GNU_SOURCE
    spCplx->R = spPolar->M * cos(spPolar->P);
    spCplx->I = spPolar->M * sin(spPolar->P);
    #else
    sincos(spPolar->P, &spCplx->I, &spCplx->R);
    spCplx->R *= spPolar->M;
    spCplx->I *= spPolar->M;
    #endif
}


INLINE void clDSPOp::Polar2Cart (utpSCoord upCoord)
{
    #ifndef _GNU_SOURCE
    upCoord->C.R = upCoord->P.M * cosf(upCoord->P.P);
    upCoord->C.I = upCoord->P.M * sinf(upCoord->P.P);
    #else
    float fReal;
    float fImag;

    sincosf(upCoord->P.P, &fImag, &fReal);
    upCoord->C.R = upCoord->P.M * fReal;
    upCoord->C.I = upCoord->P.M * fImag;
    #endif
}


INLINE void clDSPOp::Polar2Cart (utpDCoord upCoord)
{
    #ifndef _GNU_SOURCE
    upCoord->C.R = upCoord->P.M * cos(upCoord->P.P);
    upCoord->C.I = upCoord->P.M * sin(upCoord->P.P);
    #else
    double dReal;
    double dImag;

    sincos(upCoord->P.P, &dImag, &dReal);
    upCoord->C.R = upCoord->P.M * dReal;
    upCoord->C.I = upCoord->P.M * dImag;
    #endif
}


INLINE void clDSPOp::CplxAdd (stpSCplx spCplxDest, const stpSCplx spCplxSrc)
{
    spCplxDest->R += spCplxSrc->R;
    spCplxDest->I += spCplxSrc->I;
}


INLINE void clDSPOp::CplxAdd (stpDCplx spCplxDest, const stpDCplx spCplxSrc)
{
    spCplxDest->R += spCplxSrc->R;
    spCplxDest->I += spCplxSrc->I;
}


INLINE void clDSPOp::CplxAdd (stpSCplx spCplxDest, const stpSCplx spCplxSrc1,
    const stpSCplx spCplxSrc2)
{
    spCplxDest->R = spCplxSrc1->R + spCplxSrc2->R;
    spCplxDest->I = spCplxSrc1->I + spCplxSrc2->I;
}


INLINE void clDSPOp::CplxAdd (stpDCplx spCplxDest, const stpDCplx spCplxSrc1,
    const stpDCplx spCplxSrc2)
{
    spCplxDest->R = spCplxSrc1->R + spCplxSrc2->R;
    spCplxDest->I = spCplxSrc1->I + spCplxSrc2->I;
}


INLINE void clDSPOp::CplxSub (stpSCplx spCplxDest, const stpSCplx spCplxSrc)
{
    spCplxDest->R -= spCplxSrc->R;
    spCplxDest->I -= spCplxSrc->I;
}


INLINE void clDSPOp::CplxSub (stpDCplx spCplxDest, const stpDCplx spCplxSrc)
{
    spCplxDest->R -= spCplxSrc->R;
    spCplxDest->I -= spCplxSrc->I;
}


INLINE void clDSPOp::CplxSub (stpSCplx spCplxDest, const stpSCplx spCplxSrc1,
    const stpSCplx spCplxSrc2)
{
    spCplxDest->R = spCplxSrc1->R - spCplxSrc2->R;
    spCplxDest->I = spCplxSrc1->I - spCplxSrc2->I;
}


INLINE void clDSPOp::CplxSub (stpDCplx spCplxDest, const stpDCplx spCplxSrc1,
    const stpDCplx spCplxSrc2)
{
    spCplxDest->R = spCplxSrc1->R - spCplxSrc2->R;
    spCplxDest->I = spCplxSrc1->I - spCplxSrc2->I;
}


INLINE void clDSPOp::CplxMul (stpSCplx spCplxDest, float fSrc)
{
    spCplxDest->R *= fSrc;
    spCplxDest->I *= fSrc;
}


INLINE void clDSPOp::CplxMul (stpDCplx spCplxDest, double dSrc)
{
    spCplxDest->R *= dSrc;
    spCplxDest->I *= dSrc;
}


INLINE void clDSPOp::CplxMul (stpSCplx spCplxDest, const stpSCplx spCplxSrc)
{
    float fReal;
    float fImag;

    fReal = spCplxDest->R * spCplxSrc->R - spCplxDest->I * spCplxSrc->I;
    fImag = spCplxDest->R * spCplxSrc->I + spCplxDest->I * spCplxSrc->R;
    spCplxDest->R = fReal;
    spCplxDest->I = fImag;
}


INLINE void clDSPOp::CplxMul (stpDCplx spCplxDest, const stpDCplx spCplxSrc)
{
    double dReal;
    double dImag;

    dReal = spCplxDest->R * spCplxSrc->R - spCplxDest->I * spCplxSrc->I;
    dImag = spCplxDest->R * spCplxSrc->I + spCplxDest->I * spCplxSrc->R;
    spCplxDest->R = dReal;
    spCplxDest->I = dImag;
}


INLINE void clDSPOp::CplxMul (stpSCplx spCplxDest, const stpSCplx spCplxSrc1,
    const stpSCplx spCplxSrc2)
{
    spCplxDest->R = 
        spCplxSrc1->R * spCplxSrc2->R - spCplxSrc1->I * spCplxSrc2->I;
    spCplxDest->I =
        spCplxSrc1->R * spCplxSrc2->I + spCplxSrc1->I * spCplxSrc2->R;
}


INLINE void clDSPOp::CplxMul (stpDCplx spCplxDest, const stpDCplx spCplxSrc1,
    const stpDCplx spCplxSrc2)
{
    spCplxDest->R =
        spCplxSrc1->R * spCplxSrc2->R - spCplxSrc1->I * spCplxSrc2->I;
    spCplxDest->I =
        spCplxSrc1->R * spCplxSrc2->I + spCplxSrc1->I * spCplxSrc2->R;
}


INLINE void clDSPOp::CplxMulC (stpSCplx spCplxDest, const stpSCplx spCplxSrc)
{
    float fReal;
    float fImag;

    fReal = spCplxDest->R * spCplxSrc->R - spCplxDest->I * (-spCplxSrc->I);
    fImag = spCplxDest->R * (-spCplxSrc->I) + spCplxDest->I * spCplxSrc->R;
    spCplxDest->R = fReal;
    spCplxDest->I = fImag;
}


INLINE void clDSPOp::CplxMulC (stpDCplx spCplxDest, const stpDCplx spCplxSrc)
{
    double dReal;
    double dImag;

    dReal = spCplxDest->R * spCplxSrc->R - spCplxDest->I * (-spCplxSrc->I);
    dImag = spCplxDest->R * (-spCplxSrc->I) + spCplxDest->I * spCplxSrc->R;
    spCplxDest->R = dReal;
    spCplxDest->I = dImag;
}


INLINE void clDSPOp::CplxMulC (stpSCplx spCplxDest, const stpSCplx spCplxSrc1,
    const stpSCplx spCplxSrc2)
{
    spCplxDest->R =
        spCplxSrc1->R * spCplxSrc2->R - spCplxSrc1->I * (-spCplxSrc2->I);
    spCplxDest->I =
        spCplxSrc1->R * (-spCplxSrc2->I) + spCplxSrc1->I * spCplxSrc2->R;
}


INLINE void clDSPOp::CplxMulC (stpDCplx spCplxDest, const stpDCplx spCplxSrc1,
    const stpDCplx spCplxSrc2)
{
    spCplxDest->R =
        spCplxSrc1->R * spCplxSrc2->R - spCplxSrc1->I * (-spCplxSrc2->I);
    spCplxDest->I =
        spCplxSrc1->R * (-spCplxSrc2->I) + spCplxSrc1->I * spCplxSrc2->R;
}


INLINE void clDSPOp::CplxDiv (stpSCplx spCplxDest, const stpSCplx spCplxSrc)
{
    float fReal;
    float fImag;

    fReal = (spCplxDest->R * spCplxSrc->R + spCplxDest->I * spCplxSrc->I) /
        (spCplxSrc->R * spCplxSrc->R + spCplxSrc->I * spCplxSrc->I);
    fImag = (spCplxDest->I * spCplxSrc->R - spCplxDest->R * spCplxSrc->I) /
        (spCplxSrc->R * spCplxSrc->R + spCplxSrc->I * spCplxSrc->I);
    spCplxDest->R = fReal;
    spCplxDest->I = fImag;
}


INLINE void clDSPOp::CplxDiv (stpDCplx spCplxDest, const stpDCplx spCplxSrc)
{
    double dReal;
    double dImag;

    dReal = (spCplxDest->R * spCplxSrc->R + spCplxDest->I * spCplxSrc->I) /
        (spCplxSrc->R * spCplxSrc->R + spCplxSrc->I * spCplxSrc->I);
    dImag = (spCplxDest->I * spCplxSrc->R - spCplxDest->R * spCplxSrc->I) /
        (spCplxSrc->R * spCplxSrc->R + spCplxSrc->I * spCplxSrc->I);
    spCplxDest->R = dReal;
    spCplxDest->I = dImag;
}


INLINE void clDSPOp::CplxDiv (stpSCplx spCplxDest, const stpSCplx spCplxSrc1,
    const stpSCplx spCplxSrc2)
{
    spCplxDest->R =
        (spCplxSrc1->R * spCplxSrc2->R + spCplxSrc1->I * spCplxSrc2->I) /
        (spCplxSrc2->R * spCplxSrc2->R + spCplxSrc2->I * spCplxSrc2->I);
    spCplxDest->I =
        (spCplxSrc1->I * spCplxSrc2->R - spCplxSrc1->R * spCplxSrc2->I) /
        (spCplxSrc2->R * spCplxSrc2->R + spCplxSrc2->I * spCplxSrc2->I);
}


INLINE void clDSPOp::CplxDiv (stpDCplx spCplxDest, const stpDCplx spCplxSrc1,
    const stpDCplx spCplxSrc2)
{
    spCplxDest->R =
        (spCplxSrc1->R * spCplxSrc2->R + spCplxSrc1->I * spCplxSrc2->I) /
        (spCplxSrc2->R * spCplxSrc2->R + spCplxSrc2->I * spCplxSrc2->I);
    spCplxDest->I =
        (spCplxSrc1->I * spCplxSrc2->R - spCplxSrc1->R * spCplxSrc2->I) /
        (spCplxSrc2->R * spCplxSrc2->R + spCplxSrc2->I * spCplxSrc2->I);
}


INLINE void clDSPOp::CplxExp (stpSCplx spCplxDest, const stpSCplx spCplxSrc)
{
    float fRealExp;

    fRealExp = expf(spCplxSrc->R);
    #ifndef _GNU_SOURCE
    spCplxDest->R = fRealExp * cosf(spCplxSrc->I);
    spCplxDest->I = fRealExp * sinf(spCplxSrc->I);
    #else
    sincosf(spCplxSrc->I, &spCplxDest->I, &spCplxDest->R);
    spCplxDest->R *= fRealExp;
    spCplxDest->I *= fRealExp;
    #endif
}


INLINE void clDSPOp::CplxExp (stpDCplx spCplxDest, const stpDCplx spCplxSrc)
{
    double dRealExp;

    dRealExp = exp(spCplxSrc->R);
    #ifndef _GNU_SOURCE
    spCplxDest->R = dRealExp * cos(spCplxSrc->I);
    spCplxDest->I = dRealExp * sin(spCplxSrc->I);
    #else
    sincos(spCplxSrc->I, &spCplxDest->I, &spCplxDest->R);
    spCplxDest->R *= dRealExp;
    spCplxDest->I *= dRealExp;
    #endif
}


INLINE void clDSPOp::CplxLog (stpSCplx spCplxDest, const stpSCplx spCplxSrc)
{
    Cart2Polar(&spCplxDest->R, &spCplxDest->I, spCplxSrc->R, spCplxSrc->I);
    spCplxDest->R = logf(spCplxDest->R);
}


INLINE void clDSPOp::CplxLog (stpDCplx spCplxDest, const stpDCplx spCplxSrc)
{
    Cart2Polar(&spCplxDest->R, &spCplxDest->I, spCplxSrc->R, spCplxSrc->I);
    spCplxDest->R = log(spCplxDest->R);
}


INLINE void clDSPOp::CplxLog10 (stpSCplx spCplxDest, const stpSCplx spCplxSrc)
{
    Cart2Polar(&spCplxDest->R, &spCplxDest->I, spCplxSrc->R, spCplxSrc->I);
    spCplxDest->R = log10f(spCplxDest->R);
}


INLINE void clDSPOp::CplxLog10 (stpDCplx spCplxDest, const stpDCplx spCplxSrc)
{
    Cart2Polar(&spCplxDest->R, &spCplxDest->I, spCplxSrc->R, spCplxSrc->I);
    spCplxDest->R = log10(spCplxDest->R);
}


INLINE void clDSPOp::CplxPow (stpSCplx spCplxDest, const stpSCplx spCplxSrc,
    const stpSCplx spCplxExp)
{
    stSCplx sCplxTemp;

    CplxLog(&sCplxTemp, spCplxSrc);
    CplxMul(&sCplxTemp, spCplxExp);
    CplxExp(spCplxDest, &sCplxTemp);
}


INLINE void clDSPOp::CplxPow (stpDCplx spCplxDest, const stpDCplx spCplxSrc,
    const stpDCplx spCplxExp)
{
    stDCplx sCplxTemp;

    CplxLog(&sCplxTemp, spCplxSrc);
    CplxMul(&sCplxTemp, spCplxExp);
    CplxExp(spCplxDest, &sCplxTemp);
}


INLINE void clDSPOp::CplxRoot (stpSCplx spCplxDest, const stpSCplx spCplxSrc,
    const stpSCplx spCplxRoot)
{
    stSCplx sCplxExp;

    sCplxExp.R = 1.0F;
    sCplxExp.I = 0.0F;
    CplxDiv(&sCplxExp, spCplxRoot);
    CplxPow(spCplxDest, spCplxSrc, &sCplxExp);
}


INLINE void clDSPOp::CplxRoot (stpDCplx spCplxDest, const stpDCplx spCplxSrc,
    const stpDCplx spCplxRoot)
{
    stDCplx sCplxExp;

    sCplxExp.R = 1.0;
    sCplxExp.I = 0.0;
    CplxDiv(&sCplxExp, spCplxRoot);
    CplxPow(spCplxDest, spCplxSrc, &sCplxExp);
}


/*INLINE void clDSPOp::CplxConj (stpSCplx spCplx)
{
    spCplx->I = -(spCplx->I);
}


INLINE void clDSPOp::CplxConj (stpDCplx spCplx)
{
    spCplx->I = -(spCplx->I);
}*/


INLINE void clDSPOp::CplxConj (stpSCplx spCplxDest, const stpSCplx spCplxSrc)
{
    spCplxDest->R = spCplxSrc->R;
    spCplxDest->I = -(spCplxSrc->I);
}


INLINE void clDSPOp::CplxConj (stpDCplx spCplxDest, const stpDCplx spCplxSrc)
{
    spCplxDest->R = spCplxSrc->R;
    spCplxDest->I = -(spCplxSrc->I);
}


INLINE double clDSPOp::Multiple (long lValue)
{
    long lLoopCntr;
    double dMult = 1.0;

    for (lLoopCntr = 1L; lLoopCntr <= lValue; lLoopCntr++)
    {
        dMult *= (double) lLoopCntr;
    }
    return dMult;
}


INLINE float clDSPOp::ModZeroBessel (float fValue)
{
    long lLoopCntr;
    float fMZBessel = 0.0F;
    float fHalfValue;

    fHalfValue = fValue / 2.0F;
    for (lLoopCntr = 0L; lLoopCntr <= DSP_MAXBESSEL; lLoopCntr++)
    {
        fMZBessel += (float) pow(
            pow(fHalfValue, lLoopCntr) / Multiple(lLoopCntr), 2.0);
    }
    return fMZBessel;
}


INLINE double clDSPOp::ModZeroBessel (double dValue)
{
    long lLoopCntr;
    double dMZBessel = 0.0;
    double dHalfValue;

    dHalfValue = dValue / 2.0;
    for (lLoopCntr = 0L; lLoopCntr <= DSP_MAXBESSEL; lLoopCntr++)
    {
        dMZBessel += pow(
            pow(dHalfValue, lLoopCntr) / Multiple(lLoopCntr), 2.0);
    }
    return dMZBessel;
}


float clDSPOp::ChebyshevPolynom (float fOrder, float fValue)
{
    if (fabsf(fValue) <= 1.0F)
    {
        return cosf(fOrder * acosf(fValue));
    }
    else
    {
        return coshf(fOrder * acoshf(fValue));
    }
}


double clDSPOp::ChebyshevPolynom (double dOrder, double dValue)
{
    if (fabs(dValue) <= 1.0)
    {
        return cos(dOrder * acos(dValue));
    }
    else
    {
        return cosh(dOrder * acosh(dValue));
    }
}


clDSPOp::clDSPOp()
{
    #ifdef DSP_X86
    unsigned int uiNewCSR;

    bHave3DNow = (dsp_x86_have_e3dnow()) ? true : false;
    iHaveSSE = dsp_x86_have_sse();
    #ifdef __INTEL_COMPILER
    if (iHaveSSE)
    {
        _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
        _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);
    }
    #else
    if (iHaveSSE)
    {
        uiOldCSR = _mm_getcsr();
        uiNewCSR = uiOldCSR | 0x8040;
        _mm_setcsr(uiNewCSR);
    }
    #endif
    #endif
    lPrevSrcCount = 0;
    lPrevDestCount = 0;
    // This is to get maximum precision for PI
    #ifdef _ISOC9X_SOURCE
    fPI = (float) acosl(-1.0);
    dPI = (double) acosl(-1.0);
    #else
    fPI = (float) acos(-1.0);
    dPI = acos(-1.0);
    #endif
    bFFTInitialized = false;
    lpSBitRevWork = NULL;
    lpDBitRevWork = NULL;
    fpCosSinTable = NULL;
    dpCosSinTable = NULL;
}


clDSPOp::~clDSPOp()
{
    if (bFFTInitialized)
    {
        FFTUninitialize();
    }
}


void clDSPOp::Add (float *fpDest, float fSrc, long lCount)
{
    #ifdef DSP_IPP
    ippsAddC_32f_I(fSrc, fpDest, lCount);
    #else
    long lLoopCntr;

    #ifdef DSP_X86
    if (bHave3DNow)
    {
        dsp_x86_3dnow_addf(fpDest, fSrc, lCount);
    }
    else if (iHaveSSE >= 2)
    {
        dsp_x86_sse_addf(fpDest, fSrc, lCount);
    }
    else
    #endif
    {
        for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
        {
            fpDest[lLoopCntr] += fSrc;
        }
    }
    #endif
}


void clDSPOp::Add (double *dpDest, double dSrc, long lCount)
{
    #ifdef DSP_IPP
    ippsAddC_64f_I(dSrc, dpDest, lCount);
    #else
    long lLoopCntr;

    #ifdef DSP_X86
    if (iHaveSSE >= 2)
    {
        dsp_x86_sse_add(dpDest, dSrc, lCount);
    }
    else
    #endif
    {
        for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
        {
            dpDest[lLoopCntr] += dSrc;
        }
    }
    #endif
}


void clDSPOp::Add (stpSCplx spCplxDest, stSCplx sCplxSrc, long lCount)
{
    #ifdef DSP_IPP
    ippsAddC_32fc_I(*((Ipp32fc *) &sCplxSrc), (Ipp32fc *) spCplxDest, lCount);
    #else
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        CplxAdd(&spCplxDest[lLoopCntr], &sCplxSrc);
    }
    #endif
}


void clDSPOp::Add (stpDCplx spCplxDest, stDCplx sCplxSrc, long lCount)
{
    #ifdef DSP_IPP
    ippsAddC_64fc_I(*((Ipp64fc *) &sCplxSrc), (Ipp64fc *) spCplxDest, lCount);
    #else
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        CplxAdd(&spCplxDest[lLoopCntr], &sCplxSrc);
    }
    #endif
}


void clDSPOp::Add (float *fpDest, const float *fpSrc, long lCount)
{
    #ifdef DSP_IPP
    ippsAdd_32f_I(fpSrc, fpDest, lCount);
    #else
    long lLoopCntr;

    #ifdef DSP_X86
    if (bHave3DNow)
    {
        dsp_x86_3dnow_add2f(fpDest, fpSrc, lCount);
    }
    else if (iHaveSSE >= 2)
    {
        dsp_x86_sse_add2f(fpDest, fpSrc, lCount);
    }
    else
    #endif
    {
        for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
        {
            fpDest[lLoopCntr] += fpSrc[lLoopCntr];
        }
    }
    #endif
}


void clDSPOp::Add (double *dpDest, const double *dpSrc, long lCount)
{
    #ifdef DSP_IPP
    ippsAdd_64f_I(dpSrc, dpDest, lCount);
    #else
    long lLoopCntr;

    #ifdef DSP_X86
    if (iHaveSSE >= 2)
    {
        dsp_x86_sse_add2(dpDest, dpSrc, lCount);
    }
    else
    #endif
    {
        for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
        {
            dpDest[lLoopCntr] += dpSrc[lLoopCntr];
        }
    }
    #endif
}


void clDSPOp::Add (stpSCplx spCplxDest, const stpSCplx spCplxSrc, long lCount)
{
    #ifdef DSP_IPP
    ippsAdd_32fc_I((Ipp32fc *) spCplxSrc, (Ipp32fc *) spCplxDest, lCount);
    #else
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        CplxAdd(&spCplxDest[lLoopCntr], &spCplxSrc[lLoopCntr]);
    }
    #endif
}


void clDSPOp::Add (stpDCplx spCplxDest, const stpDCplx spCplxSrc, long lCount)
{
    #ifdef DSP_IPP
    ippsAdd_64fc_I((Ipp64fc *) spCplxSrc, (Ipp64fc *) spCplxDest, lCount);
    #else
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        CplxAdd(&spCplxDest[lLoopCntr], &spCplxSrc[lLoopCntr]);
    }
    #endif
}


void clDSPOp::Add (float *fpDest, const float *fpSrc1, const float *fpSrc2, 
    long lCount)
{
    #ifdef DSP_IPP
    ippsAdd_32f(fpSrc1, fpSrc2, fpDest, lCount);
    #else
    long lLoopCntr;

    #ifdef DSP_X86
    if (bHave3DNow)
    {
        dsp_x86_3dnow_add3f(fpDest, fpSrc1, fpSrc2, lCount);
    }
    else if (iHaveSSE >= 2)
    {
        dsp_x86_sse_add3f(fpDest, fpSrc1, fpSrc2, lCount);
    }
    else
    #endif
    {
        for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
        {
            fpDest[lLoopCntr] = fpSrc1[lLoopCntr] + fpSrc2[lLoopCntr];
        }
    }
    #endif
}


void clDSPOp::Add (double *dpDest, const double *dpSrc1, const double *dpSrc2,
    long lCount)
{
    #ifdef DSP_IPP
    ippsAdd_64f(dpSrc1, dpSrc2, dpDest, lCount);
    #else
    long lLoopCntr;

    #ifdef DSP_X86
    if (iHaveSSE >= 2)
    {
        dsp_x86_sse_add3(dpDest, dpSrc1, dpSrc2, lCount);
    }
    else
    #endif
    {
        for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
        {
            dpDest[lLoopCntr] = dpSrc1[lLoopCntr] + dpSrc2[lLoopCntr];
        }
    }
    #endif
}


void clDSPOp::Add (stpSCplx spCplxDest, const stpSCplx spCplxSrc1,
    const stpSCplx spCplxSrc2, long lCount)
{
    #ifdef DSP_IPP
    ippsAdd_32fc((Ipp32fc *) spCplxSrc1, (Ipp32fc *) spCplxSrc2,
        (Ipp32fc *) spCplxDest, lCount);
    #else
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        CplxAdd(&spCplxDest[lLoopCntr], &spCplxSrc1[lLoopCntr],
            &spCplxSrc2[lLoopCntr]);
    }
    #endif
}


void clDSPOp::Add (stpDCplx spCplxDest, const stpDCplx spCplxSrc1,
    const stpDCplx spCplxSrc2, long lCount)
{
    #ifdef DSP_IPP
    ippsAdd_64fc((Ipp64fc *) spCplxSrc1, (Ipp64fc *) spCplxSrc2,
        (Ipp64fc *) spCplxDest, lCount);
    #else
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        CplxAdd(&spCplxDest[lLoopCntr], &spCplxSrc1[lLoopCntr],
            &spCplxSrc2[lLoopCntr]);
    }
    #endif
}


void clDSPOp::Sub (float *fpDest, float fSrc, long lCount)
{
    #ifdef DSP_IPP
    ippsSubC_32f_I(fSrc, fpDest, lCount);
    #else
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        fpDest[lLoopCntr] -= fSrc;
    }
    #endif
}


void clDSPOp::Sub (double *dpDest, double dSrc, long lCount)
{
    #ifdef DSP_IPP
    ippsSubC_64f_I(dSrc, dpDest, lCount);
    #else
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        dpDest[lLoopCntr] -= dSrc;
    }
    #endif
}


void clDSPOp::Sub (stpSCplx spCplxDest, stSCplx sCplxSrc, long lCount)
{
    #ifdef DSP_IPP
    ippsSubC_32fc_I(*((Ipp32fc *) &sCplxSrc), (Ipp32fc *) spCplxDest, lCount);
    #else
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        CplxSub(&spCplxDest[lLoopCntr], &sCplxSrc);
    }
    #endif
}


void clDSPOp::Sub (stpDCplx spCplxDest, stDCplx sCplxSrc, long lCount)
{
    #ifdef DSP_IPP
    ippsSubC_64fc_I(*((Ipp64fc *) &sCplxSrc), (Ipp64fc *) spCplxDest, lCount);
    #else
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        CplxSub(&spCplxDest[lLoopCntr], &sCplxSrc);
    }
    #endif
}


void clDSPOp::Sub (float *fpDest, const float *fpSrc, long lCount)
{
    #ifdef DSP_IPP
    ippsSub_32f_I(fpSrc, fpDest, lCount);
    #else
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        fpDest[lLoopCntr] -= fpSrc[lLoopCntr];
    }
    #endif
}


void clDSPOp::Sub (double *dpDest, const double *dpSrc, long lCount)
{
    #ifdef DSP_IPP
    ippsSub_64f_I(dpSrc, dpDest, lCount);
    #else
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        dpDest[lLoopCntr] -= dpSrc[lLoopCntr];
    }
    #endif
}


void clDSPOp::Sub (stpSCplx spCplxDest, const stpSCplx spCplxSrc, long lCount)
{
    #ifdef DSP_IPP
    ippsSub_32fc_I((Ipp32fc *) spCplxSrc, (Ipp32fc *) spCplxDest, lCount);
    #else
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        CplxSub(&spCplxDest[lLoopCntr], &spCplxSrc[lLoopCntr]);
    }
    #endif
}


void clDSPOp::Sub (stpDCplx spCplxDest, const stpDCplx spCplxSrc, long lCount)
{
    #ifdef DSP_IPP
    ippsSub_64fc_I((Ipp64fc *) spCplxSrc, (Ipp64fc *) spCplxDest, lCount);
    #else
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        CplxSub(&spCplxDest[lLoopCntr], &spCplxSrc[lLoopCntr]);
    }
    #endif
}


void clDSPOp::Sub (float *fpDest, const float *fpSrc1, const float *fpSrc2, 
    long lCount)
{
    #ifdef DSP_IPP
    ippsSub_32f(fpSrc1, fpSrc2, fpDest, lCount);
    #else
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        fpDest[lLoopCntr] = fpSrc1[lLoopCntr] - fpSrc2[lLoopCntr];
    }
    #endif
}


void clDSPOp::Sub (double *dpDest, const double *dpSrc1, const double *dpSrc2,
    long lCount)
{
    #ifdef DSP_IPP
    ippsSub_64f(dpSrc1, dpSrc2, dpDest, lCount);
    #else
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        dpDest[lLoopCntr] = dpSrc1[lLoopCntr] - dpSrc2[lLoopCntr];
    }
    #endif
}


void clDSPOp::Sub (stpSCplx spCplxDest, const stpSCplx spCplxSrc1,
    const stpSCplx spCplxSrc2, long lCount)
{
    #ifdef DSP_IPP
    ippsSub_32fc((Ipp32fc *) spCplxSrc1, (Ipp32fc *) spCplxSrc2,
        (Ipp32fc *) spCplxDest, lCount);
    #else
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        CplxSub(&spCplxDest[lLoopCntr], &spCplxSrc1[lLoopCntr],
            &spCplxSrc2[lLoopCntr]);
    }
    #endif
}


void clDSPOp::Sub (stpDCplx spCplxDest, const stpDCplx spCplxSrc1,
    const stpDCplx spCplxSrc2, long lCount)
{
    #ifdef DSP_IPP
    ippsSub_64fc((Ipp64fc *) spCplxSrc1, (Ipp64fc *) spCplxSrc2,
        (Ipp64fc *) spCplxDest, lCount);
    #else
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        CplxSub(&spCplxDest[lLoopCntr], &spCplxSrc1[lLoopCntr],
            &spCplxSrc2[lLoopCntr]);
    }
    #endif
}



void clDSPOp::Mul (float *fpVect, float fSrc, long lCount)
{
    #ifdef DSP_IPP
    ippsMulC_32f_I(fSrc, fpVect, lCount);
    #else
    long lLoopCntr;

    #ifdef DSP_X86
    if (bHave3DNow)
    {
        dsp_x86_3dnow_mulf(fpVect, fSrc, lCount);
    }
    else if (iHaveSSE >= 2)
    {
        dsp_x86_sse_mulf(fpVect, fSrc, lCount);
    }
    else
    #endif
    {
        for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
        {
            fpVect[lLoopCntr] *= fSrc;
        }
    }
    #endif
}


void clDSPOp::Mul (double *dpVect, double dSrc, long lCount)
{
    #ifdef DSP_IPP
    ippsMulC_64f_I(dSrc, dpVect, lCount);
    #else
    long lLoopCntr;

    #ifdef DSP_X86
    if (iHaveSSE >= 2)
    {
        dsp_x86_sse_mul(dpVect, dSrc, lCount);
    }
    else
    #endif
    {
        for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
        {
            dpVect[lLoopCntr] *= dSrc;
        }
    }
    #endif
}


void clDSPOp::Mul (stpSCplx spCplxDest, float fSrc, long lCount)
{
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        CplxMul(&spCplxDest[lLoopCntr], fSrc);
    }
}


void clDSPOp::Mul (stpDCplx spCplxDest, double dSrc, long lCount)
{
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        CplxMul(&spCplxDest[lLoopCntr], dSrc);
    }
}


void clDSPOp::Mul (stpSCplx spCplxDest, stSCplx sCplxSrc, long lCount)
{
    #ifdef DSP_IPP
    ippsMulC_32fc_I(*((Ipp32fc *) &sCplxSrc), (Ipp32fc *) spCplxDest, lCount);
    #else
    long lLoopCntr;

    #ifdef DSP_X86
    if (bHave3DNow)
    {
        dsp_x86_3dnow_cmulf((float *) spCplxDest,
            (const float *) &sCplxSrc, lCount);
    }
    else if (iHaveSSE == 2)
    {
        dsp_x86_sse_cmulf((float *) spCplxDest,
            (const float *) &sCplxSrc, lCount);
    }
    else if (iHaveSSE >= 3)
    {
        dsp_x86_sse3_cmulf((float *) spCplxDest,
            (const float *) &sCplxSrc, lCount);
    }
    else
    #endif
    {
        for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
        {
            CplxMul(&spCplxDest[lLoopCntr], &sCplxSrc);
        }
    }
    #endif
}


void clDSPOp::Mul (stpDCplx spCplxDest, stDCplx sCplxSrc, long lCount)
{
    #ifdef DSP_IPP
    ippsMulC_64fc_I(*((Ipp64fc *) &sCplxSrc), (Ipp64fc *) spCplxDest, lCount);
    #else
    long lLoopCntr;

    #ifdef DSP_X86
    if (iHaveSSE == 2)
    {
        dsp_x86_sse_cmul((double *) spCplxDest,
            (const double *) &sCplxSrc, lCount);
    }
    else if (iHaveSSE >= 3)
    {
        dsp_x86_sse3_cmul((double *) spCplxDest,
            (const double *) &sCplxSrc, lCount);
    }
    else
    #endif
    {
        for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
        {
            CplxMul(&spCplxDest[lLoopCntr], &sCplxSrc);
        }
    }
    #endif
}


void clDSPOp::Mul (float *fpDest, const float *fpSrc1, float fSrc2, 
    long lCount)
{
    #ifdef DSP_IPP
    ippsMulC_32f(fpSrc1, fSrc2, fpDest, lCount);
    #else
    long lLoopCntr;

    #ifdef DSP_X86
    if (bHave3DNow)
    {
        dsp_x86_3dnow_mulf_nip(fpDest, fpSrc1, fSrc2, lCount);
    }
    else if (iHaveSSE >= 2)
    {
        dsp_x86_sse_mulf_nip(fpDest, fpSrc1, fSrc2, lCount);
    }
    else
    #endif
    {
        for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
        {
            fpDest[lLoopCntr] = fpSrc1[lLoopCntr] * fSrc2;
        }
    }
    #endif
}


void clDSPOp::Mul (double *dpDest, const double *dpSrc1, double dSrc2,
    long lCount)
{
    #ifdef DSP_IPP
    ippsMulC_64f(dpSrc1, dSrc2, dpDest, lCount);
    #else
    long lLoopCntr;

    #ifdef DSP_X86
    if (iHaveSSE >= 2)
    {
        dsp_x86_sse_mul_nip(dpDest, dpSrc1, dSrc2, lCount);
    }
    else
    #endif
    {
        for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
        {
            dpDest[lLoopCntr] = dpSrc1[lLoopCntr] * dSrc2;
        }
    }
    #endif
}


void clDSPOp::Mul (float *fpDest, const float *fpSrc, long lCount)
{
    #ifdef DSP_IPP
    ippsMul_32f_I(fpSrc, fpDest, lCount);
    #else
    long lLoopCntr;

    #ifdef DSP_X86
    if (bHave3DNow)
    {
        dsp_x86_3dnow_mul2f(fpDest, fpSrc, lCount);
    }
    else if (iHaveSSE >= 2)
    {
        dsp_x86_sse_mul2f(fpDest, fpSrc, lCount);
    }
    else
    #endif
    {
        for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
        {
            fpDest[lLoopCntr] *= fpSrc[lLoopCntr];
        }
    }
    #endif
}


void clDSPOp::Mul (double *dpDest, const double *dpSrc, long lCount)
{
    #ifdef DSP_IPP
    ippsMul_64f_I(dpSrc, dpDest, lCount);
    #else
    long lLoopCntr;

    #ifdef DSP_X86
    if (iHaveSSE >= 2)
    {
        dsp_x86_sse_mul2(dpDest, dpSrc, lCount);
    }
    else
    #endif
    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        dpDest[lLoopCntr] *= dpSrc[lLoopCntr];
    }
    #endif
}


void clDSPOp::Mul (stpSCplx spCplxDest, const float *fpSrc, long lCount)
{
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        CplxMul(&spCplxDest[lLoopCntr], fpSrc[lLoopCntr]);
    }
}


void clDSPOp::Mul (stpDCplx spCplxDest, const double *dpSrc, long lCount)
{
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        CplxMul(&spCplxDest[lLoopCntr], dpSrc[lLoopCntr]);
    }
}


void clDSPOp::Mul (stpSCplx spCplxDest, const stpSCplx spCplxSrc, long lCount)
{
    #ifdef DSP_IPP
    ippsMul_32fc_I((Ipp32fc *) spCplxSrc, (Ipp32fc *) spCplxDest, lCount);
    #else
    long lLoopCntr;

    #ifdef DSP_X86
    if (bHave3DNow)
    {
        dsp_x86_3dnow_cmul2f((float *) spCplxDest, (const float *) spCplxSrc, 
            lCount);
    }
    else if (iHaveSSE == 2)
    {
        dsp_x86_sse_cmul2f((float *) spCplxDest, (const float *) spCplxSrc,
            lCount);
    }
    else if (iHaveSSE >= 3)
    {
        dsp_x86_sse3_cmul2f((float *) spCplxDest, (const float *) spCplxSrc,
            lCount);
    }
    else
    #endif
    {
        for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
        {
            CplxMul(&spCplxDest[lLoopCntr], &spCplxSrc[lLoopCntr]);
        }
    }
    #endif
}


void clDSPOp::Mul (stpDCplx spCplxDest, const stpDCplx spCplxSrc, long lCount)
{
    #ifdef DSP_IPP
    ippsMul_64fc_I((Ipp64fc *) spCplxSrc, (Ipp64fc *) spCplxDest, lCount);
    #else
    long lLoopCntr;

    #ifdef DSP_X86
    if (iHaveSSE == 2)
    {
        dsp_x86_sse_cmul2((double *) spCplxDest, (const double *) spCplxSrc,
            lCount);
    }
    else if (iHaveSSE >= 3)
    {
        dsp_x86_sse3_cmul2((double *) spCplxDest, (const double *) spCplxSrc,
            lCount);
    }
    else
    #endif
    {
        for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
        {
            CplxMul(&spCplxDest[lLoopCntr], &spCplxSrc[lLoopCntr]);
        }
    }
    #endif
}


void clDSPOp::Mul (float *fpDest, const float *fpSrc1, 
    const float *fpSrc2, long lCount)
{
    #ifdef DSP_IPP
    ippsMul_32f(fpSrc1, fpSrc2, fpDest, lCount);
    #else
    long lLoopCntr;

    #ifdef DSP_X86
    if (bHave3DNow)
    {
        dsp_x86_3dnow_mul3f(fpDest, fpSrc1, fpSrc2, lCount);
    }
    else if (iHaveSSE >= 2)
    {
        dsp_x86_sse_mul3f(fpDest, fpSrc1, fpSrc2, lCount);
    }
    else
    #endif
    {
        for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
        {
            fpDest[lLoopCntr] = fpSrc1[lLoopCntr] * fpSrc2[lLoopCntr];
        }
    }
    #endif
}


void clDSPOp::Mul (double *dpDest, const double *dpSrc1,
    const double *dpSrc2, long lCount)
{
    #ifdef DSP_IPP
    ippsMul_64f(dpSrc1, dpSrc2, dpDest, lCount);
    #else
    long lLoopCntr;

    #ifdef DSP_X86
    if (iHaveSSE >= 2)
    {
        dsp_x86_sse_mul3(dpDest, dpSrc1, dpSrc2, lCount);
    }
    else
    #endif
    {
        for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
        {
            dpDest[lLoopCntr] = dpSrc1[lLoopCntr] * dpSrc2[lLoopCntr];
        }
    }
    #endif
}


void clDSPOp::Mul (stpSCplx spCplxDest, const stpSCplx spCplxSrc1,
    const stpSCplx spCplxSrc2, long lCount)
{
    #ifdef DSP_IPP
    ippsMul_64fc((Ipp64fc *) spCplxSrc1, (Ipp64fc *) spCplxSrc2, 
        (Ipp64fc *) spCplxDest, lCount);
    #else
    long lLoopCntr;

    #ifdef DSP_X86
    if (bHave3DNow)
    {
        dsp_x86_3dnow_cmul3f((float *) spCplxDest, 
            (const float *) spCplxSrc1,
            (const float *) spCplxSrc2,
            lCount);
    }
    else if (iHaveSSE == 2)
    {
        dsp_x86_sse_cmul3f((float *) spCplxDest,
            (const float *) spCplxSrc1,
            (const float *) spCplxSrc2,
            lCount);
    }
    else if (iHaveSSE >= 3)
    {
        dsp_x86_sse3_cmul3f((float *) spCplxDest,
            (const float *) spCplxSrc1,
            (const float *) spCplxSrc2,
            lCount);
    }
    else
    #endif
    {
        for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
        {
            CplxMul(&spCplxDest[lLoopCntr], &spCplxSrc1[lLoopCntr],
                &spCplxSrc2[lLoopCntr]);
        }
    }
    #endif
}


void clDSPOp::Mul (stpDCplx spCplxDest, const stpDCplx spCplxSrc1,
    const stpDCplx spCplxSrc2, long lCount)
{
    #ifdef DSP_IPP
    ippsMul_64fc((Ipp64fc *) spCplxSrc1, (Ipp64fc *) spCplxSrc2,
        (Ipp64fc *) spCplxDest, lCount);
    #else
    long lLoopCntr;

    #ifdef DSP_X86
    if (iHaveSSE == 2)
    {
        dsp_x86_sse_cmul3((double *) spCplxDest,
            (const double *) spCplxSrc1,
            (const double *) spCplxSrc2,
            lCount);
    }
    else if (iHaveSSE >= 3)
    {
        dsp_x86_sse3_cmul3((double *) spCplxDest,
            (const double *) spCplxSrc1,
            (const double *) spCplxSrc2,
            lCount);
    }
    else
    #endif
    {
        for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
        {
            CplxMul(&spCplxDest[lLoopCntr], &spCplxSrc1[lLoopCntr],
                &spCplxSrc2[lLoopCntr]);
        }
    }
    #endif
}


void clDSPOp::MulC (stpSCplx spCplxDest, const stpSCplx spCplxSrc, long lCount)
{
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        CplxMulC(&spCplxDest[lLoopCntr], &spCplxSrc[lLoopCntr]);
    }
}


void clDSPOp::MulC (stpDCplx spCplxDest, const stpDCplx spCplxSrc, long lCount)
{
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        CplxMulC(&spCplxDest[lLoopCntr], &spCplxSrc[lLoopCntr]);
    }
}


void clDSPOp::MulC (stpSCplx spCplxDest, const stpSCplx spCplxSrc1,
    const stpSCplx spCplxSrc2, long lCount)
{
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        CplxMulC(&spCplxDest[lLoopCntr], &spCplxSrc1[lLoopCntr],
            &spCplxSrc2[lLoopCntr]);
    }
}


void clDSPOp::MulC (stpDCplx spCplxDest, const stpDCplx spCplxSrc1,
    const stpDCplx spCplxSrc2, long lCount)
{
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        CplxMulC(&spCplxDest[lLoopCntr], &spCplxSrc1[lLoopCntr],
            &spCplxSrc2[lLoopCntr]);
    }
}


void clDSPOp::Mul2 (float *fpDst1, float *fpDst2, const float *fpSrc, 
    long lCount)
{
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        fpDst1[lLoopCntr] *= fpSrc[lLoopCntr];
        fpDst2[lLoopCntr] *= fpSrc[lLoopCntr];
    }
}


void clDSPOp::Mul2 (double *dpDst1, double *dpDst2, const double *dpSrc, 
    long lCount)
{
   long lLoopCntr;

   for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
   {
      dpDst1[lLoopCntr] *= dpSrc[lLoopCntr];
      dpDst2[lLoopCntr] *= dpSrc[lLoopCntr];
   }
}


void clDSPOp::Mul2 (float *fpDst1, float *fpDst2, const float *fpSrc1,
    const float *fpSrc2, const float *fpMul, long lCount)
{
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        fpDst1[lLoopCntr] = fpSrc1[lLoopCntr] * fpMul[lLoopCntr];
        fpDst2[lLoopCntr] = fpSrc2[lLoopCntr] * fpMul[lLoopCntr];
    }
}


void clDSPOp::Mul2 (double *dpDst1, double *dpDst2, const double *dpSrc1,
   const double *dpSrc2, const double *dpMul, long lCount)
{
   long lLoopCntr;

   for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
   {
      dpDst1[lLoopCntr] = dpSrc1[lLoopCntr] * dpMul[lLoopCntr];
      dpDst2[lLoopCntr] = dpSrc2[lLoopCntr] * dpMul[lLoopCntr];
   }
}


void clDSPOp::Div (float *fpVect, float fSrc, long lCount)
{
    #ifdef DSP_IPP
    ippsDivC_32f_I(fSrc, fpVect, lCount);
    #else
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        fpVect[lLoopCntr] /= fSrc;
    }
    #endif
}


void clDSPOp::Div (double *dpVect, double dSrc, long lCount)
{
    #ifdef DSP_IPP
    ippsDivC_64f_I(dSrc, dpVect, lCount);
    #else
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        dpVect[lLoopCntr] /= dSrc;
    }
    #endif
}


void clDSPOp::Div (stpSCplx spCplxDest, stSCplx sCplxSrc, long lCount)
{
    #ifdef DSP_IPP
    ippsDivC_32fc_I(*((Ipp32fc *) &sCplxSrc), (Ipp32fc *) spCplxDest, lCount);
    #else
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        CplxDiv(&spCplxDest[lLoopCntr], &sCplxSrc);
    }
    #endif
}


void clDSPOp::Div (stpDCplx spCplxDest, stDCplx sCplxSrc, long lCount)
{
    #ifdef DSP_IPP
    ippsDivC_64fc_I(*((Ipp64fc *) &sCplxSrc), (Ipp64fc *) spCplxDest, lCount);
    #else
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        CplxDiv(&spCplxDest[lLoopCntr], &sCplxSrc);
    }
    #endif
}


void clDSPOp::Div (float *fpDest, const float *fpSrc, long lCount)
{
    #ifdef DSP_IPP
    ippsDiv_32f_I(fpSrc, fpDest, lCount);
    #else
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        fpDest[lLoopCntr] /= fpSrc[lLoopCntr];
    }
    #endif
}


void clDSPOp::Div (double *dpDest, const double *dpSrc, long lCount)
{
    #ifdef DSP_IPP
    ippsDiv_64f_I(dpSrc, dpDest, lCount);
    #else
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        dpDest[lLoopCntr] /= dpSrc[lLoopCntr];
    }
    #endif
}


void clDSPOp::Div (stpSCplx spCplxDest, const stpSCplx spCplxSrc, long lCount)
{
    #ifdef DSP_IPP
    ippsDiv_32fc_I((Ipp32fc *) spCplxSrc, (Ipp32fc *) spCplxDest, lCount);
    #else
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        CplxDiv(&spCplxDest[lLoopCntr], &spCplxSrc[lLoopCntr]);
    }
    #endif
}


void clDSPOp::Div (stpDCplx spCplxDest, const stpDCplx spCplxSrc, long lCount)
{
    #ifdef DSP_IPP
    ippsDiv_64fc_I((Ipp64fc *) spCplxSrc, (Ipp64fc *) spCplxDest, lCount);
    #else
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        CplxDiv(&spCplxDest[lLoopCntr], &spCplxSrc[lLoopCntr]);
    }
    #endif
}


void clDSPOp::Div (float *fpDest, const float *fpSrc1, const float *fpSrc2, 
    long lCount)
{
    #ifdef DSP_IPP
    ippsDiv_32f(fpSrc1, fpSrc2, fpDest, lCount);
    #else
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        fpDest[lLoopCntr] = fpSrc1[lLoopCntr] / fpSrc2[lLoopCntr];
    }
    #endif
}


void clDSPOp::Div (double *dpDest, const double *dpSrc1, const double *dpSrc2,
    long lCount)
{
    #ifdef DSP_IPP
    ippsDiv_64f(dpSrc1, dpSrc2, dpDest, lCount);
    #else
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        dpDest[lLoopCntr] = dpSrc1[lLoopCntr] / dpSrc2[lLoopCntr];
    }
    #endif
}


void clDSPOp::Div (stpSCplx spCplxDest, const stpSCplx spCplxSrc1,
    const stpSCplx spCplxSrc2, long lCount)
{
    #ifdef DSP_IPP
    ippsDiv_32fc((Ipp32fc *) spCplxSrc1, (Ipp32fc *) spCplxSrc2,
        (Ipp32fc *) spCplxDest, lCount);
    #else
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        CplxDiv(&spCplxDest[lLoopCntr], &spCplxSrc1[lLoopCntr],
            &spCplxSrc2[lLoopCntr]);
    }
    #endif
}


void clDSPOp::Div (stpDCplx spCplxDest, const stpDCplx spCplxSrc1,
    const stpDCplx spCplxSrc2, long lCount)
{
    #ifdef DSP_IPP
    ippsDiv_64fc((Ipp64fc *) spCplxSrc1, (Ipp64fc *) spCplxSrc2,
        (Ipp64fc *) spCplxDest, lCount);
    #else
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        CplxDiv(&spCplxDest[lLoopCntr], &spCplxSrc1[lLoopCntr],
            &spCplxSrc2[lLoopCntr]);
    }
    #endif
}


void clDSPOp::Div1x (float *fpVect, long lCount)
{
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        fpVect[lLoopCntr] = 1.0F / fpVect[lLoopCntr];
    }
}


void clDSPOp::Div1x (double *dpVect, long lCount)
{
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        dpVect[lLoopCntr] = 1.0 / dpVect[lLoopCntr];
    }
}


void clDSPOp::Div1x (float *fpDest, const float *fpSrc, long lCount)
{
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        fpDest[lLoopCntr] = 1.0F / fpSrc[lLoopCntr];
    }
}


void clDSPOp::Div1x (double *dpDest, const double *dpSrc, long lCount)
{
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        dpDest[lLoopCntr] = 1.0 / dpSrc[lLoopCntr];
    }
}


void clDSPOp::MulAdd (float *fpVect, float fMul, float fAdd, long lCount)
{
    long lLoopCntr;
    
    #ifdef DSP_X86
    if (bHave3DNow)
    {
        dsp_x86_3dnow_maf(fpVect, fMul, fAdd, lCount);
    }
    else if (iHaveSSE >= 2)
    {
        dsp_x86_sse_maf(fpVect, fMul, fAdd, lCount);
    }
    else
    #endif
    {
        for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
        {
            fpVect[lLoopCntr] = fpVect[lLoopCntr] * fMul + fAdd;
        }
    }
}


void clDSPOp::MulAdd (double *dpVect, double dMul, double dAdd, long lCount)
{
    long lLoopCntr;

    #ifdef DSP_X86
    if (iHaveSSE >= 2)
    {
        dsp_x86_sse_ma(dpVect, dMul, dAdd, lCount);
    }
    else
    #endif
    {    
        for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
        {
            dpVect[lLoopCntr] = dpVect[lLoopCntr] * dMul + dAdd;
        }
    }
}


void clDSPOp::MulAdd (float *fpDest, const float *fpSrc, 
    float fMul, float fAdd, long lCount)
{
    long lLoopCntr;
    
    #ifdef DSP_X86
    if (bHave3DNow)
    {
        dsp_x86_3dnow_ma2f(fpDest, fpSrc, fMul, fAdd, lCount);
    }
    else if (iHaveSSE >= 2)
    {
        dsp_x86_sse_ma2f(fpDest, fpSrc, fMul, fAdd, lCount);
    }
    else
    #endif
    {
        for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
        {
            fpDest[lLoopCntr] = fpSrc[lLoopCntr] * fMul + fAdd;
        }
    }
}


void clDSPOp::MulAdd (double *dpDest, const double *dpSrc, 
    double dMul, double dAdd, long lCount)
{
    long lLoopCntr;
    
    #ifdef DSP_X86
    if (iHaveSSE >= 2)
    {
        dsp_x86_sse_ma2(dpDest, dpSrc, dMul, dAdd, lCount);
    }
    else
    #endif
    {
        for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
        {
            dpDest[lLoopCntr] = dpSrc[lLoopCntr] * dMul + dAdd;
        }
    }
}


void clDSPOp::MulAddC (stpSCplx spDest, const stpSCplx spSrc1, 
    const stpSCplx spSrc2, long lCount)
{
    long lLoopCntr;
    stSCplx sTemp;
    
    #ifdef DSP_X86
    if (bHave3DNow)
    {
        dsp_x86_3dnow_cmaf((float *) spDest, (const float *) spSrc1,
            (const float *) spSrc2, lCount);
    }
    else if (iHaveSSE == 2)
    {
        dsp_x86_sse_cmaf((float *) spDest, (const float *) spSrc1,
            (const float *) spSrc2, lCount);
    }
    else if (iHaveSSE >= 3)
    {
        dsp_x86_sse3_cmaf((float *) spDest, (const float *) spSrc1,
            (const float *) spSrc2, lCount);
    }
    else
    #endif
    {
        for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
        {
            CplxMul(&sTemp, &spSrc1[lLoopCntr], &spSrc2[lLoopCntr]);
            CplxAdd(&spDest[lLoopCntr], &sTemp);
        }
    }
}


void clDSPOp::MulAddC (stpDCplx spDest, const stpDCplx spSrc1, 
    const stpDCplx spSrc2, long lCount)
{
    long lLoopCntr;
    stDCplx sTemp;
    
    #ifdef DSP_X86
    if (iHaveSSE == 2)
    {
        dsp_x86_sse_cma((double *) spDest, (const double *) spSrc1,
            (const double *) spSrc2, lCount);
    }
    else if (iHaveSSE >= 3)
    {
        dsp_x86_sse3_cma((double *) spDest, (const double *) spSrc1,
            (const double *) spSrc2, lCount);
    }
    else
    #endif
    {
        for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
        {
            CplxMul(&sTemp, &spSrc1[lLoopCntr], &spSrc2[lLoopCntr]);
            CplxAdd(&spDest[lLoopCntr], &sTemp);
        }
    }
}


void clDSPOp::Abs (float *fpVect, long lCount)
{
    #ifdef DSP_IPP
    ippsAbs_32f_I(fpVect, lCount);
    #else
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        fpVect[lLoopCntr] = fabsf(fpVect[lLoopCntr]);
    }
    #endif
}


void clDSPOp::Abs (double *dpVect, long lCount)
{
    #ifdef DSP_IPP
    ippsAbs_64f_I(dpVect, lCount);
    #else
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        dpVect[lLoopCntr] = fabs(dpVect[lLoopCntr]);
    }
    #endif
}


void clDSPOp::Abs (float *fpDest, const float *fpSrc, long lCount)
{
    #ifdef DSP_IPP
    ippsAbs_32f(fpSrc, fpDest, lCount);
    #else
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        fpDest[lLoopCntr] = fabsf(fpSrc[lLoopCntr]);
    }
    #endif
}


void clDSPOp::Abs (double *dpDest, const double *dpSrc, long lCount)
{
    #ifdef DSP_IPP
    ippsAbs_64f(dpSrc, dpDest, lCount);
    #else
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        dpDest[lLoopCntr] = fabs(dpSrc[lLoopCntr]);
    }
    #endif
}


void clDSPOp::Sqrt (float *fpVect, long lCount)
{
    #ifdef DSP_IPP
    ippsSqrt_32f_I(fpVect, lCount);
    #else
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        fpVect[lLoopCntr] = sqrtf(fpVect[lLoopCntr]);
    }
    #endif
}


void clDSPOp::Sqrt (double *dpVect, long lCount)
{
    #ifdef DSP_IPP
    ippsSqrt_64f_I(dpVect, lCount);
    #else
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        dpVect[lLoopCntr] = sqrt(dpVect[lLoopCntr]);
    }
    #endif
}


void clDSPOp::Sqrt (float *fpDest, const float *fpSrc, long lCount)
{
    #ifdef DSP_IPP
    ippsSqrt_32f(fpSrc, fpDest, lCount);
    #else
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        fpDest[lLoopCntr] = sqrtf(fpSrc[lLoopCntr]);
    }
    #endif
}


void clDSPOp::Sqrt (double *dpDest, const double *dpSrc, long lCount)
{
    #ifdef DSP_IPP
    ippsSqrt_64f(dpSrc, dpDest, lCount);
    #else
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        dpDest[lLoopCntr] = sqrt(dpSrc[lLoopCntr]);
    }
    #endif
}


void clDSPOp::Zero (float *fpDest, long lCount)
{
    #ifdef DSP_IPP
    ippsZero_32f(fpDest, lCount);
    #else
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        fpDest[lLoopCntr] = 0.0F;
    }
    #endif
}


void clDSPOp::Zero (double *dpDest, long lCount)
{
    #ifdef DSP_IPP
    ippsZero_64f(dpDest, lCount);
    #else
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        dpDest[lLoopCntr] = 0.0;
    }
    #endif
}


void clDSPOp::Zero (stpSCplx spCplxDest, long lCount)
{
    #ifdef DSP_IPP
    ippsZero_32fc((Ipp32fc *) spCplxDest, lCount);
    #else
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        spCplxDest[lLoopCntr].R = 0.0F;
        spCplxDest[lLoopCntr].I = 0.0F;
    }
    #endif
}


void clDSPOp::Zero (stpDCplx spCplxDest, long lCount)
{
    #ifdef DSP_IPP
    ippsZero_64fc((Ipp64fc *) spCplxDest, lCount);
    #else
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        spCplxDest[lLoopCntr].R = 0.0;
        spCplxDest[lLoopCntr].I = 0.0;
    }
    #endif
}


void clDSPOp::Set (float *fpDest, float fSrc, long lCount)
{
    #ifdef DSP_IPP
    ippsSet_32f(fSrc, fpDest, lCount);
    #else
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        fpDest[lLoopCntr] = fSrc;
    }
    #endif
}


void clDSPOp::Set (double *dpDest, double dSrc, long lCount)
{
    #ifdef DSP_IPP
    ippsSet_64f(dSrc, dpDest, lCount);
    #else
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        dpDest[lLoopCntr] = dSrc;
    }
    #endif
}


void clDSPOp::Set (stpSCplx spCplxDest, stSCplx sCplxSrc, long lCount)
{
    #ifdef DSP_IPP
    ippsSet_32fc(*((Ipp32fc *) &sCplxSrc), (Ipp32fc *) spCplxDest, lCount);
    #else
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        spCplxDest[lLoopCntr].R = sCplxSrc.R;
        spCplxDest[lLoopCntr].I = sCplxSrc.I;
    }
    #endif
}


void clDSPOp::Set (stpDCplx spCplxDest, stDCplx sCplxSrc, long lCount)
{
    #ifdef DSP_IPP
    ippsSet_64fc(*((Ipp64fc *) &sCplxSrc), (Ipp64fc *) spCplxDest, lCount);
    #else
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        spCplxDest[lLoopCntr].R = sCplxSrc.R;
        spCplxDest[lLoopCntr].I = sCplxSrc.I;
    }
    #endif
}


void clDSPOp::Set (float *fpDest, float fSrc, long lStart, long lCount,
    long lLength)
{
    #ifdef DSP_IPP
    ippsSet_32f(fSrc, &fpDest[lStart],
        ((lStart + lCount) > lLength) ? lLength - lStart : lCount);
    #else
    long lLoopCntr;
    long lEnd;

    lEnd = ((lStart + lCount) > lLength) ? lLength : (lStart + lCount);
    for (lLoopCntr = 0L; lLoopCntr < lEnd; lLoopCntr++)
    {
        fpDest[lLoopCntr] = fSrc;
    }
    #endif
}


void clDSPOp::Set (double *dpDest, double dSrc, long lStart, long lCount,
    long lLength)
{
    #ifdef DSP_IPP
    ippsSet_64f(dSrc, &dpDest[lStart],
        ((lStart + lCount) > lLength) ? lLength - lStart : lCount);
    #else
    long lLoopCntr;
    long lEnd;

    lEnd = ((lStart + lCount) > lLength) ? lLength : (lStart + lCount);
    for (lLoopCntr = 0L; lLoopCntr < lEnd; lLoopCntr++)
    {
        dpDest[lLoopCntr] = dSrc;
    }
    #endif
}


void clDSPOp::Set (stpSCplx spCplxDest, stSCplx sCplxSrc, long lStart, 
    long lCount, long lLength)
{
    #ifdef DSP_IPP
    ippsSet_32fc(*((Ipp32fc *) &sCplxSrc), (Ipp32fc *) spCplxDest,
        ((lStart + lCount) > lLength) ? lLength - lStart : lCount);
    #else
    long lLoopCntr;
    long lEnd;

    lEnd = ((lStart + lCount) > lLength) ? lLength : (lStart + lCount);
    for (lLoopCntr = 0L; lLoopCntr < lEnd; lLoopCntr++)
    {
        spCplxDest[lLoopCntr].R = sCplxSrc.R;
        spCplxDest[lLoopCntr].I = sCplxSrc.I;
    }
    #endif
}


void clDSPOp::Set (stpDCplx spCplxDest, stDCplx sCplxSrc, long lStart,
    long lCount, long lLength)
{
    #ifdef DSP_IPP
    ippsSet_64fc(*((Ipp64fc *) &sCplxSrc), (Ipp64fc *) spCplxDest,
        ((lStart + lCount) > lLength) ? lLength - lStart : lCount);
    #else
    long lLoopCntr;
    long lEnd;

    lEnd = ((lStart + lCount) > lLength) ? lLength : (lStart + lCount);
    for (lLoopCntr = 0L; lLoopCntr < lEnd; lLoopCntr++)
    {
        spCplxDest[lLoopCntr].R = sCplxSrc.R;
        spCplxDest[lLoopCntr].I = sCplxSrc.I;
    }
    #endif
}


void clDSPOp::Clip (float *fpVect, float fValue, long lCount)
{
    #ifdef DSP_IPP
    ippsThreshold_GT_32f_I(fpVect, lCount, fValue);
    #else
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        #if (!defined(_ISOC9X_SOURCE) || ((__GNUC__ == 3) && (__GNUC_MINOR__ < 1)))
        if (fpVect[lLoopCntr] > fValue)
        {
            fpVect[lLoopCntr] = fValue;
        }
        #else
        if (isgreater(fpVect[lLoopCntr], fValue))
        {
            fpVect[lLoopCntr] = fValue;
        }
        #endif
    }
    #endif
}


void clDSPOp::Clip (double *dpVect, double dValue, long lCount)
{
    #ifdef DSP_IPP
    ippsThreshold_GT_64f_I(dpVect, lCount, dValue);
    #else
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        #if (!defined(_ISOC9X_SOURCE) || ((__GNUC__ == 3) && (__GNUC_MINOR__ < 1)))
        if (dpVect[lLoopCntr] > dValue)
        {
            dpVect[lLoopCntr] = dValue;
        }
        #else
        if (isgreater(dpVect[lLoopCntr], dValue))
        {
            dpVect[lLoopCntr] = dValue;
        }
        #endif
    }
    #endif
}


void clDSPOp::Clip (float *fpDest, const float *fpSrc, float fValue, 
    long lCount)
{
    #ifdef DSP_IPP
    ippsThreshold_GT_32f(fpSrc, fpDest, lCount, fValue);
    #else
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        #if (!defined(_ISOC9X_SOURCE) || ((__GNUC__ == 3) && (__GNUC_MINOR__ < 1)))
        fpDest[lLoopCntr] = (fpSrc[lLoopCntr] > fValue) ? 
            fValue : fpSrc[lLoopCntr];
        #else
        fpDest[lLoopCntr] = (isgreater(fpSrc[lLoopCntr], fValue)) ?
            fValue : fpSrc[lLoopCntr];
        #endif
    }
    #endif
}


void clDSPOp::Clip (double *dpDest, const double *dpSrc, double dValue,
    long lCount)
{
    #ifdef DSP_IPP
    ippsThreshold_GT_64f(dpSrc, dpDest, lCount, dValue);
    #else
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        #if (!defined(_ISOC9X_SOURCE) || ((__GNUC__ == 3) && (__GNUC_MINOR__ < 1)))
        dpDest[lLoopCntr] = (dpSrc[lLoopCntr] > dValue) ?
            dValue : dpSrc[lLoopCntr];
        #else
        dpDest[lLoopCntr] = (isgreater(dpSrc[lLoopCntr], dValue)) ?
            dValue : dpSrc[lLoopCntr];
        #endif
    }
    #endif
}


void clDSPOp::Clip (float *fpVect, float fMin, float fMax, long lCount)
{
    #ifdef DSP_IPP
    ippsThreshold_LT_32f_I(fpVect, lCount, fMin);
    ippsThreshold_GT_32f_I(fpVect, lCount, fMax);
    #else
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        #if (!defined(_ISOC9X_SOURCE) || ((__GNUC__ == 3) && (__GNUC_MINOR__ < 1)))
        if (fpVect[lLoopCntr] < fMin)
        {
            fpVect[lLoopCntr] = fMin;
        }
        else if (fpVect[lLoopCntr] > fMax)
        {
            fpVect[lLoopCntr] = fMax;
        }
        #else
        if (isless(fpVect[lLoopCntr], fMin))
        {
            fpVect[lLoopCntr] = fMin;
        }
        else if (isgreater(fpVect[lLoopCntr], fMax))
        {
            fpVect[lLoopCntr] = fMax;
        }
        #endif
    }
    #endif
}


void clDSPOp::Clip (double *dpVect, double dMin, double dMax, long lCount)
{
    #ifdef DSP_IPP
    ippsThreshold_LT_64f_I(dpVect, lCount, dMin);
    ippsThreshold_GT_64f_I(dpVect, lCount, dMax);
    #else
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        #if (!defined(_ISOC9X_SOURCE) || ((__GNUC__ == 3) && (__GNUC_MINOR__ < 1)))
        if (dpVect[lLoopCntr] < dMin)
        {
            dpVect[lLoopCntr] = dMin;
        }
        else if (dpVect[lLoopCntr] > dMax)
        {
            dpVect[lLoopCntr] = dMax;
        }
        #else
        if (isless(dpVect[lLoopCntr], dMin))
        {
            dpVect[lLoopCntr] = dMin;
        }
        else if (isgreater(dpVect[lLoopCntr], dMax))
        {
            dpVect[lLoopCntr] = dMax;
        }
        #endif
    }
    #endif
}


void clDSPOp::Clip (float *fpDest, const float *fpSrc, float fMin, 
    float fMax, long lCount)
{
    #ifdef DSP_IPP
    ippsThreshold_LT_32f(fpSrc, fpDest, lCount, fMin);
    ippsThreshold_GT_32f_I(fpDest, lCount, fMax);
    #else
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        #if (!defined(_ISOC9X_SOURCE) || ((__GNUC__ == 3) && (__GNUC_MINOR__ < 1)))
        if (fpSrc[lLoopCntr] < fMin)
        {
            fpDest[lLoopCntr] = fMin;
        }
        else if (fpSrc[lLoopCntr] > fMax)
        {
            fpDest[lLoopCntr] = fMax;
        }
        else
        {
            fpDest[lLoopCntr] = fpSrc[lLoopCntr];
        }
        #else
        if (isless(fpSrc[lLoopCntr], fMin))
        {
            fpDest[lLoopCntr] = fMin;
        }
        else if (isgreater(fpSrc[lLoopCntr], fMax))
        {
            fpDest[lLoopCntr] = fMax;
        }
        else
        {
            fpDest[lLoopCntr] = fpSrc[lLoopCntr];
        }
        #endif
    }
    #endif
}


void clDSPOp::Clip (double *dpDest, const double *dpSrc, double dMin, 
    double dMax, long lCount)
{
    #ifdef DSP_IPP
    ippsThreshold_LT_64f(dpSrc, dpDest, lCount, dMin);
    ippsThreshold_GT_64f_I(dpDest, lCount, dMin);
    #else
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        #if (!defined(_ISOC9X_SOURCE) || ((__GNUC__ == 3) && (__GNUC_MINOR__ < 1)))
        if (dpSrc[lLoopCntr] < dMin)
        {
            dpDest[lLoopCntr] = dMin;
        }
        else if (dpSrc[lLoopCntr] > dMax)
        {
            dpDest[lLoopCntr] = dMax;
        }
        else
        {
            dpDest[lLoopCntr] = dpSrc[lLoopCntr];
        }
        #else
        if (isless(dpSrc[lLoopCntr], dMin))
        {
            dpDest[lLoopCntr] = dMin;
        }
        else if (isgreater(dpSrc[lLoopCntr], dMax))
        {
            dpDest[lLoopCntr] = dMax;
        }
        else
        {
           dpDest[lLoopCntr] = dpSrc[lLoopCntr];
        }
        #endif
    }
    #endif
}


void clDSPOp::ClipZero (float *fpVect, long lCount)
{
    #ifdef DSP_IPP
    ippsThreshold_LT_32f_I(fpVect, lCount, 0.0f);
    #else
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        #if (!defined(_ISOC9X_SOURCE) || ((__GNUC__ == 3) && (__GNUC_MINOR__ < 1)))
        if (fpVect[lLoopCntr] < 0.0F)
        {
            fpVect[lLoopCntr] = 0.0F;
        }
        #else
        if (isless(fpVect[lLoopCntr], 0.0F))
        {
            fpVect[lLoopCntr] = 0.0F;
        }
        #endif
    }
    #endif
}


void clDSPOp::ClipZero (double *dpVect, long lCount)
{
    #ifdef DSP_IPP
    ippsThreshold_LT_64f_I(dpVect, lCount, 0.0);
    #else
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        #if (!defined(_ISOC9X_SOURCE) || ((__GNUC__ == 3) && (__GNUC_MINOR__ < 1)))
        if (dpVect[lLoopCntr] < 0.0)
        {
            dpVect[lLoopCntr] = 0.0;
        }
        #else
        if (isless(dpVect[lLoopCntr], 0.0))
        {
            dpVect[lLoopCntr] = 0.0;
        }
        #endif
    }
    #endif
}


void clDSPOp::ClipZero (float *fpDest, const float *fpSrc, long lCount)
{
    #ifdef DSP_IPP
    ippsThreshold_LT_32f(fpSrc, fpDest, lCount, 0.0f);
    #else
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        #if (!defined(_ISOC9X_SOURCE) || ((__GNUC__ == 3) && (__GNUC_MINOR__ < 1)))
        fpDest[lLoopCntr] = (fpSrc[lLoopCntr] < 0.0F) ? 
            0.0F : fpSrc[lLoopCntr];
        #else
        fpDest[lLoopCntr] = (isless(fpSrc[lLoopCntr], 0.0F)) ?
            0.0F : fpSrc[lLoopCntr];
        #endif
    }
    #endif
}


void clDSPOp::ClipZero (double *dpDest, const double *dpSrc, long lCount)
{
    #ifdef DSP_IPP
    ippsThreshold_LT_64f(dpSrc, dpDest, lCount, 0.0);
    #else
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        #if (!defined(_ISOC9X_SOURCE) || ((__GNUC__ == 3) && (__GNUC_MINOR__ < 1)))
        dpDest[lLoopCntr] = (dpSrc[lLoopCntr] < 0.0) ?
            0.0 : dpSrc[lLoopCntr];
        #else
        dpDest[lLoopCntr] = (isless(dpSrc[lLoopCntr], 0.0)) ?
            0.0 : dpSrc[lLoopCntr];
        #endif
    }
    #endif
}


void clDSPOp::Copy (float *fpDest, const float *fpSrc, long lCount)
{
    #ifdef DSP_IPP
    ippsMove_32f(fpSrc, fpDest, lCount);
    #else
    #ifndef USE_MEMMOVE
    long lLoopCntr;

    #ifdef DSP_X86
    if (bHave3DNow)
    {
        if (likely(fpDest <= fpSrc || abs(fpDest - fpSrc) > lCount))
        {
            dsp_x86_3dnow_copyf(fpDest, fpSrc, lCount);
        }
        else
        {
            for (lLoopCntr = lCount - 1L; lLoopCntr >= 0L; lLoopCntr--)
            {
                fpDest[lLoopCntr] = fpSrc[lLoopCntr];
            }
        }
    }
    else
    #endif
    {
        if (likely(fpDest <= fpSrc || abs(fpDest - fpSrc) > lCount))
        {
            for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
            {
                fpDest[lLoopCntr] = fpSrc[lLoopCntr];
            }
        }
        else
        {
            for (lLoopCntr = lCount - 1L; lLoopCntr >= 0L; lLoopCntr--)
            {
                fpDest[lLoopCntr] = fpSrc[lLoopCntr];
            }
        }
    }
    #else
    #ifndef DSP_X86
    memmove(fpDest, fpSrc, lCount * sizeof(float));
    #else
    memmove(fpDest, fpSrc, (lCount << 2));
    #endif
    #endif
    #endif
}


void clDSPOp::Copy (double *dpDest, const double *dpSrc, long lCount)
{
    #ifdef DSP_IPP
    ippsMove_64f(dpSrc, dpDest, lCount);
    #else
    #ifndef USE_MEMMOVE
    long lLoopCntr;

    #ifdef DSP_X86
    if (bHave3DNow)
    {
        if (likely(dpDest <= dpSrc || abs(dpDest - dpSrc) > lCount))
        {
            dsp_x86_3dnow_copyd(dpDest, dpSrc, lCount);
        }
        else
        {
            for (lLoopCntr = lCount - 1L; lLoopCntr >= 0L; lLoopCntr--)
            {
                dpDest[lLoopCntr] = dpSrc[lLoopCntr];
            }
        }
    }
    else
    #endif
    {
        if (likely(dpDest <= dpSrc || abs(dpDest - dpSrc) > lCount))
        {
            for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
            {
                dpDest[lLoopCntr] = dpSrc[lLoopCntr];
            }
        }
        else
        {
            for (lLoopCntr = lCount - 1L; lLoopCntr >= 0L; lLoopCntr--)
            {
                dpDest[lLoopCntr] = dpSrc[lLoopCntr];
            }
        }
    }
    #else
    #ifndef DSP_X86
    memmove(dpDest, dpSrc, lCount * sizeof(double));
    #else
    memmove(dpDest, dpSrc, (lCount << 3));
    #endif
    #endif
    #endif
}


void clDSPOp::Copy (stpSCplx spCplxDest, const stpSCplx spCplxSrc, long lCount)
{
    #ifdef DSP_IPP
    ippsMove_32fc((Ipp32fc *) spCplxSrc, (Ipp32fc *) spCplxDest, lCount);
    #else
    #ifndef USE_MEMMOVE
    long lLoopCntr;

    #ifdef DSP_X86
    if (bHave3DNow)
    {
        if (likely(spCplxDest <= spCplxSrc || 
            abs(spCplxDest - spCplxSrc) > lCount))
        {
            dsp_x86_3dnow_copyf((float *) spCplxDest, 
                (const float *) spCplxSrc, (lCount << 1));
        }
        else
        {
            for (lLoopCntr = lCount - 1L; lLoopCntr >= 0L; lLoopCntr--)
            {
                spCplxDest[lLoopCntr].R = spCplxSrc[lLoopCntr].R;
                spCplxDest[lLoopCntr].I = spCplxSrc[lLoopCntr].I;
            }
        }
    }
    else
    #endif
    {
        if (likely(spCplxDest <= spCplxSrc || 
            abs(spCplxDest - spCplxSrc) > lCount))
        {
            for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
            {
                spCplxDest[lLoopCntr].R = spCplxSrc[lLoopCntr].R;
                spCplxDest[lLoopCntr].I = spCplxSrc[lLoopCntr].I;
            }
        }
        else
        {
            for (lLoopCntr = lCount - 1L; lLoopCntr >= 0L; lLoopCntr--)
            {
                spCplxDest[lLoopCntr].R = spCplxSrc[lLoopCntr].R;
                spCplxDest[lLoopCntr].I = spCplxSrc[lLoopCntr].I;
            }
        }
    }
    #else
    memmove(spCplxDest, spCplxSrc, lCount * sizeof(stSCplx));
    #endif
    #endif
}


void clDSPOp::Copy (stpDCplx spCplxDest, const stpDCplx spCplxSrc, long lCount)
{
    #ifdef DSP_IPP
    ippsMove_64fc((Ipp64fc *) spCplxSrc, (Ipp64fc *) spCplxDest, lCount);
    #else
    #ifndef USE_MEMMOVE
    long lLoopCntr;

    #ifdef DSP_X86
    if (bHave3DNow)
    {
        if (likely(spCplxDest <= spCplxSrc || 
            abs(spCplxDest - spCplxSrc) > lCount))
        {
            dsp_x86_3dnow_copyd((double *) spCplxDest, 
                (const double *) spCplxSrc, (lCount << 1));
        }
        else
        {
            for (lLoopCntr = lCount - 1L; lLoopCntr >= 0L; lLoopCntr--)
            {
                spCplxDest[lLoopCntr].R = spCplxSrc[lLoopCntr].R;
                spCplxDest[lLoopCntr].I = spCplxSrc[lLoopCntr].I;
            }
        }
    }
    else
    #endif
    {
        if (likely(spCplxDest <= spCplxSrc || 
            abs(spCplxDest - spCplxSrc) > lCount))
        {
            for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
            {
                spCplxDest[lLoopCntr].R = spCplxSrc[lLoopCntr].R;
                spCplxDest[lLoopCntr].I = spCplxSrc[lLoopCntr].I;
            }
        }
        else
        {
            for (lLoopCntr = lCount - 1L; lLoopCntr >= 0L; lLoopCntr--)
            {
                spCplxDest[lLoopCntr].R = spCplxSrc[lLoopCntr].R;
                spCplxDest[lLoopCntr].I = spCplxSrc[lLoopCntr].I;
            }
        }
    }
    #else
    memmove(spCplxDest, spCplxSrc, lCount * sizeof(stDCplx));
    #endif
    #endif
}


void clDSPOp::Copy (float *fpDest1, float *fpDest2, const float *fpSrc,
    long lCount)
{
    #ifdef DSP_IPP
    ippsMove_32f(fpSrc, fpDest1, lCount);
    ippsMove_32f(fpSrc, fpDest2, lCount);
    #else
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        fpDest1[lLoopCntr] = fpDest2[lLoopCntr] = fpSrc[lLoopCntr];
    }
    #endif
}


void clDSPOp::Copy (double *dpDest1, double *dpDest2, const double *dpSrc,
    long lCount)
{
    #ifdef DSP_IPP
    ippsMove_64f(dpSrc, dpDest1, lCount);
    ippsMove_64f(dpSrc, dpDest2, lCount);
    #else
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        dpDest1[lLoopCntr] = dpDest2[lLoopCntr] = dpSrc[lLoopCntr];
    }
    #endif
}


float clDSPOp::Convolve (const float *fpSrc1, const float *fpSrc2, 
    long lCount)
{
    long lLoopCntr;
    long lMax;
    float fConv = 0.0F;

    lMax = lCount - 1L;
    for (lLoopCntr = 0L; lLoopCntr <= lMax; lLoopCntr++)
    {
        #ifndef _ISOC9X_SOURCE
        fConv += fpSrc1[lLoopCntr] * fpSrc2[lMax - lLoopCntr];
        #else
        fConv = fmaf(fpSrc1[lLoopCntr], fpSrc2[lMax - lLoopCntr], fConv);
        #endif
    }
    return fConv;
}


double clDSPOp::Convolve (const double *dpSrc1, const double *dpSrc2,
    long lCount)
{
    long lLoopCntr;
    long lMax;
    double dConv = 0.0;

    lMax = lCount - 1L;
    for (lLoopCntr = 0L; lLoopCntr <= lMax; lLoopCntr++)
    {
        #ifndef _ISOC9X_SOURCE
        dConv += dpSrc1[lLoopCntr] * dpSrc2[lMax - lLoopCntr];
        #else
        dConv = fma(dpSrc1[lLoopCntr], dpSrc2[lMax - lLoopCntr], dConv);
        #endif
    }
    return dConv;
}


void clDSPOp::Convolve (float *fpDest, const float *fpSrc1, 
    const float *fpSrc2, long lCount)
{
    long lLoopDest;
    long lLoopConv;
    long lIdx;
    long lMax;
    float fConv;

    lMax = lCount - 1L;
    for (lLoopDest = 0L; lLoopDest < lCount; lLoopDest++)
    {
        fConv = 0.0F;
        for (lLoopConv = 0L; lLoopConv <= lMax; lLoopConv++)
        {
            lIdx = ((lLoopConv - lLoopDest) < 0L) ?
                (lLoopConv - lLoopDest + lMax) : (lLoopConv - lLoopDest);
            #ifndef _ISOC9X_SOURCE
            fConv += fpSrc1[lMax - lIdx] * fpSrc2[lIdx];
            #else
            fConv = fmaf(fpSrc1[lMax - lIdx], fpSrc2[lIdx], fConv);
            #endif
        }
        fpDest[lLoopDest] = fConv;
    }
}


void clDSPOp::Convolve (double *dpDest, const double *dpSrc1,
    const double *dpSrc2, long lCount)
{
    long lLoopDest;
    long lLoopConv;
    long lIdx;
    long lMax;
    double dConv;

    lMax = lCount - 1L;
    for (lLoopDest = 0L; lLoopDest < lCount; lLoopDest++)
    {
        dConv = 0.0;
        for (lLoopConv = 0L; lLoopConv <= lMax; lLoopConv++)
        {
            lIdx = ((lLoopConv - lLoopDest) < 0L) ?
                (lLoopConv - lLoopDest + lMax) : (lLoopConv - lLoopDest);
            #ifndef _ISOC9X_SOURCE
            dConv += dpSrc1[lMax - lIdx] * dpSrc2[lIdx];
            #else
            dConv = fma(dpSrc1[lMax - lIdx], dpSrc2[lIdx], dConv);
            #endif
        }
        dpDest[lLoopDest] = dConv;
    }
}


float clDSPOp::Correlate (const float *fpSrc1, const float *fpSrc2, 
    long lCount)
{
    long lLoopCntr;
    float fCorr = 0.0F;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        #ifndef _ISOC9X_SOURCE
        fCorr += fpSrc1[lLoopCntr] * fpSrc2[lLoopCntr];
        #else
        fCorr = fmaf(fpSrc1[lLoopCntr], fpSrc2[lLoopCntr], fCorr);
        #endif
    }
    fCorr /= (float) lCount;
    return fCorr;
}


double clDSPOp::Correlate (const double *dpSrc1, const double *dpSrc2,
    long lCount)
{
    long lLoopCntr;
    double dCorr = 0.0;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        #ifndef _ISOC9X_SOURCE
        dCorr += dpSrc1[lLoopCntr] * dpSrc2[lLoopCntr];
        #else
        dCorr = fma(dpSrc1[lLoopCntr], dpSrc2[lLoopCntr], dCorr);
        #endif
    }
    dCorr /= (double) lCount;
    return dCorr;
}


void clDSPOp::Correlate (float *fpDest, const float *fpSrc1, 
    const float *fpSrc2, long lCount)
{
    long lLoopDest;
    long lLoopCorr;
    long lMax;
    long lIdx;
    float fCorr;

    lMax = lCount - 1L;
    for (lLoopDest = 0L; lLoopDest <= lMax; lLoopDest++)
    {
        fCorr = 0.0F;
        for (lLoopCorr = 0L; lLoopCorr <= lMax; lLoopCorr++)
        {
            lIdx = ((lLoopCorr + lLoopDest) > lMax) ?
                (lLoopCorr + lLoopDest - lMax) : (lLoopCorr + lLoopDest);
            #ifndef _ISOC9X_SOURCE
            fCorr += fpSrc1[lLoopCorr] * fpSrc2[lIdx];
            #else
            fCorr = fmaf(fpSrc1[lLoopCorr], fpSrc2[lIdx], fCorr);
            #endif
        }
        fpDest[lLoopDest] = fCorr / (float) lCount;
    }
}


void clDSPOp::Correlate (double *dpDest, const double *dpSrc1,
    const double *dpSrc2, long lCount)
{
    long lLoopDest;
    long lLoopCorr;
    long lMax;
    long lIdx;
    double dCorr;

    lMax = lCount - 1L;
    for (lLoopDest = 0L; lLoopDest <= lMax; lLoopDest++)
    {
        dCorr = 0.0;
        for (lLoopCorr = 0L; lLoopCorr <= lMax; lLoopCorr++)
        {
            lIdx = ((lLoopCorr + lLoopDest) > lMax) ?
                (lLoopCorr + lLoopDest - lMax) : (lLoopCorr + lLoopDest);
            #ifndef _ISOC9X_SOURCE
            dCorr += dpSrc1[lLoopCorr] * dpSrc2[lIdx];
            #else
            dCorr = fma(dpSrc1[lLoopCorr], dpSrc2[lIdx], dCorr);
            #endif
        }
        dpDest[lLoopDest] = dCorr / (double) lCount;
    }
}


float clDSPOp::AutoCorrelate (const float *fpSrc, long lCount)
{
    long lLoopCntr;
    float fAutoCorr = 0.0F;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        #ifndef _ISOC9X_SOURCE
        fAutoCorr += fpSrc[lLoopCntr] * fpSrc[lLoopCntr];
        #else
        fAutoCorr = fmaf(fpSrc[lLoopCntr], fpSrc[lLoopCntr], fAutoCorr);
        #endif
    }
    fAutoCorr /= (float) lCount;
    return fAutoCorr;
}


double clDSPOp::AutoCorrelate (const double *dpSrc, long lCount)
{
    long lLoopCntr;
    double dAutoCorr = 0.0;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        #ifndef _ISOC9X_SOURCE
        dAutoCorr += dpSrc[lLoopCntr] * dpSrc[lLoopCntr];
        #else
        dAutoCorr = fma(dpSrc[lLoopCntr], dpSrc[lLoopCntr], dAutoCorr);
        #endif
    }
    dAutoCorr /= (double) lCount;
    return dAutoCorr;
}


void clDSPOp::AutoCorrelate (float *fpDest, const float *fpSrc, long lCount)
{
    long lLoopDest;
    long lLoopCorr;
    long lMax;
    long lIdx;
    float fAutoCorr;

    lMax = lCount - 1L;
    for (lLoopDest = 0L; lLoopDest <= lMax; lLoopDest++)
    {
        fAutoCorr = 0.0F;
        for (lLoopCorr = 0L; lLoopCorr <= lMax; lLoopCorr++)
        {
            lIdx = ((lLoopCorr + lLoopDest) > lMax) ?
                (lLoopCorr + lLoopDest - lCount) : (lLoopCorr + lLoopDest);
            #ifndef _ISOC9X_SOURCE
            fAutoCorr += fpSrc[lLoopCorr] * fpSrc[lIdx];
            #else
            fAutoCorr = fmaf(fpSrc[lLoopCorr], fpSrc[lIdx], fAutoCorr);
            #endif
        }
        fpDest[lLoopDest] = fAutoCorr / (float) lCount;
    }
}


void clDSPOp::AutoCorrelate (double *dpDest, const double *dpSrc, long lCount)
{
    long lLoopDest;
    long lLoopCorr;
    long lMax;
    long lIdx;
    double dAutoCorr;

    lMax = lCount - 1L;
    for (lLoopDest = 0L; lLoopDest <= lMax; lLoopDest++)
    {
        dAutoCorr = 0.0;
        for (lLoopCorr = 0L; lLoopCorr <= lMax; lLoopCorr++)
        {
            lIdx = ((lLoopCorr + lLoopDest) > lMax) ?
                (lLoopCorr + lLoopDest - lCount) : (lLoopCorr + lLoopDest);
            #ifndef _ISOC9X_SOURCE
            dAutoCorr += dpSrc[lLoopCorr] * dpSrc[lIdx];
            #else
            dAutoCorr = fma(dpSrc[lLoopCorr], dpSrc[lIdx], dAutoCorr);
            #endif
        }
        dpDest[lLoopDest] = dAutoCorr / (double) lCount;
    }
}


float clDSPOp::DotProduct (const float *fpSrc1, const float *fpSrc2,
    long lCount)
{
    #ifdef DSP_IPP
    float fDotProd;

    ippsDotProd_32f(fpSrc1, fpSrc2, lCount, &fDotProd);
    #else
    long lLoopCntr;
    float fDotProd = 0.0F;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        #ifndef _ISOC9X_SOURCE
        fDotProd += fpSrc1[lLoopCntr] * fpSrc2[lLoopCntr];
        #else
        fDotProd = fmaf(fpSrc1[lLoopCntr], fpSrc2[lLoopCntr], fDotProd);
        #endif
    }
    #endif
    return fDotProd;
}


double clDSPOp::DotProduct (const double *dpSrc1, const double *dpSrc2,
    long lCount)
{
    #ifdef DSP_IPP
    double dDotProd;

    ippsDotProd_64f(dpSrc1, dpSrc2, lCount, &dDotProd);
    #else
    long lLoopCntr;
    double dDotProd = 0.0;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        #ifndef _ISOC9X_SOURCE
        dDotProd += dpSrc1[lLoopCntr] * dpSrc2[lLoopCntr];
        #else
        dDotProd = fma(dpSrc1[lLoopCntr], dpSrc2[lLoopCntr], dDotProd);
        #endif
    }
    #endif
    return dDotProd;
}


void clDSPOp::MinMax (float *fpMin, float *fpMax, const float *fpSrc,
    long lCount)
{
    #ifdef DSP_IPP
    ippsMin_32f(fpSrc, lCount, fpMin);
    ippsMax_32f(fpSrc, lCount, fpMax);
    #else

    #ifdef DSP_X86
    if (bHave3DNow)
    {
        dsp_x86_3dnow_minmaxf(fpMin, fpMax, fpSrc, lCount);
    }
    else if (iHaveSSE >= 2)
    {
        dsp_x86_sse_minmaxf(fpMin, fpMax, fpSrc, lCount);
    }
    else
    #endif
    #if (defined(DSP_X86_64) && (!defined(DSP_X86)))
    dsp_x86_64_minmaxf(fpMin, fpMax, fpSrc, lCount);
    #else
    {
        long lLoopCntr;
        float fTempVal;
        float fTempMin = FLT_MAX;
        float fTempMax = -FLT_MAX;

        for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
        {
            fTempVal = fpSrc[lLoopCntr];
            #ifndef _ISOC9X_SOURCE
            if (fTempVal < fTempMin)
            {
                fTempMin = fTempVal;
            }
            if (fTempVal > fTempMax)
            {
                fTempMax = fTempVal;
            }
            #else
            fTempMin = fminf(fTempVal, fTempMin);
            fTempMax = fmaxf(fTempVal, fTempMax);
            #endif
        }
        *fpMin = fTempMin;
        *fpMax = fTempMax;
    }
    #endif
    #endif
}


void clDSPOp::MinMax (double *dpMin, double *dpMax, const double *dpSrc,
    long lCount)
{
    #ifdef DSP_IPP
    ippsMin_64f(dpSrc, lCount, dpMin);
    ippsMax_64f(dpSrc, lCount, dpMax);
    #else

    #ifdef DSP_X86
    if (iHaveSSE >= 2)
    {
        dsp_x86_sse_minmax(dpMin, dpMax, dpSrc, lCount);
    }
    else
    #endif
    #if (defined(DSP_X86_64) && (!defined(DSP_X86)))
    dsp_x86_64_minmax(dpMin, dpMax, dpSrc, lCount);
    #else
    {
        long lLoopCntr;
        double dTempVal;
        double dTempMin = DBL_MAX;
        double dTempMax = -DBL_MAX;

        for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
        {
            dTempVal = dpSrc[lLoopCntr];
            #ifndef _ISOC9X_SOURCE
            if (dTempVal < dTempMin)
            {
                dTempMin = dTempVal;
            }
            if (dTempVal > dTempMax)
            {
                dTempMax = dTempVal;
            }
            #else
            dTempMin = fmin(dTempVal, dTempMin);
            dTempMax = fmax(dTempVal, dTempMax);
            #endif
        }
        *dpMin = dTempMin;
        *dpMax = dTempMax;
    }
    #endif
    #endif
}


float clDSPOp::Mean (const float *fpSrc, long lCount)
{
    #ifdef DSP_IPP
    float fMean;

    ippsMean_32f(fpSrc, lCount, &fMean, ippAlgHintFast);
    #else
    long lLoopCntr;
    float fMean = 0.0F;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        fMean += fpSrc[lLoopCntr];
    }
    fMean /= (float) lCount;
    #endif
    return fMean;
}


double clDSPOp::Mean (const double *dpSrc, long lCount)
{
    #ifdef DSP_IPP
    double dMean;

    ippsMean_64f(dpSrc, lCount, &dMean);
    #else
    long lLoopCntr;
    double dMean = 0.0;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        dMean += dpSrc[lLoopCntr];
    }
    dMean /= (double) lCount;
    #endif
    return dMean;
}


float clDSPOp::Median (const float *fpSrc, long lCount)
{
    long lMax;
    float fMedian = 0.0F;
    float *fpTemp;
    clDSPAlloc Temp(lCount * sizeof(float));

    fpTemp = Temp;
    lMax = lCount - 1L;
    if (fpTemp != NULL)
    {
        #ifdef DSP_IPP
        ippsCopy_32f(fpSrc, fpTemp, lCount);
        ippsSortAscend_32f_I(fpTemp, lCount);
        #else
        Copy(fpTemp, fpSrc, lCount);
        Sort(fpTemp, lCount);
        #endif
        fMedian = ((lCount % 2L) != 0L) ?
            fpTemp[lMax / 2L] :
            (0.5F * (fpTemp[lCount / 2L - 1L] + fpTemp[lCount / 2L]));
    }
    return fMedian;
}


double clDSPOp::Median (const double *dpSrc, long lCount)
{
    long lMax;
    double dMedian = 0.0;
    double *dpTemp;
    clDSPAlloc Temp(lCount * sizeof(double));

    dpTemp = Temp;
    lMax = lCount - 1L;
    if (dpTemp != NULL)
    {
        #ifdef DSP_IPP
        ippsCopy_64f(dpSrc, dpTemp, lCount);
        ippsSortAscend_64f_I(dpTemp, lCount);
        #else
        Copy(dpTemp, dpSrc, lCount);
        Sort(dpTemp, lCount);
        #endif
        dMedian = ((lCount % 2L) != 0L) ?
            dpTemp[lMax / 2L] :
            (0.5 * (dpTemp[lCount / 2L - 1L] + dpTemp[lCount / 2L]));
    }
    return dMedian;
}


void clDSPOp::Negate (float *fpVect, long lCount)
{
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        fpVect[lLoopCntr] = -(fpVect[lLoopCntr]);
    }
}


void clDSPOp::Negate (double *dpVect, long lCount)
{
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        dpVect[lLoopCntr] = -(dpVect[lLoopCntr]);
    }
}


void clDSPOp::Negate (float *fpDest, const float *fpSrc, long lCount)
{
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        fpDest[lLoopCntr] = -(fpSrc[lLoopCntr]);
    }
}


void clDSPOp::Negate (double *dpDest, const double *dpSrc, long lCount)
{
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        dpDest[lLoopCntr] = -(dpSrc[lLoopCntr]);
    }
}


void clDSPOp::Normalize (float *fpVect, long lCount)
{
    long lLoopCntr;
    float fMean;
    float fStdDev;

    StdDev(&fStdDev, &fMean, fpVect, lCount);
    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        fpVect[lLoopCntr] = (fpVect[lLoopCntr] - fMean) / fStdDev;
    }
}


void clDSPOp::Normalize (double *dpVect, long lCount)
{
    long lLoopCntr;
    double dMean;
    double dStdDev;

    StdDev(&dStdDev, &dMean, dpVect, lCount);
    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        dpVect[lLoopCntr] = (dpVect[lLoopCntr] - dMean) / dStdDev;
    }
}


void clDSPOp::Normalize (float *fpDest, const float *fpSrc, long lCount)
{
    long lLoopCntr;
    float fMean;
    float fStdDev;

    StdDev(&fStdDev, &fMean, fpSrc, lCount);
    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        fpDest[lLoopCntr] = (fpSrc[lLoopCntr] - fMean) / fStdDev;
    }
}


void clDSPOp::Normalize (double *dpDest, const double *dpSrc, long lCount)
{
    long lLoopCntr;
    double dMean;
    double dStdDev;

    StdDev(&dStdDev, &dMean, dpSrc, lCount);
    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        dpDest[lLoopCntr] = (dpSrc[lLoopCntr] - dMean) / dStdDev;
    }
}


float clDSPOp::Product (const float *fpSrc, long lCount)
{
    long lLoopCntr;
    float fProd = 1.0F;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        fProd *= fpSrc[lLoopCntr];
    }
    return fProd;
}


double clDSPOp::Product (const double *dpSrc, long lCount)
{
    long lLoopCntr;
    double dProd = 1.0;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        dProd *= dpSrc[lLoopCntr];
    }
    return dProd;
}


void clDSPOp::Reverse (float *fpVect, long lCount)
{
    #ifdef DSP_IPP
    ippsFlip_32f_I(fpVect, lCount);
    #else
    long lFwdIdx;
    long lRevIdx;
    long lMax;
    float fTemp1;
    float fTemp2;
    
    lMax = (lCount >> 1);
    lFwdIdx = 0;
    lRevIdx = lCount;
    while (lFwdIdx < lMax)
    {
        lRevIdx--;
        fTemp1 = fpVect[lFwdIdx];
        fTemp2 = fpVect[lRevIdx];
        fpVect[lFwdIdx] = fTemp2;
        fpVect[lRevIdx] = fTemp1;
        lFwdIdx++;
    }
    #endif
}


void clDSPOp::Reverse (double *dpVect, long lCount)
{
    #ifdef DSP_IPP
    ippsFlip_64f_I(dpVect, lCount);
    #else
    long lFwdIdx;
    long lRevIdx;
    long lMax;
    double dTemp1;
    double dTemp2;
    
    lMax = (lCount >> 1);
    lFwdIdx = 0;
    lRevIdx = lCount;
    while (lFwdIdx < lMax)
    {
        lRevIdx--;
        dTemp1 = dpVect[lFwdIdx];
        dTemp2 = dpVect[lRevIdx];
        dpVect[lFwdIdx] = dTemp2;
        dpVect[lRevIdx] = dTemp1;
        lFwdIdx++;
    }
    #endif
}


void clDSPOp::Reverse (stpSCplx spVect, long lCount)
{
    #ifdef DSP_IPP
    Ipp32fc *spTemp;
    
    spTemp = ippsMalloc_32fc(lCount);
    ippsConjFlip_32fc((Ipp32fc *) spVect, spTemp, lCount);
    ippsCopy_32fc(spTemp, (Ipp32fc *) spVect, lCount);
    ippsFree(spTemp);
    #else
    long lFwdIdx;
    long lRevIdx;
    long lMax;
    stSCplx sTemp1;
    stSCplx sTemp2;
    
    lMax = (lCount >> 1);
    lFwdIdx = 0;
    lRevIdx = lCount;
    while (lFwdIdx < lMax)
    {
        lRevIdx--;
        sTemp1.R = spVect[lFwdIdx].R;
        sTemp1.I = spVect[lFwdIdx].I;
        sTemp2.R = spVect[lRevIdx].R;
        sTemp2.I = spVect[lRevIdx].I;
        CplxConj(&sTemp1);
        CplxConj(&sTemp2);
        spVect[lFwdIdx].R = sTemp2.R;
        spVect[lFwdIdx].I = sTemp2.I;
        spVect[lRevIdx].R = sTemp1.R;
        spVect[lRevIdx].I = sTemp1.I;
        lFwdIdx++;
    }
    #endif
}


void clDSPOp::Reverse (stpDCplx spVect, long lCount)
{
    #ifdef DSP_IPP
    Ipp64fc *spTemp;
    
    spTemp = ippsMalloc_64fc(lCount);
    ippsConjFlip_64fc((Ipp64fc *) spVect, spTemp, lCount);
    ippsCopy_64fc(spTemp, (Ipp64fc *) spVect, lCount);
    ippsFree(spTemp);
    #else
    long lFwdIdx;
    long lRevIdx;
    long lMax;
    stDCplx sTemp1;
    stDCplx sTemp2;
    
    lMax = (lCount >> 1);
    lFwdIdx = 0;
    lRevIdx = lCount;
    while (lFwdIdx < lMax)
    {
        lRevIdx--;
        sTemp1.R = spVect[lFwdIdx].R;
        sTemp1.I = spVect[lFwdIdx].I;
        sTemp2.R = spVect[lRevIdx].R;
        sTemp2.I = spVect[lRevIdx].I;
        CplxConj(&sTemp1);
        CplxConj(&sTemp2);
        spVect[lFwdIdx].R = sTemp2.R;
        spVect[lFwdIdx].I = sTemp2.I;
        spVect[lRevIdx].R = sTemp1.R;
        spVect[lRevIdx].I = sTemp1.I;
        lFwdIdx++;
    }
    #endif
}


void clDSPOp::Reverse (float *fpDest, const float *fpSrc, long lCount)
{
    #ifdef DSP_IPP
    ippsFlip_32f(fpSrc, fpDest, lCount);
    #else
    long lLoopCntr;
    long lMax;

    lMax = lCount - 1L;
    for (lLoopCntr = 0L; lLoopCntr <= lMax; lLoopCntr++)
    {
        fpDest[lLoopCntr] = fpSrc[lMax - lLoopCntr];
    }
    #endif
}


void clDSPOp::Reverse (double *dpDest, const double *dpSrc, long lCount)
{
    #ifdef DSP_IPP
    ippsFlip_64f(dpSrc, dpDest, lCount);
    #else
    long lLoopCntr;
    long lMax;

    lMax = lCount - 1L;
    for (lLoopCntr = 0L; lLoopCntr <= lMax; lLoopCntr++)
    {
        dpDest[lLoopCntr] = dpSrc[lMax - lLoopCntr];
    }
    #endif
}


void clDSPOp::Reverse (stpSCplx spDest, const stpSCplx spSrc, long lCount)
{
    #ifdef DSP_IPP
    ippsConjFlip_32fc((Ipp32fc *) spSrc, (Ipp32fc *) spDest, lCount);
    #else
    long lLoopCntr;
    long lMax;
    
    lMax = lCount - 1L;
    for (lLoopCntr = 0L; lLoopCntr <= lMax; lLoopCntr++)
    {
        CplxConj(&spDest[lLoopCntr], &spSrc[lMax - lLoopCntr]);
    }
    #endif
}


void clDSPOp::Reverse (stpDCplx spDest, const stpDCplx spSrc, long lCount)
{
    #ifdef DSP_IPP
    ippsConjFlip_64fc((Ipp64fc *) spSrc, (Ipp64fc *) spDest, lCount);
    #else
    long lLoopCntr;
    long lMax;
    
    lMax = lCount - 1L;
    for (lLoopCntr = 0L; lLoopCntr <= lMax; lLoopCntr++)
    {
        CplxConj(&spDest[lLoopCntr], &spSrc[lMax - lLoopCntr]);
    }
    #endif
}


void clDSPOp::Scale (float *fpVect, long lCount)
{
    #ifdef DSP_IPP
    float fMin;
    float fMax;

    ippsMin_32f(fpVect, lCount, &fMin);
    ippsMax_32f(fpVect, lCount, &fMax);
    ippsNormalize_32f(fpVect, fpVect, lCount, fMin, (fMax - fMin) * 0.5F);
    ippsAddC_32f_I(-1.0F, fpVect, lCount);
    #else
    long lLoopCntr;
    float fMin;
    float fMax;
    float fScale;
    float fOffset;

    MinMax(&fMin, &fMax, fpVect, lCount);
    fScale = 2.0F / (fMax - fMin);
    fOffset = 1.0F - fMax * fScale;
    #ifdef DSP_X86
    if (bHave3DNow)
    {
        dsp_x86_3dnow_maf(fpVect, fScale, fOffset, lCount);
    }
    else if (iHaveSSE >= 2)
    {
        dsp_x86_sse_maf(fpVect, fScale, fOffset, lCount);
    }
    else
    #endif
    {
        for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
        {
            #ifndef _ISOC9X_SOURCE
            fpVect[lLoopCntr] = fpVect[lLoopCntr] * fScale + fOffset;
            #else
            fpVect[lLoopCntr] = 
                fmaf(fpVect[lLoopCntr], fScale, fOffset);
            #endif
        }
    }
    #endif
}


void clDSPOp::Scale (double *dpVect, long lCount)
{
    #ifdef DSP_IPP
    double dMin;
    double dMax;
    
    ippsMin_64f(dpVect, lCount, &dMin);
    ippsMax_64f(dpVect, lCount, &dMax);
    ippsNormalize_64f(dpVect, dpVect, lCount, dMin, (dMax - dMin) * 0.5);
    ippsAddC_64f_I(-1.0, dpVect, lCount);
    #else
    long lLoopCntr;
    double dMin;
    double dMax;
    double dScale;
    double dOffset;

    MinMax(&dMin, &dMax, dpVect, lCount);
    dScale = 2.0 / (dMax - dMin);
    dOffset = 1.0 - dMax * dScale;
    #ifdef DSP_X86
    if (iHaveSSE >= 2)
    {
        dsp_x86_sse_ma(dpVect, dScale, dOffset, lCount);
    }
    else
    #endif
    {
        for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
        {
            #ifndef _ISOC9X_SOURCE
            dpVect[lLoopCntr] = dpVect[lLoopCntr] * dScale + dOffset;
            #else
            dpVect[lLoopCntr] = 
                fma(dpVect[lLoopCntr], dScale, dOffset);
            #endif
        }
    }
    #endif
}


void clDSPOp::Scale (float *fpDest, const float *fpSrc, long lCount)
{
    #ifdef DSP_IPP
    float fMin;
    float fMax;

    ippsMin_32f(fpSrc, lCount, &fMin);
    ippsMax_32f(fpSrc, lCount, &fMax);
    ippsNormalize_32f(fpSrc, fpDest, lCount, fMin, (fMax - fMin) * 0.5F);
    ippsAddC_32f_I(1.0F, fpDest, lCount);
    #else
    long lLoopCntr;
    float fMin;
    float fMax;
    float fScale;
    float fOffset;

    MinMax(&fMin, &fMax, fpSrc, lCount);
    fScale = 2.0F / (fMax - fMin);
    fOffset = 1.0F - fMax * fScale;
    #ifdef DSP_X86
    if (bHave3DNow)
    {
        dsp_x86_3dnow_ma2f(fpDest, fpSrc, fScale, fOffset, lCount);
    }
    else if (iHaveSSE >= 2)
    {
        dsp_x86_sse_ma2f(fpDest, fpSrc, fScale, fOffset, lCount);
    }
    else
    #endif
    {
        for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
        {
            #ifndef _ISOC9X_SOURCE
            fpDest[lLoopCntr] = fpSrc[lLoopCntr] * fScale + fOffset;
            #else
            fpDest[lLoopCntr] = 
                fmaf(fpSrc[lLoopCntr], fScale, fOffset);
            #endif
        }
    }
    #endif
}


void clDSPOp::Scale (double *dpDest, const double *dpSrc, long lCount)
{
    #ifdef DSP_IPP
    double dMin;
    double dMax;

    ippsMin_64f(dpSrc, lCount, &dMin);
    ippsMax_64f(dpSrc, lCount, &dMax);
    ippsNormalize_64f(dpSrc, dpDest, lCount, dMin, (dMax - dMin) * 0.5);
    ippsAddC_64f_I(-1.0, dpDest, lCount);
    #else
    long lLoopCntr;
    double dMin;
    double dMax;
    double dScale;
    double dOffset;

    MinMax(&dMin, &dMax, dpSrc, lCount);
    dScale = 2.0 / (dMax - dMin);
    dOffset = 1.0 - dMax * dScale;
    #ifdef DSP_X86
    if (iHaveSSE >= 2)
    {
        dsp_x86_sse_ma2(dpDest, dpSrc, dScale, dOffset, lCount);
    }
    else
    #endif
    {
        for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
        {
            #ifndef _ISOC9X_SOURCE
            dpDest[lLoopCntr] = dpSrc[lLoopCntr] * dScale + dOffset;
            #else
            dpDest[lLoopCntr] = 
                fma(dpSrc[lLoopCntr], dScale, dOffset);
            #endif
        }
    }
    #endif
}


void clDSPOp::Scale01 (float *fpVect, long lCount)
{
    #ifdef DSP_IPP
    float fMin;
    float fMax;

    ippsMin_32f(fpVect, lCount, &fMin);
    ippsMax_32f(fpVect, lCount, &fMax);
    ippsNormalize_32f(fpVect, fpVect, lCount, fMin, fMax - fMin);
    #else
    long lLoopCntr;
    float fMin;
    float fMax;
    float fScale;
    float fOffset;

    MinMax(&fMin, &fMax, fpVect, lCount);
    fScale = 1.0F / (fMax - fMin);
    fOffset = -fMin * fScale;
    #ifdef DSP_X86
    if (bHave3DNow)
    {
        dsp_x86_3dnow_maf(fpVect, fScale, fOffset, lCount);
    }
    else if (iHaveSSE >= 2)
    {
        dsp_x86_sse_maf(fpVect, fScale, fOffset, lCount);
    }
    else
    #endif
    {
        for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
        {
            #ifndef _ISOC9X_SOURCE
            fpVect[lLoopCntr] = fpVect[lLoopCntr] * fScale + fOffset;
            #else
            fpVect[lLoopCntr] = fmaf(fpVect[lLoopCntr], fScale, fOffset);
            #endif
        }
    }
    #endif
}


void clDSPOp::Scale01 (double *dpVect, long lCount)
{
    #ifdef DSP_IPP
    double dMin;
    double dMax;

    ippsMin_64f(dpVect, lCount, &dMin);
    ippsMax_64f(dpVect, lCount, &dMax);
    ippsNormalize_64f(dpVect, dpVect, lCount, dMin, dMax - dMin);
    #else
    long lLoopCntr;
    double dMin;
    double dMax;
    double dScale;
    double dOffset;

    MinMax(&dMin, &dMax, dpVect, lCount);
    dScale = 1.0 / (dMax - dMin);
    dOffset = -dMin * dScale;
    #ifdef DSP_X86
    if (iHaveSSE >= 2)
    {
        dsp_x86_sse_ma(dpVect, dScale, dOffset, lCount);
    }
    else
    #endif
    {
        for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
        {
            #ifndef _ISOC9X_SOURCE
            dpVect[lLoopCntr] = dpVect[lLoopCntr] * dScale + dOffset;
            #else
            dpVect[lLoopCntr] = fma(dpVect[lLoopCntr], dScale, dOffset);
            #endif
        }
    }
    #endif
}


void clDSPOp::Scale01 (float *fpDest, const float *fpSrc, long lCount)
{
    #ifdef DSP_IPP
    float fMin;
    float fMax;

    ippsMin_32f(fpSrc, lCount, &fMin);
    ippsMax_32f(fpSrc, lCount, &fMax);
    ippsNormalize_32f(fpSrc, fpDest, lCount, fMin, fMax - fMin);
    #else
    long lLoopCntr;
    float fMin;
    float fMax;
    float fScale;
    float fOffset;

    MinMax(&fMin, &fMax, fpSrc, lCount);
    fScale = 1.0F / (fMax - fMin);
    fOffset = -fMin * fScale;
    #ifdef DSP_X86
    if (bHave3DNow)
    {
        dsp_x86_3dnow_ma2f(fpDest, fpSrc, fScale, fOffset, lCount);
    }
    else if (iHaveSSE >= 2)
    {
        dsp_x86_sse_ma2f(fpDest, fpSrc, fScale, fOffset, lCount);
    }
    else
    #endif
    {
        for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
        {
            #ifndef _ISOC9X_SOURCE
            fpDest[lLoopCntr] = fpSrc[lLoopCntr] * fScale + fOffset;
            #else
            fpDest[lLoopCntr] = fmaf(fpSrc[lLoopCntr], fScale, fOffset);
            #endif
        }
    }
    #endif
}


void clDSPOp::Scale01 (double *dpDest, const double *dpSrc, long lCount)
{
    #ifdef DSP_IPP
    double dMin;
    double dMax;

    ippsMin_64f(dpSrc, lCount, &dMin);
    ippsMax_64f(dpSrc, lCount, &dMax);
    ippsNormalize_64f(dpSrc, dpDest, lCount, dMin, dMax - dMin);
    #else
    long lLoopCntr;
    double dMin;
    double dMax;
    double dScale;
    double dOffset;

    MinMax(&dMin, &dMax, dpSrc, lCount);
    dScale = 1.0 / (dMax - dMin);
    dOffset = -dMin * dScale;
    #ifdef DSP_X86
    if (iHaveSSE >= 2)
    {
        dsp_x86_sse_ma2(dpDest, dpSrc, dScale, dOffset, lCount);
    }
    else
    #endif
    {
        for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
        {
            #ifndef _ISOC9X_SOURCE
            dpDest[lLoopCntr] = dpSrc[lLoopCntr] * dScale + dOffset;
            #else
            dpDest[lLoopCntr] = fma(dpSrc[lLoopCntr], dScale, dOffset);
            #endif
        }
    }
    #endif
}


void clDSPOp::Sort (float *fpVect, long lCount)
{
    #ifdef DSP_IPP
    ippsSortAscend_32f_I(fpVect, lCount);
    #else
    qsort(fpVect, (size_t) lCount, sizeof(float), FloatCompare);
    #endif
}


void clDSPOp::Sort (double *dpVect, long lCount)
{
    #ifdef DSP_IPP
    ippsSortAscend_64f_I(dpVect, lCount);
    #else
    qsort(dpVect, (size_t) lCount, sizeof(double), DoubleCompare);
    #endif
}


void clDSPOp::Sort (long *lpVect, long lCount)
{
    #ifdef DSP_IPP
    ippsSortAscend_32s_I((Ipp32s *) lpVect, lCount);
    #else
    qsort(lpVect, (size_t) lCount, sizeof(long), LongCompare);
    #endif
}


void clDSPOp::StdDev (float *fpStdDev, float *fpMean, const float *fpSrc,
    long lCount)
{
    #ifdef DSP_IPP
    ippsStdDev_32f(fpSrc, lCount, fpStdDev, ippAlgHintFast);
    ippsMean_32f(fpSrc, lCount, fpMean, ippAlgHintFast);
    #else
    long lLoopCntr;
    float fMean = 0.0F;
    float fTempVal;
    float fTemp = 0.0F;
    float fStdDev;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        fMean += fpSrc[lLoopCntr];
    }
    fMean /= (float) lCount;
    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        fTempVal = fpSrc[lLoopCntr] - fMean;
        #ifndef _ISOC9X_SOURCE
        fTemp += fTempVal * fTempVal;
        #else
        fTemp = fmaf(fTempVal, fTempVal, fTemp);
        #endif
    }
    fStdDev = sqrtf(fTemp / (float) lCount);
    *fpStdDev = fStdDev;
    *fpMean = fMean;
    #endif
}


void clDSPOp::StdDev (double *dpStdDev, double *dpMean, const double *dpSrc,
    long lCount)
{
    #ifdef DSP_IPP
    ippsStdDev_64f(dpSrc, lCount, dpStdDev);
    ippsStdDev_64f(dpSrc, lCount, dpMean);
    #else
    long lLoopCntr;
    double dMean = 0.0;
    double dTempVal;
    double dTemp = 0.0;
    double dStdDev;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        dMean += dpSrc[lLoopCntr];
    }
    dMean /= (double) lCount;
    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        dTempVal = dpSrc[lLoopCntr] - dMean;
        #ifndef _ISOC9X_SOURCE
        dTemp += dTempVal * dTempVal;
        #else
        dTemp = fma(dTempVal, dTempVal, dTemp);
        #endif
    }
    dStdDev = sqrt(dTemp / (double) lCount);
    *dpStdDev = dStdDev;
    *dpMean = dMean;
    #endif
}


float clDSPOp::Sum (const float *fpSrc, long lCount)
{
    #ifdef DSP_IPP
    float fSum;

    ippsSum_32f(fpSrc, lCount, &fSum, ippAlgHintFast);
    #else
    long lLoopCntr;
    float fSum = 0.0F;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        fSum += fpSrc[lLoopCntr];
    }
    #endif
    return fSum;
}


double clDSPOp::Sum (const double *dpSrc, long lCount)
{
    #ifdef DSP_IPP
    double dSum;

    ippsSum_64f(dpSrc, lCount, &dSum);
    #else
    long lLoopCntr;
    double dSum = 0.0;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        dSum += dpSrc[lLoopCntr];
    }
    #endif
    return dSum;
}


void clDSPOp::Square (float *fpVect, long lCount)
{
    #ifdef DSP_IPP
    ippsSqr_32f_I(fpVect, lCount);
    #else
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        fpVect[lLoopCntr] *= fpVect[lLoopCntr];
    }
    #endif
}


void clDSPOp::Square (double *dpVect, long lCount)
{
    #ifdef DSP_IPP
    ippsSqr_64f_I(dpVect, lCount);
    #else
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        dpVect[lLoopCntr] *= dpVect[lLoopCntr];
    }
    #endif
}


void clDSPOp::Square (float *fpDest, const float *fpSrc, long lCount)
{
    #ifdef DSP_IPP
    ippsSqr_32f(fpSrc, fpDest, lCount);
    #else
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        fpDest[lLoopCntr] = fpSrc[lLoopCntr] * fpSrc[lLoopCntr];
    }
    #endif
}


void clDSPOp::Square (double *dpDest, const double *dpSrc, long lCount)
{
    #ifdef DSP_IPP
    ippsSqr_64f(dpSrc, dpDest, lCount);
    #else
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        dpDest[lLoopCntr] = dpSrc[lLoopCntr] * dpSrc[lLoopCntr];
    }
    #endif
}


void clDSPOp::Convert (float *fpDest, const unsigned char *ucpSrc, long lCount)
{
    long lLoopCntr;
    float fMult;

    fMult = 1.0F / (float) UCHAR_MAX;
    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        fpDest[lLoopCntr] = (ucpSrc[lLoopCntr] * fMult - 0.5F) * 2.0F;
    }
}


void clDSPOp::Convert (float *fpDest, const signed short *sspSrc, long lCount,
    bool b12bit)
{
    #ifdef DSP_IPP
    ippsConvert_16s32f_Sfs(sspSrc, fpDest, lCount, 
        (b12bit) ? 11 : 15);
    #else
    long lLoopCntr;
    float fMult;

    #ifdef DSP_X86
    if (bHave3DNow)
    {
        dsp_x86_3dnow_i16tof(fpDest, sspSrc, lCount,
            ((b12bit) ? 0x1000 : SHRT_MAX));
    }
    else
    #endif
    {
        fMult = (b12bit) ? (1.0F / (float) 0x1000) : (1.0F / (float) SHRT_MAX);
        for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
        {
            fpDest[lLoopCntr] = sspSrc[lLoopCntr] * fMult;
        }
    }
    #endif
}


void clDSPOp::Convert (float *fpDest, const signed int *sipSrc, long lCount,
    bool b24bit)
{
    #ifdef DSP_IPP
    ippsConvert_32s32f_Sfs(sipSrc, fpDest, lCount, 31);
    ippsMulC_32f_I((float) 0x7fffffff / (float) 0x7fffff00, fpDest, lCount);
    #else
    long lLoopCntr;
    float fMult;

    #ifdef DSP_X86
    if (bHave3DNow)
    {
        dsp_x86_3dnow_i32tof(fpDest, sipSrc, lCount,
            ((b24bit) ? 0x7fffff00 : INT_MAX));
    }
    else
    #endif
    {
        fMult = (b24bit) ? 
            (1.0F / (float) 0x7fffff00) : 
            (1.0F / (float) INT_MAX);
        for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
        {
            fpDest[lLoopCntr] = sipSrc[lLoopCntr] * fMult;
        }
    }
    #endif
}


void clDSPOp::Convert (float *fpDest, const double *dpSrc, long lCount)
{
    #ifdef DSP_IPP
    ippsConvert_64f32f(dpSrc, fpDest, lCount);
    #else
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        fpDest[lLoopCntr] = (float) dpSrc[lLoopCntr];
    }
    #endif
}


void clDSPOp::Convert (double *dpDest, const unsigned char *ucpSrc, 
    long lCount)
{
    long lLoopCntr;
    double dMult;

    dMult = 1.0 / (double) UCHAR_MAX;
    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        dpDest[lLoopCntr] = (ucpSrc[lLoopCntr] * dMult - 0.5) * 2.0;
    }
}


void clDSPOp::Convert (double *dpDest, const signed short *sspSrc, long lCount,
    bool b12bit)
{
    #ifdef DSP_IPP
    ippsConvert_16s64f_Sfs(sspSrc, dpDest, lCount,
        (b12bit) ? 11 : 15);
    #else
    long lLoopCntr;
    double dMult;

    dMult = (b12bit) ? (1.0 / (double) 0x1000) : (1.0 / (double) SHRT_MAX);
    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        dpDest[lLoopCntr] = sspSrc[lLoopCntr] * dMult;
    }
    #endif
}


void clDSPOp::Convert (double *dpDest, const signed int *sipSrc, long lCount,
    bool b24bit)
{
    #ifdef DSP_IPP
    ippsConvert_32s64f_Sfs(sipSrc, dpDest, lCount, 31);
    ippsMulC_64f_I((double) 0x7fffffff / (double) 0x7fffff00, dpDest, lCount);
    #else
    long lLoopCntr;
    double dMult;

    dMult = (b24bit) ? 
        (1.0 / (double) 0x7fffff00) : 
        (1.0 / (double) INT_MAX);
    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        dpDest[lLoopCntr] = sipSrc[lLoopCntr] * dMult;
    }
    #endif
}


void clDSPOp::Convert (double *dpDest, const float *fpSrc, long lCount)
{
    #ifdef DSP_IPP
    ippsConvert_32f64f(fpSrc, dpDest, lCount);
    #else
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        dpDest[lLoopCntr] = (double) fpSrc[lLoopCntr];
    }
    #endif
}


void clDSPOp::Convert (unsigned char *ucpDest, const float *fpSrc, long lCount)
{
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        ucpDest[lLoopCntr] = (unsigned char)
            Round((fpSrc[lLoopCntr] + 1.0F) / 2.0F * (float) UCHAR_MAX);
    }
}


void clDSPOp::Convert (unsigned char *ucpDest, const double *dpSrc, 
    long lCount)
{
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        ucpDest[lLoopCntr] = (unsigned char)
            Round((dpSrc[lLoopCntr] + 1.0) / 2.0 * (double) UCHAR_MAX);
    }
}


void clDSPOp::Convert (signed short *sspDest, const float *fpSrc, long lCount,
    bool b12bit)
{
    #ifdef DSP_IPP
    ippsConvert_32f16s_Sfs(fpSrc, sspDest, lCount, ippRndNear,
        (b12bit) ? -11 : -15);
    #else
    long lLoopCntr;
    float fMult;

    fMult = (b12bit) ? (float) 0x1000 : (float) SHRT_MAX;
    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        sspDest[lLoopCntr] = (signed short)
            Round(fpSrc[lLoopCntr] * fMult);
    }
    #endif
}


void clDSPOp::Convert (signed short *sspDest, const double *dpSrc, long lCount,
    bool b12bit)
{
    long lLoopCntr;
    double dMult;

    dMult = (b12bit) ? (double) 0x1000 : (double) SHRT_MAX;
    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        sspDest[lLoopCntr] = (signed short) 
            Round(dpSrc[lLoopCntr] * dMult);
    }
}


void clDSPOp::Convert (signed int *sipDest, const float *fpSrc, long lCount,
    bool b24bit)
{
    #ifdef DSP_IPP
    ippsConvert_32f32s_Sfs(fpSrc, sipDest, lCount, ippRndNear, -31);
    #else
    long lLoopCntr;
    float fMult;

    fMult = (b24bit) ? ((float) 0x7fffff00) : ((float) INT_MAX);
    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        sipDest[lLoopCntr] = (signed int)
            Round(fpSrc[lLoopCntr] * fMult);
    }
    #endif
}


void clDSPOp::Convert (signed int *sipDest, const double *dpSrc, long lCount,
    bool b24bit)
{
    #ifdef DSP_IPP
    ippsConvert_64f32s_Sfs(dpSrc, sipDest, lCount, ippRndNear, -31);
    #else
    long lLoopCntr;
    double dMult;

    dMult = (b24bit) ? ((double) 0x7fffff00) : ((double) INT_MAX);
    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        sipDest[lLoopCntr] = (signed int)
            Round(dpSrc[lLoopCntr] * dMult);
    }
    #endif
}


void clDSPOp::CartToPolar (float *fpMagn, float *fpPhase, 
    const float *fpReal, const float *fpImag, long lCount)
{
    #ifdef DSP_IPP
    ippsCartToPolar_32f(fpReal, fpImag, fpMagn, fpPhase, lCount);
    #else
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        Cart2Polar(&fpMagn[lLoopCntr], &fpPhase[lLoopCntr], 
            fpReal[lLoopCntr], fpImag[lLoopCntr]);
    }
    #endif
}


void clDSPOp::CartToPolar (double *dpMagn, double *dpPhase, 
    const double *dpReal, const double *dpImag, long lCount)
{
    #ifdef DSP_IPP
    ippsCartToPolar_64f(dpReal, dpImag, dpMagn, dpPhase, lCount);
    #else
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        Cart2Polar(&dpMagn[lLoopCntr], &dpPhase[lLoopCntr],
            dpReal[lLoopCntr], dpImag[lLoopCntr]);
    }
    #endif
}


void clDSPOp::CartToPolar (float *fpMagn, float *fpPhase,
    const stpSCplx spCart, long lCount)
{
    #ifdef DSP_IPP
    ippsCartToPolar_32fc((Ipp32fc *) spCart, fpMagn, fpPhase, lCount);
    #else
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        Cart2Polar(&fpMagn[lLoopCntr], &fpPhase[lLoopCntr],
            &spCart[lLoopCntr]);
    }
    #endif
}


void clDSPOp::CartToPolar (double *dpMagn, double *dpPhase,
    const stpDCplx spCart, long lCount)
{
    #ifdef DSP_IPP
    ippsCartToPolar_64fc((Ipp64fc *) spCart, dpMagn, dpPhase, lCount);
    #else
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        Cart2Polar(&dpMagn[lLoopCntr], &dpPhase[lLoopCntr],
            &spCart[lLoopCntr]);
    }
    #endif
}


void clDSPOp::CartToPolar (stpSPolar spPolar, const stpSCplx spCart,
    long lCount)
{
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        Cart2Polar(&spPolar[lLoopCntr], &spCart[lLoopCntr]);
    }
}


void clDSPOp::CartToPolar (stpDPolar spPolar, const stpDCplx spCart,
    long lCount)
{
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        Cart2Polar(&spPolar[lLoopCntr], &spCart[lLoopCntr]);
    }
}


void clDSPOp::CartToPolar (utpSCoord upCoord, long lCount)
{
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        Cart2Polar(&upCoord[lLoopCntr]);
    }
}


void clDSPOp::CartToPolar (utpDCoord upCoord, long lCount)
{
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        Cart2Polar(&upCoord[lLoopCntr]);
    }
}


void clDSPOp::PolarToCart (float *fpReal, float *fpImag, 
    const float *fpMagn, const float *fpPhase, long lCount)
{
    #ifdef DSP_IPP
    ippsPolarToCart_32f(fpMagn, fpPhase, fpReal, fpImag, lCount);
    #else
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        Polar2Cart(&fpReal[lLoopCntr], &fpImag[lLoopCntr],
            fpMagn[lLoopCntr], fpPhase[lLoopCntr]);
    }
    #endif
}


void clDSPOp::PolarToCart (double *dpReal, double *dpImag,
    const double *dpMagn, const double *dpPhase, long lCount)
{
    #ifdef DSP_IPP
    ippsPolarToCart_64f(dpMagn, dpPhase, dpReal, dpImag, lCount);
    #else
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        Polar2Cart(&dpReal[lLoopCntr], &dpImag[lLoopCntr],
            dpMagn[lLoopCntr], dpPhase[lLoopCntr]);
    }
    #endif
}


void clDSPOp::PolarToCart (stpSCplx spCart, const float *fpMagn,
    const float *fpPhase, long lCount)
{
    #ifdef DSP_IPP
    ippsPolarToCart_32fc(fpMagn, fpPhase, (Ipp32fc *) spCart, lCount);
    #else
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        Polar2Cart(&spCart[lLoopCntr], fpMagn[lLoopCntr], fpPhase[lLoopCntr]);
    }
    #endif
}


void clDSPOp::PolarToCart (stpDCplx spCart, const double *dpMagn,
    const double *dpPhase, long lCount)
{
    #ifdef DSP_IPP
    ippsPolarToCart_64fc(dpMagn, dpPhase, (Ipp64fc *) spCart, lCount);
    #else
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        Polar2Cart(&spCart[lLoopCntr], dpMagn[lLoopCntr], dpPhase[lLoopCntr]);
    }
    #endif
}


void clDSPOp::PolarToCart (stpSCplx spCart, const stpSPolar spPolar,
    long lCount)
{
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        Polar2Cart(&spCart[lLoopCntr], &spPolar[lLoopCntr]);
    }
}


void clDSPOp::PolarToCart (stpDCplx spCart, const stpDPolar spPolar,
    long lCount)
{
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        Polar2Cart(&spCart[lLoopCntr], &spPolar[lLoopCntr]);
    }
}


void clDSPOp::PolarToCart (utpSCoord upCoord, long lCount)
{
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        Polar2Cart(&upCoord[lLoopCntr]);
    }
}


void clDSPOp::PolarToCart (utpDCoord upCoord, long lCount)
{
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        Polar2Cart(&upCoord[lLoopCntr]);
    }
}


float clDSPOp::CrossCorr (const float *fpSrc1, const float *fpSrc2, 
    long lCount)
{
    #ifdef DSP_X86
    if (bHave3DNow)
    {
        return dsp_x86_3dnow_crosscorrf(fpSrc1, fpSrc2, lCount);
    }
    #if ((__GNUC__ > 1) && (__GNUC__ < 4))
    else if (iHaveSSE >= 2)
    {
        return dsp_x86_sse_crosscorrf(fpSrc1, fpSrc2, lCount);
    }
    #endif
    else
    #endif
    #if (defined(DSP_X86_64) && (!defined(DSP_X86)))
    return dsp_x86_64_crosscorrf(fpSrc1, fpSrc2, lCount);
    #else
    {
        long lLoopCntr;
        #ifndef DSP_EXTPREC
        float fScale;
        float fNormFact;
        float fProdSum = 0.0F;
        float fSqSum1 = 0.0F;
        float fSqSum2 = 0.0F;
        #else
        double fScale;
        double fNormFact;
        double fProdSum = 0.0;
        double fSqSum1 = 0.0;
        double fSqSum2 = 0.0;
        #endif

        for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
        {
            #ifndef _ISOC9X_SOURCE
            fProdSum += fpSrc1[lLoopCntr] * fpSrc2[lLoopCntr];
            fSqSum1 += fpSrc1[lLoopCntr] * fpSrc1[lLoopCntr];
            fSqSum2 += fpSrc2[lLoopCntr] * fpSrc2[lLoopCntr];
            #else
            #ifndef DSP_EXTPREC
            fProdSum = fmaf(fpSrc1[lLoopCntr], fpSrc2[lLoopCntr], fProdSum);
            fSqSum1 = fmaf(fpSrc1[lLoopCntr], fpSrc1[lLoopCntr], fSqSum1);
            fSqSum2 = fmaf(fpSrc2[lLoopCntr], fpSrc2[lLoopCntr], fSqSum2);
            #else
            fProdSum = fma(fpSrc1[lLoopCntr], fpSrc2[lLoopCntr], fProdSum);
            fSqSum1 = fma(fpSrc1[lLoopCntr], fpSrc1[lLoopCntr], fSqSum1);
            fSqSum2 = fma(fpSrc2[lLoopCntr], fpSrc2[lLoopCntr], fSqSum2);
            #endif
            #endif
        }
        #ifndef DSP_EXTPREC
        fScale = 1.0F / lCount;
        fNormFact = sqrtf(fSqSum1 * fSqSum2) * fScale;
        #else
        fScale = 1.0 / lCount;
        fNormFact = sqrt(fSqSum1 * fSqSum2) * fScale;
        #endif
        return ((fProdSum * fScale) / fNormFact);
    }
    #endif
}


double clDSPOp::CrossCorr (const double *dpSrc1, const double *dpSrc2, 
    long lCount)
{
    #ifdef DSP_X86
    #if ((__GNUC__ > 1) && (__GNUC__ < 4))
    if (iHaveSSE >= 2)
    {
        return dsp_x86_sse_crosscorr(dpSrc1, dpSrc2, lCount);
    }
    else
    #endif
    #endif
    #if (defined(DSP_X86_64) && (!defined(DSP_X86)))
    return dsp_x86_64_crosscorr(dpSrc1, dpSrc2, lCount);
    #else
    {
        long lLoopCntr;
        #ifndef DSP_EXTPREC
        double dScale;
        double dNormFact;
        double dProdSum = 0.0;
        double dSqSum1 = 0.0;
        double dSqSum2 = 0.0;
        #else
        long double dScale;
        long double dNormFact;
        long double dProdSum = 0.0L;
        long double dSqSum1 = 0.0L;
        long double dSqSum2 = 0.0L;
        #endif

        for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
        {
            #ifndef _ISOC9X_SOURCE
            dProdSum += dpSrc1[lLoopCntr] * dpSrc2[lLoopCntr];
            dSqSum1 += dpSrc1[lLoopCntr] * dpSrc1[lLoopCntr];
            dSqSum2 += dpSrc2[lLoopCntr] * dpSrc2[lLoopCntr];
            #else
            #ifndef DSP_EXTPREC
            dProdSum = fma(dpSrc1[lLoopCntr], dpSrc2[lLoopCntr], dProdSum);
            dSqSum1 = fma(dpSrc1[lLoopCntr], dpSrc1[lLoopCntr], dSqSum1);
            dSqSum2 = fma(dpSrc2[lLoopCntr], dpSrc2[lLoopCntr], dSqSum2);
            #else
            dProdSum = fmal(dpSrc1[lLoopCntr], dpSrc2[lLoopCntr], dProdSum);
            dSqSum1 = fmal(dpSrc1[lLoopCntr], dpSrc1[lLoopCntr], dSqSum1);
            dSqSum2 = fmal(dpSrc2[lLoopCntr], dpSrc2[lLoopCntr], dSqSum2);
            #endif
            #endif
        }
        #ifndef DSP_EXTPREC
        dScale = 1.0 / lCount;
        dNormFact = sqrt(dSqSum1 * dSqSum2) * dScale;
        #else
        dScale = 1.0L / lCount;
        dNormFact = sqrtl(dSqSum1 * dSqSum2) * dScale;
        #endif
        return ((dProdSum * dScale) / dNormFact);
    }
    #endif
}


float clDSPOp::DelCrossCorr (const float *fpSrc1, const float *fpSrc2,
    long lDelay, long lCount)
{
    long lMax;

    lMax = lCount - lDelay;
    #ifdef DSP_X86
    if (bHave3DNow)
    {
        return dsp_x86_3dnow_crosscorrf(fpSrc1, &fpSrc2[lDelay], lMax);
    }
    #if ((__GNUC__ > 1) && (__GNUC__ < 4))
    else  if (iHaveSSE >= 2)
    {
        return dsp_x86_sse_crosscorrf(fpSrc1, &fpSrc2[lDelay], lMax);
    }
    #endif
    else
    #endif
    #if (defined(DSP_X86_64) && (!defined(DSP_X86)))
    return dsp_x86_64_crosscorrf(fpSrc1, &fpSrc2[lDelay], lMax);
    #else
    {
        long lLoopCntr;
        long lIdx2;
        float fScale;
        float fNormFact;
        float fProdSum = 0.0F;
        float fSqSum1 = 0.0F;
        float fSqSum2 = 0.0F;

        for (lLoopCntr = 0L; lLoopCntr < lMax; lLoopCntr++)
        {
            lIdx2 = lLoopCntr + lDelay;
            #ifndef _ISOC9X_SOURCE
            fProdSum += fpSrc1[lLoopCntr] * fpSrc2[lIdx2];
            fSqSum1 += fpSrc1[lLoopCntr] * fpSrc1[lLoopCntr];
            fSqSum2 += fpSrc2[lIdx2] * fpSrc2[lIdx2];
            #else
            fProdSum = fmaf(fpSrc1[lLoopCntr], fpSrc2[lIdx2], fProdSum);
            fSqSum1 = fmaf(fpSrc1[lLoopCntr], fpSrc1[lLoopCntr], fSqSum1);
            fSqSum2 = fmaf(fpSrc2[lIdx2], fpSrc2[lIdx2], fSqSum2);
            #endif
        }
        fScale = 1.0F / (float) lMax;
        fNormFact = sqrtf(fSqSum1 * fSqSum2) * fScale;
        return ((fProdSum * fScale) / fNormFact);
    }
    #endif
}


double clDSPOp::DelCrossCorr (const double *dpSrc1, const double *dpSrc2,
    long lDelay, long lCount)
{
    long lMax;

    lMax = lCount - lDelay;
    #ifdef DSP_X86
    #if ((__GNUC__ > 1) && (__GNUC__ < 4))
    if (iHaveSSE >= 2)
    {
        return dsp_x86_sse_crosscorr(dpSrc1, &dpSrc2[lDelay], lMax);
    }
    else
    #endif
    #endif
    #if (defined(DSP_X86_64) && (!defined(DSP_X86)))
    return dsp_x86_64_crosscorr(dpSrc1, &dpSrc2[lDelay], lMax);
    #else
    {
        long lLoopCntr;
        long lIdx2;
        double dScale;
        double dNormFact;
        double dProdSum = 0.0;
        double dSqSum1 = 0.0;
        double dSqSum2 = 0.0;

        for (lLoopCntr = 0L; lLoopCntr < lMax; lLoopCntr++)
        {
            lIdx2 = lLoopCntr + lDelay;
            #ifndef _ISOC9X_SOURCE
            dProdSum += dpSrc1[lLoopCntr] * dpSrc2[lIdx2];
            dSqSum1 += dpSrc1[lLoopCntr] * dpSrc1[lLoopCntr];
            dSqSum2 += dpSrc2[lIdx2] * dpSrc2[lIdx2];
            #else
            dProdSum = fma(dpSrc1[lLoopCntr], dpSrc2[lIdx2], dProdSum);
            dSqSum1 = fma(dpSrc1[lLoopCntr], dpSrc1[lLoopCntr], dSqSum1);
            dSqSum2 = fma(dpSrc2[lIdx2], dpSrc2[lIdx2], dSqSum2);
            #endif
        }
        dScale = 1.0 / (double) lMax;
        dNormFact = sqrt(dSqSum1 * dSqSum2) * dScale;
        return ((dProdSum * dScale) / dNormFact);
    }
    #endif
}


void clDSPOp::DelCrossCorr (float *fpDest, const float *fpSrc1, 
    const float *fpSrc2, long lCount, const long *lpDelay, long lDelayCount)
{
    long lLoopDelay;
    long lLoopCorr;
    long lDelay;
    long lMax;
    long lIdx2;
    float fScale;
    float fNormFact;
    float fProdSum;
    float fSqSum1;
    float fSqSum2;

    for (lLoopDelay = 0L; lLoopDelay < lDelayCount; lLoopDelay++)
    {
        fProdSum = 0.0F;
        fSqSum1 = 0.0F;
        fSqSum2 = 0.0F;
        lDelay = lpDelay[lLoopDelay];
        lMax = lCount - lDelay;
        for (lLoopCorr = 0L; lLoopCorr < lMax; lLoopCorr++)
        {
            lIdx2 = lLoopCorr + lDelay;
            #ifndef _ISOC9X_SOURCE
            fProdSum += fpSrc1[lLoopCorr] * fpSrc2[lIdx2];
            fSqSum1 += fpSrc1[lLoopCorr] * fpSrc1[lLoopCorr];
            fSqSum2 += fpSrc2[lIdx2] * fpSrc2[lIdx2];
            #else
            fProdSum = fmaf(fpSrc1[lLoopCorr], fpSrc2[lIdx2], fProdSum);
            fSqSum1 = fmaf(fpSrc1[lLoopCorr], fpSrc1[lLoopCorr], fSqSum1);
            fSqSum2 = fmaf(fpSrc2[lIdx2], fpSrc2[lIdx2], fSqSum2);
            #endif
        }
        fScale = 1.0F / (float) lMax;
        fNormFact = sqrtf(fSqSum1 * fSqSum2) * fScale;
        fpDest[lLoopDelay] = (fProdSum * fScale) / fNormFact;
    }
}


void clDSPOp::DelCrossCorr (double *dpDest, const double *dpSrc1,
    const double *dpSrc2, long lCount, const long *lpDelay, long lDelayCount)
{
    long lLoopDelay;
    long lLoopCorr;
    long lDelay;
    long lMax;
    long lIdx2;
    double dScale;
    double dNormFact;
    double dProdSum;
    double dSqSum1;
    double dSqSum2;

    for (lLoopDelay = 0L; lLoopDelay < lDelayCount; lLoopDelay++)
    {
        dProdSum = 0.0;
        dSqSum1 = 0.0;
        dSqSum2 = 0.0;
        lDelay = lpDelay[lLoopDelay];
        lMax = lCount - lDelay;
        for (lLoopCorr = 0L; lLoopCorr < lMax; lLoopCorr++)
        {
            lIdx2 = lLoopCorr + lDelay;
            #ifndef _ISOC9X_SOURCE
            dProdSum += dpSrc1[lLoopCorr] * dpSrc2[lIdx2];
            dSqSum1 += dpSrc1[lLoopCorr] * dpSrc1[lLoopCorr];
            dSqSum2 += dpSrc2[lIdx2] * dpSrc2[lIdx2];
            #else
            dProdSum = fma(dpSrc1[lLoopCorr], dpSrc2[lIdx2], dProdSum);
            dSqSum1 = fma(dpSrc1[lLoopCorr], dpSrc1[lLoopCorr], dSqSum1);
            dSqSum2 = fma(dpSrc2[lIdx2], dpSrc2[lIdx2], dSqSum2);
            #endif
        }
        dScale = 1.0 / (double) lMax;
        dNormFact = sqrt(dSqSum1 * dSqSum2) * dScale;
        dpDest[lLoopDelay] = (dProdSum * dScale) / dNormFact;
    }
}


float clDSPOp::Energy (const float *fpSrc, long lCount)
{
    long lLoopCntr;
    float fEnergy = 0.0F;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        #ifndef _ISOC9X_SOURCE
        fEnergy += fpSrc[lLoopCntr] * fpSrc[lLoopCntr];
        #else
        fEnergy = fmaf(fpSrc[lLoopCntr], fpSrc[lLoopCntr], fEnergy);
        #endif
    }
    return fEnergy;
}


double clDSPOp::Energy (const double *dpSrc, long lCount)
{
    long lLoopCntr;
    double dEnergy = 0.0;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        #ifndef _ISOC9X_SOURCE
        dEnergy += dpSrc[lLoopCntr] * dpSrc[lLoopCntr];
        #else
        dEnergy = fma(dpSrc[lLoopCntr], dpSrc[lLoopCntr], dEnergy);
        #endif
    }
    return dEnergy;
}


void clDSPOp::Magnitude (float *fpMagn, const stpSCplx spCplx, long lCount)
{
    #ifdef DSP_IPP
    ippsMagnitude_32fc((Ipp32fc *) spCplx, fpMagn, lCount);
    #else
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        #ifndef _ISOC9X_SOURCE
        fpMagn[lLoopCntr] = sqrtf(
            spCplx[lLoopCntr].R * spCplx[lLoopCntr].R +
            spCplx[lLoopCntr].I * spCplx[lLoopCntr].I);
        #else
        fpMagn[lLoopCntr] = 
            hypotf(spCplx[lLoopCntr].R, spCplx[lLoopCntr].I);
        #endif
    }
    #endif
}


void clDSPOp::Magnitude (double *dpMagn, const stpDCplx spCplx, long lCount)
{
    #ifdef DSP_IPP
    ippsMagnitude_64fc((Ipp64fc *) spCplx, dpMagn, lCount);
    #else
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        #ifndef _ISOC9X_SOURCE
        dpMagn[lLoopCntr] = sqrt(
            spCplx[lLoopCntr].R * spCplx[lLoopCntr].R +
            spCplx[lLoopCntr].I * spCplx[lLoopCntr].I);
        #else
        dpMagn[lLoopCntr] =
            hypot(spCplx[lLoopCntr].R, spCplx[lLoopCntr].I);
        #endif
    }
    #endif
}


void clDSPOp::Power (float *fpPower, const stpSCplx spCplx, long lCount)
{
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        #ifndef _ISOC9X_SOURCE
        fpPower[lLoopCntr] = 20.0F * log10f(sqrtf(
            spCplx[lLoopCntr].R * spCplx[lLoopCntr].R +
            spCplx[lLoopCntr].I * spCplx[lLoopCntr].I));
        #else
        fpPower[lLoopCntr] = 20.0F * log10f(
            hypotf(spCplx[lLoopCntr].R, spCplx[lLoopCntr].I));
        #endif
    }
}


void clDSPOp::Power (double *dpPower, const stpDCplx spCplx, long lCount)
{
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        #ifndef _ISOC9X_SOURCE
        dpPower[lLoopCntr] = 20.0 * log10(sqrt(
            spCplx[lLoopCntr].R * spCplx[lLoopCntr].R +
            spCplx[lLoopCntr].I * spCplx[lLoopCntr].I));
        #else
        dpPower[lLoopCntr] = 20.0 * log10(
            hypot(spCplx[lLoopCntr].R, spCplx[lLoopCntr].I));
        #endif
    }
}


void clDSPOp::Phase (float *fpPhase, const stpSCplx spCplx, long lCount)
{
    #ifdef DSP_IPP
    ippsPhase_32fc((Ipp32fc *) spCplx, fpPhase, lCount);
    #else
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        fpPhase[lLoopCntr] =
            atan2f(spCplx[lLoopCntr].I, spCplx[lLoopCntr].R);
    }
    #endif
}


void clDSPOp::Phase (double *dpPhase, const stpDCplx spCplx, long lCount)
{
    #ifdef DSP_IPP
    ippsPhase_64fc((Ipp64fc *) spCplx, dpPhase, lCount);
    #else
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        dpPhase[lLoopCntr] = 
            atan2(spCplx[lLoopCntr].I, spCplx[lLoopCntr].R);
    }
    #endif
}



void clDSPOp::PowerPhase (float *fpPower, float *fpPhase, 
    const stpSCplx spCplx, long lCount)
{
    long lLoopCntr;
    float fReal;
    float fImag;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        fReal = spCplx[lLoopCntr].R;
        fImag = spCplx[lLoopCntr].I;
        #ifndef _ISOC9X_SOURCE
        fpPower[lLoopCntr] = 
            20.0F * log10f(sqrtf(fReal * fReal + fImag * fImag));
        #else
        fpPower[lLoopCntr] =
            20.0F * log10f(hypotf(fReal, fImag));
        #endif
        fpPhase[lLoopCntr] = atan2f(fImag, fReal);
    }
}


void clDSPOp::PowerPhase (double *dpPower, double *dpPhase,
    const stpDCplx spCplx, long lCount)
{
    long lLoopCntr;
    double dReal;
    double dImag;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        dReal = spCplx[lLoopCntr].R;
        dImag = spCplx[lLoopCntr].I;
        #ifndef _ISOC9X_SOURCE
        dpPower[lLoopCntr] = 
            20.0 * log10(sqrt(dReal * dReal + dImag * dImag));
        #else
        dpPower[lLoopCntr] =
            20.0 * log10(hypot(dReal, dImag));
        #endif
        dpPhase[lLoopCntr] = atan2(dImag, dReal);
    }
}


void clDSPOp::Decimate (float *fpDest, const float *fpSrc, long lFactor,
    long lCount)
{
    #ifdef DSP_IPP
    int iDestCount = lCount / lFactor;
    int iPhase = 0;

    ippsSampleDown_32f(fpSrc, lCount, fpDest, &iDestCount, lFactor, &iPhase);
    #else
    long lLoopCntr;
    long lMax;

    lMax = lCount / lFactor;
    for (lLoopCntr = 0L; lLoopCntr < lMax; lLoopCntr++)
    {
        fpDest[lLoopCntr] = fpSrc[lLoopCntr * lFactor];
    }
    #endif
}


void clDSPOp::Decimate (double *dpDest, const double *dpSrc, long lFactor,
    long lCount)
{
    #ifdef DSP_IPP
    int iDestCount = lCount / lFactor;
    int iPhase = 0;

    ippsSampleDown_64f(dpSrc, lCount, dpDest, &iDestCount, lFactor, &iPhase);
    #else
    long lLoopCntr;
    long lMax;

    lMax = lCount / lFactor;
    for (lLoopCntr = 0L; lLoopCntr < lMax; lLoopCntr++)
    {
        dpDest[lLoopCntr] = dpSrc[lLoopCntr * lFactor];
    }
    #endif
}


void clDSPOp::DecimateAvg (float *fpDest, const float *fpSrc, long lFactor,
    long lCount)
{
    long lLoopDest;
    long lLoopAvg;
    long lMax;
    float fAvg;

    lMax = lCount / lFactor;
    for (lLoopDest = 0L; lLoopDest < lMax; lLoopDest++)
    {
        fAvg = 0.0F;
        for (lLoopAvg = 0L; lLoopAvg < lFactor; lLoopAvg++)
        {
            fAvg += fpSrc[lLoopDest * lFactor + lLoopAvg];
        }
        fpDest[lLoopDest] = fAvg / (float) lFactor;
    }
}


void clDSPOp::DecimateAvg (double *dpDest, const double *dpSrc, long lFactor,
    long lCount)
{
    long lLoopDest;
    long lLoopAvg;
    long lMax;
    double dAvg;

    lMax = lCount / lFactor;
    for (lLoopDest = 0L; lLoopDest < lMax; lLoopDest++)
    {
        dAvg = 0.0;
        for (lLoopAvg = 0L; lLoopAvg < lFactor; lLoopAvg++)
        {
            dAvg += dpSrc[lLoopDest * lFactor + lLoopAvg];
        }
        dpDest[lLoopDest] = dAvg / (double) lFactor;
    }
}


void clDSPOp::Interpolate (float *fpDest, const float *fpSrc, long lFactor,
    long lCount)
{
    #ifdef DSP_IPP
    int iDestCount = lCount * lFactor;
    int iPhase = 0;

    ippsSampleUp_32f(fpSrc, lCount, fpDest, &iDestCount, lFactor, &iPhase);
    #else
    long lLoopSrc;
    long lLoopInt;

    for (lLoopSrc = 0L; lLoopSrc < lCount; lLoopSrc++)
    {
        fpDest[lLoopSrc * lFactor] = fpSrc[lLoopSrc];
        for (lLoopInt = 1L; lLoopInt < lFactor; lLoopInt++)
        {
            fpDest[lLoopSrc * lFactor + lLoopInt] = 0.0F;
        }
    }
    #endif
}


void clDSPOp::Interpolate (double *dpDest, const double *dpSrc, long lFactor,
    long lCount)
{
    #ifdef DSP_IPP
    int iDestCount = lCount * lFactor;
    int iPhase = 0;

    ippsSampleUp_64f(dpSrc, lCount, dpDest, &iDestCount, lFactor, &iPhase);
    #else
    long lLoopSrc;
    long lLoopInt;

    for (lLoopSrc = 0L; lLoopSrc < lCount; lLoopSrc++)
    {
        dpDest[lLoopSrc * lFactor] = dpSrc[lLoopSrc];
        for (lLoopInt = 1L; lLoopInt < lFactor; lLoopInt++)
        {
            dpDest[lLoopSrc * lFactor + lLoopInt] = 0.0;
        }
    }
    #endif
}


void clDSPOp::InterpolateAvg (float *fpDest, const float *fpSrc, 
    long lFactor, long lCount)
{
    long lSrcCntr;
    long lIntCntr;
    long lX;
    long lX0;
    long lX1;
    float fF0;
    float fF1;
    float fL0;
    float fL1;
    float fP1;

    for (lSrcCntr = 0L; lSrcCntr < lCount; lSrcCntr++)
    {
        lX0 = lSrcCntr * lFactor;
        lX1 = (lSrcCntr + 1L) * lFactor;
        if (lSrcCntr != (lCount - 1L))
        {
            fF0 = fpSrc[lSrcCntr];
            fF1 = fpSrc[lSrcCntr + 1L];
        }
        else
        {
            fF0 = fpSrc[lSrcCntr];
            fF1 = fpSrc[lSrcCntr] - (fpSrc[lSrcCntr - 1L] - fpSrc[lSrcCntr]);
        }
        fpDest[lSrcCntr * lFactor] = fF0;
        for (lIntCntr = 1L; lIntCntr < lFactor; lIntCntr++)
        {
            lX = lSrcCntr * lFactor + lIntCntr;
            fL0 = (float) (lX - lX1) / (float) (lX0 - lX1);
            fL1 = (float) (lX - lX0) / (float) (lX1 - lX0);
            fP1 = fL0 * fF0 + fL1 * fF1;
            fpDest[lX] = fP1;
        }
    }
}


void clDSPOp::InterpolateAvg (double *dpDest, const double *dpSrc,
    long lFactor, long lCount)
{
    long lSrcCntr;
    long lIntCntr;
    long lX;
    long lX0;
    long lX1;
    double dF0;
    double dF1;
    double dL0;
    double dL1;
    double dP1;

    for (lSrcCntr = 0L; lSrcCntr < lCount; lSrcCntr++)
    {
        lX0 = lSrcCntr * lFactor;
        lX1 = (lSrcCntr + 1L) * lFactor;
        if (lSrcCntr != (lCount - 1L))
        {
            dF0 = dpSrc[lSrcCntr];
            dF1 = dpSrc[lSrcCntr + 1L];
        }
        else
        {
            dF0 = dpSrc[lSrcCntr];
            dF1 = dpSrc[lSrcCntr] - (dpSrc[lSrcCntr - 1L] - dpSrc[lSrcCntr]);
        }
        dpDest[lSrcCntr * lFactor] = dF0;
        for (lIntCntr = 1L; lIntCntr < lFactor; lIntCntr++)
        {
            lX = lSrcCntr * lFactor + lIntCntr;
            dL0 = (double) (lX - lX1) / (double) (lX0 - lX1);
            dL1 = (double) (lX - lX0) / (double) (lX1 - lX0);
            dP1 = dL0 * dF0 + dL1 * dF1;
            dpDest[lX] = dP1;
        }
    }
}


void clDSPOp::Resample (float *fpDest, long lDestCount, 
    const float *fpSrc, long lSrcCount)
{
    long lDestCntr;
    long lSrcIdx;
    float fScale;

    fScale = (float) lSrcCount / (float) lDestCount;
    for (lDestCntr = 0; lDestCntr < lDestCount; lDestCntr++)
    {
        lSrcIdx = (long) ((float) lDestCntr * fScale + 0.5F);
        if (unlikely(lSrcIdx >= lSrcCount))
            lSrcIdx = lSrcCount - 1;
        fpDest[lDestCntr] = fpSrc[lSrcIdx];
    }
}


void clDSPOp::Resample (double *dpDest, long lDestCount,
    const double *dpSrc, long lSrcCount)
{
    long lDestCntr;
    long lSrcIdx;
    double dScale;

    dScale = (double) lSrcCount / (double) lDestCount;
    for (lDestCntr = 0; lDestCntr < lDestCount; lDestCntr++)
    {
        lSrcIdx = (long) ((double) lDestCntr * dScale + 0.5);
        if (unlikely(lSrcIdx >= lSrcCount))
            lSrcIdx = lSrcCount - 1;
        dpDest[lDestCntr] = dpSrc[lSrcIdx];
    }
}


void clDSPOp::ResampleAvg (float *fpDest, long lDestCount,
    const float *fpSrc, long lSrcCount)
{
    long lDestCntr;
    long lSrcCntr;
    long lSrcIdx1;
    long lSrcIdx2;
    long lDist;
    float fScale;
    float fValue;

    fScale = (float) lSrcCount / (float) lDestCount;
    if (lDestCount < lSrcCount)
    {
        for (lDestCntr = 0; lDestCntr < lDestCount; lDestCntr++)
        {
            lSrcIdx1 = (long) ((float) lDestCntr * fScale + 0.5F);
            lSrcIdx2 = (long) ((float) (lDestCntr + 1) * fScale + 0.5F);
            if (unlikely(lSrcIdx1 >= lSrcCount))
                lSrcIdx1 = lSrcCount - 1;
            if (unlikely(lSrcIdx2 >= lSrcCount))
                lSrcIdx2 = lSrcCount - 1;
            lDist = lSrcIdx2 - lSrcIdx1;
            if (lDist > 0)
            {
                fValue = 0.0f;
                for (lSrcCntr = lSrcIdx1; lSrcCntr < lSrcIdx2; lSrcCntr++)
                    fValue += fpSrc[lSrcCntr];
                fValue /= (float) lDist;
            }
            else fValue = fpSrc[lSrcIdx1];
            fpDest[lDestCntr] = fValue;
        }
    }
    else if (lDestCount > lSrcCount)
    {
        fpDest[0] = fpSrc[0];
        for (lDestCntr = 1; lDestCntr < lDestCount; lDestCntr++)
        {
            lSrcIdx1 = (long) ((float) lDestCntr * fScale + 0.5F);
            lSrcIdx2 = (long) ((float) (lDestCntr + 1) * fScale + 0.5F);
            if (unlikely(lSrcIdx1 >= lSrcCount))
                lSrcIdx1 = lSrcCount - 1;
            if (unlikely(lSrcIdx2 >= lSrcCount))
                lSrcIdx2 = lSrcCount - 1;
            fpDest[lDestCntr] = fpSrc[lSrcIdx1] + 
                ((fpSrc[lSrcIdx2] - fpSrc[lSrcIdx1]) * fScale);
        }
    }
    else Copy(fpDest, fpSrc, lDestCount);
}


void clDSPOp::ResampleAvg (double *dpDest, long lDestCount,
    const double *dpSrc, long lSrcCount)
{
    long lDestCntr;
    long lSrcCntr;
    long lSrcIdx1;
    long lSrcIdx2;
    long lDist;
    double dScale;
    double dValue;

    dScale = (double) lSrcCount / (double) lDestCount;
    if (lDestCount < lSrcCount)
    {
        for (lDestCntr = 0; lDestCntr < lDestCount; lDestCntr++)
        {
            lSrcIdx1 = (long) ((double) lDestCntr * dScale + 0.5);
            lSrcIdx2 = (long) ((double) (lDestCntr + 1) * dScale + 0.5);
            if (unlikely(lSrcIdx1 >= lSrcCount))
                lSrcIdx1 = lSrcCount - 1;
            if (unlikely(lSrcIdx2 >= lSrcCount))
                lSrcIdx2 = lSrcCount - 1;
            lDist = lSrcIdx2 - lSrcIdx1;
            if (lDist > 0)
            {
                dValue = 0.0;
                for (lSrcCntr = lSrcIdx1; lSrcCntr < lSrcIdx2; lSrcCntr++)
                    dValue += dpSrc[lSrcCntr];
                dValue /= (double) lDist;
            }
            else dValue = dpSrc[lSrcIdx1];
            dpDest[lDestCntr] = dValue;
        }
    }
    else if (lDestCount > lSrcCount)
    {
        dpDest[0] = dpSrc[0];
        for (lDestCntr = 1; lDestCntr < lDestCount; lDestCntr++)
        {
            lSrcIdx1 = (long) ((double) lDestCntr * dScale + 0.5);
            lSrcIdx2 = (long) ((double) (lDestCntr + 1) * dScale + 0.5);
            if (unlikely(lSrcIdx1 >= lSrcCount))
                lSrcIdx1 = lSrcCount - 1;
            if (unlikely(lSrcIdx2 >= lSrcCount))
                lSrcIdx2 = lSrcCount - 1;
            dpDest[lDestCntr] = dpSrc[lSrcIdx1] +
                ((dpSrc[lSrcIdx2] - dpSrc[lSrcIdx1]) * dScale);
        }
    }
    else Copy(dpDest, dpSrc, lDestCount);
}


float clDSPOp::RMS (const float *fpSrc, long lCount)
{
    long lLoopCntr;
    float fSqSum = 0.0F;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        #ifndef _ISOC9X_SOURCE
        fSqSum += fpSrc[lLoopCntr] * fpSrc[lLoopCntr];
        #else
        fSqSum = fmaf(fpSrc[lLoopCntr], fpSrc[lLoopCntr], fSqSum);
        #endif
    }
    return ((float) sqrtf(fSqSum / (float) lCount));
}


double clDSPOp::RMS (const double *dpSrc, long lCount)
{
    long lLoopCntr;
    double dSqSum = 0.0;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        #ifndef _ISOC9X_SOURCE
        dSqSum += dpSrc[lLoopCntr] * dpSrc[lLoopCntr];
        #else
        dSqSum = fma(dpSrc[lLoopCntr], dpSrc[lLoopCntr], dSqSum);
        #endif
    }
    return (sqrt(dSqSum / (double) lCount));
}


float clDSPOp::Variance (float *fpVariance, float *fpMean, 
    const float *fpSrc, long lCount)
{
    long lLoopCntr;
    float fMean = 0.0F;
    float fVariance = 0.0F;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        fMean += fpSrc[lLoopCntr];
    }
    fMean /= (float) lCount;
    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        fVariance += (float) powf(fpSrc[lLoopCntr] - fMean, 2.0F);
    }
    fVariance /= (float) lCount;
    if (fpVariance != NULL)
    {
        *fpVariance = fVariance;
    }
    if (fpMean != NULL)
    {
        *fpMean = fMean;
    }
    return fVariance;
}


double clDSPOp::Variance (double *dpVariance, double *dpMean,
    const double *dpSrc, long lCount)
{
    long lLoopCntr;
    double dMean = 0.0;
    double dVariance = 0.0;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        dMean += dpSrc[lLoopCntr];
    }
    dMean /= (double) lCount;
    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        dVariance += (double) pow(dpSrc[lLoopCntr] - dMean, 2.0);
    }
    dVariance /= (double) lCount;
    if (dpVariance != NULL)
    {
        *dpVariance = dVariance;
    }
    if (dpMean != NULL)
    {
        *dpMean = dMean;
    }
    return dVariance;
}


float clDSPOp::PeakLevel (const float *fpSrc, long lSize)
{
    float fMin;
    float fMax;
    float fLinMax;

    MinMax(&fMin, &fMax, fpSrc, lSize);
    fMin = fabsf(fMin);
    fMax = fabsf(fMax);
    #ifndef _ISOC9X_SOURCE
    fLinMax = (fMax >= fMin) ? fMax : fMin;
    #else
    fLinMax = fmaxf(fMin, fMax);
    #endif
    return (20.0F * log10f(fLinMax));
}


double clDSPOp::PeakLevel (const double *dpSrc, long lSize)
{
    double dMin;
    double dMax;
    double dLinMax;

    MinMax(&dMin, &dMax, dpSrc, lSize);
    dMin = fabs(dMin);
    dMax = fabs(dMax);
    #ifndef _ISOC9X_SOURCE
    dLinMax = (dMax >= dMin) ? dMax : dMin;
    #else
    dLinMax = fmax(dMin, dMax);
    #endif
    return (20.0 * log10(dLinMax));
}


void clDSPOp::WinBartlett (float *fpDest, long lSize)
{
    #ifdef DSP_IPP
    ippsSet_32f(1.0F, fpDest, lSize);
    ippsWinBartlett_32f_I(fpDest, lSize);
    #else
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr <= ((lSize - 1L) / 2L); lLoopCntr++)
    {
        fpDest[lLoopCntr] = 2.0F * (float) lLoopCntr / (float) (lSize - 1L);
    }
    for (lLoopCntr = ((lSize - 1L) / 2L + 1L); lLoopCntr < lSize; lLoopCntr++)
    {
        fpDest[lLoopCntr] = 
            2.0F - 2.0F * (float) lLoopCntr / (float) (lSize - 1L);
    }
    #endif
}


void clDSPOp::WinBartlett (double *dpDest, long lSize)
{
    #ifdef DSP_IPP
    ippsSet_64f(1.0, dpDest, lSize);
    ippsWinBartlett_64f_I(dpDest, lSize);
    #else
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr <= ((lSize - 1L) / 2L); lLoopCntr++)
    {
        dpDest[lLoopCntr] = 2.0 * (double) lLoopCntr / (double) (lSize - 1L);
    }
    for (lLoopCntr = ((lSize - 1L) / 2L + 1L); lLoopCntr < lSize; lLoopCntr++)
    {
        dpDest[lLoopCntr] =
            2.0 - 2.0 * (double) lLoopCntr / (double) (lSize - 1L);
    }
    #endif
}


void clDSPOp::WinBlackman (float *fpDest, long lSize)
{
    #ifdef DSP_IPP
    ippsSet_32f(1.0F, fpDest, lSize);
    ippsWinBlackmanStd_32f_I(fpDest, lSize);
    #else
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lSize; lLoopCntr++)
    {
        fpDest[lLoopCntr] = 0.42F - 
            0.5F * cos((2.0F * fPI * (float) lLoopCntr) / (float) lSize) +
            0.08F * cos((4.0F * fPI * (float) lLoopCntr) / (float) lSize);
    }
    #endif
}


void clDSPOp::WinBlackman (double *dpDest, long lSize)
{
    #ifdef DSP_IPP
    ippsSet_64f(1.0, dpDest, lSize);
    ippsWinBlackmanStd_64f_I(dpDest, lSize);
    #else
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lSize; lLoopCntr++)
    {
        dpDest[lLoopCntr] = 0.42 - 
            0.5 * cos((2.0 * dPI * (double) lLoopCntr) / (double) lSize) +
            0.08 * cos((4.0 * dPI * (double) lLoopCntr) / (double) lSize);
    }
    #endif
}


void clDSPOp::WinBlackman (float *fpDest, long lSize, float fAlphaP)
{
    #ifdef DSP_IPP
    ippsSet_32f(1.0F, fpDest, lSize);
    if (fAlphaP != 0.0F)
        ippsWinBlackman_32f_I(fpDest, lSize, fAlphaP);
    else
        ippsWinBlackmanOpt_32f_I(fpDest, lSize);
    #else
    long lLoopCntr;
    float fAlpha;

    fAlpha = (fAlphaP != 0.0F) ?
        fAlphaP :
        (0.5F / (1.0F + cos((2.0F * fPI) / (float) (lSize - 1))));
    for (lLoopCntr = 0L; lLoopCntr < lSize; lLoopCntr++)
    {
        fpDest[lLoopCntr] = ((fAlpha + 1.0F) / 2.0F) - 0.5F * 
            cos((2.0F * fPI * (float) lLoopCntr) / ((float) (lSize - 1))) -
            (fAlpha / 2.0F) *
            cos((4.0F * fPI * (float) lLoopCntr) / ((float) (lSize - 1)));
    }
    #endif
}


void clDSPOp::WinBlackman (double *dpDest, long lSize, double dAlphaP)
{
    #ifdef DSP_IPP
    ippsSet_64f(1.0, dpDest, lSize);
    if (dAlphaP != 0.0)
        ippsWinBlackman_64f_I(dpDest, lSize, dAlphaP);
    else
        ippsWinBlackmanOpt_64f_I(dpDest, lSize);
    #else
    long lLoopCntr;
    double dAlpha;

    dAlpha = (dAlphaP != 0.0) ?
        dAlphaP :
        (0.5 / (1.0 + cos((2.0 * dPI) / (double) (lSize - 1))));
    for (lLoopCntr = 0L; lLoopCntr < lSize; lLoopCntr++)
    {
        dpDest[lLoopCntr] = ((dAlpha + 1.0) / 2.0) - 0.5 * 
            cos((2.0 * dPI * (double) lLoopCntr) / ((double) (lSize - 1))) -
            (dAlpha / 2.0) *
            cos((4.0 * dPI * (double) lLoopCntr) / ((double) (lSize - 1)));
    }
    #endif
}


void clDSPOp::WinBlackmanHarris (float *fpDest, long lSize)
{
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lSize; lLoopCntr++)
    {
        fpDest[lLoopCntr] = 0.42323F -
            0.49855F * 
            cos((2.0F * fPI * (float) lLoopCntr) / (float) lSize) +
            0.07922F * 
            cos((4.0F * fPI * (float) lLoopCntr) / (float) lSize);
    }
}


void clDSPOp::WinBlackmanHarris (double *dpDest, long lSize)
{
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lSize; lLoopCntr++)
    {
        dpDest[lLoopCntr] = 0.42323 -
            0.49855 *
            cos((2.0 * dPI * (double) lLoopCntr) / (double) lSize) +
            0.07922 *
            cos((4.0 * dPI * (double) lLoopCntr) / (double) lSize);
    }
}


void clDSPOp::WinCosTapered (float *fpDest, long lSize)
{
    long lLoopCntr;
    long lM;

    lM = Round((float) lSize / 10.0F);
    for (lLoopCntr = 0L; lLoopCntr < lSize; lLoopCntr++)
    {
        if ((lLoopCntr < lM) || (lLoopCntr > (lSize - lM - 1L)))
        {
            fpDest[lLoopCntr] = 0.5F * (1.0F -
                cos((2.0F * fPI * (float) lLoopCntr) / (float) lSize));
        }
        else
        {
            fpDest[lLoopCntr] = 1.0F;
        }
    }
}


void clDSPOp::WinCosTapered (double *dpDest, long lSize)
{
    long lLoopCntr;
    long lM;

    lM = Round((double) lSize / 10.0);
    for (lLoopCntr = 0L; lLoopCntr < lSize; lLoopCntr++)
    {
        if ((lLoopCntr < lM) || (lLoopCntr > (lSize - lM - 1L)))
        {
            dpDest[lLoopCntr] = 0.5 * (1.0 - 
                cos((2.0 * dPI * (double) lLoopCntr) / (double) lSize));
        }
        else
        {
            dpDest[lLoopCntr] = 1.0;
        }
    }
}


void clDSPOp::WinCosTaperedA (float *fpVect, long lSize)
{
    long lLoopCntr;
    long lM;

    lM = Round((float) lSize / 10.0F);
    for (lLoopCntr = 0L; lLoopCntr < lM; lLoopCntr++)
    {
        fpVect[lLoopCntr] = 0.5F * fpVect[lLoopCntr] *
            (1.0F - cos((2.0F * fPI * (float) lLoopCntr) / (float) lSize));
    }
    for (lLoopCntr = (lSize - lM); lLoopCntr < lSize; lLoopCntr++)
    {
        fpVect[lLoopCntr] = 0.5F * fpVect[lLoopCntr] *
            (1.0F - cos((2.0F * fPI * (float) lLoopCntr) / (float) lSize));
    }
}


void clDSPOp::WinCosTaperedA (double *dpVect, long lSize)
{
    long lLoopCntr;
    long lM;

    lM = Round((double) lSize / 10.0);
    for (lLoopCntr = 0L; lLoopCntr < lM; lLoopCntr++)
    {
        dpVect[lLoopCntr] = 0.5 * dpVect[lLoopCntr] *
            (1.0 - cos((2.0 * dPI * (double) lLoopCntr) / (double) lSize));
    }
    for (lLoopCntr = (lSize - lM); lLoopCntr < lSize; lLoopCntr++)
    {
        dpVect[lLoopCntr] = 0.5 * dpVect[lLoopCntr] *
            (1.0 - cos((2.0 * dPI * (double) lLoopCntr) / (double) lSize));
    }
}


void clDSPOp::WinCosTaperedA (float *fpDest, const float *fpSrc, long lSize)
{
    long lLoopCntr;
    long lM;

    lM = Round((float) lSize / 10.0F);
    for (lLoopCntr = 0L; lLoopCntr < lSize; lLoopCntr++)
    {
        if ((lLoopCntr < lM) || (lLoopCntr > (lSize - lM)))
        {
            fpDest[lLoopCntr] = 0.5F * fpSrc[lLoopCntr] * (1.0F - 
                cos((2.0F * fPI * (float) lLoopCntr) / (float) lSize));
        }
        else
        {
            fpDest[lLoopCntr] = fpSrc[lLoopCntr];
        }
    }
}


void clDSPOp::WinCosTaperedA (double *dpDest, const double *dpSrc, long lSize)
{
    long lLoopCntr;
    long lM;

    lM = Round((double) lSize / 10.0);
    for (lLoopCntr = 0L; lLoopCntr < lSize; lLoopCntr++)
    {
        if ((lLoopCntr < lM) || (lLoopCntr > (lSize - lM)))
        {
            dpDest[lLoopCntr] = 0.5 * dpSrc[lLoopCntr] * (1.0 -
                cos((2.0 * dPI * (double) lLoopCntr) / (double) lSize));
        }
        else
        {
            dpDest[lLoopCntr] = dpSrc[lLoopCntr];
        }
    }
}


void clDSPOp::WinExactBlackman (float *fpDest, long lSize)
{
    long lLoopCntr;
    float fA0;
    float fA1;
    float fA2;

    fA0 = 7938.0F / 18608.0F;
    fA1 = 9240.0F / 18608.0F;
    fA2 = 1430.0F / 18608.0F;
    for (lLoopCntr = 0L; lLoopCntr < lSize; lLoopCntr++)
    {
        fpDest[lLoopCntr] = fA0 - 
            fA1 * cos((2.0F * fPI * (float) lLoopCntr) / (float) lSize) +
            fA2 * cos((4.0F * fPI * (float) lLoopCntr) / (float) lSize);
    }
}


void clDSPOp::WinExactBlackman (double *dpDest, long lSize)
{
    long lLoopCntr;
    double dA0;
    double dA1;
    double dA2;

    dA0 = 7938.0 / 18608.0;
    dA1 = 9240.0 / 18608.0;
    dA2 = 1430.0 / 18608.0;
    for (lLoopCntr = 0L; lLoopCntr < lSize; lLoopCntr++)
    {
        dpDest[lLoopCntr] = dA0 -
            dA1 * cos((2.0 * dPI * (double) lLoopCntr) / (double) lSize) +
            dA2 * cos((4.0 * dPI * (double) lLoopCntr) / (double) lSize);
    }
}


void clDSPOp::WinExp (float *fpDest, float fFinal, long lSize)
{
    long lLoopCntr;
    float fA;

    fA = log(fFinal + 1.0F) / ((float) lSize / 2.0F);
    for (lLoopCntr = 0L; lLoopCntr < (lSize / 2L + 1L); lLoopCntr++)
    {
        fpDest[lLoopCntr] = fpDest[lSize - lLoopCntr - 1L] = 
            exp(fA * (float) lLoopCntr) - 1.0F;
    }
}


void clDSPOp::WinExp (double *dpDest, double dFinal, long lSize)
{
    long lLoopCntr;
    double dA;

    dA = log(dFinal + 1.0) / ((double) lSize / 2.0);
    for (lLoopCntr = 0L; lLoopCntr < (lSize / 2L + 1L); lLoopCntr++)
    {
        dpDest[lLoopCntr] = dpDest[lSize - lLoopCntr - 1L] = 
            exp(dA * (double) lLoopCntr) - 1.0;
    }
}


void clDSPOp::WinFlatTop (float *fpDest, long lSize)
{
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lSize; lLoopCntr++)
    {
        fpDest[lLoopCntr] = 0.2810639F - 
            0.5208972F * 
            cos((2.0F * fPI * (float) lLoopCntr) / (float) lSize) +
            0.1980399F * 
            cos((4.0F * fPI * (float) lLoopCntr) / (float) lSize);
    }
}


void clDSPOp::WinFlatTop (double *dpDest, long lSize)
{
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lSize; lLoopCntr++)
    {
        dpDest[lLoopCntr] = 0.2810639 -
            0.5208972 *
            cos((2.0 * dPI * (double) lLoopCntr) / (double) lSize) +
            0.1980399 *
            cos((4.0 * dPI * (double) lLoopCntr) / (double) lSize);
    }
}


void clDSPOp::WinGenericCos (float *fpDest, long lSize, 
    const float *fpCoeff, long lCoeffCount)
{
    long lLoopWin;
    long lLoopCoeff;
    float fWinCoeff;

    for (lLoopWin = 0L; lLoopWin < lSize; lLoopWin++)
    {
        fWinCoeff = 0.0F;
        for (lLoopCoeff = 0L; lLoopCoeff < lCoeffCount; lLoopCoeff++)
        {
            fWinCoeff += pow(-1.0F, (float) lLoopCoeff) * 
                fpCoeff[lLoopCoeff] * cos(
                (2.0F * fPI * (float) lLoopCoeff * (float) lLoopWin) /
                (float) lSize);
        }
        fpDest[lLoopWin] = fWinCoeff;
    }
}


void clDSPOp::WinGenericCos (double *dpDest, long lSize,
    const double *dpCoeff, long lCoeffCount)
{
    long lLoopWin;
    long lLoopCoeff;
    double dWinCoeff;

    for (lLoopWin = 0L; lLoopWin < lSize; lLoopWin++)
    {
        dWinCoeff = 0.0;
        for (lLoopCoeff = 0L; lLoopCoeff < lCoeffCount; lLoopCoeff++)
        {
            dWinCoeff += pow(-1.0, (double) lLoopCoeff) *
                dpCoeff[lLoopCoeff] * cos(
                (2.0 * dPI * (double) lLoopCoeff * (double) lLoopWin) /
                (double) lSize);
        }
        dpDest[lLoopWin] = dWinCoeff;
    }
}


void clDSPOp::WinHamming (float *fpDest, long lSize)
{
    #ifdef DSP_IPP
    ippsSet_32f(1.0F, fpDest, lSize);
    ippsWinHamming_32f_I(fpDest, lSize);
    #else
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lSize; lLoopCntr++)
    {
        fpDest[lLoopCntr] = 0.54F -
            0.46F * cos((2.0F * fPI * (float) lLoopCntr) / (float) lSize);
    }
    #endif
}


void clDSPOp::WinHamming (double *dpDest, long lSize)
{
    #ifdef DSP_IPP
    ippsSet_64f(1.0, dpDest, lSize);
    ippsWinHamming_64f_I(dpDest, lSize);
    #else
    long lLoopCntr;
   
    for (lLoopCntr = 0L; lLoopCntr < lSize; lLoopCntr++)
    {
        dpDest[lLoopCntr] = 0.54 - 
            0.46 * cos((2.0 * dPI * (double) lLoopCntr) / (double) lSize);
    }
    #endif
}


void clDSPOp::WinHanning (float *fpDest, long lSize)
{
    #ifdef DSP_IPP
    ippsSet_32f(1.0F, fpDest, lSize);
    ippsWinHann_32f_I(fpDest, lSize);
    #else
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lSize; lLoopCntr++)
    {
        fpDest[lLoopCntr] = 0.5F -
            0.5F * cos((2.0F * fPI * (float) lLoopCntr) / (float) lSize);
    }
    #endif
}


void clDSPOp::WinHanning (double *dpDest, long lSize)
{
    #ifdef DPS_IPP
    ippsSet_64f(1.0, dpDest, lSize);
    ippsWinHann_64f_I(dpDest, lSize);
    #else
    long lLoopCntr;
   
    for (lLoopCntr = 0L; lLoopCntr < lSize; lLoopCntr++)
    {
        dpDest[lLoopCntr] = 0.5 - 
            0.5 * cos((2.0 * dPI * (double) lLoopCntr) / (double) lSize);
    }
    #endif
}


void clDSPOp::WinKaiser (float *fpDest, float fBeta, long lSize)
{
    #ifdef DSP_IPP
    ippsSet_32f(1.0F, fpDest, lSize);
    ippsWinKaiser_32f_I(fpDest, lSize, 
        fBeta * (1.0f / ((float) (lSize - 1) / 2.0f)));
    #else
    long lLoopCntr;
    float fA;

    for (lLoopCntr = 0L; lLoopCntr < lSize; lLoopCntr++)
    {
        fA = fabs(1.0F - (2.0F * (float) lLoopCntr) / (float) lSize);
        fpDest[lLoopCntr] =
            ModZeroBessel(fBeta * sqrt(1.0F - fA * fA)) /
            ModZeroBessel(fBeta);
    }
    #endif
}


void clDSPOp::WinKaiser (double *dpDest, double dBeta, long lSize)
{
    #ifdef DSP_IPP
    ippsSet_64f(1.0, dpDest, lSize);
    ippsWinKaiser_64f_I(dpDest, lSize, 
        dBeta * (1.0 / ((double) (lSize - 1) / 2.0)));
    #else
    long lLoopCntr;
    double dA;

    for (lLoopCntr = 0L; lLoopCntr < lSize; lLoopCntr++)
    {
        dA = fabs(1.0 - (2.0 * (double) lLoopCntr) / (double) lSize);
        dpDest[lLoopCntr] =
            ModZeroBessel(dBeta * sqrt(1.0 - dA * dA)) /
            ModZeroBessel(dBeta);
    }
    #endif
}


void clDSPOp::WinKaiserBessel (float *fpDest, float fAlpha, long lSize)
{
    long lLoopCntr;
    float fHalfN;
    float fPiAlpha;

    fHalfN = (float) lSize / 2.0F;
    fPiAlpha = fPI * fAlpha;
    for (lLoopCntr = 0L; lLoopCntr < lSize; lLoopCntr++)
    {
        fpDest[lLoopCntr] =
            ModZeroBessel(fPiAlpha * 
            sqrt(1.0F - pow(((float) lLoopCntr - fHalfN) / fHalfN, 2.0F))) /
            ModZeroBessel(fPiAlpha);
    }
}


void clDSPOp::WinKaiserBessel (double *dpDest, double dAlpha, long lSize)
{
    long lLoopCntr;
    double dHalfN;
    double dPiAlpha;

    dHalfN = (double) lSize / 2.0;
    dPiAlpha = dPI * dAlpha;
    for (lLoopCntr = 0L; lLoopCntr < lSize; lLoopCntr++)
    {
        dpDest[lLoopCntr] =
            ModZeroBessel(dPiAlpha *
            sqrt(1.0 - pow(((double) lLoopCntr - dHalfN) / dHalfN, 2.0))) /
            ModZeroBessel(dPiAlpha);
    }
}


void clDSPOp::WinTukey (float *fpDest, long lSize)
{
    long lLoopCntr;
    long lHalfSize;

    lHalfSize = lSize / 2L;
    for (lLoopCntr = 0L; lLoopCntr < lSize; lLoopCntr++)
    {
        fpDest[lLoopCntr] = 0.5F * 
            (1.0F + cos((float) (lLoopCntr - lHalfSize) * fPI / lHalfSize));
    }
}


void clDSPOp::WinTukey (double *dpDest, long lSize)
{
    long lLoopCntr;
    long lHalfSize;

    lHalfSize = lSize / 2L;
    for (lLoopCntr = 0L; lLoopCntr < lSize; lLoopCntr++)
    {
        dpDest[lLoopCntr] = 0.5 *
            (1.0 + cos((float) (lLoopCntr - lHalfSize) * dPI / lHalfSize));
    }
}


void clDSPOp::WinDolphChebyshev (float *fpDest, float fGamma, long lSize)
{
    long lHalfSize;
    long lLoopCntr;
    long lSumCntr;
    float fScale;
    float fInvGamma;
    float fX0;
    float fThetaS;
    float fChebyshevSum;
    float fTheta;
    float fTheta2;
    float fMin;
    float fMax;
    
    lHalfSize = lSize / 2L;
    fScale = 1.0F / (float) lSize;
    fInvGamma = 1.0F / fGamma;
    fX0 = coshf((1.0F / (float) (lSize - 1)) * acoshf(fInvGamma));
    fThetaS = 2.0F * acosf(1.0F / fX0);
    for (lLoopCntr = 0L; lLoopCntr < lSize; lLoopCntr++)
    {
        fChebyshevSum = 0.0F;
        for (lSumCntr = 1L; lSumCntr <= lHalfSize; lSumCntr++)
        {
            fTheta = (2.0F * fPI * (float) lSumCntr) / (float) lSize;
            fTheta2 = (fPI * (float) (2L * lLoopCntr - lSize + 1L)) /
                (float) lSize;
            fChebyshevSum += 
                ChebyshevPolynom((float) (lSize - 1L),
                (float) (fX0 * cosf(fTheta / 2.0F))) * 
                cosf((float) lSumCntr * fTheta2);
        }
        fpDest[lLoopCntr] = fScale * 
            (1.0F + 2.0F * fGamma * fChebyshevSum);
    }

    MinMax(&fMin, &fMax, fpDest, lSize);
    fScale = 1.0F / fMax;
    Mul(fpDest, fScale, lSize);
}


void clDSPOp::WinDolphChebyshev (double *dpDest, double dGamma, long lSize)
{
    long lHalfSize;
    long lLoopCntr;
    long lSumCntr;
    double dScale;
    double dInvGamma;
    double dX0;
    double dThetaS;
    double dChebyshevSum;
    double dTheta;
    double dTheta2;
    double dMin;
    double dMax;
    
    lHalfSize = lSize / 2L;
    dScale = 1.0F / (double) lSize;
    dInvGamma = 1.0F / dGamma;
    dX0 = cosh((1.0 / (double) (lSize - 1)) * acosh(dInvGamma));
    dThetaS = 2.0 * acos(1.0 / dX0);
    for (lLoopCntr = 0L; lLoopCntr < lSize; lLoopCntr++)
    {
        dChebyshevSum = 0.0;
        for (lSumCntr = 1L; lSumCntr <= lHalfSize; lSumCntr++)
        {
            dTheta = (2.0 * dPI * (double) lSumCntr) / (double) lSize;
            dTheta2 = (dPI * (double) (2L * lLoopCntr - lSize + 1L)) /
                (double) lSize;
            dChebyshevSum += 
                ChebyshevPolynom((double) (lSize - 1L),
                dX0 * cos(dTheta / 2.0)) * 
                cos((double) lSumCntr * dTheta2);
        }
        dpDest[lLoopCntr] = dScale * 
            (1.0 + 2.0 * dGamma * dChebyshevSum);
    }

    MinMax(&dMin, &dMax, dpDest, lSize);
    dScale = 1.0 / dMax;
    Mul(dpDest, dScale, lSize);
}


void clDSPOp::Mix (float *fpDest, const float *fpSrc, long lCount)
{
    long lLoopCntr;
    long lSrcIdx;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        lSrcIdx = lLoopCntr << 1;
        fpDest[lLoopCntr] = (fpSrc[lSrcIdx] + fpSrc[lSrcIdx + 1L]) * 0.5F;
    }
}


void clDSPOp::Mix (double *dpDest, const double *dpSrc, long lCount)
{
    long lLoopCntr;
    long lSrcIdx;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        lSrcIdx = lLoopCntr << 1;
        dpDest[lLoopCntr] = (dpSrc[lSrcIdx] + dpSrc[lSrcIdx + 1L]) * 0.5;
    }
}


void clDSPOp::Mix (float *fpDest, const float *fpSrc1, const float *fpSrc2,
    long lCount)
{
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        fpDest[lLoopCntr] = (fpSrc1[lLoopCntr] + fpSrc2[lLoopCntr]) * 0.5F;
    }
}


void clDSPOp::Mix (double *dpDest, const double *dpSrc1, const double *dpSrc2,
    long lCount)
{
    long lLoopCntr;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        dpDest[lLoopCntr] = (dpSrc1[lLoopCntr] + dpSrc2[lLoopCntr]) * 0.5;
    }
}


void clDSPOp::Mix (float *fpDest, const float *fpSrc, long lChCount,
    long lDestCount)
{
    long lDestCntr;
    long lSrcCntr;
    long lSrcMax;
    float fMix;
    float fScaler;

    fScaler = 1.0F / (float) lChCount;
    for (lDestCntr = 0L; lDestCntr < lDestCount; lDestCntr++)
    {
        fMix = 0.0F;
        lSrcMax = lDestCntr * lChCount + lChCount;
        for (lSrcCntr = lDestCntr * lChCount; lSrcCntr < lSrcMax; lSrcCntr++)
        {
            fMix += fpSrc[lSrcCntr];
        }
        fpDest[lDestCntr] = fMix * fScaler;
    }
}


void clDSPOp::Mix (double *dpDest, const double *dpSrc, long lChCount,
    long lDestCount)
{
    long lDestCntr;
    long lSrcCntr;
    long lSrcMax;
    double dMix;
    double dScaler;

    dScaler = 1.0 / (double) lChCount;
    for (lDestCntr = 0L; lDestCntr < lDestCount; lDestCntr++)
    {
        dMix = 0.0;
        lSrcMax = lDestCntr * lChCount + lChCount;
        for (lSrcCntr = lDestCntr * lChCount; lSrcCntr < lSrcMax; lSrcCntr++)
        {
            dMix += dpSrc[lSrcCntr];
        }
        dpDest[lDestCntr] = dMix * dScaler;
    }
}


void clDSPOp::Spatialize (float *fpDest1, float *fpDest2, 
    const float *fpSrc, long lCount)
{
    long lLoopCntr;
    float fTemp;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        fTemp = fpSrc[lLoopCntr];
        fpDest1[lLoopCntr] = fTemp;
        fpDest2[lLoopCntr] = -fTemp;
    }
}


void clDSPOp::Spatialize (double *dpDest1, double *dpDest2,
    const double *dpSrc, long lCount)
{
    long lLoopCntr;
    double dTemp;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        dTemp = dpSrc[lLoopCntr];
        dpDest1[lLoopCntr] = dTemp;
        dpDest2[lLoopCntr] = -dTemp;
    }
}


void clDSPOp::Spatialize (float *fpDest, const float *fpSrc, long lCount)
{
    long lSrcCntr;
    long lDestIdx;
    float fTemp;

    for (lSrcCntr = 0L; lSrcCntr < lCount; lSrcCntr++)
    {
        fTemp = fpSrc[lSrcCntr];
        lDestIdx = lSrcCntr << 1;
        fpDest[lDestIdx] = fTemp;
        fpDest[lDestIdx + 1L] = -fTemp;
    }
}


void clDSPOp::Spatialize (double *dpDest, const double *dpSrc, long lCount)
{
    long lSrcCntr;
    long lDestIdx;
    double dTemp;

    for (lSrcCntr = 0L; lSrcCntr < lCount; lSrcCntr++)
    {
        dTemp = dpSrc[lSrcCntr];
        lDestIdx = lSrcCntr << 1;
        dpDest[lDestIdx] = dTemp;
        dpDest[lDestIdx + 1L] = -dTemp;
    }
}


void clDSPOp::Extract (float *fpDest, const float *fpSrc, long lChannel,
    long lChCount, long lSrcLength)
{
    long lLoopCntr;
    long lMax;

    lMax = lSrcLength / lChCount;
    for (lLoopCntr = 0; lLoopCntr < lMax; lLoopCntr++)
    {
        fpDest[lLoopCntr] = fpSrc[lLoopCntr * lChCount + lChannel];
    }
}


void clDSPOp::Extract (double *dpDest, const double *dpSrc, long lChannel,
    long lChCount, long lSrcLength)
{
    long lLoopCntr;
    long lMax;

    lMax = lSrcLength / lChCount;
    for (lLoopCntr = 0; lLoopCntr < lMax; lLoopCntr++)
    {
        dpDest[lLoopCntr] = dpSrc[lLoopCntr * lChCount + lChannel];
    }
}


void clDSPOp::Pack (float *fpDest, const float *fpSrc, long lChannel,
    long lChCount, long lSrcLength)
{
    long lLoopCntr;

    for (lLoopCntr = 0; lLoopCntr < lSrcLength; lLoopCntr++)
    {
        fpDest[lLoopCntr * lChCount + lChannel] = fpSrc[lLoopCntr];
    }
}


void clDSPOp::Pack (double *dpDest, const double *dpSrc, long lChannel,
    long lChCount, long lSrcLength)
{
    long lLoopCntr;

    for (lLoopCntr = 0; lLoopCntr < lSrcLength; lLoopCntr++)
    {
        dpDest[lLoopCntr * lChCount + lChannel] = dpSrc[lLoopCntr];
    }
}


long clDSPOp::ReBuffer (float *fpDest, const float *fpSrc, long lDestSize,
    long lSrcSize)
{
    long lDestCount;
    long lSrcCount;
    long lCopyCount;

    lDestCount = lDestSize - lPrevDestCount;
    lSrcCount = lSrcSize - lPrevSrcCount;
    lCopyCount = (lDestCount < lSrcCount) ? lDestCount : lSrcCount;
    Copy(&fpDest[lPrevDestCount], &fpSrc[lPrevSrcCount], lCopyCount);
    lPrevDestCount += lCopyCount;
    lPrevSrcCount += lCopyCount;
    if ((lPrevDestCount == lDestSize) && (lPrevSrcCount == lSrcSize))
    {
        lPrevDestCount = 0L;
        lPrevSrcCount = 0L;
        return 1;
    }
    else if (lPrevDestCount == lDestSize)
    {
        lPrevDestCount = 0L;
        return 2;
    }
    else if (lPrevSrcCount == lSrcSize)
    {
        lPrevSrcCount = 0L;
        return 0;
    }
    else
    {
        fprintf(stderr, "clDSPOp::ReBuffer(): Fatal error; bug found\n");
    }
    return 0;
}


long clDSPOp::ReBuffer (double *dpDest, const double *dpSrc, long lDestSize,
    long lSrcSize)
{
    long lDestCount;
    long lSrcCount;
    long lCopyCount;

    lDestCount = lDestSize - lPrevDestCount;
    lSrcCount = lSrcSize - lPrevSrcCount;
    lCopyCount = (lDestCount < lSrcCount) ? lDestCount : lSrcCount;
    Copy(&dpDest[lPrevDestCount], &dpSrc[lPrevSrcCount], lCopyCount);
    lPrevDestCount += lCopyCount;
    lPrevSrcCount += lCopyCount;
    if (lPrevDestCount == lDestSize && lPrevSrcCount == lSrcSize)
    {
        lPrevDestCount = 0L;
        lPrevSrcCount = 0L;
        return 1;
    }
    else if (lPrevDestCount == lDestSize)
    {
        lPrevDestCount = 0L;
        return 2;
    }
    else if (lPrevSrcCount == lSrcSize)
    {
        lPrevSrcCount = 0L;
        return 0;
    }
    else
    {
        fprintf(stderr, "clDSPOp::ReBuffer(): Fatal error; bug found\n");
    }
    return 0;
}


/* I know, here we actually play with one unneccessary value in FIRBuf,
   but it makes code a lot cleaner, small price for one extra value! */

void clDSPOp::FIRAllocate (const float *fpCoeff, long lCount)
{
    #ifdef DSP_IPP
    iFIRDlyIdx = 0;
    lFIRLength = lCount;
    FIRCoeff.Size(lCount * sizeof(float));
    FIRBuf.Size(lCount * 2 * sizeof(float));
    Copy((float *) FIRCoeff, fpCoeff, lCount);
    Zero((float *) FIRBuf, lCount * 2);
    #else
    lFIRLength = lCount;
    FIRCoeff.Size(lCount * sizeof(float));
    FIRBuf.Size(lCount * sizeof(float));
    Copy((float *) FIRCoeff, fpCoeff, lCount);
    Zero((float *) FIRBuf, lCount);
    #endif
}


void clDSPOp::FIRAllocate (const double *dpCoeff, long lCount)
{
    #ifdef DSP_IPP
    iFIRDlyIdx = 0;
    lFIRLength = lCount;
    FIRCoeff.Size(lCount * sizeof(double));
    FIRBuf.Size(lCount * 2 * sizeof(double));
    Copy((double *) FIRCoeff, dpCoeff, lCount);
    Zero((double *) FIRBuf, lCount * 2);
    #else
    lFIRLength = lCount;
    FIRCoeff.Size(lCount * sizeof(double));
    FIRBuf.Size(lCount * sizeof(double));
    Copy((double *) FIRCoeff, dpCoeff, lCount);
    Zero((double *) FIRBuf, lCount);
    #endif
}


void clDSPOp::FIRFilter (float *fpVect, long lCount)
{
    #ifdef DSP_IPP
    ippsFIR_Direct_32f_I(fpVect, lCount, FIRCoeff, lFIRLength,
        FIRBuf, &iFIRDlyIdx);
    #else
    long lMax;
    float *fpFIRCoeff = FIRCoeff;
    float *fpFIRBuf = FIRBuf;
    float *fpFIRWork;

    fpFIRWork = (float *) FIRWork.Size((lCount + lFIRLength) * sizeof(float));
    Copy(fpFIRWork, fpFIRBuf, lFIRLength);
    Copy(&fpFIRWork[lFIRLength], fpVect, lCount);
    lMax = lCount + lFIRLength;
    #ifdef DSP_X86
    if (bHave3DNow)
    {
        dsp_x86_3dnow_firf(fpVect, fpFIRWork, lCount,
            fpFIRCoeff, lFIRLength);
    }
    else if (iHaveSSE >= 2)
    {
        dsp_x86_sse_firf(fpVect, fpFIRWork, lCount,
            fpFIRCoeff, lFIRLength);
    }
    else
    #endif
    #if (defined(DSP_X86_64) && (!defined(DSP_X86)))
    dsp_x86_64_firf(fpVect, fpFIRWork, lCount, fpFIRCoeff, lFIRLength);
    #else
    {
        long lSrcCntr;
        long lConvCntr;
        long lDestCntr;
        float fTempVal;

        lDestCntr = 0L;
        for (lSrcCntr = lFIRLength; lSrcCntr < lMax; lSrcCntr++)
        {
            fTempVal = 0.0F;
            for (lConvCntr = 0L; lConvCntr < lFIRLength; lConvCntr++)
            {
                #ifndef _ISOC9X_SOURCE
                fTempVal += fpFIRCoeff[lConvCntr] * 
                    fpFIRWork[lSrcCntr - lConvCntr];
                #else
                fTempVal = fmaf(fpFIRCoeff[lConvCntr], 
                    fpFIRWork[lSrcCntr - lConvCntr], fTempVal);
                #endif
            }
            fpVect[lDestCntr++] = fTempVal;
        }
    }
    #endif
    Copy(fpFIRBuf, &fpFIRWork[lMax - lFIRLength], lFIRLength);
    #endif
}


void clDSPOp::FIRFilter (double *dpVect, long lCount)
{
    #ifdef DSP_IPP
    ippsFIR_Direct_64f_I(dpVect, lCount, FIRCoeff, lFIRLength,
        FIRBuf, &iFIRDlyIdx);
    #else
    long lMax;
    double *dpFIRCoeff = FIRCoeff;
    double *dpFIRBuf = FIRBuf;
    double *dpFIRWork;

    lMax = lCount + lFIRLength;
    dpFIRWork = (double *) FIRWork.Size((lCount + lFIRLength) * sizeof(double));
    Copy(dpFIRWork, dpFIRBuf, lFIRLength);
    Copy(&dpFIRWork[lFIRLength], dpVect, lCount);
    #ifdef DSP_X86
    if (iHaveSSE >= 2)
    {
        dsp_x86_sse_fir(dpVect, dpFIRWork, lCount,
            dpFIRCoeff, lFIRLength);
    }
    else
    #endif
    #if (defined(DSP_X86_64) && (!defined(DSP_X86)))
    dsp_x86_64_fir(dpVect, dpFIRWork, lCount, dpFIRCoeff, lFIRLength);
    #else
    {
        long lSrcCntr;
        long lConvCntr;
        long lDestCntr;
        double dTempVal;

        lDestCntr = 0L;
        for (lSrcCntr = lFIRLength; lSrcCntr < lMax; lSrcCntr++)
        {
            dTempVal = 0.0;
            for (lConvCntr = 0L; lConvCntr < lFIRLength; lConvCntr++)
            {
                #ifndef _ISOC9X_SOURCE
                dTempVal += dpFIRCoeff[lConvCntr] *
                    dpFIRWork[lSrcCntr - lConvCntr];
                #else
                dTempVal = fma(dpFIRCoeff[lConvCntr],
                    dpFIRWork[lSrcCntr - lConvCntr], dTempVal);
                #endif
            }
            dpVect[lDestCntr++] = dTempVal;
        }
    }
    #endif
    Copy(dpFIRBuf, &dpFIRWork[lMax - lFIRLength], lFIRLength);
    #endif
}


void clDSPOp::FIRFilter (float *fpDest, const float *fpSrc, long lCount)
{
    #ifdef DSP_IPP
    ippsFIR_Direct_32f(fpSrc, fpDest, lCount, FIRCoeff, lFIRLength,
        FIRBuf, &iFIRDlyIdx);
    #else
    long lMax;
    float *fpFIRCoeff = FIRCoeff;
    float *fpFIRBuf = FIRBuf;
    float *fpFIRWork;

    fpFIRWork = (float *) FIRWork.Size((lCount + lFIRLength) * sizeof(float));
    Copy(fpFIRWork, fpFIRBuf, lFIRLength);
    Copy(&fpFIRWork[lFIRLength], fpSrc, lCount);
    lMax = lCount + lFIRLength;
    #ifdef DSP_X86
    if (bHave3DNow)
    {
        dsp_x86_3dnow_firf(fpDest, fpFIRWork, lCount, 
            fpFIRCoeff, lFIRLength);
    }
    else if (iHaveSSE >= 2)
    {
        dsp_x86_sse_firf(fpDest, fpFIRWork, lCount,
            fpFIRCoeff, lFIRLength);
    }
    else
    #endif
    #if (defined(DSP_X86_64) && (!defined(DSP_X86)))
    dsp_x86_64_firf(fpDest, fpFIRWork, lCount, fpFIRCoeff, lFIRLength);
    #else
    {
        long lSrcCntr;
        long lConvCntr;
        long lDestCntr;
        float fTempVal;

        lDestCntr = 0L;
        for (lSrcCntr = lFIRLength; lSrcCntr < lMax; lSrcCntr++)
        {
            fTempVal = 0.0F;
            for (lConvCntr = 0L; lConvCntr < lFIRLength; lConvCntr++)
            {
                #ifndef _ISOC9X_SOURCE
                fTempVal += fpFIRCoeff[lConvCntr] * 
                    fpFIRWork[lSrcCntr - lConvCntr];
                #else
                fTempVal = fmaf(fpFIRCoeff[lConvCntr],
                    fpFIRWork[lSrcCntr - lConvCntr], fTempVal);
                #endif
            }
            fpDest[lDestCntr++] = fTempVal;
        }
    }
    #endif
    Copy(fpFIRBuf, &fpFIRWork[lMax - lFIRLength], lFIRLength);
    #endif
}


void clDSPOp::FIRFilter (double *dpDest, const double *dpSrc, long lCount)
{
    #ifdef DSP_IPP
    ippsFIR_Direct_64f(dpSrc, dpDest, lCount, FIRCoeff, lFIRLength,
        FIRBuf, &iFIRDlyIdx);
    #else
    long lMax;
    double *dpFIRCoeff = FIRCoeff;
    double *dpFIRBuf = FIRBuf;
    double *dpFIRWork;

    lMax = lCount + lFIRLength;
    dpFIRWork = (double *) FIRWork.Size((lCount + lFIRLength) * sizeof(double));
    Copy(dpFIRWork, dpFIRBuf, lFIRLength);
    Copy(&dpFIRWork[lFIRLength], dpSrc, lCount);
    #ifdef DSP_X86
    if (iHaveSSE >= 2)
    {
        dsp_x86_sse_fir(dpDest, dpFIRWork, lCount,
            dpFIRCoeff, lFIRLength);
    }
    else
    #endif
    #if (defined(DSP_X86_64) && (!defined(DSP_X86)))
    dsp_x86_64_fir(dpDest, dpFIRWork, lCount, dpFIRCoeff, lFIRLength);
    #else
    {
        long lSrcCntr;
        long lConvCntr;
        long lDestCntr;
        double dTempVal;

        lDestCntr = 0L;
        for (lSrcCntr = lFIRLength; lSrcCntr < lMax; lSrcCntr++)
        {
            dTempVal = 0.0;
            for (lConvCntr = 0L; lConvCntr < lFIRLength; lConvCntr++)
            {
                #ifndef _ISOC9X_SOURCE
                dTempVal += dpFIRCoeff[lConvCntr] *
                    dpFIRWork[lSrcCntr - lConvCntr];
                #else
                dTempVal = fma(dpFIRCoeff[lConvCntr],
                    dpFIRWork[lSrcCntr - lConvCntr], dTempVal);
                #endif
            }
            dpDest[lDestCntr++] = dTempVal;
        }
    }
    #endif
    Copy(dpFIRBuf, &dpFIRWork[lMax - lFIRLength], lFIRLength);
    #endif
}


void clDSPOp::FIRFilterF (float *fpDest, float *fpSrc, long lCount)
{
    long lSrcCntr;
    long lConvCntr;
    long lDestCntr;
    long lMax;
    float fTempVal;
    float *fpFIRCoeff = FIRCoeff;

    lDestCntr = 0L;
    lMax = lCount + lFIRLength;
    for (lSrcCntr = lFIRLength; lSrcCntr < lMax; lSrcCntr++)
    {
        fTempVal = 0.0F;
        for (lConvCntr = 0L; lConvCntr < lFIRLength; lConvCntr++)
        {
            #ifndef _ISOC9X_SOURCE
            fTempVal += fpFIRCoeff[lConvCntr] *
                fpSrc[lSrcCntr - lConvCntr];
            #else
            fTempVal = fmaf(fpFIRCoeff[lConvCntr],
                fpSrc[lSrcCntr - lConvCntr], fTempVal);
            #endif
        }
        fpDest[lDestCntr++] = fTempVal;
    }
    Copy(fpSrc, &fpSrc[lMax - lFIRLength], lFIRLength);
}


void clDSPOp::FIRFilterF (double *dpDest, double *dpSrc, long lCount)
{
    long lSrcCntr;
    long lConvCntr;
    long lDestCntr;
    long lMax;
    double dTempVal;
    double *dpFIRCoeff = FIRCoeff;

    lDestCntr = 0L;
    lMax = lCount + lFIRLength;
    for (lSrcCntr = lFIRLength; lSrcCntr < lMax; lSrcCntr++)
    {
        dTempVal = 0.0;
        for (lConvCntr = 0L; lConvCntr < lFIRLength; lConvCntr++)
        {
            #ifndef _ISOC9X_SOURCE
            dTempVal += dpFIRCoeff[lConvCntr] *
                dpSrc[lSrcCntr - lConvCntr];
            #else
            dTempVal = fma(dpFIRCoeff[lConvCntr],
                dpSrc[lSrcCntr - lConvCntr], dTempVal);
            #endif
        }
        dpDest[lDestCntr++] = dTempVal;
    }
    Copy(dpSrc, &dpSrc[lMax - lFIRLength], lFIRLength);
}


void clDSPOp::FIRFree ()
{
    FIRCoeff.Free();
    FIRBuf.Free();
}


void clDSPOp::IIRInitialize (const float *fpCoeffs)
{
    Copy(fpIIR_C, fpCoeffs, 5);
}


void clDSPOp::IIRInitialize (const double *dpCoeffs)
{
    Copy(dpIIR_C, dpCoeffs, 5);
}


void clDSPOp::IIRFilter (float *fpVect, long lCount)
{
    #ifdef DSP_X86
    if (bHave3DNow)
    {
        dsp_x86_3dnow_iirf(fpVect, lCount, fpIIR_C, fpIIR_X, fpIIR_Y);
    }
    else if (iHaveSSE >= 2)
    {
        dsp_x86_sse_iirf(fpVect, lCount, fpIIR_C, fpIIR_X, fpIIR_Y);
    }
    else
    #endif
    #if (defined(DSP_X86_64) && (!defined(DSP_X86)))
    dsp_x86_64_iirf(fpVect, lCount, fpIIR_C, fpIIR_X, fpIIR_Y);
    #else
    {
        long lLoopCntr;
        float fBXk;
        float fAYk;

        for (lLoopCntr = 0; lLoopCntr < lCount; lLoopCntr++)
        {
            fpIIR_X[0] = fpIIR_X[1];
            fpIIR_X[1] = fpIIR_X[2];
            fpIIR_X[2] = fpVect[lLoopCntr];
            fBXk =
                fpIIR_C[0] * fpIIR_X[2] + 
                fpIIR_C[1] * fpIIR_X[1] +
                fpIIR_C[2] * fpIIR_X[0];
            fAYk = 
                fpIIR_C[3] * fpIIR_Y[1] +
                fpIIR_C[4] * fpIIR_Y[0];
            fpVect[lLoopCntr] = fAYk + fBXk;
            fpIIR_Y[0] = fpIIR_Y[1];
            fpIIR_Y[1] = fpVect[lLoopCntr];
        }
    }
    #endif
}


void clDSPOp::IIRFilter (double *dpVect, long lCount)
{
    #ifdef DSP_X86
    if (iHaveSSE >= 2)
    {
        dsp_x86_sse_iir(dpVect, lCount, dpIIR_C, dpIIR_X, dpIIR_Y);
    }
    else
    #endif
    #if (defined(DSP_X86_64) && (!defined(DSP_X86)))
    dsp_x86_64_iir(dpVect, lCount, dpIIR_C, dpIIR_X, dpIIR_Y);
    #else
    {
        long lLoopCntr;
        double dBXk;
        double dAYk;

        for (lLoopCntr = 0; lLoopCntr < lCount; lLoopCntr++)
        {
            dpIIR_X[0] = dpIIR_X[1];
            dpIIR_X[1] = dpIIR_X[2];
            dpIIR_X[2] = dpVect[lLoopCntr];
            dBXk =
                dpIIR_C[0] * dpIIR_X[2] + 
                dpIIR_C[1] * dpIIR_X[1] +
                dpIIR_C[2] * dpIIR_X[0];
            dAYk = 
                dpIIR_C[3] * dpIIR_Y[1] +
                dpIIR_C[4] * dpIIR_Y[0];
            dpVect[lLoopCntr] = dAYk + dBXk;
            dpIIR_Y[0] = dpIIR_Y[1];
            dpIIR_Y[1] = dpVect[lLoopCntr];
        }
    }
    #endif
}


void clDSPOp::IIRFilter (float *fpDest, const float *fpSrc, long lCount)
{
    #ifdef DSP_X86
    if (bHave3DNow)
    {
        dsp_x86_3dnow_iirf_nip(fpDest, fpSrc, lCount, 
            fpIIR_C, fpIIR_X, fpIIR_Y);
    }
    else if (iHaveSSE >= 2)
    {
        dsp_x86_sse_iirf_nip(fpDest, fpSrc, lCount, 
            fpIIR_C, fpIIR_X, fpIIR_Y);
    }
    else
    #endif
    #if (defined(DSP_X86_64) && (!defined(DSP_X86)))
    dsp_x86_64_iirf_nip(fpDest, fpSrc, lCount, fpIIR_C, fpIIR_X, fpIIR_Y);
    #else
    {
        long lLoopCntr;
        float fBXk;
        float fAYk;

        for (lLoopCntr = 0; lLoopCntr < lCount; lLoopCntr++)
        {
            fpIIR_X[0] = fpIIR_X[1];
            fpIIR_X[1] = fpIIR_X[2];
            fpIIR_X[2] = fpSrc[lLoopCntr];
            fBXk =
                fpIIR_C[0] * fpIIR_X[2] + 
                fpIIR_C[1] * fpIIR_X[1] +
                fpIIR_C[2] * fpIIR_X[0];
            fAYk = 
                fpIIR_C[3] * fpIIR_Y[1] +
                fpIIR_C[4] * fpIIR_Y[0];
            fpDest[lLoopCntr] = fAYk + fBXk;
            fpIIR_Y[0] = fpIIR_Y[1];
            fpIIR_Y[1] = fpDest[lLoopCntr];
        }
    }
    #endif
}


void clDSPOp::IIRFilter (double *dpDest, const double *dpSrc, long lCount)
{
    #ifdef DSP_X86
    if (iHaveSSE >= 2)
    {
        dsp_x86_sse_iir_nip(dpDest, dpSrc, lCount, dpIIR_C, dpIIR_X, dpIIR_Y);
    }
    else
    #endif
    #if (defined(DSP_X86_64) && (!defined(DSP_X86)))
    dsp_x86_64_iir_nip(dpDest, dpSrc, lCount, dpIIR_C, dpIIR_X, dpIIR_Y);
    #else
    {
        long lLoopCntr;
        double dBXk;
        double dAYk;

        for (lLoopCntr = 0; lLoopCntr < lCount; lLoopCntr++)
        {
            dpIIR_X[0] = dpIIR_X[1];
            dpIIR_X[1] = dpIIR_X[2];
            dpIIR_X[2] = dpSrc[lLoopCntr];
            dBXk =
                dpIIR_C[0] * dpIIR_X[2] + 
                dpIIR_C[1] * dpIIR_X[1] +
                dpIIR_C[2] * dpIIR_X[0];
            dAYk = 
                dpIIR_C[3] * dpIIR_Y[1] +
                dpIIR_C[4] * dpIIR_Y[0];
            dpDest[lLoopCntr] = dAYk + dBXk;
            dpIIR_Y[0] = dpIIR_Y[1];
            dpIIR_Y[1] = dpDest[lLoopCntr];
        }
    }
    #endif
}


void clDSPOp::IIRClear ()
{
    Zero(fpIIR_X, 3);
    Zero(fpIIR_Y, 2);
    Zero(dpIIR_X, 3);
    Zero(dpIIR_Y, 2);
}


void clDSPOp::FFTInitialize(long lSize, bool bIsReal)
{
    #if defined(DSP_IPP)
    int iSWorkSize;
    int iDWorkSize;
    int iOrder;

    iOrder = (int) (log((double) lSize) / log(2.0) + 0.5);
    if (bIsReal)
    {
        IppsFFTSpec_R_32f *pSFFTSpec;
        IppsFFTSpec_R_64f *pDFFTSpec;
        ippsFFTInitAlloc_R_32f(&pSFFTSpec, iOrder, IPP_FFT_DIV_FWD_BY_N,
            ippAlgHintFast);
        ippsFFTInitAlloc_R_64f(&pDFFTSpec, iOrder, IPP_FFT_DIV_FWD_BY_N,
            ippAlgHintFast);
        ippsFFTGetBufSize_R_32f(pSFFTSpec, &iSWorkSize);
        ippsFFTGetBufSize_R_64f(pDFFTSpec, &iDWorkSize);
        SBitRevWork.Size(iSWorkSize);
        DBitRevWork.Size(iDWorkSize);
        vpSTfrm = (void *) pSFFTSpec;
        vpDTfrm = (void *) pDFFTSpec;
    }
    else
    {
        IppsFFTSpec_C_32fc *pSFFTSpec;
        IppsFFTSpec_C_64fc *pDFFTSpec;
        ippsFFTInitAlloc_C_32fc(&pSFFTSpec, iOrder, IPP_FFT_DIV_FWD_BY_N,
            ippAlgHintFast);
        ippsFFTInitAlloc_C_64fc(&pDFFTSpec, iOrder, IPP_FFT_DIV_FWD_BY_N,
            ippAlgHintFast);
        ippsFFTGetBufSize_C_32fc(pSFFTSpec, &iSWorkSize);
        ippsFFTGetBufSize_C_64fc(pDFFTSpec, &iDWorkSize);
        SBitRevWork.Size(iSWorkSize);
        DBitRevWork.Size(iDWorkSize);
        vpSTfrm = (void *) pSFFTSpec;
        vpDTfrm = (void *) pDFFTSpec;
    }
    #elif defined(DSP_USE_FFTW)
    if (bIsReal)
    {
        clDSPAlloc FFTBufS(lSize * sizeof(float));
        clDSPAlloc FFTBufD(lSize * sizeof(double));
        clDSPAlloc FFTBufSC((lSize / 2 + 1) * sizeof(fftwf_complex));
        clDSPAlloc FFTBufDC((lSize / 2 + 1) * sizeof(fftw_complex));
        fftwpSPlan = fftwf_plan_dft_r2c_1d(lSize, FFTBufS, FFTBufSC,
            FFTW_ESTIMATE|FFTW_PRESERVE_INPUT|FFTW_UNALIGNED);
        fftwpSIPlan = fftwf_plan_dft_c2r_1d(lSize, FFTBufSC, FFTBufS,
            FFTW_ESTIMATE|FFTW_PRESERVE_INPUT|FFTW_UNALIGNED);
        fftwpDPlan = fftw_plan_dft_r2c_1d(lSize, FFTBufD, FFTBufDC,
            FFTW_ESTIMATE|FFTW_PRESERVE_INPUT|FFTW_UNALIGNED);
        fftwpDIPlan = fftw_plan_dft_c2r_1d(lSize, FFTBufDC, FFTBufD,
            FFTW_ESTIMATE|FFTW_PRESERVE_INPUT|FFTW_UNALIGNED);
        FFTBuf.Size(lSize * sizeof(double));
    }
    else
    {
        clDSPAlloc FFTBufSI(lSize * sizeof(fftwf_complex));
        clDSPAlloc FFTBufSO(lSize * sizeof(fftwf_complex));
        clDSPAlloc FFTBufDI(lSize * sizeof(fftw_complex));
        clDSPAlloc FFTBufDO(lSize * sizeof(fftw_complex));

        fftwpSPlan = fftwf_plan_dft_1d(lSize, FFTBufSI, FFTBufSO,
            FFTW_FORWARD, FFTW_ESTIMATE|FFTW_PRESERVE_INPUT|FFTW_UNALIGNED);
        fftwpSIPlan = fftwf_plan_dft_1d(lSize, FFTBufSO, FFTBufSI,
            FFTW_BACKWARD, FFTW_ESTIMATE|FFTW_PRESERVE_INPUT|FFTW_UNALIGNED);
        fftwpDPlan = fftw_plan_dft_1d(lSize, FFTBufDI, FFTBufDO,
            FFTW_FORWARD, FFTW_ESTIMATE|FFTW_PRESERVE_INPUT|FFTW_UNALIGNED);
        fftwpDIPlan = fftw_plan_dft_1d(lSize, FFTBufDO, FFTBufDI,
            FFTW_BACKWARD, FFTW_ESTIMATE|FFTW_PRESERVE_INPUT|FFTW_UNALIGNED);
        FFTBuf.Size(lSize * 2 * sizeof(double));
    }
    fFFTScale = 1.0F / (float) lSize;
    dFFTScale = 1.0 / (double) lSize;
    #else
    float *fpFFTBuf;
    double *dpFFTBuf;
    clDSPAlloc FFTBufS(lSize * 2 * sizeof(float));
    clDSPAlloc FFTBufD(lSize * 2 * sizeof(double));

    fpFFTBuf = FFTBufS;
    dpFFTBuf = FFTBufD;
    if (bIsReal)
    {
        fFFTScale = 2.0F / (float) lSize;
        dFFTScale = 2.0 / (double) lSize;
    }
    else
    {
        fFFTScale = 1.0F / (float) lSize;
        dFFTScale = 1.0 / (double) lSize;
    }
    lpSBitRevWork = (long *)
        SBitRevWork.Size(((size_t) ceil(2.0 + sqrt(lSize))) * sizeof(long));
    lpDBitRevWork = (long *)
        DBitRevWork.Size(((size_t) ceil(2.0 + sqrt(lSize))) * sizeof(long));
    fpCosSinTable = (float *) 
        SCosSinTable.Size((lSize / 2 + 1) * sizeof(float));
    dpCosSinTable = (double *) 
        DCosSinTable.Size((lSize / 2 + 1) * sizeof(double));
    lpSBitRevWork[0] = 0;
    lpSBitRevWork[1] = 0;
    lpDBitRevWork[0] = 0;
    lpDBitRevWork[1] = 0;
    if (bIsReal)
    {
        Tfrm.rdft(lSize, 1, fpFFTBuf, lpSBitRevWork, fpCosSinTable);
        Tfrm.rdft(lSize, 1, dpFFTBuf, lpDBitRevWork, dpCosSinTable);
        FFTBuf.Size(lSize * sizeof(double));
    }
    else
    {
        Tfrm.cdft(lSize * 2, 1, fpFFTBuf, lpSBitRevWork, fpCosSinTable);
        Tfrm.cdft(lSize * 2, 1, dpFFTBuf, lpDBitRevWork, dpCosSinTable);
        FFTBuf.Size(lSize * 2 * sizeof(double));
    }
    #endif
    bFFTInitialized = true;
    bRealTransform = bIsReal;
    lFFTLength = lSize;
}


void clDSPOp::FFTUninitialize()
{
    #if defined(DSP_IPP)
    if (bFFTInitialized)
    {
        if (bRealTransform)
        {
            ippsFFTFree_R_32f((IppsFFTSpec_R_32f *) vpSTfrm);
            ippsFFTFree_R_64f((IppsFFTSpec_R_64f *) vpDTfrm);
        }
        else
        {
            ippsFFTFree_C_32f((IppsFFTSpec_C_32f *) vpSTfrm);
            ippsFFTFree_C_64f((IppsFFTSpec_C_64f *) vpDTfrm);
        }
    }
    #elif defined(DSP_USE_FFTW)
    if (bFFTInitialized)
    {
        fftwf_destroy_plan(fftwpSPlan);
        fftwf_destroy_plan(fftwpSIPlan);
        fftw_destroy_plan(fftwpDPlan);
        fftw_destroy_plan(fftwpDIPlan);
    }
    #else
    // NOP
    #endif
    SBitRevWork.Free();
    DBitRevWork.Free();
    SCosSinTable.Free();
    DCosSinTable.Free();
    FFTBuf.Free();
    bFFTInitialized = false;
}


void clDSPOp::FFTi(stpSCplx spDest, float *fpSrc)
{
    #if defined(DSP_IPP)
    ippsFFTFwd_RToCCS_32f(fpSrc, (Ipp32f *) spDest, 
        (IppsFFTSpec_R_32f *) vpSTfrm, SBitRevWork);
    #elif defined(DSP_USE_FFTW)
    Mul(fpSrc, fFFTScale, lFFTLength);
    fftwf_execute_dft_r2c(fftwpSPlan, fpSrc, (fftwf_complex *) spDest);
    #else
    long lLoopCntr;
    long lMax;
    
    Mul(fpSrc, fFFTScale, lFFTLength);
    Tfrm.rdft(lFFTLength, 1, fpSrc, lpSBitRevWork, fpCosSinTable);
    lMax = lFFTLength / 2 - 1;
    spDest[0].R = fpSrc[0];
    spDest[0].I = 0.0F;
    for (lLoopCntr = 1; lLoopCntr <= lMax; lLoopCntr++)
    {
        spDest[lLoopCntr].R = fpSrc[lLoopCntr * 2];
        spDest[lLoopCntr].I = fpSrc[lLoopCntr * 2 + 1];
    }
    spDest[lMax + 1].R = fpSrc[1];
    spDest[lMax + 1].I = 0.0F;
    #endif
}


void clDSPOp::FFTi(stpDCplx spDest, double *dpSrc)
{
    #if defined(DSP_IPP)
    ippsFFTFwd_RToCCS_64f(dpSrc, (Ipp64f *) spDest, 
        (IppsFFTSpec_R_64f *) vpDTfrm, DBitRevWork);
    #elif defined(DSP_USE_FFTW)
    Mul(dpSrc, dFFTScale, lFFTLength);
    fftw_execute_dft_r2c(fftwpDPlan, dpSrc, (fftw_complex *) spDest);
    #else
    long lLoopCntr;
    long lMax;

    Mul(dpSrc, dFFTScale, lFFTLength);
    Tfrm.rdft(lFFTLength, 1, dpSrc, lpDBitRevWork, dpCosSinTable);
    lMax = lFFTLength / 2 - 1;
    spDest[0].R = dpSrc[0];
    spDest[0].I = 0.0;
    for (lLoopCntr = 1; lLoopCntr <= lMax; lLoopCntr++)
    {
        spDest[lLoopCntr].R = dpSrc[lLoopCntr * 2];
        spDest[lLoopCntr].I = dpSrc[lLoopCntr * 2 + 1];
    }
    spDest[lMax + 1].R = dpSrc[1];
    spDest[lMax + 1].I = 0.0;
    #endif
}


void clDSPOp::FFTo(stpSCplx spDest, const float *fpSrc)
{
    #if defined(DSP_IPP)
    ippsFFTFwd_RToCCS_32f(fpSrc, (Ipp32f *) spDest, 
        (IppsFFTSpec_R_32f *) vpSTfrm, SBitRevWork);
    #elif defined(DSP_USE_FFTW)
    float *fpFFTBuf;

    fpFFTBuf = FFTBuf;
    Mul(fpFFTBuf, fpSrc, fFFTScale, lFFTLength);
    fftwf_execute_dft_r2c(fftwpSPlan, fpFFTBuf, (fftwf_complex *) spDest);
    #else
    long lLoopCntr;
    long lMax;
    float *fpFFTBuf;

    fpFFTBuf = FFTBuf;
    Mul(fpFFTBuf, fpSrc, fFFTScale, lFFTLength);
    Tfrm.rdft(lFFTLength, 1, fpFFTBuf, lpSBitRevWork, fpCosSinTable);
    lMax = lFFTLength / 2 - 1;
    spDest[0].R = fpFFTBuf[0];
    spDest[0].I = 0.0F;
    for (lLoopCntr = 1; lLoopCntr <= lMax; lLoopCntr++)
    {
        spDest[lLoopCntr].R = fpFFTBuf[lLoopCntr * 2];
        spDest[lLoopCntr].I = fpFFTBuf[lLoopCntr * 2 + 1];
    }
    spDest[lMax + 1].R = fpFFTBuf[1];
    spDest[lMax + 1].I = 0.0F;
    #endif
}


void clDSPOp::FFTo(stpDCplx spDest, const double *dpSrc)
{
    #if defined(DSP_IPP)
    ippsFFTFwd_RToCCS_64f(dpSrc, (Ipp64f *) spDest, 
        (IppsFFTSpec_R_64f *) vpDTfrm, DBitRevWork);
    #elif defined(DSP_USE_FFTW)
    double *dpFFTBuf;

    dpFFTBuf = FFTBuf;
    Mul(dpFFTBuf, dpSrc, dFFTScale, lFFTLength);
    fftw_execute_dft_r2c(fftwpDPlan, dpFFTBuf, (fftw_complex *) spDest);
    #else
    long lLoopCntr;
    long lMax;
    double *dpFFTBuf;

    dpFFTBuf = FFTBuf;
    Mul(dpFFTBuf, dpSrc, dFFTScale, lFFTLength);
    Tfrm.rdft(lFFTLength, 1, dpFFTBuf, lpDBitRevWork, dpCosSinTable);
    lMax = lFFTLength / 2 - 1;
    spDest[0].R = dpFFTBuf[0];
    spDest[0].I = 0.0;
    for (lLoopCntr = 1; lLoopCntr <= lMax; lLoopCntr++)
    {
        spDest[lLoopCntr].R = dpFFTBuf[lLoopCntr * 2];
        spDest[lLoopCntr].I = dpFFTBuf[lLoopCntr * 2 + 1];
    }
    spDest[lMax + 1].R = dpFFTBuf[1];
    spDest[lMax + 1].I = 0.0;
    #endif
}


void clDSPOp::FFTo(stpSCplx spDest, const stpSCplx spSrc)
{
    #if defined(DSP_IPP)
    ippsFFTFwd_CToC_32fc((Ipp32fc *) spSrc, (Ipp32fc *) spDest,
        (IppsFFTSpec_C_32fc *) vpSTfrm, SBitRevWork);
    #elif defined(DSP_USE_FFTW)
    Mul((float *) FFTBuf, (const float *) spSrc, fFFTScale, lFFTLength * 2);
    fftwf_execute_dft(fftwpSPlan, FFTBuf, (fftwf_complex *) spDest);
    #else
    long lLoopCntr;
    float *fpFFTBuf;

    fpFFTBuf = FFTBuf;
    for (lLoopCntr = 0; lLoopCntr < lFFTLength; lLoopCntr++)
    {
        fpFFTBuf[lLoopCntr * 2] = spSrc[lLoopCntr].R;
        fpFFTBuf[lLoopCntr * 2 + 1] = spSrc[lLoopCntr].I;
    }
    Mul(fpFFTBuf, fFFTScale, lFFTLength * 2);
    Tfrm.cdft(lFFTLength * 2, 1, fpFFTBuf, lpSBitRevWork, fpCosSinTable);
    for (lLoopCntr = 0; lLoopCntr < lFFTLength; lLoopCntr++)
    {
        spDest[lLoopCntr].R = fpFFTBuf[lLoopCntr * 2];
        spDest[lLoopCntr].I = fpFFTBuf[lLoopCntr * 2 + 1];
    }
    #endif
}


void clDSPOp::FFTo(stpDCplx spDest, const stpDCplx spSrc)
{
    #if defined(DSP_IPP)
    ippsFFTFwd_CToC_64fc((Ipp64fc *) spSrc, (Ipp64fc *) spDest,
        (IppsFFTSpec_C_64fc *) vpDTfrm, DBitRevWork);
    #elif defined(DSP_USE_FFTW)
    Mul((double *) FFTBuf, (const double *) spSrc, dFFTScale, lFFTLength * 2);
    fftw_execute_dft(fftwpDPlan, FFTBuf, (fftw_complex *) spDest);
    #else
    long lLoopCntr;
    double *dpFFTBuf;

    dpFFTBuf = FFTBuf;
    for (lLoopCntr = 0; lLoopCntr < lFFTLength; lLoopCntr++)
    {
        dpFFTBuf[lLoopCntr * 2] = spSrc[lLoopCntr].R;
        dpFFTBuf[lLoopCntr * 2 + 1] = spSrc[lLoopCntr].I;
    }
    Mul(dpFFTBuf, dFFTScale, lFFTLength * 2);
    Tfrm.cdft(lFFTLength * 2, 1, dpFFTBuf, lpDBitRevWork, dpCosSinTable);
    for (lLoopCntr = 0; lLoopCntr < lFFTLength; lLoopCntr++)
    {
        spDest[lLoopCntr].R = dpFFTBuf[lLoopCntr * 2];
        spDest[lLoopCntr].I = dpFFTBuf[lLoopCntr * 2 + 1];
    }
    #endif
}


void clDSPOp::IFFTo(float *fpDest, const stpSCplx spSrc)
{
    #if defined(DSP_IPP)
    ippsFFTInv_CCSToR_32f((Ipp32f *) spSrc, fpDest, 
        (IppsFFTSpec_R_32f *) vpSTfrm, SBitRevWork);
    #elif defined(DSP_USE_FFTW)
    fftwf_execute_dft_c2r(fftwpSIPlan, (fftwf_complex *) spSrc, fpDest);
    #else
    long lLoopCntr;
    long lMax;

    lMax = lFFTLength / 2 - 1;
    fpDest[0] = spSrc[0].R;
    for (lLoopCntr = 1; lLoopCntr <= lMax; lLoopCntr++)
    {
        fpDest[lLoopCntr * 2] = spSrc[lLoopCntr].R;
        fpDest[lLoopCntr * 2 + 1] = spSrc[lLoopCntr].I;
    }
    fpDest[1] = spSrc[lMax + 1].R;
    Tfrm.rdft(lFFTLength, -1, fpDest, lpSBitRevWork, fpCosSinTable);
    #endif
}


void clDSPOp::IFFTo(double *dpDest, const stpDCplx spSrc)
{
    #if defined(DSP_IPP)
    ippsFFTInv_CCSToR_64f((Ipp64f *) spSrc, dpDest, 
        (IppsFFTSpec_R_64f *) vpDTfrm, DBitRevWork);
    #elif defined(DSP_USE_FFTW)
    fftw_execute_dft_c2r(fftwpDIPlan, (fftw_complex *) spSrc, dpDest);
    #else
    long lLoopCntr;
    long lMax;

    lMax = lFFTLength / 2 - 1;
    dpDest[0] = spSrc[0].R;
    for (lLoopCntr = 1; lLoopCntr <= lMax; lLoopCntr++)
    {
        dpDest[lLoopCntr * 2] = spSrc[lLoopCntr].R;
        dpDest[lLoopCntr * 2 + 1] = spSrc[lLoopCntr].I;
    }
    dpDest[1] = spSrc[lMax + 1].R;
    Tfrm.rdft(lFFTLength, -1, dpDest, lpDBitRevWork, dpCosSinTable);
    #endif
}


void clDSPOp::IFFTo(stpSCplx spDest, const stpSCplx spSrc)
{
    #if defined(DSP_IPP)
    ippsFFTInv_CToC_32fc((Ipp32fc *) spSrc, (Ipp32fc *) spDest,
        (IppsFFTSpec_C_32fc *) vpSTfrm, SBitRevWork);
    #elif defined(DSP_USE_FFTW)
    fftwf_execute_dft(fftwpSIPlan, (fftwf_complex *) spSrc,
        (fftwf_complex *) spDest);
    #else
    long lLoopCntr;
    float fScale;
    float *fpFFTBuf;

    fpFFTBuf = FFTBuf;
    for (lLoopCntr = 0; lLoopCntr < lFFTLength; lLoopCntr++)
    {
        fpFFTBuf[lLoopCntr * 2] = spSrc[lLoopCntr].R;
        fpFFTBuf[lLoopCntr * 2 + 1] = spSrc[lLoopCntr].I;
    }
    Tfrm.cdft(lFFTLength * 2, -1, fpFFTBuf, lpSBitRevWork, fpCosSinTable);
    fScale = 1.0F / (float) lFFTLength;
    for (lLoopCntr = 0; lLoopCntr < lFFTLength; lLoopCntr++)
    {
        spDest[lLoopCntr].R = fpFFTBuf[lLoopCntr * 2] * fScale;
        spDest[lLoopCntr].I = fpFFTBuf[lLoopCntr * 2 + 1] * fScale;
    }
    #endif
}


void clDSPOp::IFFTo(stpDCplx spDest, const stpDCplx spSrc)
{
    #if defined(DSP_IPP)
    ippsFFTInv_CToC_64fc((Ipp64fc *) spSrc, (Ipp64fc *) spDest,
        (IppsFFTSpec_C_64fc *) vpDTfrm, DBitRevWork);
    #elif defined(DSP_USE_FFTW)
    fftw_execute_dft(fftwpDIPlan, (fftw_complex *) spSrc,
        (fftw_complex *) spDest);
    #else
    long lLoopCntr;
    double dScale;
    double *dpFFTBuf;

    dpFFTBuf = FFTBuf;
    for (lLoopCntr = 0; lLoopCntr < lFFTLength; lLoopCntr++)
    {
        dpFFTBuf[lLoopCntr * 2] = spSrc[lLoopCntr].R;
        dpFFTBuf[lLoopCntr * 2 + 1] = spSrc[lLoopCntr].I;
    }
    Tfrm.cdft(lFFTLength * 2, -1, dpFFTBuf, lpDBitRevWork, dpCosSinTable);
    dScale = 1.0 / (double) lFFTLength;
    for (lLoopCntr = 0; lLoopCntr < lFFTLength; lLoopCntr++)
    {
        spDest[lLoopCntr].R = dpFFTBuf[lLoopCntr * 2] * dScale;
        spDest[lLoopCntr].I = dpFFTBuf[lLoopCntr * 2 + 1] * dScale;
    }
    #endif
}


void clDSPOp::DCTInitialize (long lSize)
{
    float *fpFFTBuf;
    double *dpFFTBuf;
    clDSPAlloc FFTBufS(lSize * sizeof(float));
    clDSPAlloc FFTBufD(lSize * sizeof(double));

    fpFFTBuf = FFTBufS;
    dpFFTBuf = FFTBufD;
    fFFTScale = 2.0F / (float) lSize;
    dFFTScale = 2.0 / (double) lSize;
    lpSBitRevWork = (long *)
        SBitRevWork.Size(((size_t) ceil(2.0 + sqrt(lSize / 2))) * sizeof(long));
    lpDBitRevWork = (long *)
        DBitRevWork.Size(((size_t) ceil(2.0 + sqrt(lSize / 2))) * sizeof(long));
    fpCosSinTable = (float *) 
        SCosSinTable.Size((lSize * 5 / 4) * sizeof(float));
    dpCosSinTable = (double *) 
        DCosSinTable.Size((lSize * 5 / 4) * sizeof(double));
    lpSBitRevWork[0] = 0;
    lpSBitRevWork[1] = 0;
    lpDBitRevWork[0] = 0;
    lpDBitRevWork[1] = 0;
    Tfrm.ddct(lSize, 1, fpFFTBuf, lpSBitRevWork, fpCosSinTable);
    Tfrm.ddct(lSize, 1, dpFFTBuf, lpDBitRevWork, dpCosSinTable);
    FFTBuf.Size(lSize * sizeof(double));

    bFFTInitialized = true;
    lFFTLength = lSize;
}


void clDSPOp::DCTUninitialize ()
{
    SBitRevWork.Free();
    DBitRevWork.Free();
    SCosSinTable.Free();
    DCosSinTable.Free();
    FFTBuf.Free();
    bFFTInitialized = false;
}


void clDSPOp::DCTi (float *fpVect)
{
    Tfrm.ddct(lFFTLength, 1, fpVect, lpSBitRevWork, fpCosSinTable);
}


void clDSPOp::DCTi (double *dpVect)
{
    Tfrm.ddct(lFFTLength, 1, dpVect, lpDBitRevWork, dpCosSinTable);
}


void clDSPOp::DCTo (float *fpDest, const float *fpSrc)
{
    Copy(fpDest, fpSrc, lFFTLength);
    Tfrm.ddct(lFFTLength, 1, fpDest, lpSBitRevWork, fpCosSinTable);
}


void clDSPOp::DCTo (double *dpDest, const double *dpSrc)
{
    Copy(dpDest, dpSrc, lFFTLength);
    Tfrm.ddct(lFFTLength, 1, dpDest, lpDBitRevWork, dpCosSinTable);
}


void clDSPOp::IDCTi (float *fpVect)
{
    fpVect[0] *= 0.5F;
    Mul(&fpVect[1], fFFTScale, lFFTLength - 1);
    Tfrm.ddct(lFFTLength, -1, fpVect, lpSBitRevWork, fpCosSinTable);
}


void clDSPOp::IDCTi (double *dpVect)
{
    dpVect[0] *= 0.5;
    Mul(&dpVect[1], dFFTScale, lFFTLength - 1);
    Tfrm.ddct(lFFTLength, -1, dpVect, lpDBitRevWork, dpCosSinTable);
}


void clDSPOp::IDCTo (float *fpDest, const float *fpSrc)
{
    fpDest[0] = fpSrc[0] * 0.5F;
    Mul(&fpDest[1], &fpSrc[1], fFFTScale, lFFTLength - 1);
    Tfrm.ddct(lFFTLength, -1, fpDest, lpSBitRevWork, fpCosSinTable);
}


void clDSPOp::IDCTo (double *dpDest, const double *dpSrc)
{
    dpDest[0] = dpSrc[0] * 0.5;
    Mul(&dpDest[1], &dpSrc[1], dFFTScale, lFFTLength - 1);
    Tfrm.ddct(lFFTLength, -1, dpDest, lpDBitRevWork, dpCosSinTable);
}
        

void clDSPOp::FFTWConvert(stpSCplx spDest, const float *fpSrc, long lLength)
{
    long lLoopCntr;
    long lMax;

    lMax = lLength / 2;
    spDest[0].R = fpSrc[0];
    spDest[0].I = 0.0F;
    for (lLoopCntr = 1; lLoopCntr < lMax; lLoopCntr++)
    {
        spDest[lLoopCntr].R = fpSrc[lLoopCntr];
        spDest[lLoopCntr].I = fpSrc[lLength - lLoopCntr];
    }
    spDest[lMax].R = fpSrc[lMax];
    spDest[lMax].I = 0.0F;
}


void clDSPOp::FFTWConvert(stpDCplx spDest, const float *fpSrc, long lLength)
{
    long lLoopCntr;
    long lMax;

    lMax = lLength / 2;
    spDest[0].R = fpSrc[0];
    spDest[0].I = 0.0;
    for (lLoopCntr = 1; lLoopCntr < lMax; lLoopCntr++)
    {
        spDest[lLoopCntr].R = fpSrc[lLoopCntr];
        spDest[lLoopCntr].I = fpSrc[lLength - lLoopCntr];
    }
    spDest[lMax].R = fpSrc[lMax];
    spDest[lMax].I = 0.0;
}


void clDSPOp::FFTWConvert(stpSCplx spDest, const double *dpSrc, long lLength)
{
    long lLoopCntr;
    long lMax;

    lMax = lLength / 2;
    spDest[0].R = dpSrc[0];
    spDest[0].I = 0.0F;
    for (lLoopCntr = 1; lLoopCntr < lMax; lLoopCntr++)
    {
        spDest[lLoopCntr].R = dpSrc[lLoopCntr];
        spDest[lLoopCntr].I = dpSrc[lLength - lLoopCntr];
    }
    spDest[lMax].R = dpSrc[lMax];
    spDest[lMax].I = 0.0F;
}


void clDSPOp::FFTWConvert(stpDCplx spDest, const double *dpSrc, long lLength)
{
    long lLoopCntr;
    long lMax;

    lMax = lLength / 2;
    spDest[0].R = dpSrc[0];
    spDest[0].I = 0.0;
    for (lLoopCntr = 1; lLoopCntr < lMax; lLoopCntr++)
    {
        spDest[lLoopCntr].R = dpSrc[lLoopCntr];
        spDest[lLoopCntr].I = dpSrc[lLength - lLoopCntr];
    }
    spDest[lMax].R = dpSrc[lMax];
    spDest[lMax].I = 0.0;
}


void clDSPOp::FFTWConvert(float *fpDest, const stpSCplx spSrc, long lLength)
{
    long lLoopCntr;
    long lMax;

    lMax = lLength / 2;
    fpDest[0] = spSrc[0].R;
    for (lLoopCntr = 1; lLoopCntr < lMax; lLoopCntr++)
    {
        fpDest[lLoopCntr] = spSrc[lLoopCntr].R;
        fpDest[lLength - lLoopCntr] = spSrc[lLoopCntr].I;
    }
    fpDest[lMax] = spSrc[lMax].R;
}


void clDSPOp::FFTWConvert(float *fpDest, const stpDCplx spSrc, long lLength)
{
    long lLoopCntr;
    long lMax;

    lMax = lLength / 2;
    fpDest[0] = spSrc[0].R;
    for (lLoopCntr = 1; lLoopCntr < lMax; lLoopCntr++)
    {
        fpDest[lLoopCntr] = spSrc[lLoopCntr].R;
        fpDest[lLength - lLoopCntr] = spSrc[lLoopCntr].I;
    }
    fpDest[lMax] = spSrc[lMax].R;
}


void clDSPOp::FFTWConvert(double *dpDest, const stpSCplx spSrc, long lLength)
{
    long lLoopCntr;
    long lMax;

    lMax = lLength / 2;
    dpDest[0] = spSrc[0].R;
    for (lLoopCntr = 1; lLoopCntr < lMax; lLoopCntr++)
    {
        dpDest[lLoopCntr] = spSrc[lLoopCntr].R;
        dpDest[lLength - lLoopCntr] = spSrc[lLoopCntr].I;
    }
    dpDest[lMax] = spSrc[lMax].R;
}


void clDSPOp::FFTWConvert(double *dpDest, const stpDCplx spSrc, long lLength)
{
    long lLoopCntr;
    long lMax;

    lMax = lLength / 2;
    dpDest[0] = spSrc[0].R;
    for (lLoopCntr = 1; lLoopCntr < lMax; lLoopCntr++)
    {
        dpDest[lLoopCntr] = spSrc[lLoopCntr].R;
        dpDest[lLength - lLoopCntr] = spSrc[lLoopCntr].I;
    }
    dpDest[lMax] = spSrc[lMax].R;
}


/*INLINE double DegToRad (double dSrc)
{
   return (M_PI / 180.0) * dSrc;
}*/


/*INLINE double RadToDeg (double dSrc)
{
   return (180.0 / M_PI) * dSrc;
}*/

