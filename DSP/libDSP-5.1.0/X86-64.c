/*

    x86-64 specific optimized assembler dsp routines
    Copyright (C) 2003-2004 Jussi Laako

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


#if (defined(DSP_X86_64) && (!defined(DSP_X86)))

#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <float.h>

#include "dsp/X86-64.h"


#ifdef __cplusplus
extern "C"
{
#endif


void dsp_x86_64_minmaxf (float *fpMin, float *fpMax, const float *fpSrc, 
    long lDataLength)
{
    long lDataCntr;

    *fpMin = FLT_MAX;
    *fpMax = -FLT_MAX;
    X86_64_ASM (
        "movss %0, %%xmm0\n\t" \
        "movss %1, %%xmm1\n\t"
        :
        : "m" (*fpMin),
          "m" (*fpMax)
        : "xmm0", "xmm1", "memory");
    for (lDataCntr = 0; lDataCntr < lDataLength; lDataCntr++)
    {
        X86_64_ASM (
            "movss %0, %%xmm2\n\t" \
            "minss %%xmm2, %%xmm0\n\t" \
            "maxss %%xmm2, %%xmm1\n\t"
            :
            : "m" (fpSrc[lDataCntr])
            : "xmm0", "xmm1", "xmm2", "memory");
    }
    X86_64_ASM (
        "movss %%xmm0, %0\n\t" \
        "movss %%xmm1, %1\n\t"
        : "=m" (*fpMin),
          "=m" (*fpMax)
        :
        : "xmm0", "xmm1", "memory");
}


void dsp_x86_64_minmax (double *dpMin, double *dpMax, const double *dpSrc, 
    long lDataLength)
{
    long lDataCntr;

    *dpMin = FLT_MAX;
    *dpMax = -FLT_MAX;
    X86_64_ASM (
        "movsd %0, %%xmm0\n\t" \
        "movsd %1, %%xmm1\n\t"
        :
        : "m" (*dpMin),
          "m" (*dpMax)
        : "xmm0", "xmm1", "memory");
    for (lDataCntr = 0; lDataCntr < lDataLength; lDataCntr++)
    {
        X86_64_ASM (
            "movsd %0, %%xmm2\n\t" \
            "minsd %%xmm2, %%xmm0\n\t" \
            "maxsd %%xmm2, %%xmm1\n\t"
            :
            : "m" (dpSrc[lDataCntr])
            : "xmm0", "xmm1", "xmm2", "memory");
    }
    X86_64_ASM (
        "movss %%xmm0, %0\n\t" \
        "movss %%xmm1, %1\n\t"
        : "=m" (*dpMin),
          "=m" (*dpMax)
        :
        : "xmm0", "xmm1", "memory");
}


float dsp_x86_64_crosscorrf (const float *fpSrc1, const float *fpSrc2,
    long lDataLength)
{
    long lDataCntr;
    float fScale;
    float fNormFact;
    float fProdSum;
    float fSqSum1;
    float fSqSum2;
    float fRes;
    
    X86_64_ASM (
        "xorps %%xmm0, %%xmm0\n\t" \
        "xorps %%xmm1, %%xmm1\n\t" \
        "xorps %%xmm2, %%xmm2\n\t"
        :
        :
        : "xmm0", "xmm1", "xmm2");
    for (lDataCntr = 0; lDataCntr < lDataLength; lDataCntr++)
    {
        X86_64_ASM (
            "movss %3, %%xmm3\n\t" \
            "movss %4, %%xmm4\n\t" \
            \
            "movss %%xmm4, %%xmm5\n\t" \
            "mulss %%xmm3, %%xmm5\n\t" \
            "addss %%xmm5, %%xmm0\n\t" \
            \
            "movss %%xmm3, %%xmm5\n\t" \
            "mulss %%xmm3, %%xmm5\n\t" \
            "addss %%xmm5, %%xmm1\n\t" \
            \
            "movss %%xmm4, %%xmm5\n\t" \
            "mulss %%xmm4, %%xmm5\n\t" \
            "addss %%xmm5, %%xmm2\n\t" \
            \
            "movss %%xmm0, %0\n\t" \
            "movss %%xmm1, %1\n\t" \
            "movss %%xmm2, %2\n\t"
            : "=m" (fProdSum),
              "=m" (fSqSum1),
              "=m" (fSqSum2)
            : "m" (fpSrc1[lDataCntr]),
              "m" (fpSrc2[lDataCntr])
            : "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5");
    }
    fScale = 1.0F / lDataLength;
    fNormFact = sqrtf(fSqSum1 * fSqSum2) * fScale;
    fRes = (fProdSum * fScale) / fNormFact;

    return fRes;
}


double dsp_x86_64_crosscorr (const double *dpSrc1, const double *dpSrc2,
    long lDataLength)
{
    long lDataCntr;
    double dScale;
    double dNormFact;
    double dProdSum;
    double dSqSum1;
    double dSqSum2;
    double dRes;
    
    X86_64_ASM (
        "xorpd %%xmm0, %%xmm0\n\t" \
        "xorpd %%xmm1, %%xmm1\n\t" \
        "xorpd %%xmm2, %%xmm2\n\t"
        :
        :
        : "xmm0", "xmm1", "xmm2");
    for (lDataCntr = 0; lDataCntr < lDataLength; lDataCntr++)
    {
        X86_64_ASM (
            "movsd %3, %%xmm3\n\t" \
            "movsd %4, %%xmm4\n\t" \
            \
            "movsd %%xmm4, %%xmm5\n\t" \
            "mulsd %%xmm3, %%xmm5\n\t" \
            "addsd %%xmm5, %%xmm0\n\t" \
            \
            "movsd %%xmm3, %%xmm5\n\t" \
            "mulsd %%xmm3, %%xmm5\n\t" \
            "addsd %%xmm5, %%xmm1\n\t" \
            \
            "movsd %%xmm4, %%xmm5\n\t" \
            "mulsd %%xmm4, %%xmm5\n\t" \
            "addsd %%xmm5, %%xmm2\n\t" \
            \
            "movsd %%xmm0, %0\n\t" \
            "movsd %%xmm1, %1\n\t" \
            "movsd %%xmm2, %2\n\t"
            : "=m" (dProdSum),
              "=m" (dSqSum1),
              "=m" (dSqSum2)
            : "m" (dpSrc1[lDataCntr]),
              "m" (dpSrc2[lDataCntr])
            : "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5");
    }
    dScale = 1.0 / lDataLength;
    dNormFact = sqrt(dSqSum1 * dSqSum2) * dScale;
    dRes = (dProdSum * dScale) / dNormFact;

    return dRes;
}


void dsp_x86_64_firf (float *fpDest, const float *fpSrc, long lDataLength, 
    const float *fpCoeff, long lCoeffLength)
{
    long lDestCntr;
    long lSrcCntr;
    long lCoeffCntr;
    long lSrcCount;

    lDestCntr = 0;
    lSrcCount = lDataLength + lCoeffLength;
    for (lSrcCntr = lCoeffLength;
        lSrcCntr < lSrcCount;
        lSrcCntr++)
    {
        X86_64_ASM (
            "xorps %%xmm0, %%xmm0\n\t"
            :
            :
            : "xmm0");
        for (lCoeffCntr = 0;
            lCoeffCntr < lCoeffLength;
            lCoeffCntr++)
        {
            X86_64_ASM (
                "movss %0, %%xmm1\n\t"
                "mulss %1, %%xmm1\n\t"
                "addss %%xmm1, %%xmm0\n\t"
                :
                : "m" (fpSrc[lSrcCntr - lCoeffCntr]),
                  "m" (fpCoeff[lCoeffCntr])
                : "xmm0", "xmm1", "memory");
        }
        X86_64_ASM (
            "movss %%xmm0, %0\n\t"
            : "=m" (fpDest[lDestCntr++])
            :
            : "xmm0", "memory");
    }
}


void dsp_x86_64_fir (double *dpDest, const double *dpSrc, long lDataLength, 
    const double *dpCoeff, long lCoeffLength)
{
    long lDestCntr;
    long lSrcCntr;
    long lCoeffCntr;
    long lSrcCount;

    lDestCntr = 0;
    lSrcCount = lDataLength + lCoeffLength;
    for (lSrcCntr = lCoeffLength;
        lSrcCntr < lSrcCount;
        lSrcCntr++)
    {
        X86_64_ASM (
            "xorpd %%xmm0, %%xmm0\n\t"
            :
            :
            : "xmm0");
        for (lCoeffCntr = 0;
            lCoeffCntr < lCoeffLength;
            lCoeffCntr++)
        {
            X86_64_ASM (
                "movsd %0, %%xmm1\n\t"
                "mulsd %1, %%xmm1\n\t"
                "addsd %%xmm1, %%xmm0\n\t"
                :
                : "m" (dpSrc[lSrcCntr - lCoeffCntr]),
                  "m" (dpCoeff[lCoeffCntr])
                : "xmm0", "xmm1", "memory");
        }
        X86_64_ASM (
            "movsd %%xmm0, %0\n\t"
            : "=m" (dpDest[lDestCntr++])
            :
            : "xmm0", "memory");
    }
}


void dsp_x86_64_iirf (float *fpVect, long lDataLength, const float *fpCoeff, 
    float *fpX, float *fpY)
{
    long lDataCntr;
    
    X86_64_ASM (
        "movss %0, %%xmm1\n\t" \
        "movss %1, %%xmm2\n\t" \
        "movss %2, %%xmm3\n\t" \
        "movss %3, %%xmm4\n\t" \
        "prefetchnta %4\n\t"
        :
        : "m" (fpX[1]),
          "m" (fpX[2]),
          "m" (fpY[0]),
          "m" (fpY[1]),
          "m" (fpCoeff[0])
        : "xmm1", "xmm2", "xmm3", "xmm4", "memory");
    for (lDataCntr = 0; lDataCntr < lDataLength; lDataCntr++)
    {
        X86_64_ASM (
            "movss %%xmm1, %%xmm0\n\t" \
            "movss %%xmm2, %%xmm1\n\t" \
            "movss %1, %%xmm2\n\t" \
            \
            "movss %2, %%xmm5\n\t" \
            "mulss %%xmm2, %%xmm5\n\t" \
            "movss %3, %%xmm6\n\t" \
            "mulss %%xmm1, %%xmm6\n\t" \
            "addss %%xmm6, %%xmm5\n\t" \
            "movss %4, %%xmm6\n\t" \
            "mulss %%xmm0, %%xmm6\n\t" \
            "addss %%xmm6, %%xmm5\n\t" \
            \
            "movss %5, %%xmm6\n\t" \
            "mulss %%xmm4, %%xmm6\n\t" \
            "movss %6, %%xmm7\n\t" \
            "mulss %%xmm3, %%xmm7\n\t" \
            "addss %%xmm7, %%xmm6\n\t" \
            \
            "addss %%xmm5, %%xmm6\n\t" \
            "movss %%xmm4, %%xmm3\n\t" \
            "movss %%xmm6, %%xmm4\n\t" \
            \
            "movss %%xmm6, %0\n\t"
            : "=m" (fpVect[lDataCntr])
            : "m0" (fpVect[lDataCntr]),
              "m" (fpCoeff[0]),
              "m" (fpCoeff[1]),
              "m" (fpCoeff[2]),
              "m" (fpCoeff[3]),
              "m" (fpCoeff[4])
            : "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7",
              "memory");
    }
    X86_64_ASM (
        "movss %%xmm0, %0\n\t" \
        "movss %%xmm1, %1\n\t" \
        "movss %%xmm2, %2\n\t" \
        "movss %%xmm3, %3\n\t" \
        "movss %%xmm4, %4\n\t"
        : "=m" (fpX[0]),
          "=m" (fpX[1]),
          "=m" (fpX[2]),
          "=m" (fpY[0]),
          "=m" (fpY[1])
        :
        : "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "memory");
}


void dsp_x86_64_iir (double *dpVect, long lDataLength, const double *dpCoeff, 
    double *dpX, double *dpY)
{
    long lDataCntr;
    
    X86_64_ASM (
        "movsd %0, %%xmm1\n\t" \
        "movsd %1, %%xmm2\n\t" \
        "movsd %2, %%xmm3\n\t" \
        "movsd %3, %%xmm4\n\t" \
        "prefetchnta %4\n\t" \
        "prefetchnta %5\n\t"
        :
        : "m" (dpX[1]),
          "m" (dpX[2]),
          "m" (dpY[0]),
          "m" (dpY[1]),
          "m" (dpCoeff[0]),
          "m" (dpCoeff[3])
        : "xmm1", "xmm2", "xmm3", "xmm4", "memory");
    for (lDataCntr = 0; lDataCntr < lDataLength; lDataCntr++)
    {
        X86_64_ASM (
            "movsd %%xmm1, %%xmm0\n\t" \
            "movsd %%xmm2, %%xmm1\n\t" \
            "movsd %1, %%xmm2\n\t" \
            \
            "movsd %2, %%xmm5\n\t" \
            "mulsd %%xmm2, %%xmm5\n\t" \
            "movsd %3, %%xmm6\n\t" \
            "mulsd %%xmm1, %%xmm6\n\t" \
            "addsd %%xmm6, %%xmm5\n\t" \
            "movsd %4, %%xmm6\n\t" \
            "mulsd %%xmm0, %%xmm6\n\t" \
            "addsd %%xmm6, %%xmm5\n\t" \
            \
            "movsd %5, %%xmm6\n\t" \
            "mulsd %%xmm4, %%xmm6\n\t" \
            "movsd %6, %%xmm7\n\t" \
            "mulsd %%xmm3, %%xmm7\n\t" \
            "addsd %%xmm7, %%xmm6\n\t" \
            \
            "addsd %%xmm5, %%xmm6\n\t" \
            "movsd %%xmm4, %%xmm3\n\t" \
            "movsd %%xmm6, %%xmm4\n\t" \
            \
            "movsd %%xmm6, %0\n\t"
            : "=m" (dpVect[lDataCntr])
            : "m0" (dpVect[lDataCntr]),
              "m" (dpCoeff[0]),
              "m" (dpCoeff[1]),
              "m" (dpCoeff[2]),
              "m" (dpCoeff[3]),
              "m" (dpCoeff[4])
            : "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7",
              "memory");
    }
    X86_64_ASM (
        "movsd %%xmm0, %0\n\t" \
        "movsd %%xmm1, %1\n\t" \
        "movsd %%xmm2, %2\n\t" \
        "movsd %%xmm3, %3\n\t" \
        "movsd %%xmm4, %4\n\t"
        : "=m" (dpX[0]),
          "=m" (dpX[1]),
          "=m" (dpX[2]),
          "=m" (dpY[0]),
          "=m" (dpY[1])
        :
        : "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "memory");
}


void dsp_x86_64_iirf_nip (float *fpDest, const float *fpSrc, long lDataLength, 
    const float *fpCoeff, float *fpX, float *fpY)
{
    long lDataCntr;
    
    X86_64_ASM (
        "movss %0, %%xmm1\n\t" \
        "movss %1, %%xmm2\n\t" \
        "movss %2, %%xmm3\n\t" \
        "movss %3, %%xmm4\n\t" \
        "prefetchnta %4\n\t"
        :
        : "m" (fpX[1]),
          "m" (fpX[2]),
          "m" (fpY[0]),
          "m" (fpY[1]),
          "m" (fpCoeff[0])
        : "xmm1", "xmm2", "xmm3", "xmm4", "memory");
    for (lDataCntr = 0; lDataCntr < lDataLength; lDataCntr++)
    {
        X86_64_ASM (
            "movss %%xmm1, %%xmm0\n\t" \
            "movss %%xmm2, %%xmm1\n\t" \
            "movss %1, %%xmm2\n\t" \
            \
            "movss %2, %%xmm5\n\t" \
            "mulss %%xmm2, %%xmm5\n\t" \
            "movss %3, %%xmm6\n\t" \
            "mulss %%xmm1, %%xmm6\n\t" \
            "addss %%xmm6, %%xmm5\n\t" \
            "movss %4, %%xmm6\n\t" \
            "mulss %%xmm0, %%xmm6\n\t" \
            "addss %%xmm6, %%xmm5\n\t" \
            \
            "movss %5, %%xmm6\n\t" \
            "mulss %%xmm4, %%xmm6\n\t" \
            "movss %6, %%xmm7\n\t" \
            "mulss %%xmm3, %%xmm7\n\t" \
            "addss %%xmm7, %%xmm6\n\t" \
            \
            "addss %%xmm5, %%xmm6\n\t" \
            "movss %%xmm4, %%xmm3\n\t" \
            "movss %%xmm6, %%xmm4\n\t" \
            \
            "movss %%xmm6, %0\n\t"
            : "=m" (fpDest[lDataCntr])
            : "m" (fpSrc[lDataCntr]),
              "m" (fpCoeff[0]),
              "m" (fpCoeff[1]),
              "m" (fpCoeff[2]),
              "m" (fpCoeff[3]),
              "m" (fpCoeff[4])
            : "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7",
              "memory");
    }
    X86_64_ASM (
        "movss %%xmm0, %0\n\t" \
        "movss %%xmm1, %1\n\t" \
        "movss %%xmm2, %2\n\t" \
        "movss %%xmm3, %3\n\t" \
        "movss %%xmm4, %4\n\t"
        : "=m" (fpX[0]),
          "=m" (fpX[1]),
          "=m" (fpX[2]),
          "=m" (fpY[0]),
          "=m" (fpY[1])
        :
        : "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "memory");
}


void dsp_x86_64_iir_nip (double *dpDest, const double *dpSrc, long lDataLength, 
    const double *dpCoeff, double *dpX, double *dpY)
{
    long lDataCntr;
    
    X86_64_ASM (
        "movsd %0, %%xmm1\n\t" \
        "movsd %1, %%xmm2\n\t" \
        "movsd %2, %%xmm3\n\t" \
        "movsd %3, %%xmm4\n\t" \
        "prefetchnta %4\n\t" \
        "prefetchnta %5\n\t"
        :
        : "m" (dpX[1]),
          "m" (dpX[2]),
          "m" (dpY[0]),
          "m" (dpY[1]),
          "m" (dpCoeff[0]),
          "m" (dpCoeff[3])
        : "xmm1", "xmm2", "xmm3", "xmm4", "memory");
    for (lDataCntr = 0; lDataCntr < lDataLength; lDataCntr++)
    {
        X86_64_ASM (
            "movsd %%xmm1, %%xmm0\n\t" \
            "movsd %%xmm2, %%xmm1\n\t" \
            "movsd %1, %%xmm2\n\t" \
            \
            "movsd %2, %%xmm5\n\t" \
            "mulsd %%xmm2, %%xmm5\n\t" \
            "movsd %3, %%xmm6\n\t" \
            "mulsd %%xmm1, %%xmm6\n\t" \
            "addsd %%xmm6, %%xmm5\n\t" \
            "movsd %4, %%xmm6\n\t" \
            "mulsd %%xmm0, %%xmm6\n\t" \
            "addsd %%xmm6, %%xmm5\n\t" \
            \
            "movsd %5, %%xmm6\n\t" \
            "mulsd %%xmm4, %%xmm6\n\t" \
            "movsd %6, %%xmm7\n\t" \
            "mulsd %%xmm3, %%xmm7\n\t" \
            "addsd %%xmm7, %%xmm6\n\t" \
            \
            "addsd %%xmm5, %%xmm6\n\t" \
            "movsd %%xmm4, %%xmm3\n\t" \
            "movsd %%xmm6, %%xmm4\n\t" \
            \
            "movsd %%xmm6, %0\n\t"
            : "=m" (dpDest[lDataCntr])
            : "m" (dpSrc[lDataCntr]),
              "m" (dpCoeff[0]),
              "m" (dpCoeff[1]),
              "m" (dpCoeff[2]),
              "m" (dpCoeff[3]),
              "m" (dpCoeff[4])
            : "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7",
              "memory");
    }
    X86_64_ASM (
        "movsd %%xmm0, %0\n\t" \
        "movsd %%xmm1, %1\n\t" \
        "movsd %%xmm2, %2\n\t" \
        "movsd %%xmm3, %3\n\t" \
        "movsd %%xmm4, %4\n\t"
        : "=m" (dpX[0]),
          "=m" (dpX[1]),
          "=m" (dpX[2]),
          "=m" (dpY[0]),
          "=m" (dpY[1])
        :
        : "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "memory");
}


#ifdef __cplusplus
}
#endif

#endif
