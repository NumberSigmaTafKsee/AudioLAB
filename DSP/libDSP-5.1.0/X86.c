/*

    x86 specific optimized assembler dsp routines
    Copyright (C) 2001-2005 Jussi Laako

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


#ifdef DSP_X86


#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <float.h>

#include "dsp/X86.h"


#ifndef DSP_X86_64
static char cpCPUid[13];
#endif


#ifdef __cplusplus
extern "C"
{
#endif


#ifndef DSP_X86_64
const char *dsp_x86_cpuid ()
{
    unsigned int *ipCPUid = (unsigned int *) cpCPUid;
    
    X86_ASM (
        "pushl %%ebx\n\t" \
        "xorl %%eax, %%eax\n\t" \
        "cpuid\n\t" \
        "movl %%ebx, %0\n\t" \
        "movl %%ecx, %2\n\t" \
        "movl %%edx, %1\n\t" \
        "popl %%ebx\n\t"
        : "=m" (ipCPUid[0]),
          "=m" (ipCPUid[1]),
          "=m" (ipCPUid[2])
        :
        : "eax", "ecx", "edx", "memory");
    cpCPUid[12] = '\0';

    return cpCPUid;
}


void dsp_x86_features (unsigned int *uipFeatures1, unsigned int *uipFeatures2)
{
    unsigned int uiF1;
    unsigned int uiF2;

    X86_ASM (
        "pushl %%ebx\n\t" \
        "movl $1, %%eax\n\t" \
        "cpuid\n\t" \
        "movl %%edx, %0\n\t" \
        "movl %%ecx, %1\n\t" \
        "popl %%ebx\n\t"
        : "=m" (uiF1),
          "=m" (uiF2)
        :
        : "eax", "ecx", "edx", "memory");
    *uipFeatures1 = uiF1;
    *uipFeatures2 = uiF2;
}


unsigned int dsp_x86_amd_features ()
{
    unsigned int uiFunction = 0x80000001;
    unsigned int uiFeatures = 0;
    
    X86_ASM (
        "pushl %%ebx\n\t" \
        "movl %1, %%eax\n\t" \
        "cpuid\n\t" \
        "movl %%edx, %0\n\t" \
        "popl %%ebx\n\t"
        : "=m" (uiFeatures)
        : "m" (uiFunction)
        : "eax", "ecx", "edx", "memory");
    
    return uiFeatures;
}
#endif


extern int dsp_x86_have_e3dnow ()
{
    #ifndef DSP_X86_64
    unsigned int uiExtSup = 0;
    unsigned int uiFeatures;

    X86_ASM (
        "pushl %%ebx\n\t" \
        "movl $0x80000000, %%eax\n\t" \
        "cpuid\n\t" \
        "cmpl $0x80000001, %%eax\n\t" \
        "jl have3dnowxit\n\t" \
        "movl $1, %0\n\t" \
        "have3dnowxit:\n\t" \
        "popl %%ebx\n\t"
        : "=m" (uiExtSup)
        :
        : "eax", "ecx", "edx", "memory");
    if (uiExtSup)
    {
        uiFeatures = dsp_x86_amd_features();
        if ((uiFeatures & (1 << 31)) && (uiFeatures & (1 << 30)))
            return 1;
    }
    return 0;
    #else
    return 1;
    #endif
}


extern int dsp_x86_have_sse ()
{
    #ifndef DSP_X86_64
    int iSSELevel = 0;
    unsigned int uiFeatures1;
    unsigned int uiFeatures2;
    
    dsp_x86_features(&uiFeatures1, &uiFeatures2);
    if (uiFeatures1 & (1 << 25))
    {
        iSSELevel++;
        if (uiFeatures1 & (1 << 26))
        {
            iSSELevel++;
            if (uiFeatures2 & 0x1)
                iSSELevel++;
        }
    }
    return iSSELevel;
    #else
    return 2;
    #endif
}


// --- inline code snippets


inline void dsp_x86_prefetchntf_init (const float *fpSrc)
{
    pv2sf m64pSrc = (pv2sf) fpSrc;

    X86_ASM (
        "prefetchnta %0\n\t" \
        "prefetchnta %1\n\t" \
        "prefetchnta %2\n\t" \
        "prefetchnta %3\n\t"
        :
        : "m" (m64pSrc[0]),
          "m" (m64pSrc[8]),
          "m" (m64pSrc[16]),
          "m" (m64pSrc[24]));
}


inline void dsp_x86_prefetchnt_init (const double *dpSrc)
{
    pv2sf m64pSrc = (pv2sf) dpSrc;

    X86_ASM (
        "prefetchnta %0\n\t" \
        "prefetchnta %1\n\t" \
        "prefetchnta %2\n\t" \
        "prefetchnta %3\n\t"
        :
        : "m" (m64pSrc[0]),
          "m" (m64pSrc[8]),
          "m" (m64pSrc[16]),
          "m" (m64pSrc[24]));
}


inline void dsp_x86_prefetchtf_init (const float *fpSrc)
{
    pv2sf m64pSrc = (pv2sf) fpSrc;

    X86_ASM (
        "prefetcht0 %0\n\t" \
        "prefetcht0 %1\n\t" \
        "prefetcht0 %2\n\t" \
        "prefetcht0 %3\n\t"
        :
        : "m" (m64pSrc[0]),
          "m" (m64pSrc[8]),
          "m" (m64pSrc[16]),
          "m" (m64pSrc[24]));
}


inline void dsp_x86_prefetcht_init (const double *dpSrc)
{
    pv2sf m64pSrc = (pv2sf) dpSrc;

    X86_ASM (
        "prefetcht0 %0\n\t" \
        "prefetcht0 %1\n\t" \
        "prefetcht0 %2\n\t" \
        "prefetcht0 %3\n\t"
        :
        : "m" (m64pSrc[0]),
          "m" (m64pSrc[8]),
          "m" (m64pSrc[16]),
          "m" (m64pSrc[24]));
}


inline void dsp_x86_prefetchntf_next (const float *fpSrc)
{
    pv2sf m64pSrc = (pv2sf) fpSrc;

    X86_ASM (
        "prefetchnta %0\n\t"
        :
        : "m" (m64pSrc[32]));
}


inline void dsp_x86_prefetchnt_next (const double *dpSrc)
{
    pv2sf m64pSrc = (pv2sf) dpSrc;

    X86_ASM (
        "prefetchnta %0\n\t"
        :
        : "m" (m64pSrc[32]));
}


inline void dsp_x86_prefetchtf_next (const float *fpSrc)
{
    pv2sf m64pSrc = (pv2sf) fpSrc;

    X86_ASM (
        "prefetcht0 %0\n\t"
        :
        : "m" (m64pSrc[32]));
}


inline void dsp_x86_prefetcht_next (const double *dpSrc)
{
    pv2sf m64pSrc = (pv2sf) dpSrc;

    X86_ASM (
        "prefetcht0 %0\n\t"
        :
        : "m" (m64pSrc[32]));
}


// ---


void dsp_x86_3dnow_copyf (float *fpDest, const float *fpSrc, int iDataLength)
{
    int iStartIdx;
    int iDataCntr;
    int iDataCount;
    pv2sf m64pDest = (pv2sf) fpDest;
    pv2sf m64pSrc = (pv2sf) fpSrc;
    
    iStartIdx = 0;
    X86_ASM (
        "prefetchnta %0\n\t" \
        "prefetchnta %1\n\t" \
        "prefetchnta %2\n\t" \
        "prefetchnta %3\n\t"
        :
        : "m" (m64pSrc[0]),
          "m" (m64pSrc[8]),
          "m" (m64pSrc[16]),
          "m" (m64pSrc[24]));
    iDataCount = ((iDataLength & 0xfffffff0) >> 1);
    for (iDataCntr = iStartIdx; iDataCntr < iDataCount; iDataCntr += 8)
    {
        X86_ASM (
            "prefetchnta %16\n\t" \
            "movq %8, %%mm0\n\t" \
            "movq %9, %%mm1\n\t" \
            "movq %10, %%mm2\n\t" \
            "movq %11, %%mm3\n\t" \
            "movq %12, %%mm4\n\t" \
            "movq %13, %%mm5\n\t" \
            "movq %14, %%mm6\n\t" \
            "movq %15, %%mm7\n\t" \
            "movntq %%mm0, %0\n\t" \
            "movntq %%mm1, %1\n\t" \
            "movntq %%mm2, %2\n\t" \
            "movntq %%mm3, %3\n\t" \
            "movntq %%mm4, %4\n\t" \
            "movntq %%mm5, %5\n\t" \
            "movntq %%mm6, %6\n\t" \
            "movntq %%mm7, %7\n\t"
            : "=m" (m64pDest[iDataCntr]),
              "=m" (m64pDest[iDataCntr + 1]),
              "=m" (m64pDest[iDataCntr + 2]),
              "=m" (m64pDest[iDataCntr + 3]),
              "=m" (m64pDest[iDataCntr + 4]),
              "=m" (m64pDest[iDataCntr + 5]),
              "=m" (m64pDest[iDataCntr + 6]),
              "=m" (m64pDest[iDataCntr + 7])
            : "m" (m64pSrc[iDataCntr]),
              "m" (m64pSrc[iDataCntr + 1]),
              "m" (m64pSrc[iDataCntr + 2]),
              "m" (m64pSrc[iDataCntr + 3]),
              "m" (m64pSrc[iDataCntr + 4]),
              "m" (m64pSrc[iDataCntr + 5]),
              "m" (m64pSrc[iDataCntr + 6]),
              "m" (m64pSrc[iDataCntr + 7]),
              "m" (m64pSrc[iDataCntr + 32])
            : "mm0", "mm1", "mm2", "mm3", "mm4", "mm5", "mm6", "mm7", "memory");
    }
    iStartIdx = iDataCount;
    iDataCount = ((iDataLength & 0xfffffffe) >> 1);
    for (iDataCntr = iStartIdx; iDataCntr < iDataCount; iDataCntr++)
    {
        X86_ASM (
            "prefetchnta %2\n\t" \
            "movq %1, %%mm0\n\t" \
            "movntq %%mm0, %0\n\t"
            : "=m" (m64pDest[iDataCntr])
            : "m" (m64pSrc[iDataCntr]),
              "m" (m64pSrc[iDataCntr + 32])
            : "mm0", "memory");
    }
    if (iDataLength & 0x1)
    {
        X86_ASM (
            "movd %1, %%mm0\n\t" \
            "movd %%mm0, %0\n\t"
            : "=m" (fpDest[iDataLength - 1])
            : "m" (fpSrc[iDataLength - 1])
            : "mm0", "memory");
    }
    X86_ASM (
        "femms\n\t" \
        "sfence\n\t");
}


void dsp_x86_3dnow_copyd (double *dpDest, const double *dpSrc, int iDataLength)
{
    int iStartIdx;
    int iDataCntr;
    int iDataCount;
    
    iStartIdx = 0;
    X86_ASM (
        "prefetchnta %0\n\t" \
        "prefetchnta %1\n\t" \
        "prefetchnta %2\n\t" \
        "prefetchnta %3\n\t"
        :
        : "m" (dpSrc[0]),
          "m" (dpSrc[8]),
          "m" (dpSrc[16]),
          "m" (dpSrc[24]));
    iDataCount = (iDataLength & 0xfffffff8);
    for (iDataCntr = iStartIdx; iDataCntr < iDataCount; iDataCntr += 8)
    {
        X86_ASM (
            "prefetchnta %16\n\t" \
            "movq %8, %%mm0\n\t" \
            "movq %9, %%mm1\n\t" \
            "movq %10, %%mm2\n\t" \
            "movq %11, %%mm3\n\t" \
            "movq %12, %%mm4\n\t" \
            "movq %13, %%mm5\n\t" \
            "movq %14, %%mm6\n\t" \
            "movq %15, %%mm7\n\t" \
            "movntq %%mm0, %0\n\t" \
            "movntq %%mm1, %1\n\t" \
            "movntq %%mm2, %2\n\t" \
            "movntq %%mm3, %3\n\t" \
            "movntq %%mm4, %4\n\t" \
            "movntq %%mm5, %5\n\t" \
            "movntq %%mm6, %6\n\t" \
            "movntq %%mm7, %7\n\t"
            : "=m" (dpDest[iDataCntr]),
              "=m" (dpDest[iDataCntr + 1]),
              "=m" (dpDest[iDataCntr + 2]),
              "=m" (dpDest[iDataCntr + 3]),
              "=m" (dpDest[iDataCntr + 4]),
              "=m" (dpDest[iDataCntr + 5]),
              "=m" (dpDest[iDataCntr + 6]),
              "=m" (dpDest[iDataCntr + 7])
            : "m" (dpSrc[iDataCntr]),
              "m" (dpSrc[iDataCntr + 1]),
              "m" (dpSrc[iDataCntr + 2]),
              "m" (dpSrc[iDataCntr + 3]),
              "m" (dpSrc[iDataCntr + 4]),
              "m" (dpSrc[iDataCntr + 5]),
              "m" (dpSrc[iDataCntr + 6]),
              "m" (dpSrc[iDataCntr + 7]),
              "m" (dpSrc[iDataCntr + 32])
            : "mm0", "mm1", "mm2", "mm3", "mm4", "mm5", "mm6", "mm7", "memory");
    }
    iStartIdx = iDataCount;
    iDataCount = iDataLength;
    for (iDataCntr = iStartIdx; iDataCntr < iDataCount; iDataCntr++)
    {
        X86_ASM (
            "prefetchnta %2\n\t" \
            "movq %1, %%mm0\n\t" \
            "movntq %%mm0, %0\n\t"
            : "=m" (dpDest[iDataCntr])
            : "m" (dpSrc[iDataCntr]),
              "m" (dpSrc[iDataCntr + 32])
            : "mm0", "memory");
    }
    X86_ASM (
        "femms\n\t" \
        "sfence\n\t");
}


void dsp_x86_3dnow_addf (float *fpVect, float fSrc, int iDataLength)
{
    int iDataCntr;
    int iDataCount;
    pv2sf m64pVect = (pv2sf) fpVect;
    stm64 m64Src;

    m64Src.f[0] = m64Src.f[1] = fSrc;
    iDataCount = (iDataLength >> 1);
    X86_ASM (
        "movq %0, %%mm1\n\t"
        :
        : "m" (m64Src)
        : "mm1", "memory");
    for (iDataCntr = 0; iDataCntr < iDataCount; iDataCntr++)
    {
        X86_ASM (
            "movq %1, %%mm0\n\t" \
            "pfadd %%mm1, %%mm0\n\t" \
            "movntq %%mm0, %0\n\t"
            : "=m" (m64pVect[iDataCntr])
            : "m0" (m64pVect[iDataCntr])
            : "mm0", "mm1", "memory");
    }
    if (iDataLength & 0x1)
    {
        X86_ASM (
            "movd %1, %%mm0\n\t" \
            "pfadd %%mm1, %%mm0\n\t" \
            "movd %%mm0, %0\n\t"
            : "=m" (fpVect[iDataLength - 1])
            : "m0" (fpVect[iDataLength - 1])
            : "mm0", "mm1", "memory");
    }
    X86_ASM (
        "femms\n\t" \
        "sfence\n\t");
}


void dsp_x86_sse_addf (float *fpVect, float fSrc, int iDataLength)
{
    int iDataCntr;
    
    X86_ASM (
        "movss %0, %%xmm1\n\t"
        :
        : "m" (fSrc)
        : "xmm1", "memory");
    for (iDataCntr = 0; iDataCntr < iDataLength; iDataCntr++)
    {
        X86_ASM (
            "movss %1, %%xmm0\n\t" \
            "addss %%xmm1, %%xmm0\n\t" \
            "movss %%xmm0, %0\n\t"
            : "=m" (fpVect[iDataCntr])
            : "m0" (fpVect[iDataCntr])
            : "xmm0", "xmm1", "memory");
    }
}


void dsp_x86_sse_add (double *dpVect, double dSrc, int iDataLength)
{
    int iDataCntr;
    
    X86_ASM (
        "movsd %0, %%xmm1\n\t"
        :
        : "m" (dSrc)
        : "xmm1", "memory");
    for (iDataCntr = 0; iDataCntr < iDataLength; iDataCntr++)
    {
        X86_ASM (
            "movsd %1, %%xmm0\n\t" \
            "addsd %%xmm1, %%xmm0\n\t" \
            "movsd %%xmm0, %0\n\t"
            : "=m" (dpVect[iDataCntr])
            : "m0" (dpVect[iDataCntr])
            : "xmm0", "xmm1", "memory");
    }
}


void dsp_x86_3dnow_mulf (float *fpVect, float fSrc, int iDataLength)
{
    int iDataCntr;
    int iDataCount;
    pv2sf m64pVect = (pv2sf) fpVect;
    stm64 m64Src;

    m64Src.f[0] = m64Src.f[1] = fSrc;
    iDataCount = (iDataLength >> 1);
    X86_ASM (
        "movq %0, %%mm1\n\t"
        :
        : "m" (m64Src)
        : "mm1", "memory");
    for (iDataCntr = 0; iDataCntr < iDataCount; iDataCntr++)
    {
        X86_ASM (
            "movq %1, %%mm0\n\t" \
            "pfmul %%mm1, %%mm0\n\t" \
            "movntq %%mm0, %0\n\t"
            : "=m" (m64pVect[iDataCntr])
            : "m0" (m64pVect[iDataCntr])
            : "mm0", "mm1", "memory");
    }
    if (iDataLength & 0x1)
    {
        X86_ASM (
            "movd %1, %%mm0\n\t" \
            "pfmul %%mm1, %%mm0\n\t" \
            "movd %%mm0, %0\n\t"
            : "=m" (fpVect[iDataLength - 1])
            : "m0" (fpVect[iDataLength - 1])
            : "mm0", "mm1", "memory");
    }
    X86_ASM (
        "femms\n\t" \
        "sfence\n\t");
}


void dsp_x86_sse_mulf (float *fpVect, float fSrc, int iDataLength)
{
    int iDataCntr;

    X86_ASM (
        "movss %0, %%xmm1\n\t"
        :
        : "m" (fSrc)
        : "xmm1", "memory");
    for (iDataCntr = 0; iDataCntr < iDataLength; iDataCntr++)
    {
        X86_ASM (
            "movss %1, %%xmm0\n\t" \
            "mulss %%xmm1, %%xmm0\n\t" \
            "movss %%xmm0, %0\n\t"
            : "=m" (fpVect[iDataCntr])
            : "m0" (fpVect[iDataCntr])
            : "xmm0", "xmm1", "memory");
    }
}


void dsp_x86_sse_mul (double *dpVect, double dSrc, int iDataLength)
{
    int iDataCntr;
    
    X86_ASM (
        "movsd %0, %%xmm1\n\t"
        :
        : "m" (dSrc)
        : "xmm1", "memory");
    for (iDataCntr = 0; iDataCntr <iDataLength; iDataCntr++)
    {
        X86_ASM (
            "movsd %1, %%xmm0\n\t" \
            "mulsd %%xmm1, %%xmm0\n\t" \
            "movsd %%xmm0, %0\n\t"
            : "=m" (dpVect[iDataCntr])
            : "m0" (dpVect[iDataCntr])
            : "xmm0", "xmm1", "memory");
    }
}


void dsp_x86_3dnow_mulf_nip (float *fpDest, const float *fpSrc1, float fSrc2, 
    int iDataLength)
{
    int iDataCntr;
    int iDataCount;
    pv2sf m64pDest = (pv2sf) fpDest;
    pv2sf m64pSrc1 = (pv2sf) fpSrc1;
    stm64 m64Src2;

    m64Src2.f[0] = m64Src2.f[1] = fSrc2;
    iDataCount = (iDataLength >> 1);
    X86_ASM (
        "movq %0, %%mm1\n\t"
        :
        : "m" (m64Src2)
        : "mm1", "memory");
    for (iDataCntr = 0; iDataCntr < iDataCount; iDataCntr++)
    {
        X86_ASM (
            "movq %1, %%mm0\n\t" \
            "pfmul %%mm1, %%mm0\n\t" \
            "movntq %%mm0, %0\n\t"
            : "=m" (m64pDest[iDataCntr])
            : "m" (m64pSrc1[iDataCntr])
            : "mm0", "mm1", "memory");
    }
    if (iDataLength & 0x1)
    {
        X86_ASM (
            "movd %1, %%mm0\n\t" \
            "pfmul %%mm1, %%mm0\n\t" \
            "movd %%mm0, %0\n\t"
            : "=m" (fpDest[iDataLength - 1])
            : "m" (fpSrc1[iDataLength - 1])
            : "mm0", "mm1", "memory");
    }
    X86_ASM (
        "femms\n\t" \
        "sfence\n\t");
}


void dsp_x86_sse_mulf_nip (float *fpDest, const float *fpSrc1, float fSrc2,
    int iDataLength)
{
    int iDataCntr;

    X86_ASM (
        "movss %0, %%xmm1\n\t"
        :
        : "m" (fSrc2)
        : "xmm1", "memory");
    for (iDataCntr = 0; iDataCntr < iDataLength; iDataCntr++)
    {
        X86_ASM (
            "movss %1, %%xmm0\n\t" \
            "mulss %%xmm1, %%xmm0\n\t" \
            "movss %%xmm0, %0\n\t"
            : "=m" (fpDest[iDataCntr])
            : "m" (fpSrc1[iDataCntr])
            : "xmm0", "xmm1", "memory");
    }
}


void dsp_x86_sse_mul_nip (double *dpDest, const double *dpSrc1, double dSrc2,
    int iDataLength)
{
    int iDataCntr;
    
    X86_ASM (
        "movsd %0, %%xmm1\n\t"
        :
        : "m" (dSrc2)
        : "xmm1", "memory");
    for (iDataCntr = 0; iDataCntr < iDataLength; iDataCntr++)
    {
        X86_ASM (
            "movsd %1, %%xmm0\n\t" \
            "mulsd %%xmm1, %%xmm0\n\t" \
            "movsd %%xmm0, %0\n\t"
            : "=m" (dpDest[iDataCntr])
            : "m" (dpSrc1[iDataCntr])
            : "xmm0", "xmm1", "memory");
    }
}


void dsp_x86_3dnow_add2f (float *fpDest, const float *fpSrc, int iDataLength)
{
    int iDataCntr;
    int iDataCount;
    pv2sf m64pDest = (pv2sf) fpDest;
    pv2sf m64pSrc = (pv2sf) fpSrc;

    iDataCount = (iDataLength >> 1);
    for (iDataCntr = 0; iDataCntr < iDataCount; iDataCntr++)
    {
        X86_ASM (
            "movq %1, %%mm0\n\t" \
            "movq %2, %%mm1\n\t" \
            "pfadd %%mm1, %%mm0\n\t" \
            "movntq %%mm0, %0\n\t"
            : "=m" (m64pDest[iDataCntr])
            : "m0" (m64pDest[iDataCntr]),
              "m" (m64pSrc[iDataCntr])
            : "mm0", "mm1", "memory");
    }
    if (iDataLength & 0x1)
    {
        X86_ASM (
            "movd %1, %%mm0\n\t" \
            "movd %2, %%mm1\n\t" \
            "pfadd %%mm1, %%mm0\n\t" \
            "movd %%mm0, %0\n\t"
            : "=m" (fpDest[iDataLength - 1])
            : "m0" (fpDest[iDataLength - 1]),
              "m" (fpSrc[iDataLength - 1])
            : "mm0", "mm1", "memory");
    }
    X86_ASM (
        "femms\n\t" \
        "sfence\n\t");
}


void dsp_x86_sse_add2f (float *fpDest, const float *fpSrc, int iDataLength)
{
    int iDataCntr;
    
    for (iDataCntr = 0; iDataCntr < iDataLength; iDataCntr++)
    {
        X86_ASM (
            "movss %1, %%xmm0\n\t" \
            "addss %2, %%xmm0\n\t" \
            "movss %%xmm0, %0\n\t"
            : "=m" (fpDest[iDataCntr])
            : "m0" (fpDest[iDataCntr]),
              "m" (fpSrc[iDataCntr])
            : "xmm0", "memory");
    }
}


void dsp_x86_sse_add2 (double *dpDest, const double *dpSrc, int iDataLength)
{
    int iDataCntr;
    
    for (iDataCntr = 0; iDataCntr < iDataLength; iDataCntr++)
    {
        X86_ASM (
            "movsd %1, %%xmm0\n\t" \
            "addsd %2, %%xmm0\n\t" \
            "movsd %%xmm0, %0\n\t"
            : "=m" (dpDest[iDataCntr])
            : "m0" (dpDest[iDataCntr]),
              "m" (dpSrc[iDataCntr])
            : "xmm0", "memory");
    }
}


void dsp_x86_3dnow_mul2f (float *fpDest, const float *fpSrc, int iDataLength)
{
    int iDataCntr;
    int iDataCount;
    pv2sf m64pDest = (pv2sf) fpDest;
    pv2sf m64pSrc = (pv2sf) fpSrc;

    iDataCount = (iDataLength >> 1);
    for (iDataCntr = 0; iDataCntr < iDataCount; iDataCntr++)
    {
        X86_ASM (
            "movq %1, %%mm0\n\t" \
            "movq %2, %%mm1\n\t" \
            "pfmul %%mm1, %%mm0\n\t" \
            "movntq %%mm0, %0\n\t"
            : "=m" (m64pDest[iDataCntr])
            : "m0" (m64pDest[iDataCntr]),
              "m" (m64pSrc[iDataCntr])
            : "mm0", "mm1", "memory");
    }
    if (iDataLength & 0x1)
    {
        X86_ASM (
            "movd %1, %%mm0\n\t" \
            "movd %2, %%mm1\n\t" \
            "pfmul %%mm1, %%mm0\n\t" \
            "movd %%mm0, %0\n\t"
            : "=m" (fpDest[iDataLength - 1])
            : "m0" (fpDest[iDataLength - 1]),
              "m" (fpSrc[iDataLength - 1])
            : "mm0", "mm1", "memory");
    }
    X86_ASM (
        "femms\n\t" \
        "sfence\n\t");
}


void dsp_x86_sse_mul2f (float *fpDest, const float *fpSrc, int iDataLength)
{
    int iDataCntr;
    
    for (iDataCntr = 0; iDataCntr < iDataLength; iDataCntr++)
    {
        X86_ASM (
            "movss %1, %%xmm0\n\t" \
            "mulss %2, %%xmm0\n\t" \
            "movss %%xmm0, %0\n\t"
            : "=m" (fpDest[iDataCntr])
            : "m0" (fpDest[iDataCntr]),
              "m" (fpSrc[iDataCntr])
            : "xmm0", "memory");
    }
}


void dsp_x86_sse_mul2 (double *dpDest, const double *dpSrc, int iDataLength)
{
    int iDataCntr;
    
    for (iDataCntr = 0; iDataCntr < iDataLength; iDataCntr++)
    {
        X86_ASM (
            "movsd %1, %%xmm0\n\t" \
            "mulsd %2, %%xmm0\n\t" \
            "movsd %%xmm0, %0\n\t"
            : "=m" (dpDest[iDataCntr])
            : "m0" (dpDest[iDataCntr]),
              "m" (dpSrc[iDataCntr])
            : "xmm0", "memory");
    }
}


void dsp_x86_3dnow_add3f (float *fpDest, const float *fpSrc1, 
    const float *fpSrc2, int iDataLength)
{
    int iDataCntr;
    int iDataCount;
    pv2sf m64pDest = (pv2sf) fpDest;
    pv2sf m64pSrc1 = (pv2sf) fpSrc1;
    pv2sf m64pSrc2 = (pv2sf) fpSrc2;

    iDataCount = (iDataLength >> 1);
    for (iDataCntr = 0; iDataCntr < iDataCount; iDataCntr++)
    {
        X86_ASM (
            "movq %1, %%mm0\n\t" \
            "movq %2, %%mm1\n\t" \
            "pfadd %%mm1, %%mm0\n\t" \
            "movntq %%mm0, %0\n\t"
            : "=m" (m64pDest[iDataCntr])
            : "m" (m64pSrc1[iDataCntr]),
              "m" (m64pSrc2[iDataCntr])
            : "mm0", "mm1", "memory");
    }
    if (iDataLength & 0x1)
    {
        X86_ASM (
            "movd %1, %%mm0\n\t" \
            "movd %2, %%mm1\n\t" \
            "pfadd %%mm1, %%mm0\n\t" \
            "movd %%mm0, %0\n\t"
            : "=m" (fpDest[iDataLength - 1])
            : "m" (fpSrc1[iDataLength - 1]),
              "m" (fpSrc2[iDataLength - 1])
            : "mm0", "mm1", "memory");
    }
    X86_ASM (
        "femms\n\t" \
        "sfence\n\t");
}


void dsp_x86_sse_add3f (float *fpDest, const float *fpSrc1, 
    const float *fpSrc2, int iDataLength)
{
    int iDataCntr;
    
    for (iDataCntr = 0; iDataCntr < iDataLength; iDataCntr++)
    {
        X86_ASM (
            "movss %1, %%xmm0\n\t" \
            "addss %2, %%xmm0\n\t" \
            "movss %%xmm0, %0\n\t"
            : "=m" (fpDest[iDataCntr])
            : "m" (fpSrc1[iDataCntr]),
              "m" (fpSrc2[iDataCntr])
            : "xmm0", "memory");
    }
}


void dsp_x86_sse_add3 (double *dpDest, const double *dpSrc1, 
    const double *dpSrc2, int iDataLength)
{
    int iDataCntr;
    
    for (iDataCntr = 0; iDataCntr < iDataLength; iDataCntr++)
    {
        X86_ASM (
            "movsd %1, %%xmm0\n\t" \
            "addsd %2, %%xmm0\n\t" \
            "movsd %%xmm0, %0\n\t"
            : "=m" (dpDest[iDataCntr])
            : "m" (dpSrc1[iDataCntr]),
              "m" (dpSrc2[iDataCntr])
            : "xmm0", "memory");
    }
}


void dsp_x86_3dnow_mul3f (float *fpDest, const float *fpSrc1, 
    const float *fpSrc2, int iDataLength)
{
    int iDataCntr;
    int iDataCount;
    pv2sf m64pDest = (pv2sf) fpDest;
    pv2sf m64pSrc1 = (pv2sf) fpSrc1;
    pv2sf m64pSrc2 = (pv2sf) fpSrc2;

    iDataCount = (iDataLength >> 1);
    for (iDataCntr = 0; iDataCntr < iDataCount; iDataCntr++)
    {
        X86_ASM (
            "movq %1, %%mm0\n\t" \
            "movq %2, %%mm1\n\t" \
            "pfmul %%mm1, %%mm0\n\t" \
            "movntq %%mm0, %0\n\t"
            : "=m" (m64pDest[iDataCntr])
            : "m" (m64pSrc1[iDataCntr]),
              "m" (m64pSrc2[iDataCntr])
            : "mm0", "mm1", "memory");
    }
    if (iDataLength & 0x1)
    {
        X86_ASM (
            "movd %1, %%mm0\n\t" \
            "movd %2, %%mm1\n\t" \
            "pfmul %%mm1, %%mm0\n\t" \
            "movd %%mm0, %0\n\t"
            : "=m" (fpDest[iDataLength - 1])
            : "m" (fpSrc1[iDataLength - 1]),
              "m" (fpSrc2[iDataLength - 1])
            : "mm0", "mm1", "memory");
    }
    X86_ASM (
        "femms\n\t" \
        "sfence\n\t");
}


void dsp_x86_sse_mul3f (float *fpDest, const float *fpSrc1, 
    const float *fpSrc2, int iDataLength)
{
    int iDataCntr;
    
    for (iDataCntr = 0; iDataCntr < iDataLength; iDataCntr++)
    {
        X86_ASM (
            "movss %1, %%xmm0\n\t" \
            "mulss %2, %%xmm0\n\t" \
            "movss %%xmm0, %0\n\t"
            : "=m" (fpDest[iDataCntr])
            : "m" (fpSrc1[iDataCntr]),
              "m" (fpSrc2[iDataCntr])
            : "xmm0", "memory");
    }
}


void dsp_x86_sse_mul3 (double *dpDest, const double *dpSrc1, 
    const double *dpSrc2, int iDataLength)
{
    int iDataCntr;
    
    for (iDataCntr = 0; iDataCntr < iDataLength; iDataCntr++)
    {
        X86_ASM (
            "movsd %1, %%xmm0\n\t" \
            "mulsd %2, %%xmm0\n\t" \
            "movsd %%xmm0, %0\n\t"
            : "=m" (dpDest[iDataCntr])
            : "m" (dpSrc1[iDataCntr]),
              "m" (dpSrc2[iDataCntr])
            : "xmm0", "memory");
    }
}


void dsp_x86_3dnow_cmulf (float *fpDest, const float *fpSrc, int iDataLength)
{
    int iDataCntr;
    pv2sf m64pDest = (pv2sf) fpDest;
    
    X86_ASM (
        "movq %0, %%mm3\n\t"
        :
        : "m" (fpSrc[0])
        : "mm3", "memory");
    for (iDataCntr = 0; iDataCntr < iDataLength; iDataCntr++)
    {
        X86_ASM (
            "movq %1, %%mm0\n\t" \
            "movq %%mm3, %%mm1\n\t" \
            "pswapd %%mm1, %%mm2\n\t" \
            "pfmul %%mm0, %%mm1\n\t" \
            "pfmul %%mm0, %%mm2\n\t" \
            "pfpnacc %%mm2, %%mm1\n\t" \
            "movntq %%mm1, %0\n\t"
            : "=m" (m64pDest[iDataCntr])
            : "m0" (m64pDest[iDataCntr])
            : "mm0", "mm1", "mm2", "mm3", "memory");
    }
    X86_ASM (
        "femms\n\t" \
        "sfence\n\t");
}


void dsp_x86_sse_cmulf (float *fpDest, const float *fpSrc, int iDataLength)
{
    int iDataCntr;
    int iDataCount;
    
    X86_ASM (
        "movss %0, %%xmm2\n\t" \
        "movss %1, %%xmm3\n\t"
        :
        : "m" (fpSrc[0]),
          "m" (fpSrc[1])
        : "xmm2", "xmm3", "memory");
    iDataCount = (iDataLength << 1);
    for (iDataCntr = 0; iDataCntr < iDataCount; iDataCntr += 2)
    {
        X86_ASM (
            "movss %2, %%xmm0\n\t" \
            "movss %%xmm0, %%xmm1\n\t" \
            "movss %3, %%xmm4\n\t" \
            \
            "mulss %%xmm2, %%xmm0\n\t" \
            "movss %%xmm4, %%xmm5\n\t" \
            "mulss %%xmm3, %%xmm5\n\t" \
            "subss %%xmm5, %%xmm0\n\t" \
            \
            "mulss %%xmm3, %%xmm1\n\t" \
            "movss %%xmm4, %%xmm5\n\t" \
            "mulss %%xmm2, %%xmm5\n\t" \
            "addss %%xmm5, %%xmm1\n\t" \
            \
            "movss %%xmm0, %0\n\t" \
            "movss %%xmm1, %1\n\t"
            : "=m" (fpDest[iDataCntr]),
              "=m" (fpDest[iDataCntr + 1])
            : "m0" (fpDest[iDataCntr]),
              "m1" (fpDest[iDataCntr + 1])
            : "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "memory");
    }
}


void dsp_x86_sse3_cmulf (float *fpDest, const float *fpSrc, int iDataLength)
{
    int iDataCntr;
    int iDataCount;
    unsigned long ulDestAddr;
    unsigned long ulSrcAddr;
    pv4sf m128pDest = (pv4sf) fpDest;
    
    ulDestAddr = (unsigned long) fpDest;
    ulSrcAddr = (unsigned long) fpSrc;
    if (unlikely((ulDestAddr & 0xf) || (ulSrcAddr & 0x8)))
    {
        dsp_x86_sse_cmulf(fpDest, fpSrc, iDataLength);
        return;
    }
    X86_ASM (
        "movlps %0, %%xmm5\n\t" \
        "movhps %0, %%xmm5\n\t" \
        "movsldup %%xmm5, %%xmm0\n\t" \
        "movshdup %%xmm5, %%xmm2\n\t"
        :
        : "m" (*fpSrc)
        : "xmm0", "xmm2", "xmm5", "memory");
    iDataCount = (iDataLength >> 1);
    for (iDataCntr = 0; iDataCntr < iDataCount; iDataCntr++)
    {
        X86_ASM (
            "movaps %1, %%xmm1\n\t" \
            "movaps %%xmm1, %%xmm3\n\t" \
            "shufps $0x11, %%xmm1, %%xmm3\n\t" \
            \
            "mulps %%xmm0, %%xmm1\n\t" \
            "mulps %%xmm2, %%xmm3\n\t" \
            "addsubps %%xmm3, %%xmm1\n\t" \
            \
            "movaps %%xmm1, %0\n\t"
            : "=m" (m128pDest[iDataCntr])
            : "m0" (m128pDest[iDataCntr])
            : "xmm0", "xmm1", "xmm2", "xmm3", "memory");
    }
    if (iDataLength & 1)
    {
        dsp_x86_sse_cmulf(fpDest + (iDataCount << 2), fpSrc, 1);
    }
}


void dsp_x86_sse_cmul (double *dpDest, const double *dpSrc, int iDataLength)
{
    int iDataCntr;
    int iDataCount;
    
    X86_ASM (
        "movsd %0, %%xmm2\n\t" \
        "movsd %1, %%xmm3\n\t"
        :
        : "m" (dpSrc[0]),
          "m" (dpSrc[1])
        : "xmm2", "xmm3", "memory");
    iDataCount = (iDataLength << 1);
    for (iDataCntr = 0; iDataCntr < iDataCount; iDataCntr += 2)
    {
        X86_ASM (
            "movsd %2, %%xmm0\n\t" \
            "movsd %%xmm0, %%xmm1\n\t" \
            "movsd %3, %%xmm4\n\t" \
            \
            "mulsd %%xmm2, %%xmm0\n\t" \
            "movsd %%xmm4, %%xmm5\n\t" \
            "mulsd %%xmm3, %%xmm5\n\t" \
            "subsd %%xmm5, %%xmm0\n\t" \
            \
            "mulsd %%xmm3, %%xmm1\n\t" \
            "movsd %%xmm4, %%xmm5\n\t" \
            "mulsd %%xmm2, %%xmm5\n\t" \
            "addsd %%xmm5, %%xmm1\n\t" \
            \
            "movsd %%xmm0, %0\n\t" \
            "movsd %%xmm1, %1\n\t"
            : "=m" (dpDest[iDataCntr]),
              "=m" (dpDest[iDataCntr + 1])
            : "m0" (dpDest[iDataCntr]),
              "m1" (dpDest[iDataCntr + 1])
            : "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "memory");
    }
}


void dsp_x86_sse3_cmul (double *dpDest, const double *dpSrc, int iDataLength)
{
    int iDataCntr;
    unsigned long ulDestAddr;
    unsigned long ulSrcAddr;
    pv2df m128pDest = (pv2df) dpDest;
    
    ulDestAddr = (unsigned long) dpDest;
    ulSrcAddr = (unsigned long) dpSrc;
    if (unlikely((ulDestAddr & 0xf) || (ulSrcAddr & 0xf)))
    {
        dsp_x86_sse_cmul(dpDest, dpSrc, iDataLength);
        return;
    }
    X86_ASM (
        "movddup %0, %%xmm0\n\t" \
        "movddup %1, %%xmm2\n\t"
        :
        : "m" (*dpSrc),
          "m" (*(dpSrc + 1))
        : "xmm0", "xmm2", "memory");
    for (iDataCntr = 0; iDataCntr < iDataLength; iDataCntr++)
    {
        X86_ASM (
            "movapd %1, %%xmm1\n\t" \
            "movapd %%xmm1, %%xmm3\n\t" \
            "shufpd $0x01, %%xmm1, %%xmm3\n\t" \
            \
            "mulpd %%xmm0, %%xmm1\n\t" \
            "mulpd %%xmm2, %%xmm3\n\t" \
            "addsubpd %%xmm3, %%xmm1\n\t" \
            \
            "movapd %%xmm1, %0\n\t"
            : "=m" (m128pDest[iDataCntr])
            : "m0" (m128pDest[iDataCntr])
            : "xmm0", "xmm1", "xmm2", "xmm3", "memory");
    }
}


void dsp_x86_3dnow_cmul2f (float *fpDest, const float *fpSrc, int iDataLength)
{
    int iDataCntr;
    pv2sf m64pDest = (pv2sf) fpDest;
    pv2sf m64pSrc = (pv2sf) fpSrc;
    
    for (iDataCntr = 0; iDataCntr < iDataLength; iDataCntr++)
    {
        X86_ASM (
            "movq %1, %%mm0\n\t" \
            "movq %2, %%mm1\n\t" \
            "pswapd %%mm1, %%mm2\n\t" \
            "pfmul %%mm0, %%mm1\n\t" \
            "pfmul %%mm0, %%mm2\n\t" \
            "pfpnacc %%mm2, %%mm1\n\t" \
            "movntq %%mm1, %0\n\t"
            : "=m" (m64pDest[iDataCntr])
            : "m0" (m64pDest[iDataCntr]),
              "m" (m64pSrc[iDataCntr])
            : "mm0", "mm1", "mm2", "memory");
    }
    X86_ASM (
        "femms\n\t" \
        "sfence\n\t");
}


void dsp_x86_sse_cmul2f (float *fpDest, const float *fpSrc, int iDataLength)
{
    int iDataCntr;
    int iDataCount;
    
    iDataCount = (iDataLength << 1);
    for (iDataCntr = 0; iDataCntr < iDataCount; iDataCntr += 2)
    {
        X86_ASM (
            "movss %4, %%xmm2\n\t" \
            "movss %5, %%xmm3\n\t" \
            \
            "movss %2, %%xmm0\n\t" \
            "movss %%xmm0, %%xmm1\n\t" \
            "movss %3, %%xmm4\n\t" \
            \
            "mulss %%xmm2, %%xmm0\n\t" \
            "movss %%xmm4, %%xmm5\n\t" \
            "mulss %%xmm3, %%xmm5\n\t" \
            "subss %%xmm5, %%xmm0\n\t" \
            \
            "mulss %%xmm3, %%xmm1\n\t" \
            "movss %%xmm4, %%xmm5\n\t" \
            "mulss %%xmm2, %%xmm5\n\t" \
            "addss %%xmm5, %%xmm1\n\t" \
            \
            "movss %%xmm0, %0\n\t" \
            "movss %%xmm1, %1\n\t"
            : "=m" (fpDest[iDataCntr]),
              "=m" (fpDest[iDataCntr + 1])
            : "m0" (fpDest[iDataCntr]),
              "m1" (fpDest[iDataCntr + 1]),
              "m" (fpSrc[iDataCntr]),
              "m" (fpSrc[iDataCntr + 1])
            : "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "memory");
    }
}


void dsp_x86_sse3_cmul2f (float *fpDest, const float *fpSrc, int iDataLength)
{
    int iDataCntr;
    int iDataCount;
    unsigned long ulDestAddr;
    unsigned long ulSrcAddr;
    pv4sf m128pDest = (pv4sf) fpDest;
    pv4sf m128pSrc = (pv4sf) fpSrc;
    
    ulDestAddr = (unsigned long) fpDest;
    ulSrcAddr = (unsigned long) fpSrc;
    if (unlikely((ulDestAddr & 0xf) || (ulSrcAddr & 0xf)))
    {
        dsp_x86_sse_cmul2f(fpDest, fpSrc, iDataLength);
        return;
    }
    iDataCount = (iDataLength >> 1);
    for (iDataCntr = 0; iDataCntr < iDataCount; iDataCntr++)
    {
        X86_ASM (
            "movsldup %1, %%xmm0\n\t" \
            "movshdup %1, %%xmm2\n\t" \
            "movaps %2, %%xmm1\n\t" \
            "movaps %%xmm1, %%xmm3\n\t" \
            "shufps $0x11, %%xmm1, %%xmm3\n\t" \
            \
            "mulps %%xmm1, %%xmm0\n\t" \
            "mulps %%xmm3, %%xmm2\n\t" \
            "addsubps %%xmm2, %%xmm0\n\t" \
            \
            "movaps %%xmm0, %0\n\t"
            : "=m" (m128pDest[iDataCntr])
            : "m0" (m128pDest[iDataCntr]),
              "m" (m128pSrc[iDataCntr])
            : "xmm0", "xmm1", "xmm2", "xmm3", "memory");
    }
    if (iDataLength & 1)
    {
        dsp_x86_sse_cmul2f(fpDest + (iDataCount << 2),
            fpSrc + (iDataCount << 2), 1);
    }
}


void dsp_x86_sse_cmul2 (double *dpDest, const double *dpSrc, int iDataLength)
{
    int iDataCntr;
    int iDataCount;
    
    iDataCount = (iDataLength << 1);
    for (iDataCntr = 0; iDataCntr < iDataCount; iDataCntr += 2)
    {
        X86_ASM (
            "movsd %4, %%xmm2\n\t" \
            "movsd %5, %%xmm3\n\t" \
            \
            "movsd %2, %%xmm0\n\t" \
            "movsd %%xmm0, %%xmm1\n\t" \
            "movsd %3, %%xmm4\n\t" \
            \
            "mulsd %%xmm2, %%xmm0\n\t" \
            "movsd %%xmm4, %%xmm5\n\t" \
            "mulsd %%xmm3, %%xmm5\n\t" \
            "subsd %%xmm5, %%xmm0\n\t" \
            \
            "mulsd %%xmm3, %%xmm1\n\t" \
            "movsd %%xmm4, %%xmm5\n\t" \
            "mulsd %%xmm2, %%xmm5\n\t" \
            "addsd %%xmm5, %%xmm1\n\t" \
            \
            "movsd %%xmm0, %0\n\t" \
            "movsd %%xmm1, %1\n\t"
            : "=m" (dpDest[iDataCntr]),
              "=m" (dpDest[iDataCntr + 1])
            : "m0" (dpDest[iDataCntr]),
              "m1" (dpDest[iDataCntr + 1]),
              "m" (dpSrc[iDataCntr]),
              "m" (dpSrc[iDataCntr + 1])
            : "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "memory");
    }
}


void dsp_x86_sse3_cmul2 (double *dpDest, const double *dpSrc, int iDataLength)
{
    int iDataCntr;
    int iDataIdx2 = 0;
    unsigned long ulDestAddr;
    unsigned long ulSrcAddr;
    pv2df m128pDest = (pv2df) dpDest;
    
    ulDestAddr = (unsigned long) dpDest;
    ulSrcAddr = (unsigned long) dpSrc;
    if (unlikely((ulDestAddr & 0xf) || (ulSrcAddr & 0xf)))
    {
        dsp_x86_sse_cmul2(dpDest, dpSrc, iDataLength);
        return;
    }
    for (iDataCntr = 0; iDataCntr < iDataLength; iDataCntr++)
    {
        X86_ASM (
            "movddup %2, %%xmm0\n\t" \
            "movddup %3, %%xmm2\n\t" \
            \
            "movapd %1, %%xmm1\n\t" \
            "movapd %%xmm1, %%xmm3\n\t" \
            "shufpd $0x01, %%xmm1, %%xmm3\n\t" \
            \
            "mulpd %%xmm0, %%xmm1\n\t" \
            "mulpd %%xmm2, %%xmm3\n\t" \
            "addsubpd %%xmm3, %%xmm1\n\t" \
            \
            "movapd %%xmm1, %0\n\t"
            : "=m" (m128pDest[iDataCntr])
            : "m0" (m128pDest[iDataCntr]),
              "m" (dpSrc[iDataIdx2]),
              "m" (dpSrc[iDataIdx2 + 1])
            : "xmm0", "xmm1", "xmm2", "xmm3", "memory");
        iDataIdx2 += 2;
    }
}


void dsp_x86_3dnow_cmul3f (float *fpDest, const float *fpSrc1, 
    const float *fpSrc2, int iDataLength)
{
    int iDataCntr;
    pv2sf m64pDest = (pv2sf) fpDest;
    pv2sf m64pSrc1 = (pv2sf) fpSrc1;
    pv2sf m64pSrc2 = (pv2sf) fpSrc2;
    
    for (iDataCntr = 0; iDataCntr < iDataLength; iDataCntr++)
    {
        X86_ASM (
            "movq %1, %%mm0\n\t" \
            "movq %2, %%mm1\n\t" \
            "pswapd %%mm1, %%mm2\n\t" \
            "pfmul %%mm0, %%mm1\n\t" \
            "pfmul %%mm0, %%mm2\n\t" \
            "pfpnacc %%mm2, %%mm1\n\t" \
            "movntq %%mm1, %0\n\t"
            : "=m" (m64pDest[iDataCntr])
            : "m" (m64pSrc1[iDataCntr]),
              "m" (m64pSrc2[iDataCntr])
            : "mm0", "mm1", "mm2", "memory");
    }
    X86_ASM (
        "femms\n\t" \
        "sfence\n\t");
}


void dsp_x86_sse_cmul3f (float *fpDest, const float *fpSrc1, 
    const float *fpSrc2, int iDataLength)
{
    int iDataCntr;
    int iDataCount;
    
    iDataCount = (iDataLength << 1);
    for (iDataCntr = 0; iDataCntr < iDataCount; iDataCntr += 2)
    {
        X86_ASM (
            "movss %4, %%xmm2\n\t" \
            "movss %5, %%xmm3\n\t" \
            \
            "movss %2, %%xmm0\n\t" \
            "movss %%xmm0, %%xmm1\n\t" \
            "movss %3, %%xmm4\n\t" \
            \
            "mulss %%xmm2, %%xmm0\n\t" \
            "movss %%xmm4, %%xmm5\n\t" \
            "mulss %%xmm3, %%xmm5\n\t" \
            "subss %%xmm5, %%xmm0\n\t" \
            \
            "mulss %%xmm3, %%xmm1\n\t" \
            "movss %%xmm4, %%xmm5\n\t" \
            "mulss %%xmm2, %%xmm5\n\t" \
            "addss %%xmm5, %%xmm1\n\t" \
            \
            "movss %%xmm0, %0\n\t" \
            "movss %%xmm1, %1\n\t"
            : "=m" (fpDest[iDataCntr]),
              "=m" (fpDest[iDataCntr + 1])
            : "m" (fpSrc1[iDataCntr]),
              "m" (fpSrc1[iDataCntr + 1]),
              "m" (fpSrc2[iDataCntr]),
              "m" (fpSrc2[iDataCntr + 1])
            : "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "memory");
    }
}


void dsp_x86_sse3_cmul3f (float *fpDest, const float *fpSrc1,
    const float *fpSrc2, int iDataLength)
{
    int iDataCntr;
    int iDataCount;
    unsigned long ulDestAddr;
    unsigned long ulSrc1Addr;
    unsigned long ulSrc2Addr;
    pv4sf m128pDest = (pv4sf) fpDest;
    pv4sf m128pSrc1 = (pv4sf) fpSrc1;
    pv4sf m128pSrc2 = (pv4sf) fpSrc2;
    
    ulDestAddr = (unsigned long) fpDest;
    ulSrc1Addr = (unsigned long) fpSrc1;
    ulSrc2Addr = (unsigned long) fpSrc2;
    if (unlikely((ulDestAddr & 0xf) ||
        (ulSrc1Addr & 0xf) || (ulSrc2Addr & 0xf)))
    {
        dsp_x86_sse_cmul3f(fpDest, fpSrc1, fpSrc2, iDataLength);
        return;
    }
    iDataCount = (iDataLength >> 1);
    for (iDataCntr = 0; iDataCntr < iDataCount; iDataCntr++)
    {
        X86_ASM (
            "movsldup %1, %%xmm0\n\t" \
            "movshdup %1, %%xmm2\n\t" \
            "movaps %2, %%xmm1\n\t" \
            "movaps %%xmm1, %%xmm3\n\t" \
            "shufps $0x11, %%xmm1, %%xmm3\n\t" \
            \
            "mulps %%xmm1, %%xmm0\n\t" \
            "mulps %%xmm3, %%xmm2\n\t" \
            "addsubps %%xmm2, %%xmm0\n\t" \
            \
            "movaps %%xmm0, %0\n\t"
            : "=m" (m128pDest[iDataCntr])
            : "m" (m128pSrc1[iDataCntr]),
              "m" (m128pSrc2[iDataCntr])
            : "xmm0", "xmm1", "xmm2", "xmm3", "memory");
    }
    if (iDataLength & 1)
    {
        dsp_x86_sse_cmul3f(fpDest + (iDataCount << 2),
            fpSrc1 + (iDataCount << 2), fpSrc2 + (iDataCount << 2), 1);
    }
}


void dsp_x86_sse_cmul3 (double *dpDest, const double *dpSrc1, 
    const double *dpSrc2, int iDataLength)
{
    int iDataCntr;
    int iDataCount;
    
    iDataCount = (iDataLength << 1);
    for (iDataCntr = 0; iDataCntr < iDataCount; iDataCntr += 2)
    {
        X86_ASM (
            "movsd %4, %%xmm2\n\t" \
            "movsd %5, %%xmm3\n\t" \
            \
            "movsd %2, %%xmm0\n\t" \
            "movsd %%xmm0, %%xmm1\n\t" \
            "movsd %3, %%xmm4\n\t" \
            \
            "mulsd %%xmm2, %%xmm0\n\t" \
            "movsd %%xmm4, %%xmm5\n\t" \
            "mulsd %%xmm3, %%xmm5\n\t" \
            "subsd %%xmm5, %%xmm0\n\t" \
            \
            "mulsd %%xmm3, %%xmm1\n\t" \
            "movsd %%xmm4, %%xmm5\n\t" \
            "mulsd %%xmm2, %%xmm5\n\t" \
            "addsd %%xmm5, %%xmm1\n\t" \
            \
            "movsd %%xmm0, %0\n\t" \
            "movsd %%xmm1, %1\n\t"
            : "=m" (dpDest[iDataCntr]),
              "=m" (dpDest[iDataCntr + 1])
            : "m" (dpSrc1[iDataCntr]),
              "m" (dpSrc1[iDataCntr + 1]),
              "m" (dpSrc2[iDataCntr]),
              "m" (dpSrc2[iDataCntr + 1])
            : "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "memory");
    }
}


void dsp_x86_sse3_cmul3 (double *dpDest, const double *dpSrc1, 
    const double *dpSrc2, int iDataLength)
{
    int iDataCntr;
    int iDataIdx2 = 0;
    unsigned long ulDestAddr;
    unsigned long ulSrc1Addr;
    unsigned long ulSrc2Addr;
    pv2df m128pDest = (pv2df) dpDest;
    pv2df m128pSrc1 = (pv2df) dpSrc1;
    
    ulDestAddr = (unsigned long) dpDest;
    ulSrc1Addr = (unsigned long) dpSrc1;
    ulSrc2Addr = (unsigned long) dpSrc2;
    if (unlikely((ulDestAddr & 0xf) ||
        (ulSrc1Addr & 0xf) || (ulSrc2Addr & 0xf)))
    {
        dsp_x86_sse_cmul3(dpDest, dpSrc1, dpSrc2, iDataLength);
        return;
    }
    for (iDataCntr = 0; iDataCntr < iDataLength; iDataCntr++)
    {
        X86_ASM (
            "movddup %2, %%xmm0\n\t" \
            "movddup %3, %%xmm2\n\t" \
            \
            "movapd %1, %%xmm1\n\t" \
            "movapd %%xmm1, %%xmm3\n\t" \
            "shufpd $0x01, %%xmm1, %%xmm3\n\t" \
            \
            "mulpd %%xmm0, %%xmm1\n\t" \
            "mulpd %%xmm2, %%xmm3\n\t" \
            "addsubpd %%xmm3, %%xmm1\n\t" \
            \
            "movapd %%xmm1, %0\n\t"
            : "=m" (m128pDest[iDataCntr])
            : "m" (m128pSrc1[iDataCntr]),
              "m" (dpSrc2[iDataIdx2]),
              "m" (dpSrc2[iDataIdx2 + 1])
            : "xmm0", "xmm1", "xmm2", "xmm3", "memory");
        iDataIdx2 += 2;
    }
}


void dsp_x86_3dnow_maf (float *fpVect, float fMul, float fAdd, int iDataLength)
{
    int iDataCntr;
    int iDataCount;
    pv2sf m64pVect = (pv2sf) fpVect;
    stm64 m64Mul;
    stm64 m64Add;

    m64Mul.f[0] = m64Mul.f[1] = fMul;
    m64Add.f[0] = m64Add.f[1] = fAdd;
    iDataCount = (iDataLength >> 1);
    X86_ASM (
        "movq %0, %%mm1\n\t" \
        "movq %1, %%mm2\n\t"
        :
        : "m" (m64Mul),
          "m" (m64Add)
        : "mm1", "mm2", "memory");
    for (iDataCntr = 0; iDataCntr < iDataCount; iDataCntr++)
    {
        X86_ASM (
            "movq %1, %%mm0\n\t" \
            "pfmul %%mm1, %%mm0\n\t" \
            "pfadd %%mm2, %%mm0\n\t" \
            "movntq %%mm0, %0\n\t"
            : "=m" (m64pVect[iDataCntr])
            : "m0" (m64pVect[iDataCntr])
            : "mm0", "mm1", "mm2", "memory");
    }
    if (iDataLength & 0x1)
    {
        X86_ASM (
            "movd %1, %%mm0\n\t" \
            "pfmul %%mm1, %%mm0\n\t" \
            "pfadd %%mm2, %%mm0\n\t" \
            "movd %%mm0, %0\n\t"
            : "=m" (fpVect[iDataLength - 1])
            : "m0" (fpVect[iDataLength - 1])
            : "mm0", "mm1", "mm2", "memory");
    }
    X86_ASM (
        "femms\n\t" \
        "sfence\n\t");
}


void dsp_x86_sse_maf (float *fpVect, float fMul, float fAdd, int iDataLength)
{
    int iDataCntr;
    
    X86_ASM (
        "movss %0, %%xmm1\n\t" \
        "movss %1, %%xmm2\n\t"
        :
        : "m" (fMul),
          "m" (fAdd)
        : "xmm1", "xmm2", "memory");
    for (iDataCntr = 0; iDataCntr < iDataLength; iDataCntr++)
    {
        X86_ASM (
            "movss %1, %%xmm0\n\t" \
            "mulss %%xmm1, %%xmm0\n\t" \
            "addss %%xmm2, %%xmm0\n\t" \
            "movss %%xmm0, %0\n\t"
            : "=m" (fpVect[iDataCntr])
            : "m0" (fpVect[iDataCntr])
            : "xmm0", "xmm1", "xmm2", "memory");
    }
}


void dsp_x86_sse_ma (double *dpVect, double dMul, double dAdd, int iDataLength)
{
    int iDataCntr;
    
    X86_ASM (
        "movsd %0, %%xmm1\n\t" \
        "movsd %1, %%xmm2\n\t"
        :
        : "m" (dMul),
          "m" (dAdd)
        : "xmm1", "xmm2", "memory");
    for (iDataCntr = 0; iDataCntr < iDataLength; iDataCntr++)
    {
        X86_ASM (
            "movsd %1, %%xmm0\n\t" \
            "mulsd %%xmm1, %%xmm0\n\t" \
            "addsd %%xmm2, %%xmm0\n\t" \
            "movsd %%xmm0, %0\n\t"
            : "=m" (dpVect[iDataCntr])
            : "m0" (dpVect[iDataCntr])
            : "xmm0", "xmm1", "xmm2", "memory");
    }
}


void dsp_x86_3dnow_ma2f (float *fpDest, const float *fpSrc,
    float fMul, float fAdd, int iDataLength)
{
    int iDataCntr;
    int iDataCount;
    pv2sf m64pDest = (pv2sf) fpDest;
    pv2sf m64pSrc = (pv2sf) fpSrc;
    stm64 m64Mul;
    stm64 m64Add;

    m64Mul.f[0] = m64Mul.f[1] = fMul;
    m64Add.f[0] = m64Add.f[1] = fAdd;
    iDataCount = (iDataLength >> 1);
    X86_ASM (
        "movq %0, %%mm1\n\t" \
        "movq %1, %%mm2\n\t"
        :
        : "m" (m64Mul),
          "m" (m64Add)
        : "mm1", "mm2", "memory");
    for (iDataCntr = 0; iDataCntr < iDataCount; iDataCntr++)
    {
        X86_ASM (
            "movq %1, %%mm0\n\t" \
            "pfmul %%mm1, %%mm0\n\t" \
            "pfadd %%mm2, %%mm0\n\t" \
            "movntq %%mm0, %0\n\t"
            : "=m" (m64pDest[iDataCntr])
            : "m" (m64pSrc[iDataCntr])
            : "mm0", "mm1", "mm2", "memory");
    }
    if (iDataLength & 0x1)
    {
        X86_ASM (
            "movd %1, %%mm0\n\t" \
            "pfmul %%mm1, %%mm0\n\t" \
            "pfadd %%mm2, %%mm0\n\t" \
            "movd %%mm0, %0\n\t"
            : "=m" (fpDest[iDataLength - 1])
            : "m" (fpSrc[iDataLength - 1])
            : "mm0", "mm1", "mm2", "memory");
    }
    X86_ASM (
        "femms\n\t" \
        "sfence\n\t");
}


void dsp_x86_sse_ma2f (float *fpDest, const float *fpSrc, 
    float fMul, float fAdd, int iDataLength)
{
    int iDataCntr;
    
    X86_ASM (
        "movss %0, %%xmm1\n\t" \
        "movss %1, %%xmm2\n\t"
        :
        : "m" (fMul),
          "m" (fAdd)
        : "xmm1", "xmm2", "memory");
    for (iDataCntr = 0; iDataCntr < iDataLength; iDataCntr++)
    {
        X86_ASM (
            "movss %1, %%xmm0\n\t" \
            "mulss %%xmm1, %%xmm0\n\t" \
            "addss %%xmm2, %%xmm0\n\t" \
            "movss %%xmm0, %0\n\t"
            : "=m" (fpDest[iDataCntr])
            : "m" (fpSrc[iDataCntr])
            : "xmm0", "xmm1", "xmm2", "memory");
    }
}


void dsp_x86_sse_ma2 (double *dpDest, const double *dpSrc, 
    double dMul, double dAdd, int iDataLength)
{
    int iDataCntr;
    
    X86_ASM (
        "movsd %0, %%xmm1\n\t" \
        "movsd %1, %%xmm2\n\t"
        :
        : "m" (dMul),
          "m" (dAdd)
        : "xmm1", "xmm2", "memory");
    for (iDataCntr = 0; iDataCntr < iDataLength; iDataCntr++)
    {
        X86_ASM (
            "movsd %1, %%xmm0\n\t" \
            "mulsd %%xmm1, %%xmm0\n\t" \
            "addsd %%xmm2, %%xmm0\n\t" \
            "movsd %%xmm0, %0\n\t"
            : "=m" (dpDest[iDataCntr])
            : "m" (dpSrc[iDataCntr])
            : "xmm0", "xmm1", "xmm2", "memory");
    }
}


void dsp_x86_3dnow_cmaf (float *fpDest, const float *fpSrc1, 
    const float *fpSrc2, int iDataLength)
{
    int iDataCntr;
    pv2sf m64pDest = (pv2sf) fpDest;
    pv2sf m64pSrc1 = (pv2sf) fpSrc1;
    pv2sf m64pSrc2 = (pv2sf) fpSrc2;
    
    for (iDataCntr = 0; iDataCntr < iDataLength; iDataCntr++)
    {
        X86_ASM (
            "movq %2, %%mm0\n\t" \
            "movq %3, %%mm1\n\t" \
            "movq %1, %%mm3\n\t" \
            "pswapd %%mm1, %%mm2\n\t" \
            "pfmul %%mm0, %%mm1\n\t" \
            "pfmul %%mm0, %%mm2\n\t" \
            "pfpnacc %%mm2, %%mm1\n\t" \
            "pfadd %%mm1, %%mm3\n\t" \
            "movntq %%mm3, %0\n\t"
            : "=m" (m64pDest[iDataCntr])
            : "m0" (m64pDest[iDataCntr]),
              "m" (m64pSrc1[iDataCntr]),
              "m" (m64pSrc2[iDataCntr])
            : "mm0", "mm1", "mm2", "mm3", "memory");
    }
    X86_ASM (
        "femms\n\t" \
        "sfence\n\t");
}


void dsp_x86_sse_cmaf (float *fpDest, const float *fpSrc1, 
    const float *fpSrc2, int iDataLength)
{
    int iDataCntr;
    int iDataCount;
    
    iDataCount = (iDataLength << 1);
    for (iDataCntr = 0; iDataCntr < iDataCount; iDataCntr += 2)
    {
        X86_ASM (
            "movss %6, %%xmm2\n\t" \
            "movss %7, %%xmm3\n\t" \
            \
            "movss %4, %%xmm0\n\t" \
            "movss %%xmm0, %%xmm1\n\t" \
            "movss %5, %%xmm4\n\t" \
            \
            "movss %2, %%xmm6\n\t" \
            "movss %3, %%xmm7\n\t" \
            \
            "mulss %%xmm2, %%xmm0\n\t" \
            "movss %%xmm4, %%xmm5\n\t" \
            "mulss %%xmm3, %%xmm5\n\t" \
            "subss %%xmm5, %%xmm0\n\t" \
            \
            "mulss %%xmm3, %%xmm1\n\t" \
            "movss %%xmm4, %%xmm5\n\t" \
            "mulss %%xmm2, %%xmm5\n\t" \
            "addss %%xmm5, %%xmm1\n\t" \
            \
            "addss %%xmm0, %%xmm6\n\t" \
            "addss %%xmm1, %%xmm7\n\t" \
            \
            "movss %%xmm6, %0\n\t" \
            "movss %%xmm7, %1\n\t"
            : "=m" (fpDest[iDataCntr]),
              "=m" (fpDest[iDataCntr + 1])
            : "m0" (fpDest[iDataCntr]),
              "m1" (fpDest[iDataCntr + 1]),
              "m" (fpSrc1[iDataCntr]),
              "m" (fpSrc1[iDataCntr + 1]),
              "m" (fpSrc2[iDataCntr]),
              "m" (fpSrc2[iDataCntr + 1])
            : "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7",
              "memory");
    }
}


void dsp_x86_sse3_cmaf (float *fpDest, const float *fpSrc1,
    const float *fpSrc2, int iDataLength)
{
    int iDataCntr;
    int iDataCount;
    unsigned long ulDestAddr;
    unsigned long ulSrc1Addr;
    unsigned long ulSrc2Addr;
    pv4sf m128pDest = (pv4sf) fpDest;
    pv4sf m128pSrc1 = (pv4sf) fpSrc1;
    pv4sf m128pSrc2 = (pv4sf) fpSrc2;
    
    ulDestAddr = (unsigned long) fpDest;
    ulSrc1Addr = (unsigned long) fpSrc1;
    ulSrc2Addr = (unsigned long) fpSrc2;
    if (unlikely((ulDestAddr & 0xf) ||
        (ulSrc1Addr & 0xf) || (ulSrc2Addr & 0xf)))
    {
        dsp_x86_sse_cmaf(fpDest, fpSrc1, fpSrc2, iDataLength);
        return;
    }
    iDataCount = (iDataLength >> 1);
    for (iDataCntr = 0; iDataCntr < iDataCount; iDataCntr++)
    {
        X86_ASM (
            "movsldup %2, %%xmm0\n\t" \
            "movshdup %2, %%xmm2\n\t" \
            "movaps %3, %%xmm1\n\t" \
            "movaps %%xmm1, %%xmm3\n\t" \
            "shufps $0x11, %%xmm1, %%xmm3\n\t" \
            \
            "mulps %%xmm1, %%xmm0\n\t" \
            "mulps %%xmm3, %%xmm2\n\t" \
            "addsubps %%xmm2, %%xmm0\n\t" \
            \
            "movaps %1, %%xmm4\n\t" \
            "addps %%xmm0, %%xmm4\n\t" \
            "movaps %%xmm4, %0\n\t"
            : "=m" (m128pDest[iDataCntr])
            : "m0" (m128pDest[iDataCntr]),
              "m" (m128pSrc1[iDataCntr]),
              "m" (m128pSrc2[iDataCntr])
            : "xmm0", "xmm1", "xmm2", "xmm3", "memory");
    }
    if (iDataLength & 1)
    {
        dsp_x86_sse_cmaf(fpDest + (iDataCount << 2),
            fpSrc1 + (iDataCount << 2), fpSrc2 + (iDataCount << 2), 1);
    }
}


void dsp_x86_sse_cma (double *dpDest, const double *dpSrc1, 
    const double *dpSrc2, int iDataLength)
{
    int iDataCntr;
    int iDataCount;
    
    iDataCount = (iDataLength << 1);
    for (iDataCntr = 0; iDataCntr < iDataCount; iDataCntr += 2)
    {
        X86_ASM (
            "movsd %6, %%xmm2\n\t" \
            "movsd %7, %%xmm3\n\t" \
            \
            "movsd %4, %%xmm0\n\t" \
            "movsd %%xmm0, %%xmm1\n\t" \
            "movsd %5, %%xmm4\n\t" \
            \
            "movsd %2, %%xmm6\n\t" \
            "movsd %3, %%xmm7\n\t" \
            \
            "mulsd %%xmm2, %%xmm0\n\t" \
            "movsd %%xmm4, %%xmm5\n\t" \
            "mulsd %%xmm3, %%xmm5\n\t" \
            "subsd %%xmm5, %%xmm0\n\t" \
            \
            "mulsd %%xmm3, %%xmm1\n\t" \
            "movsd %%xmm4, %%xmm5\n\t" \
            "mulsd %%xmm2, %%xmm5\n\t" \
            "addsd %%xmm5, %%xmm1\n\t" \
            \
            "addsd %%xmm0, %%xmm6\n\t" \
            "addsd %%xmm1, %%xmm7\n\t" \
            \
            "movsd %%xmm6, %0\n\t" \
            "movsd %%xmm7, %1\n\t"
            : "=m" (dpDest[iDataCntr]),
              "=m" (dpDest[iDataCntr + 1])
            : "m0" (dpDest[iDataCntr]),
              "m1" (dpDest[iDataCntr + 1]),
              "m" (dpSrc1[iDataCntr]),
              "m" (dpSrc1[iDataCntr + 1]),
              "m" (dpSrc2[iDataCntr]),
              "m" (dpSrc2[iDataCntr + 1])
            : "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7",
              "memory");
    }
}


void dsp_x86_sse3_cma (double *dpDest, const double *dpSrc1, 
    const double *dpSrc2, int iDataLength)
{
    int iDataCntr;
    int iDataIdx2 = 0;
    unsigned long ulDestAddr;
    unsigned long ulSrc1Addr;
    unsigned long ulSrc2Addr;
    pv2df m128pDest = (pv2df) dpDest;
    pv2df m128pSrc1 = (pv2df) dpSrc1;
    
    ulDestAddr = (unsigned long) dpDest;
    ulSrc1Addr = (unsigned long) dpSrc1;
    ulSrc2Addr = (unsigned long) dpSrc2;
    if (unlikely((ulDestAddr & 0xf) || (ulSrc1Addr & 0xf) ||
        (ulSrc2Addr & 0xf)))
    {
        dsp_x86_sse_cma(dpDest, dpSrc1, dpSrc2, iDataLength);
        return;
    }
    for (iDataCntr = 0; iDataCntr < iDataLength; iDataCntr++)
    {
        X86_ASM (
            "movddup %3, %%xmm0\n\t" \
            "movddup %4, %%xmm2\n\t" \
            \
            "movapd %2, %%xmm1\n\t" \
            "movapd %%xmm1, %%xmm3\n\t" \
            "shufpd $0x01, %%xmm1, %%xmm3\n\t" \
            \
            "mulpd %%xmm0, %%xmm1\n\t" \
            "mulpd %%xmm2, %%xmm3\n\t" \
            "addsubpd %%xmm3, %%xmm1\n\t" \
            \
            "movapd %1, %%xmm4\n\t" \
            "addpd %%xmm1, %%xmm4\n\t" \
            "movapd %%xmm4, %0\n\t"
            : "=m" (m128pDest[iDataCntr])
            : "m0" (m128pDest[iDataCntr]),
              "m" (m128pSrc1[iDataCntr]),
              "m" (dpSrc2[iDataIdx2]),
              "m" (dpSrc2[iDataIdx2 + 1])
            : "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "memory");
        iDataIdx2 += 2;
    }
}


void dsp_x86_3dnow_amf (float *fpVect, float fAdd, float fMul, int iDataLength)
{
    int iDataCntr;
    int iDataCount;
    pv2sf m64pVect = (pv2sf) fpVect;
    stm64 m64Add;
    stm64 m64Mul;

    m64Add.f[0] = m64Add.f[1] = fAdd;
    m64Mul.f[0] = m64Mul.f[1] = fMul;
    iDataCount = (iDataLength >> 1);
    X86_ASM (
        "movq %0, %%mm1\n\t" \
        "movq %1, %%mm2\n\t"
        :
        : "m" (m64Add),
          "m" (m64Mul)
        : "mm1", "mm2", "memory");
    for (iDataCntr = 0; iDataCntr < iDataCount; iDataCntr++)
    {
        X86_ASM (
            "movq %1, %%mm0\n\t" \
            "pfadd %%mm1, %%mm0\n\t" \
            "pfmul %%mm2, %%mm0\n\t" \
            "movntq %%mm0, %0\n\t"
            : "=m" (m64pVect[iDataCntr])
            : "m0" (m64pVect[iDataCntr])
            : "mm0", "mm1", "mm2", "memory");
    }
    if (iDataLength & 0x1)
    {
        X86_ASM (
            "movd %1, %%mm0\n\t" \
            "pfadd %%mm1, %%mm0\n\t" \
            "pfmul %%mm2, %%mm0\n\t" \
            "movd %%mm0, %0\n\t"
            : "=m" (fpVect[iDataLength - 1])
            : "m0" (fpVect[iDataLength - 1])
            : "mm0", "mm1", "mm2", "memory");
    }
    X86_ASM (
        "femms\n\t" \
        "sfence\n\t");
}


void dsp_x86_sse_amf (float *fpVect, float fAdd, float fMul, int iDataLength)
{
    int iDataCntr;

    X86_ASM (
        "movss %0, %%xmm1\n\t" \
        "movss %1, %%xmm2\n\t"
        :
        : "m" (fAdd),
          "m" (fMul)
        : "xmm1", "xmm2", "memory");
    for (iDataCntr = 0; iDataCntr < iDataLength; iDataCntr++)
    {
        X86_ASM (
            "movss %1, %%xmm0\n\t" \
            "addss %%xmm1, %%xmm0\n\t" \
            "mulss %%xmm2, %%xmm0\n\t" \
            "movss %%xmm0, %0\n\t"
            : "=m" (fpVect[iDataCntr])
            : "m0" (fpVect[iDataCntr])
            : "xmm0", "xmm1", "xmm2", "memory");
    }
}


float dsp_x86_3dnow_macf (const float *fpSrc1, const float *fpSrc2, 
    int iDataLength)
{
    int iDataCntr;
    int iDataCount;
    float fRes;
    pv2sf m64pSrc1 = (pv2sf) fpSrc1;
    pv2sf m64pSrc2 = (pv2sf) fpSrc2;

    iDataCount = (iDataLength >> 1);
    X86_ASM (
        "pxor %%mm0, %%mm0\n\t"
        :
        :
        : "mm0");
    for (iDataCntr = 0; iDataCntr < iDataCount; iDataCntr++)
    {
        X86_ASM (
            "movq %0, %%mm1\n\t" \
            "movq %1, %%mm2\n\t" \
            "pfmul %%mm2, %%mm1\n\t" \
            "pfacc %%mm1, %%mm0\n\t"
            :
            : "m" (m64pSrc1[iDataCntr]),
              "m" (m64pSrc2[iDataCntr])
            : "mm0", "mm1", "mm2", "memory");
    }
    if (iDataLength & 0x1)
    {
        X86_ASM (
            "movd %0, %%mm1\n\t" \
            "movd %1, %%mm2\n\t" \
            "pfmul %%mm2, %%mm1\n\t" \
            "pfacc %%mm1, %%mm0\n\t"
            :
            : "m" (fpSrc1[iDataLength - 1]),
              "m" (fpSrc2[iDataLength - 1])
            : "mm0", "mm1", "mm2", "memory");
    }
    X86_ASM (
        "pfacc %%mm0, %%mm0\n\t" \
        "movd %%mm0, %0\n\t"
        : "=m" (fRes)
        :
        : "mm0", "memory");
    X86_ASM ("femms\n\t");

    return fRes;
}


float dsp_x86_sse_macf (const float *fpSrc1, const float *fpSrc2,
    int iDataLength)
{
    int iDataCntr;
    float fRes;
    
    X86_ASM (
        "xorps %%xmm0, %%xmm0\n\t"
        :
        :
        : "xmm0");
    for (iDataCntr = 0; iDataCntr < iDataLength; iDataCntr++)
    {
        X86_ASM (
            "movss %0, %%xmm1\n\t" \
            "mulss %1, %%xmm1\n\t" \
            "addss %%xmm1, %%xmm0\n\t"
            :
            : "m" (fpSrc1[iDataCntr]),
              "m" (fpSrc2[iDataCntr])
            : "xmm0", "xmm1", "xmm2", "memory");
    }
    X86_ASM (
        "movss %%xmm0, %0\n\t"
        : "=m" (fRes)
        :
        : "xmm0");

    return fRes;
}


double dsp_x86_sse_mac (const double *dpSrc1, const double *dpSrc2,
    int iDataLength)
{
    int iDataCntr;
    double dRes;
    
    X86_ASM (
        "xorpd %%xmm0, %%xmm0\n\t"
        :
        :
        : "xmm0");
    for (iDataCntr = 0; iDataCntr < iDataLength; iDataCntr++)
    {
        X86_ASM (
            "movsd %0, %%xmm1\n\t" \
            "mulsd %1, %%xmm1\n\t" \
            "addsd %%xmm1, %%xmm0\n\t"
            :
            : "m" (dpSrc1[iDataCntr]),
              "m" (dpSrc2[iDataCntr])
            : "xmm0", "xmm1", "xmm2", "memory");
    }
    X86_ASM (
        "movsd %%xmm0, %0\n\t"
        : "=m" (dRes)
        :
        : "xmm0");

    return dRes;
}


void dsp_x86_3dnow_minmaxf (float *fpMin, float *fpMax, const float *fpSrc, 
    int iDataLength)
{
    int iDataCntr;
    int iDataCount;
    stm64 m64Min;
    stm64 m64Max;
    pv2sf m64pSrc = (pv2sf) fpSrc;
    
    m64Min.f[0] = m64Min.f[1] = FLT_MAX;
    m64Max.f[0] = m64Max.f[1] = -FLT_MAX;
    iDataCount = (iDataLength >> 1);
    X86_ASM (
        "movq %0, %%mm1\n\t" \
        "movq %1, %%mm2\n\t"
        :
        : "m" (m64Min),
          "m" (m64Max)
        : "mm1", "mm2", "memory");
    for (iDataCntr = 0; iDataCntr < iDataCount; iDataCntr++)
    {
        X86_ASM (
            "movq %0, %%mm0\n\t" \
            "pfmin %%mm0, %%mm1\n\t" \
            "pfmax %%mm0, %%mm2\n\t"
            :
            : "m" (m64pSrc[iDataCntr])
            : "mm0", "mm1", "mm2", "memory");
    }
    if (iDataLength & 0x1)
    {
        X86_ASM (
            "movd %0, %%mm0\n\t" \
            "pfmin %%mm0, %%mm1\n\t" \
            "pfmax %%mm0, %%mm2\n\t"
            :
            : "m" (fpSrc[iDataLength - 1])
            : "mm0", "mm1", "mm2", "memory");
    }
    X86_ASM (
        "pswapd %%mm1, %%mm3\n\t" \
        "pfmin %%mm3, %%mm1\n\t" \
        "pswapd %%mm2, %%mm3\n\t" \
        "pfmax %%mm3, %%mm2\n\t" \
        "movd %%mm1, %0\n\t" \
        "movd %%mm2, %1\n\t"
        : "=m" (*fpMin),
          "=m" (*fpMax)
        :
        : "mm1", "mm2", "mm3", "memory");
    X86_ASM ("femms\n\t");
}


void dsp_x86_sse_minmaxf (float *fpMin, float *fpMax, const float *fpSrc, 
    int iDataLength)
{
    int iDataCntr;

    *fpMin = FLT_MAX;
    *fpMax = -FLT_MAX;
    X86_ASM (
        "movss %0, %%xmm0\n\t" \
        "movss %1, %%xmm1\n\t"
        :
        : "m" (*fpMin),
          "m" (*fpMax)
        : "xmm0", "xmm1", "memory");
    for (iDataCntr = 0; iDataCntr < iDataLength; iDataCntr++)
    {
        X86_ASM (
            "movss %0, %%xmm2\n\t" \
            "minss %%xmm2, %%xmm0\n\t" \
            "maxss %%xmm2, %%xmm1\n\t"
            :
            : "m" (fpSrc[iDataCntr])
            : "xmm0", "xmm1", "xmm2", "memory");
    }
    X86_ASM (
        "movss %%xmm0, %0\n\t" \
        "movss %%xmm1, %1\n\t"
        : "=m" (*fpMin),
          "=m" (*fpMax)
        :
        : "xmm0", "xmm1", "memory");
}


void dsp_x86_sse_minmax (double *dpMin, double *dpMax, const double *dpSrc, 
    int iDataLength)
{
    int iDataCntr;

    *dpMin = FLT_MAX;
    *dpMax = -FLT_MAX;
    X86_ASM (
        "movsd %0, %%xmm0\n\t" \
        "movsd %1, %%xmm1\n\t"
        :
        : "m" (*dpMin),
          "m" (*dpMax)
        : "xmm0", "xmm1", "memory");
    for (iDataCntr = 0; iDataCntr < iDataLength; iDataCntr++)
    {
        X86_ASM (
            "movsd %0, %%xmm2\n\t" \
            "minsd %%xmm2, %%xmm0\n\t" \
            "maxsd %%xmm2, %%xmm1\n\t"
            :
            : "m" (dpSrc[iDataCntr])
            : "xmm0", "xmm1", "xmm2", "memory");
    }
    X86_ASM (
        "movss %%xmm0, %0\n\t" \
        "movss %%xmm1, %1\n\t"
        : "=m" (*dpMin),
          "=m" (*dpMax)
        :
        : "xmm0", "xmm1", "memory");
}


float dsp_x86_3dnow_crosscorrf (const float *fpSrc1, const float *fpSrc2,
    int iDataLength)
{
    int iDataCntr;
    int iDataCount;
    float fRes;
    pv2sf m64pSrc1 = (pv2sf) fpSrc1;
    pv2sf m64pSrc2 = (pv2sf) fpSrc2;
    
    iDataCount = (iDataLength >> 1);
    X86_ASM (
        "pxor %%mm3, %%mm3\n\t" \
        "pxor %%mm4, %%mm4\n\t" \
        "pxor %%mm5, %%mm5\n\t"
        :
        :
        : "mm3", "mm4", "mm5");
    for (iDataCntr = 0; iDataCntr < iDataCount; iDataCntr++)
    {
        X86_ASM (
            "movq %0, %%mm0\n\t" \
            "movq %1, %%mm1\n\t" \
            "movq %%mm1, %%mm2\n\t" \
            "pfmul %%mm0, %%mm2\n\t" \
            "pfacc %%mm2, %%mm5\n\t" \
            "pfmul %%mm0, %%mm0\n\t" \
            "pfacc %%mm0, %%mm3\n\t" \
            "pfmul %%mm1, %%mm1\n\t" \
            "pfacc %%mm1, %%mm4\n\t"
            :
            : "m" (m64pSrc1[iDataCntr]),
              "m" (m64pSrc2[iDataCntr])
            : "mm0", "mm1", "mm2", "mm3", "mm4", "mm5", "memory");
    }
    if (iDataLength & 0x1)
    {
        X86_ASM (
            "movd %0, %%mm0\n\t" \
            "movd %1, %%mm1\n\t" \
            "movq %%mm1, %%mm2\n\t" \
            "pfmul %%mm0, %%mm2\n\t" \
            "pfacc %%mm2, %%mm5\n\t" \
            "pfmul %%mm0, %%mm0\n\t" \
            "pfacc %%mm0, %%mm3\n\t" \
            "pfmul %%mm1, %%mm1\n\t" \
            "pfacc %%mm1, %%mm4\n\t"
            :
            : "m" (fpSrc1[iDataLength - 1]),
              "m" (fpSrc2[iDataLength - 1])
            : "mm0", "mm1", "mm2", "mm3", "mm4", "mm5", "memory");
    }
    X86_ASM (
        "pfacc %%mm3, %%mm3\n\t" \
        "pfacc %%mm4, %%mm4\n\t" \
        "pfacc %%mm5, %%mm5\n\t" \
        \
        "movd %1, %%mm6\n\t" \
        "pswapd %%mm6, %%mm7\n\t" \
        "paddd %%mm7, %%mm6\n\t" \
        "pi2fd %%mm6, %%mm7\n\t" \
        \
        "pfrcp %%mm7, %%mm6\n\t" \
        "pfrcpit1 %%mm6, %%mm7\n\t" \
        "pfrcpit2 %%mm6, %%mm7\n\t" \
        \
        "pfmul %%mm3, %%mm4\n\t" \
        \
        "movq %%mm4, %%mm0\n\t" \
        "pfrsqrt %%mm4, %%mm1\n\t" \
        "movq %%mm1, %%mm2\n\t" \
        "pfmul %%mm1, %%mm1\n\t" \
        "pfrsqit1 %%mm4, %%mm1\n\t" \
        "pfrcpit2 %%mm2, %%mm1\n\t" \
        "pfmul %%mm1, %%mm4\n\t" \
        \
        "pfmul %%mm6, %%mm4\n\t" \
        \
        "pfrcp %%mm4, %%mm0\n\t" \
        "pfrcpit1 %%mm0, %%mm4\n\t" \
        "pfrcpit2 %%mm0, %%mm4\n\t" \
        \
        "pfmul %%mm6, %%mm5\n\t" \
        "pfmul %%mm4, %%mm5\n\t" \
        "movd %%mm5, %0\n\t"
        : "=m" (fRes)
        : "m" (iDataLength)
        : "mm0", "mm1", "mm2", "mm3", "mm4", "mm5", "mm6", "mm7", "memory");
    X86_ASM ("femms\n\t");

    return fRes;
}


float dsp_x86_sse_crosscorrf (const float *fpSrc1, const float *fpSrc2,
    int iDataLength)
{
    int iDataCntr;
    float fScale;
    float fNormFact;
    float fProdSum;
    float fSqSum1;
    float fSqSum2;
    float fRes;
    
    X86_ASM (
        "xorps %%xmm0, %%xmm0\n\t" \
        "xorps %%xmm1, %%xmm1\n\t" \
        "xorps %%xmm2, %%xmm2\n\t"
        :
        :
        : "xmm0", "xmm1", "xmm2");
    for (iDataCntr = 0; iDataCntr < iDataLength; iDataCntr++)
    {
        X86_ASM (
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
            : "m" (fpSrc1[iDataCntr]),
              "m" (fpSrc2[iDataCntr])
            : "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5");
    }
    fScale = 1.0F / iDataLength;
    fNormFact = sqrtf(fSqSum1 * fSqSum2) * fScale;
    fRes = (fProdSum * fScale) / fNormFact;

    return fRes;
}


double dsp_x86_sse_crosscorr (const double *dpSrc1, const double *dpSrc2,
    int iDataLength)
{
    int iDataCntr;
    double dScale;
    double dNormFact;
    double dProdSum;
    double dSqSum1;
    double dSqSum2;
    double dRes;
    
    X86_ASM (
        "xorpd %%xmm0, %%xmm0\n\t" \
        "xorpd %%xmm1, %%xmm1\n\t" \
        "xorpd %%xmm2, %%xmm2\n\t"
        :
        :
        : "xmm0", "xmm1", "xmm2");
    for (iDataCntr = 0; iDataCntr < iDataLength; iDataCntr++)
    {
        X86_ASM (
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
            : "m" (dpSrc1[iDataCntr]),
              "m" (dpSrc2[iDataCntr])
            : "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5");
    }
    dScale = 1.0 / iDataLength;
    dNormFact = sqrt(dSqSum1 * dSqSum2) * dScale;
    dRes = (dProdSum * dScale) / dNormFact;

    return dRes;
}


void dsp_x86_3dnow_i16tof (float *fpDest, const short *ipSrc, int iDataLength,
    int iIntMax)
{
    int iDataCntr;
    float fScale;
    
    X86_ASM (
        "movd %1, %%mm1\n\t" \
        "pswapd %%mm1, %%mm2\n\t" \
        "paddd %%mm2, %%mm1\n\t" \
        "pi2fd %%mm1, %%mm1\n\t" \
        "pfrcp %%mm1, %%mm2\n\t" \
        "pfrcpit1 %%mm2, %%mm1\n\t" \
        "pfrcpit2 %%mm2, %%mm1\n\t" \
        "movd %%mm1, %0\n\t"
        : "=m" (fScale)
        : "m" (iIntMax)
        : "mm1", "mm2", "memory");
    for (iDataCntr = 0; iDataCntr < iDataLength; iDataCntr += 2)
    {
        X86_ASM (
            "movd %1, %%mm0\n\t" \
            "punpcklwd %%mm0, %%mm0\n\t" \
            "pi2fw %%mm0, %%mm0\n\t" \
            "pfmul %%mm1, %%mm0\n\t" \
            "movntq %%mm0, %0\n\t"
            : "=m" (fpDest[iDataCntr])
            : "m" (ipSrc[iDataCntr])
            : "mm0", "mm1", "memory");
    }
    X86_ASM (
        "femms\n\t" \
        "sfence\n\t");
    if ((iDataLength % 2) != 0)
    {
        fpDest[iDataLength - 1] = ((float) ipSrc[iDataLength - 1]) * fScale;
    }
}


void dsp_x86_3dnow_i32tof (float *fpDest, const int *ipSrc, int iDataLength,
    int iIntMax)
{
    int iDataCntr;
    float fScale;
    
    X86_ASM (
        "movd %1, %%mm1\n\t" \
        "pswapd %%mm1, %%mm2\n\t" \
        "paddd %%mm2, %%mm1\n\t" \
        "pi2fd %%mm1, %%mm1\n\t" \
        "pfrcp %%mm1, %%mm2\n\t" \
        "pfrcpit1 %%mm2, %%mm1\n\t" \
        "pfrcpit2 %%mm2, %%mm1\n\t" \
        "movd %%mm1, %0\n\t"
        : "=m" (fScale)
        : "m" (iIntMax)
        : "mm1", "mm2", "memory");
    for (iDataCntr = 0; iDataCntr < iDataLength; iDataCntr += 2)
    {
        X86_ASM (
            "movq %1, %%mm0\n\t" \
            "pi2fd %%mm0, %%mm0\n\t" \
            "pfmul %%mm1, %%mm0\n\t" \
            "movntq %%mm0, %0\n\t"
            : "=m" (fpDest[iDataCntr])
            : "m" (ipSrc[iDataCntr])
            : "mm0", "mm1", "memory");
    }
    X86_ASM (
        "femms\n\t" \
        "sfence\n\t");
    if ((iDataLength % 2) != 0)
    {
        fpDest[iDataLength - 1] = ((float) ipSrc[iDataLength - 1]) * fScale;
    }
}


void dsp_x86_3dnow_firf (float *fpDest, const float *fpSrc, int iDataLength, 
    const float *fpCoeff, int iCoeffLength)
{
    int iSrcCntr;
    int iDestCntr;
    int iCoeffCntr;
    int iSrcCount;
    pv2sf m64pDest = (pv2sf) fpDest;

    iDestCntr = 0;
    iSrcCount = iDataLength + iCoeffLength;
    for (iSrcCntr = iCoeffLength; 
        iSrcCntr < iSrcCount; 
        iSrcCntr += 2)
    {
        X86_ASM (
            "pxor %%mm0, %%mm0\n\t" 
            :
            :
            : "mm0");
        for (iCoeffCntr = 0; 
            iCoeffCntr < iCoeffLength;
            iCoeffCntr++)
        {
            X86_ASM (
                "movq %0, %%mm1\n\t" \
                "movd %1, %%mm2\n\t" \
                "pswapd %%mm2, %%mm3\n\t" \
                "pfadd %%mm3, %%mm2\n\t" \
                "pfmul %%mm2, %%mm1\n\t" \
                "pfadd %%mm1, %%mm0\n\t" 
                :
                : "m" (fpSrc[iSrcCntr - iCoeffCntr]),
                  "m" (fpCoeff[iCoeffCntr])
                : "mm0", "mm1", "mm2", "mm3", "memory");
        }
        X86_ASM (
            "movntq %%mm0, %0\n\t"
            : "=m" (m64pDest[iDestCntr++])
            :
            : "mm0", "memory");
    }
    if (iDataLength & 0x1)
    {
        X86_ASM (
            "pxor %%mm0, %%mm0\n\t" 
            :
            :
            : "mm0");
        for (iCoeffCntr = 0; 
            iCoeffCntr < iCoeffLength;
            iCoeffCntr++)
        {
            X86_ASM (
                "movd %0, %%mm1\n\t" \
                "movd %1, %%mm2\n\t" \
                "pfmul %%mm2, %%mm1\n\t" \
                "pfadd %%mm1, %%mm0\n\t" 
                :
                : "m" (fpSrc[iDataLength - 1 - iCoeffCntr]),
                  "m" (fpCoeff[iCoeffCntr])
                : "mm0", "mm1", "mm2", "memory");
        }
        X86_ASM (
            "movd %%mm0, %0\n\t"
            : "=m" (fpDest[iDataLength - 1])
            :
            : "mm0", "memory");
    }
    X86_ASM (
        "femms\n\t" \
        "sfence\n\t");
}


void dsp_x86_sse_firf (float *fpDest, const float *fpSrc, int iDataLength, 
    const float *fpCoeff, int iCoeffLength)
{
    int iDestCntr;
    int iSrcCntr;
    int iCoeffCntr;
    int iSrcCount;

    iDestCntr = 0;
    iSrcCount = iDataLength + iCoeffLength;
    for (iSrcCntr = iCoeffLength;
        iSrcCntr < iSrcCount;
        iSrcCntr++)
    {
        X86_ASM (
            "xorps %%xmm0, %%xmm0\n\t"
            :
            :
            : "xmm0");
        for (iCoeffCntr = 0;
            iCoeffCntr < iCoeffLength;
            iCoeffCntr++)
        {
            X86_ASM (
                "movss %0, %%xmm1\n\t"
                "mulss %1, %%xmm1\n\t"
                "addss %%xmm1, %%xmm0\n\t"
                :
                : "m" (fpSrc[iSrcCntr - iCoeffCntr]),
                  "m" (fpCoeff[iCoeffCntr])
                : "xmm0", "xmm1", "memory");
        }
        X86_ASM (
            "movss %%xmm0, %0\n\t"
            : "=m" (fpDest[iDestCntr++])
            :
            : "xmm0", "memory");
    }
}


void dsp_x86_sse_fir (double *dpDest, const double *dpSrc, int iDataLength, 
    const double *dpCoeff, int iCoeffLength)
{
    int iDestCntr;
    int iSrcCntr;
    int iCoeffCntr;
    int iSrcCount;

    iDestCntr = 0;
    iSrcCount = iDataLength + iCoeffLength;
    for (iSrcCntr = iCoeffLength;
        iSrcCntr < iSrcCount;
        iSrcCntr++)
    {
        X86_ASM (
            "xorpd %%xmm0, %%xmm0\n\t"
            :
            :
            : "xmm0");
        for (iCoeffCntr = 0;
            iCoeffCntr < iCoeffLength;
            iCoeffCntr++)
        {
            X86_ASM (
                "movsd %0, %%xmm1\n\t"
                "mulsd %1, %%xmm1\n\t"
                "addsd %%xmm1, %%xmm0\n\t"
                :
                : "m" (dpSrc[iSrcCntr - iCoeffCntr]),
                  "m" (dpCoeff[iCoeffCntr])
                : "xmm0", "xmm1", "memory");
        }
        X86_ASM (
            "movsd %%xmm0, %0\n\t"
            : "=m" (dpDest[iDestCntr++])
            :
            : "xmm0", "memory");
    }
}


void dsp_x86_3dnow_iirf (float *fpVect, int iDataLength, const float *fpCoeff, 
    float *fpX, float *fpY)
{
    int iDataCntr;
    pv2sf m64pCoeff = (pv2sf) &fpCoeff[1];
    pv2sf m64pCoeff2 = (pv2sf) &fpCoeff[3];
    pv2sf m64pX = (pv2sf) fpX;
    pv2sf m64pY = (pv2sf) fpY;

    X86_ASM (
        "movq %0, %%mm0\n\t" \
        "pswapd %%mm0, %%mm2\n\t" \
        "movd %1, %%mm3\n\t" \
        "movq %2, %%mm0\n\t" \
        "pswapd %%mm0, %%mm4\n\t" \
        "movq %3, %%mm5\n\t" \
        "movq %4, %%mm7\n\t" \
        :
        : "m" (*m64pCoeff),
          "m" (fpCoeff[0]),
          "m" (*m64pCoeff2),
          "m" (*m64pX),
          "m" (*m64pY)
        : "mm0", "mm2", "mm3", "mm4", "mm5", "mm7", "memory");
    for (iDataCntr = 0; 
        iDataCntr < iDataLength; 
        iDataCntr++)
    {
        X86_ASM (
            "pxor %%mm0, %%mm0\n\t" \
            "movd %1, %%mm6\n\t" \
            "movq %%mm5, %%mm1\n\t" \
            "pfmul %%mm2, %%mm1\n\t" \
            "pfacc %%mm1, %%mm0\n\t" \
            "movq %%mm6, %%mm1\n\t" \
            "pfmul %%mm3, %%mm1\n\t" \
            "pfacc %%mm1, %%mm0\n\t" \
            "movq %%mm7, %%mm1\n\t" \
            "pfmul %%mm4, %%mm1\n\t" \
            "pfacc %%mm1, %%mm0\n\t" \
            "pfacc %%mm0, %%mm0\n\t" \
            \
            "pswapd %%mm7, %%mm1\n\t" \
            "movq %%mm1, %%mm7\n\t" \
            "punpckldq %%mm0, %%mm7\n\t" \
            \
            "pswapd %%mm5, %%mm1\n\t" \
            "movq %%mm1, %%mm5\n\t" \
            "movq %%mm6, %%mm1\n\t" \
            "punpckldq %%mm1, %%mm5\n\t" \
            \
            "movd %%mm0, %0\n\t"
            : "=m" (fpVect[iDataCntr])
            : "m0" (fpVect[iDataCntr])
            : "mm0", "mm1", "mm2", "mm3", "mm4", "mm5", "mm6", "mm7", "memory");
    }
    X86_ASM (
        "movq %%mm5, %0\n\t" \
        "movd %%mm6, %1\n\t" \
        "movq %%mm7, %2\n\t"
        : "=m" (*m64pX),
          "=m" (fpX[2]),
          "=m" (*m64pY)
        :
        : "mm5", "mm6", "mm7", "memory");
    X86_ASM ("femms\n\t");
}


void dsp_x86_sse_iirf (float *fpVect, int iDataLength, const float *fpCoeff, 
    float *fpX, float *fpY)
{
    int iDataCntr;
    
    X86_ASM (
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
    for (iDataCntr = 0; iDataCntr < iDataLength; iDataCntr++)
    {
        X86_ASM (
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
            : "=m" (fpVect[iDataCntr])
            : "m0" (fpVect[iDataCntr]),
              "m" (fpCoeff[0]),
              "m" (fpCoeff[1]),
              "m" (fpCoeff[2]),
              "m" (fpCoeff[3]),
              "m" (fpCoeff[4])
            : "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7",
              "memory");
    }
    X86_ASM (
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


void dsp_x86_sse_iir (double *dpVect, int iDataLength, const double *dpCoeff, 
    double *dpX, double *dpY)
{
    int iDataCntr;
    
    X86_ASM (
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
    for (iDataCntr = 0; iDataCntr < iDataLength; iDataCntr++)
    {
        X86_ASM (
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
            : "=m" (dpVect[iDataCntr])
            : "m0" (dpVect[iDataCntr]),
              "m" (dpCoeff[0]),
              "m" (dpCoeff[1]),
              "m" (dpCoeff[2]),
              "m" (dpCoeff[3]),
              "m" (dpCoeff[4])
            : "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7",
              "memory");
    }
    X86_ASM (
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


void dsp_x86_3dnow_iirf_nip (float *fpDest, const float *fpSrc, int iDataLength, 
    const float *fpCoeff, float *fpX, float *fpY)
{
    int iDataCntr;
    pv2sf m64pCoeff = (pv2sf) &fpCoeff[1];
    pv2sf m64pCoeff2 = (pv2sf) &fpCoeff[3];
    pv2sf m64pX = (pv2sf) fpX;
    pv2sf m64pY = (pv2sf) fpY;

    X86_ASM (
        "movq %0, %%mm0\n\t" \
        "pswapd %%mm0, %%mm2\n\t" \
        "movd %1, %%mm3\n\t" \
        "movq %2, %%mm0\n\t" \
        "pswapd %%mm0, %%mm4\n\t" \
        "movq %3, %%mm5\n\t" \
        "movq %4, %%mm7\n\t" \
        :
        : "m" (*m64pCoeff),
          "m" (fpCoeff[0]),
          "m" (*m64pCoeff2),
          "m" (*m64pX),
          "m" (*m64pY)
        : "mm0", "mm2", "mm3", "mm4", "mm5", "mm7", "memory");
    for (iDataCntr = 0; 
        iDataCntr < iDataLength; 
        iDataCntr++)
    {
        X86_ASM (
            "pxor %%mm0, %%mm0\n\t" \
            "movd %1, %%mm6\n\t" \
            "movq %%mm5, %%mm1\n\t" \
            "pfmul %%mm2, %%mm1\n\t" \
            "pfacc %%mm1, %%mm0\n\t" \
            "movq %%mm6, %%mm1\n\t" \
            "pfmul %%mm3, %%mm1\n\t" \
            "pfacc %%mm1, %%mm0\n\t" \
            "movq %%mm7, %%mm1\n\t" \
            "pfmul %%mm4, %%mm1\n\t" \
            "pfacc %%mm1, %%mm0\n\t" \
            "pfacc %%mm0, %%mm0\n\t" \
            \
            "pswapd %%mm7, %%mm1\n\t" \
            "movq %%mm1, %%mm7\n\t" \
            "punpckldq %%mm0, %%mm7\n\t" \
            \
            "pswapd %%mm5, %%mm1\n\t" \
            "movq %%mm1, %%mm5\n\t" \
            "movq %%mm6, %%mm1\n\t" \
            "punpckldq %%mm1, %%mm5\n\t" \
            \
            "movd %%mm0, %0\n\t"
            : "=m" (fpDest[iDataCntr])
            : "m" (fpSrc[iDataCntr])
            : "mm0", "mm1", "mm2", "mm3", "mm4", "mm5", "mm6", "mm7", "memory");
    }
    X86_ASM (
        "movq %%mm5, %0\n\t" \
        "movd %%mm6, %1\n\t" \
        "movq %%mm7, %2\n\t"
        : "=m" (*m64pX),
          "=m" (fpX[2]),
          "=m" (*m64pY)
        :
        : "mm5", "mm6", "mm7", "memory");
    X86_ASM ("femms\n\t");
}


void dsp_x86_sse_iirf_nip (float *fpDest, const float *fpSrc, int iDataLength, 
    const float *fpCoeff, float *fpX, float *fpY)
{
    int iDataCntr;
    
    X86_ASM (
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
    for (iDataCntr = 0; iDataCntr < iDataLength; iDataCntr++)
    {
        X86_ASM (
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
            : "=m" (fpDest[iDataCntr])
            : "m" (fpSrc[iDataCntr]),
              "m" (fpCoeff[0]),
              "m" (fpCoeff[1]),
              "m" (fpCoeff[2]),
              "m" (fpCoeff[3]),
              "m" (fpCoeff[4])
            : "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7",
              "memory");
    }
    X86_ASM (
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


void dsp_x86_sse_iir_nip (double *dpDest, const double *dpSrc, int iDataLength, 
    const double *dpCoeff, double *dpX, double *dpY)
{
    int iDataCntr;
    
    X86_ASM (
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
    for (iDataCntr = 0; iDataCntr < iDataLength; iDataCntr++)
    {
        X86_ASM (
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
            : "=m" (dpDest[iDataCntr])
            : "m" (dpSrc[iDataCntr]),
              "m" (dpCoeff[0]),
              "m" (dpCoeff[1]),
              "m" (dpCoeff[2]),
              "m" (dpCoeff[3]),
              "m" (dpCoeff[4])
            : "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7",
              "memory");
    }
    X86_ASM (
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
