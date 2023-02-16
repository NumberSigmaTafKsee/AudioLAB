/*

    DSP operations C-wrapper
    Copyright (C) 2002-2005 Jussi Laako

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


#include "dsp/dspop.h"
#include "dsp/DSPOp.hh"
#include "dsp/IIRCascade.hh"
#include "dsp/FFTDecimator.hh"
#include "dsp/FIRDecimator.hh"
#include "dsp/IIRDecimator.hh"
#include "dsp/RecDecimator.hh"
#include "dsp/FFTInterpolator.hh"
#include "dsp/FIRInterpolator.hh"
#include "dsp/IIRInterpolator.hh"
#include "dsp/RecInterpolator.hh"
#include "dsp/Filter.hh"
#include "dsp/FlipBand.hh"
#include "dsp/ReBuffer.hh"


#ifdef DSP_X86
#include <xmmintrin.h>
#include "dsp/X86.h"

extern bool bHave3DNow;
extern int iHaveSSE;
extern unsigned int uiOldCSR;
#endif


extern "C"
{

/* clDSPOp */


void dsp_init ()
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
}


signed long dsp_roundf (float fSrc)
{
    return clDSPOp::Round(fSrc);
}


signed long dsp_round (double dSrc)
{
    return clDSPOp::Round(dSrc);
}


void dsp_addf (float *fpDest, float fSrc, long lCount)
{
    clDSPOp::Add(fpDest, fSrc, lCount);
}


void dsp_add (double *dpDest, double dSrc, long lCount)
{
    clDSPOp::Add(dpDest, dSrc, lCount);
}


void dsp_caddf (stpSCplx spDest, stSCplx sSrc, long lCount)
{
    clDSPOp::Add(spDest, sSrc, lCount);
}


void dsp_cadd (stpDCplx spDest, stDCplx sSrc, long lCount)
{
    clDSPOp::Add(spDest, sSrc, lCount);
}


void dsp_add2f (float *fpDest, const float *fpSrc, long lCount)
{
    clDSPOp::Add(fpDest, fpSrc, lCount);
}


void dsp_add2 (double *dpDest, const double *dpSrc, long lCount)
{
    clDSPOp::Add(dpDest, dpSrc, lCount);
}


void dsp_cadd2f (stpSCplx spDest, const stpSCplx spSrc, long lCount)
{
    clDSPOp::Add(spDest, spSrc, lCount);
}


void dsp_cadd2 (stpDCplx spDest, const stpDCplx spSrc, long lCount)
{
    clDSPOp::Add(spDest, spSrc, lCount);
}


void dsp_add3f (float *fpDest, const float *fpSrc1, const float *fpSrc2,
    long lCount)
{
    clDSPOp::Add(fpDest, fpSrc1, fpSrc2, lCount);
}


void dsp_add3 (double *dpDest, const double *dpSrc1, const double *dpSrc2,
    long lCount)
{
    clDSPOp::Add(dpDest, dpSrc1, dpSrc2, lCount);
}


void dsp_cadd3f (stpSCplx spDest, const stpSCplx spSrc1, const stpSCplx spSrc2,
    long lCount)
{
    clDSPOp::Add(spDest, spSrc1, spSrc2, lCount);
}


void dsp_cadd3 (stpDCplx spDest, const stpDCplx spSrc1, const stpDCplx spSrc2,
    long lCount)
{
    clDSPOp::Add(spDest, spSrc1, spSrc2, lCount);
}


void dsp_subf (float *fpDest, float fSrc, long lCount)
{
    clDSPOp::Sub(fpDest, fSrc, lCount);
}


void dsp_sub (double *dpDest, double dSrc, long lCount)
{
    clDSPOp::Sub(dpDest, dSrc, lCount);
}


void dsp_csubf (stpSCplx spDest, stSCplx sSrc, long lCount)
{
    clDSPOp::Sub(spDest, sSrc, lCount);
}


void dsp_csub (stpDCplx spDest, stDCplx sSrc, long lCount)
{
    clDSPOp::Sub(spDest, sSrc, lCount);
}


void dsp_sub2f (float *fpDest, const float *fpSrc, long lCount)
{
    clDSPOp::Sub(fpDest, fpSrc, lCount);
}


void dsp_sub2 (double *dpDest, const double *dpSrc, long lCount)
{
    clDSPOp::Sub(dpDest, dpSrc, lCount);
}


void dsp_csub2f (stpSCplx spDest, const stpSCplx spSrc, long lCount)
{
    clDSPOp::Sub(spDest, spSrc, lCount);
}


void dsp_csub2 (stpDCplx spDest, const stpDCplx spSrc, long lCount)
{
    clDSPOp::Sub(spDest, spSrc, lCount);
}


void dsp_sub3f (float *fpDest, const float *fpSrc1, const float *fpSrc2,
    long lCount)
{
    clDSPOp::Sub(fpDest, fpSrc1, fpSrc2, lCount);
}


void dsp_sub3 (double *dpDest, const double *dpSrc1, const double *dpSrc2,
    long lCount)
{
    clDSPOp::Sub(dpDest, dpSrc1, dpSrc2, lCount);
}


void dsp_csub3f (stpSCplx spDest, const stpSCplx spSrc1, const stpSCplx spSrc2,
    long lCount)
{
    clDSPOp::Sub(spDest, spSrc1, spSrc2, lCount);
}


void dsp_csub3 (stpDCplx spDest, const stpDCplx spSrc1, const stpDCplx spSrc2,
    long lCount)
{
    clDSPOp::Sub(spDest, spSrc1, spSrc2, lCount);
}


void dsp_mulf (float *fpDest, float fSrc, long lCount)
{
    clDSPOp::Mul(fpDest, fSrc, lCount);
}


void dsp_mul (double *dpDest, double dSrc, long lCount)
{
    clDSPOp::Mul(dpDest, dSrc, lCount);
}


void dsp_chmulf (stpSCplx spDest, float fSrc, long lCount)
{
    clDSPOp::Mul(spDest, fSrc, lCount);
}


void dsp_chmul (stpDCplx spDest, double dSrc, long lCount)
{
    clDSPOp::Mul(spDest, dSrc, lCount);
}


void dsp_cmulf (stpSCplx spDest, stSCplx sSrc, long lCount)
{
    clDSPOp::Mul(spDest, sSrc, lCount);
}


void dsp_cmul (stpDCplx spDest, stDCplx sSrc, long lCount)
{
    clDSPOp::Mul(spDest, sSrc, lCount);
}


void dsp_mulf_nip (float *fpDest, const float *fpSrc1, float fSrc2,
    long lCount)
{
    clDSPOp::Mul(fpDest, fpSrc1, fSrc2, lCount);
}


void dsp_mul_nip (double *dpDest, const double *dpSrc1, double dSrc2,
    long lCount)
{
    clDSPOp::Mul(dpDest, dpSrc1, dSrc2, lCount);
}


void dsp_mul2f (float *fpDest, const float *fpSrc, long lCount)
{
    clDSPOp::Mul(fpDest, fpSrc, lCount);
}


void dsp_mul2 (double *dpDest, const double *dpSrc, long lCount)
{
    clDSPOp::Mul(dpDest, dpSrc, lCount);
}


void dsp_chmul2f (stpSCplx spDest, const float *fpSrc, long lCount)
{
    clDSPOp::Mul(spDest, fpSrc, lCount);
}


void dsp_chmul2 (stpDCplx spDest, const double *dpSrc, long lCount)
{
    clDSPOp::Mul(spDest, dpSrc, lCount);
}


void dsp_cmul2f (stpSCplx spDest, const stpSCplx spSrc, long lCount)
{
    clDSPOp::Mul(spDest, spSrc, lCount);
}


void dsp_cmul2 (stpDCplx spDest, const stpDCplx spSrc, long lCount)
{
    clDSPOp::Mul(spDest, spSrc, lCount);
}


void dsp_mul3f (float *fpDest, const float *fpSrc1, const float *fpSrc2,
    long lCount)
{
    clDSPOp::Mul(fpDest, fpSrc1, fpSrc2, lCount);
}


void dsp_mul3 (double *dpDest, const double *dpSrc1, const double *dpSrc2,
    long lCount)
{
    clDSPOp::Mul(dpDest, dpSrc1, dpSrc2, lCount);
}


void dsp_cmul3f (stpSCplx spDest, const stpSCplx spSrc1, const stpSCplx spSrc2,
    long lCount)
{
    clDSPOp::Mul(spDest, spSrc1, spSrc2, lCount);
}


void dsp_cmul3 (stpDCplx spDest, const stpDCplx spSrc1, const stpDCplx spSrc2,
    long lCount)
{
    clDSPOp::Mul(spDest, spSrc1, spSrc2, lCount);
}


void dsp_ccmulf (stpSCplx spDest, const stpSCplx spSrc, long lCount)
{
    clDSPOp::MulC(spDest, spSrc, lCount);
}


void dsp_ccmul (stpDCplx spDest, const stpDCplx spSrc, long lCount)
{
    clDSPOp::MulC(spDest, spSrc, lCount);
}


void dsp_ccmulf_nip (stpSCplx spDest, const stpSCplx spSrc1,
    const stpSCplx spSrc2, long lCount)
{
    clDSPOp::MulC(spDest, spSrc1, spSrc2, lCount);
}


void dsp_ccmul_nip (stpDCplx spDest, const stpDCplx spSrc1,
    const stpDCplx spSrc2, long lCount)
{
    clDSPOp::MulC(spDest, spSrc1, spSrc2, lCount);
}


void dsp_divf (float *fpDest, float fSrc, long lCount)
{
    clDSPOp::Div(fpDest, fSrc, lCount);
}


void dsp_div (double *dpDest, double dSrc, long lCount)
{
    clDSPOp::Div(dpDest, dSrc, lCount);
}


void dsp_cdivf (stpSCplx spDest, stSCplx sSrc, long lCount)
{
    clDSPOp::Div(spDest, sSrc, lCount);
}


void dsp_cdiv (stpDCplx spDest, stDCplx sSrc, long lCount)
{
    clDSPOp::Div(spDest, sSrc, lCount);
}


void dsp_div2f (float *fpDest, const float *fpSrc, long lCount)
{
    clDSPOp::Div(fpDest, fpSrc, lCount);
}


void dsp_div2 (double *dpDest, const double *dpSrc, long lCount)
{
    clDSPOp::Div(dpDest, dpSrc, lCount);
}


void dsp_cdiv2f (stpSCplx spDest, const stpSCplx spSrc, long lCount)
{
    clDSPOp::Div(spDest, spSrc, lCount);
}


void dsp_cdiv2 (stpDCplx spDest, const stpDCplx spSrc, long lCount)
{
    clDSPOp::Div(spDest, spSrc, lCount);
}


void dsp_div3f (float *fpDest, const float *fpSrc1, const float *fpSrc2,
    long lCount)
{
    clDSPOp::Div(fpDest, fpSrc1, fpSrc2, lCount);
}


void dsp_div3 (double *dpDest, const double *dpSrc1, const double *dpSrc2,
    long lCount)
{
    clDSPOp::Div(dpDest, dpSrc1, dpSrc2, lCount);
}


void dsp_cdiv3f (stpSCplx spDest, const stpSCplx spSrc1, const stpSCplx spSrc2,
    long lCount)
{
    clDSPOp::Div(spDest, spSrc1, spSrc2, lCount);
}


void dsp_cdiv3 (stpDCplx spDest, const stpDCplx spSrc1, const stpDCplx spSrc2,
    long lCount)
{
    clDSPOp::Div(spDest, spSrc1, spSrc2, lCount);
}


void dsp_div1xf (float *fpVect, long lCount)
{
    clDSPOp::Div1x(fpVect, lCount);
}


void dsp_div1x (double *dpVect, long lCount)
{
    clDSPOp::Div1x(dpVect, lCount);
}


void dsp_div1xf_nip (float *fpDest, const float *fpSrc, long lCount)
{
    clDSPOp::Div1x(fpDest, fpSrc, lCount);
}


void dsp_div1x_nip (double *dpDest, const double *dpSrc, long lCount)
{
    clDSPOp::Div1x(dpDest, dpSrc, lCount);
}


void dsp_muladdf (float *fpVect, float fMul, float fAdd, long lCount)
{
    clDSPOp::MulAdd(fpVect, fMul, fAdd, lCount);
}


void dsp_muladd (double *dpVect, double dMul, double dAdd, long lCount)
{
    clDSPOp::MulAdd(dpVect, dMul, dAdd, lCount);
}


void dsp_muladdf_nip (float *fpDest, const float *fpSrc, 
    float fMul, float fAdd, long lCount)
{
    clDSPOp::MulAdd(fpDest, fpSrc, fMul, fAdd, lCount);
}


void dsp_muladd_nip (double *dpDest, const double *dpSrc,
    double dMul, double dAdd, long lCount)
{
    clDSPOp::MulAdd(dpDest, dpSrc, dMul, dAdd, lCount);
}


void dsp_cmuladdf (stpSCplx spDest, const stpSCplx spSrc1,
    const stpSCplx spSrc2, long lCount)
{
    clDSPOp::MulAddC(spDest, spSrc1, spSrc2, lCount);
}


void dsp_cmuladd (stpDCplx spDest, const stpDCplx spSrc1,
    const stpDCplx spSrc2, long lCount)
{
    clDSPOp::MulAddC(spDest, spSrc1, spSrc2, lCount);
}


void dsp_absf (float *fpVect, long lCount)
{
    clDSPOp::Abs(fpVect, lCount);
}


void dsp_abs (double *dpVect, long lCount)
{
    clDSPOp::Abs(dpVect, lCount);
}


void dsp_absf_nip (float *fpDest, const float *fpSrc, long lCount)
{
    clDSPOp::Abs(fpDest, fpSrc, lCount);
}


void dsp_abs_nip (double *dpDest, const double *dpSrc, long lCount)
{
    clDSPOp::Abs(dpDest, dpSrc, lCount);
}


void dsp_sqrtf (float *fpVect, long lCount)
{
    clDSPOp::Sqrt(fpVect, lCount);
}


void dsp_sqrt (double *dpVect, long lCount)
{
    clDSPOp::Sqrt(dpVect, lCount);
}


void dsp_sqrtf_nip (float *fpDest, const float *fpSrc, long lCount)
{
    clDSPOp::Sqrt(fpDest, fpSrc, lCount);
}


void dsp_sqrt_nip (double *dpDest, const double *dpSrc, long lCount)
{
    clDSPOp::Sqrt(dpDest, dpSrc, lCount);
}


void dsp_zerof (float *fpVect, long lCount)
{
    clDSPOp::Zero(fpVect, lCount);
}


void dsp_zero (double *dpVect, long lCount)
{
    clDSPOp::Zero(dpVect, lCount);
}


void dsp_czerof (stpSCplx spVect, long lCount)
{
    clDSPOp::Zero(spVect, lCount);
}


void dsp_czero (stpDCplx spVect, long lCount)
{
    clDSPOp::Zero(spVect, lCount);
}


void dsp_setf (float *fpDest, float fSrc, long lCount)
{
    clDSPOp::Set(fpDest, fSrc, lCount);
}


void dsp_set (double *dpDest, double dSrc, long lCount)
{
    clDSPOp::Set(dpDest, dSrc, lCount);
}


void dsp_csetf (stpSCplx spDest, stSCplx sSrc, long lCount)
{
    clDSPOp::Set(spDest, sSrc, lCount);
}


void dsp_cset (stpDCplx spDest, stDCplx sSrc, long lCount)
{
    clDSPOp::Set(spDest, sSrc, lCount);
}


void dsp_set2f (float *fpDest, float fSrc, long lStart, long lCount,
    long lLength)
{
    clDSPOp::Set(fpDest, fSrc, lStart, lCount, lLength);
}


void dsp_set2 (double *dpDest, double dSrc, long lStart, long lCount,
    long lLength)
{
    clDSPOp::Set(dpDest, dSrc, lStart, lCount, lLength);
}


void dsp_cset2f (stpSCplx spDest, stSCplx sSrc, long lStart, long lCount,
    long lLength)
{
    clDSPOp::Set(spDest, sSrc, lStart, lCount, lLength);
}


void dsp_cset2 (stpDCplx spDest, stDCplx sSrc, long lStart, long lCount,
    long lLength)
{
    clDSPOp::Set(spDest, sSrc, lStart, lCount, lLength);
}


void dsp_clipf (float *fpVect, float fValue, long lCount)
{
    clDSPOp::Clip(fpVect, fValue, lCount);
}


void dsp_clip (double *dpVect, double dValue, long lCount)
{
    clDSPOp::Clip(dpVect, dValue, lCount);
}


void dsp_clipf_nip (float *fpDest, const float *fpSrc, float fValue, 
    long lCount)
{
    clDSPOp::Clip(fpDest, fpSrc, fValue, lCount);
}


void dsp_clip_nip (double *dpDest, const double *dpSrc, double dValue,
    long lCount)
{
    clDSPOp::Clip(dpDest, dpSrc, dValue, lCount);
}


void dsp_clip2f (float *fpVect, float fMin, float fMax, long lCount)
{
    clDSPOp::Clip(fpVect, fMin, fMax, lCount);
}


void dsp_clip2 (double *dpVect, double dMin, double dMax, long lCount)
{
    clDSPOp::Clip(dpVect, dMin, dMax, lCount);
}


void dsp_clip2f_nip (float *fpDest, const float *fpSrc, 
    float fMin, float fMax, long lCount)
{
    clDSPOp::Clip(fpDest, fpSrc, fMin, fMax, lCount);
}


void dsp_clip2_nip (double *dpDest, const double *dpSrc,
    double dMin, double dMax, long lCount)
{
    clDSPOp::Clip(dpDest, dpSrc, dMin, dMax, lCount);
}


void dsp_clipzerof (float *fpVect, long lCount)
{
    clDSPOp::ClipZero(fpVect, lCount);
}


void dsp_clipzero (double *dpVect, long lCount)
{
    clDSPOp::ClipZero(dpVect, lCount);
}


void dsp_clipzerof_nip (float *fpDest, const float *fpSrc, long lCount)
{
    clDSPOp::ClipZero(fpDest, fpSrc, lCount);
}


void dsp_clipzero_nip (double *dpDest, const double *dpSrc, long lCount)
{
    clDSPOp::ClipZero(dpDest, dpSrc, lCount);
}


void dsp_copyf (float *fpDest, const float *fpSrc, long lCount)
{
    clDSPOp::Copy(fpDest, fpSrc, lCount);
}


void dsp_copy (double *dpDest, const double *dpSrc, long lCount)
{
    clDSPOp::Copy(dpDest, dpSrc, lCount);
}


void dsp_ccopyf (stpSCplx spDest, const stpSCplx spSrc, long lCount)
{
    clDSPOp::Copy(spDest, spSrc, lCount);
}


void dsp_ccopy (stpDCplx spDest, const stpDCplx spSrc, long lCount)
{
    clDSPOp::Copy(spDest, spSrc, lCount);
}


float dsp_convolvef (const float *fpSrc1, const float *fpSrc2, long lCount)
{
    return clDSPOp::Convolve(fpSrc1, fpSrc2, lCount);
}


double dsp_convolve (const double *dpSrc1, const double *dpSrc2, long lCount)
{
    return clDSPOp::Convolve(dpSrc1, dpSrc2, lCount);
}


void dsp_convolve2f (float *fpDest, const float *fpSrc1, const float *fpSrc2,
    long lCount)
{
    clDSPOp::Convolve(fpDest, fpSrc1, fpSrc2, lCount);
}


void dsp_convolve2 (double *dpDest, const double *dpSrc1, const double *dpSrc2,
    long lCount)
{
    clDSPOp::Convolve(dpDest, dpSrc1, dpSrc2, lCount);
}


float dsp_correlatef (const float *fpSrc1, const float *fpSrc2, long lCount)
{
    return clDSPOp::Correlate(fpSrc1, fpSrc2, lCount);
}


double dsp_correlate (const double *dpSrc1, const double *dpSrc2, long lCount)
{
    return clDSPOp::Correlate(dpSrc1, dpSrc2, lCount);
}


void dsp_correlate2f (float *fpDest, const float *fpSrc1, const float *fpSrc2,
    long lCount)
{
    clDSPOp::Correlate(fpDest, fpSrc1, fpSrc2, lCount);
}


void dsp_correlate2 (double *dpDest, const double *dpSrc1, const double *dpSrc2,
    long lCount)
{
    clDSPOp::Correlate(dpDest, dpSrc1, dpSrc2, lCount);
}


float dsp_autocorrf (const float *fpSrc, long lCount)
{
    return clDSPOp::AutoCorrelate(fpSrc, lCount);
}


double dsp_autocorr (const double *dpSrc, long lCount)
{
    return clDSPOp::AutoCorrelate(dpSrc, lCount);
}


void dsp_autocorr2f (float *fpDest, const float *fpSrc, long lCount)
{
    clDSPOp::AutoCorrelate(fpDest, fpSrc, lCount);
}


void dsp_autocorr2 (double *dpDest, const double *dpSrc, long lCount)
{
    clDSPOp::AutoCorrelate(dpDest, dpSrc, lCount);
}


float dsp_dotproductf (const float *fpSrc1, const float *fpSrc2, long lCount)
{
    return clDSPOp::DotProduct(fpSrc1, fpSrc2, lCount);
}


double dsp_dotproduct (const double *dpSrc1, const double *dpSrc2, long lCount)
{
    return clDSPOp::DotProduct(dpSrc1, dpSrc2, lCount);
}


void dsp_minmaxf (float *fpMin, float *fpMax, const float *fpSrc, long lCount)
{
    clDSPOp::MinMax(fpMin, fpMax, fpSrc, lCount);
}


void dsp_minmax (double *dpMin, double *dpMax, const double *dpSrc, long lCount)
{
    clDSPOp::MinMax(dpMin, dpMax, dpSrc, lCount);
}


float dsp_meanf (const float *fpSrc, long lCount)
{
    return clDSPOp::Mean(fpSrc, lCount);
}


double dsp_mean (const double *dpSrc, long lCount)
{
    return clDSPOp::Mean(dpSrc, lCount);
}


float dsp_medianf (const float *fpSrc, long lCount)
{
    return clDSPOp::Median(fpSrc, lCount);
}


double dsp_median (const double *dpSrc, long lCount)
{
    return clDSPOp::Median(dpSrc, lCount);
}


void dsp_negatef (float *fpVect, long lCount)
{
    clDSPOp::Negate(fpVect, lCount);
}


void dsp_negate (double *dpVect, long lCount)
{
    clDSPOp::Negate(dpVect, lCount);
}


void dsp_negatef_nip (float *fpDest, const float *fpSrc, long lCount)
{
    clDSPOp::Negate(fpDest, fpSrc, lCount);
}


void dsp_negate_nip (double *dpDest, const double *dpSrc, long lCount)
{
    clDSPOp::Negate(dpDest, dpSrc, lCount);
}


void dsp_normalizef (float *fpVect, long lCount)
{
    clDSPOp::Normalize(fpVect, lCount);
}


void dsp_normalize (double *dpVect, long lCount)
{
    clDSPOp::Normalize(dpVect, lCount);
}


void dsp_normalizef_nip (float *fpDest, const float *fpSrc, long lCount)
{
    clDSPOp::Normalize(fpDest, fpSrc, lCount);
}


void dsp_normalize_nip (double *dpDest, const double *dpSrc, long lCount)
{
    clDSPOp::Normalize(dpDest, dpSrc, lCount);
}


float dsp_productf (const float *fpSrc, long lCount)
{
    return clDSPOp::Product(fpSrc, lCount);
}


double dsp_product (const double *dpSrc, long lCount)
{
    return clDSPOp::Product(dpSrc, lCount);
}


void dsp_reversef (float *fpVect, long lCount)
{
    clDSPOp::Reverse(fpVect, lCount);
}


void dsp_reverse (double *dpVect, long lCount)
{
    clDSPOp::Reverse(dpVect, lCount);
}


void dsp_reversef_nip (float *fpDest, const float *fpSrc, long lCount)
{
    clDSPOp::Reverse(fpDest, fpSrc, lCount);
}


void dsp_reverse_nip (double *dpDest, const double *dpSrc, long lCount)
{
    clDSPOp::Reverse(dpDest, dpSrc, lCount);
}


void dsp_scalef (float *fpVect, long lCount)
{
    clDSPOp::Scale(fpVect, lCount);
}


void dsp_scale (double *dpVect, long lCount)
{
    clDSPOp::Scale(dpVect, lCount);
}


void dsp_scalef_nip (float *fpDest, const float *fpSrc, long lCount)
{
    clDSPOp::Scale(fpDest, fpSrc, lCount);
}


void dsp_scale_nip (double *dpDest, const double *dpSrc, long lCount)
{
    clDSPOp::Scale(dpDest, dpSrc, lCount);
}


void dsp_scale01f (float *fpVect, long lCount)
{
    clDSPOp::Scale01(fpVect, lCount);
}


void dsp_scale01 (double *dpVect, long lCount)
{
    clDSPOp::Scale01(dpVect, lCount);
}


void dsp_scale01f_nip (float *fpDest, const float *fpSrc, long lCount)
{
    clDSPOp::Scale01(fpDest, fpSrc, lCount);
}


void dsp_scale01_nip (double *dpDest, const double *dpSrc, long lCount)
{
    clDSPOp::Scale01(dpDest, dpSrc, lCount);
}


void dsp_sortf (float *fpVect, long lCount)
{
    clDSPOp::Sort(fpVect, lCount);
}


void dsp_sort (double *dpVect, long lCount)
{
    clDSPOp::Sort(dpVect, lCount);
}


void dsp_sortl (long *lpVect, long lCount)
{
    clDSPOp::Sort(lpVect, lCount);
}


void dsp_stddevf (float *fpStdDev, float *fpMean, const float *fpSrc,
    long lCount)
{
    clDSPOp::StdDev(fpStdDev, fpMean, fpSrc, lCount);
}


void dsp_stddev (double *dpStdDev, double *dpMean, const double *dpSrc,
    long lCount)
{
    clDSPOp::StdDev(dpStdDev, dpMean, dpSrc, lCount);
}


float dsp_sumf (const float *fpSrc, long lCount)
{
    return clDSPOp::Sum(fpSrc, lCount);
}


double dsp_sum (const double *dpSrc, long lCount)
{
    return clDSPOp::Sum(dpSrc, lCount);
}


void dsp_squaref (float *fpVect, long lCount)
{
    clDSPOp::Square(fpVect, lCount);
}


void dsp_square (double *dpVect, long lCount)
{
    clDSPOp::Square(dpVect, lCount);
}


void dsp_squaref_nip (float *fpDest, const float *fpSrc, long lCount)
{
    clDSPOp::Square(fpDest, fpSrc, lCount);
}


void dsp_square_nip (double *dpDest, const double *dpSrc, long lCount)
{
    clDSPOp::Square(dpDest, dpSrc, lCount);
}


void dsp_convertu8f (float *fpDest, const unsigned char *ucpSrc, long lCount)
{
    clDSPOp::Convert(fpDest, ucpSrc, lCount);
}


void dsp_convertu8 (double *dpDest, const unsigned char *ucpSrc, long lCount)
{
    clDSPOp::Convert(dpDest, ucpSrc, lCount);
}


void dsp_converts16f (float *fpDest, const signed short *ipSrc, long lCount,
    int i12bit)
{
    clDSPOp::Convert(fpDest, ipSrc, lCount, (i12bit) ? true : false);
}


void dsp_converts16 (double *dpDest, const signed short *ipSrc, long lCount,
    int i12bit)
{
    clDSPOp::Convert(dpDest, ipSrc, lCount, (i12bit) ? true : false);
}


void dsp_converts32f (float *fpDest, const signed int *ipSrc, long lCount,
    int i24bit)
{
    clDSPOp::Convert(fpDest, ipSrc, lCount, (i24bit) ? true : false);
}


void dsp_converts32 (double *dpDest, const signed int *ipSrc, long lCount,
    int i24bit)
{
    clDSPOp::Convert(dpDest, ipSrc, lCount, (i24bit) ? true : false);
}


void dsp_convertd64f (float *fpDest, const double *dpSrc, long lCount)
{
    clDSPOp::Convert(fpDest, dpSrc, lCount);
}


void dsp_convertf32 (double *dpDest, const float *fpSrc, long lCount)
{
    clDSPOp::Convert(dpDest, fpSrc, lCount);
}


void dsp_convertf32c (unsigned char *ucpDest, const float *fpSrc, long lCount)
{
    clDSPOp::Convert(ucpDest, fpSrc, lCount);
}


void dsp_convertd64c (unsigned char *ucpDest, const double *dpSrc, long lCount)
{
    clDSPOp::Convert(ucpDest, dpSrc, lCount);
}


void dsp_convertf32s (signed short *ipDest, const float *fpSrc, long lCount,
    int i12bit)
{
    clDSPOp::Convert(ipDest, fpSrc, lCount, (i12bit) ? true : false);
}


void dsp_convertd64s (signed short *ipDest, const double *dpSrc, long lCount,
    int i12bit)
{
    clDSPOp::Convert(ipDest, dpSrc, lCount, (i12bit) ? true : false);
}


void dsp_convertf32i (signed int *ipDest, const float *fpSrc, long lCount,
    int i24bit)
{
    clDSPOp::Convert(ipDest, fpSrc, lCount, (i24bit) ? true : false);
}


void dsp_convertd64i (signed int *ipDest, const double *dpSrc, long lCount,
    int i24bit)
{
    clDSPOp::Convert(ipDest, dpSrc, lCount, (i24bit) ? true : false);
}


void dsp_cart2polarf (float *fpMagn, float *fpPhase, const float *fpReal,
    const float *fpImag, long lCount)
{
    clDSPOp::CartToPolar(fpMagn, fpPhase, fpReal, fpImag, lCount);
}


void dsp_cart2polar (double *dpMagn, double *dpPhase, const double *dpReal,
    const double *dpImag, long lCount)
{
    clDSPOp::CartToPolar(dpMagn, dpPhase, dpReal, dpImag, lCount);
}


void dsp_cart2polar2f (float *fpMagn, float *fpPhase, const stpSCplx spCart,
    long lCount)
{
    clDSPOp::CartToPolar(fpMagn, fpPhase, spCart, lCount);
}


void dsp_cart2polar2 (double *dpMagn, double *dpPhase, const stpDCplx spCart,
    long lCount)
{
    clDSPOp::CartToPolar(dpMagn, dpPhase, spCart, lCount);
}


void dsp_cart2polar3f (stpSPolar spPolar, const stpSCplx spCart, long lCount)
{
    clDSPOp::CartToPolar(spPolar, spCart, lCount);
}


void dsp_cart2polar3 (stpDPolar spPolar, const stpDCplx spCart, long lCount)
{
    clDSPOp::CartToPolar(spPolar, spCart, lCount);
}


void dsp_cart2polar4f (utpSCoord upCoord, long lCount)
{
    clDSPOp::CartToPolar(upCoord, lCount);
}


void dsp_cart2polar4 (utpDCoord upCoord, long lCount)
{
    clDSPOp::CartToPolar(upCoord, lCount);
}


void dsp_polar2cartf (float *fpReal, float *fpImag, const float *fpMagn,
    const float *fpPhase, long lCount)
{
    clDSPOp::PolarToCart(fpReal, fpImag, fpMagn, fpPhase, lCount);
}


void dsp_polar2cart (double *dpReal, double *dpImag, const double *dpMagn,
    const double *dpPhase, long lCount)
{
    clDSPOp::PolarToCart(dpReal, dpImag, dpMagn, dpPhase, lCount);
}


void dsp_polar2cart2f (stpSCplx spCart, const float *fpMagn, 
    const float *fpPhase, long lCount)
{
    clDSPOp::PolarToCart(spCart, fpMagn, fpPhase, lCount);
}


void dsp_polar2cart2 (stpDCplx spCart, const double *dpMagn,
    const double *dpPhase, long lCount)
{
    clDSPOp::PolarToCart(spCart, dpMagn, dpPhase, lCount);
}


void dsp_polar2cart3f (stpSCplx spCart, const stpSPolar spPolar, long lCount)
{
    clDSPOp::PolarToCart(spCart, spPolar, lCount);
}


void dsp_polar2cart3 (stpDCplx spCart, const stpDPolar spPolar, long lCount)
{
    clDSPOp::PolarToCart(spCart, spPolar, lCount);
}


void dsp_polar2cart4f (utpSCoord upCoord, long lCount)
{
    clDSPOp::PolarToCart(upCoord, lCount);
}


void dsp_polar2cart4 (utpDCoord upCoord, long lCount)
{
    clDSPOp::PolarToCart(upCoord, lCount);
}


float dsp_crosscorrf (const float *fpSrc1, const float *fpSrc2, long lCount)
{
    return clDSPOp::CrossCorr(fpSrc1, fpSrc2, lCount);
}


double dsp_crosscorr (const double *dpSrc1, const double *dpSrc2, long lCount)
{
    return clDSPOp::CrossCorr(dpSrc1, dpSrc2, lCount);
}


float dsp_crosscorr2f (const float *fpSrc1, const float *fpSrc2,
    long lDelay, long lCount)
{
    return clDSPOp::DelCrossCorr(fpSrc1, fpSrc2, lDelay, lCount);
}


double dsp_crosscorr2 (const double *dpSrc1, const double *dpSrc2,
    long lDelay, long lCount)
{
    return clDSPOp::DelCrossCorr(dpSrc1, dpSrc2, lDelay, lCount);
}


void dsp_crosscorr3f (float *fpDest, const float *fpSrc1, 
    const float *fpSrc2, long lCount, const long *lpDelays, long lDelayCount)
{
    clDSPOp::DelCrossCorr(fpDest, fpSrc1, fpSrc2, lCount, 
        lpDelays, lDelayCount);
}


void dsp_crosscorr3 (double *dpDest, const double *dpSrc1,
    const double *dpSrc2, long lCount, const long *lpDelays, long lDelayCount)
{
    clDSPOp::DelCrossCorr(dpDest, dpSrc1, dpSrc2, lCount,
        lpDelays, lDelayCount);
}


float dsp_energyf (const float *fpSrc, long lCount)
{
    return clDSPOp::Energy(fpSrc, lCount);
}


double dsp_energy (const double *dpSrc, long lCount)
{
    return clDSPOp::Energy(dpSrc, lCount);
}


void dsp_magnitudef (float *fpDest, const stpSCplx spSrc, long lCount)
{
    clDSPOp::Magnitude(fpDest, spSrc, lCount);
}


void dsp_magnitude (double *dpDest, const stpDCplx spSrc, long lCount)
{
    clDSPOp::Magnitude(dpDest, spSrc, lCount);
}


void dsp_powerf (float *fpDest, const stpSCplx spSrc, long lCount)
{
    clDSPOp::Power(fpDest, spSrc, lCount);
}


void dsp_power (double *dpDest, const stpDCplx spSrc, long lCount)
{
    clDSPOp::Power(dpDest, spSrc, lCount);
}


void dsp_phasef (float *fpDest, const stpSCplx spSrc, long lCount)
{
    clDSPOp::Phase(fpDest, spSrc, lCount);
}


void dsp_phase (double *dpDest, const stpDCplx spSrc, long lCount)
{
    clDSPOp::Phase(dpDest, spSrc, lCount);
}


void dsp_powerphasef (float *fpPower, float *fpPhase, const stpSCplx spSrc,
    long lCount)
{
    clDSPOp::PowerPhase(fpPower, fpPhase, spSrc, lCount);
}


void dsp_powerphase (double *dpPower, double *dpPhase, const stpDCplx spSrc,
    long lCount)
{
    clDSPOp::PowerPhase(dpPower, dpPhase, spSrc, lCount);
}


void dsp_decimatef (float *fpDest, const float *fpSrc, long lFactor,
    long lCount)
{
    clDSPOp::Decimate(fpDest, fpSrc, lFactor, lCount);
}


void dsp_decimate (double *dpDest, const double *dpSrc, long lFactor,
    long lCount)
{
    clDSPOp::Decimate(dpDest, dpSrc, lFactor, lCount);
}


void dsp_decimateavgf (float *fpDest, const float *fpSrc, long lFactor,
    long lCount)
{
    clDSPOp::DecimateAvg(fpDest, fpSrc, lFactor, lCount);
}


void dsp_decimateavg (double *dpDest, const double *dpSrc, long lFactor,
    long lCount)
{
    clDSPOp::DecimateAvg(dpDest, dpSrc, lFactor, lCount);
}


void dsp_interpolatef (float *fpDest, const float *fpSrc, long lFactor,
    long lCount)
{
    clDSPOp::Interpolate(fpDest, fpSrc, lFactor, lCount);
}


void dsp_interpolate (double *dpDest, const double *dpSrc, long lFactor,
    long lCount)
{
    clDSPOp::Interpolate(dpDest, dpSrc, lFactor, lCount);
}


void dsp_resamplef (float *fpDest, long lDestCount, const float *fpSrc,
    long lSrcCount)
{
    clDSPOp::Resample(fpDest, lDestCount, fpSrc, lSrcCount);
}


void dsp_resample (double *dpDest, long lDestCount, const double *dpSrc,
    long lSrcCount)
{
    clDSPOp::Resample(dpDest, lDestCount, dpSrc, lSrcCount);
}


void dsp_resampleavgf (float *fpDest, long lDestCount, const float *fpSrc,
    long lSrcCount)
{
    clDSPOp::ResampleAvg(fpDest, lDestCount, fpSrc, lSrcCount);
}


void dsp_resampleavg (double *dpDest, long lDestCount, const double *dpSrc,
    long lSrcCount)
{
    clDSPOp::ResampleAvg(dpDest, lDestCount, dpSrc, lSrcCount);
}


void dsp_interpolateavgf (float *fpDest, const float *fpSrc, long lFactor,
    long lCount)
{
    clDSPOp::InterpolateAvg(fpDest, fpSrc, lFactor, lCount);
}


void dsp_interpolateavg (double *dpDest, const double *dpSrc, long lFactor,
    long lCount)
{
    clDSPOp::InterpolateAvg(dpDest, dpSrc, lFactor, lCount);
}


float dsp_rmsf (const float *fpSrc, long lCount)
{
    return clDSPOp::RMS(fpSrc, lCount);
}


double dsp_rms (const double *dpSrc, long lCount)
{
    return clDSPOp::RMS(dpSrc, lCount);
}


float dsp_variancef (float *fpVariance, float *fpMean, const float *fpSrc,
    long lCount)
{
    return clDSPOp::Variance(fpVariance, fpMean, fpSrc, lCount);
}


double dsp_variance (double *dpVariance, double *dpMean, const double *dpSrc,
    long lCount)
{
    return clDSPOp::Variance(dpVariance, dpMean, dpSrc, lCount);
}


float dsp_peaklevelf (const float *fpSrc, long lCount)
{
    return clDSPOp::PeakLevel(fpSrc, lCount);
}


double dsp_peaklevel (const double *dpSrc, long lCount)
{
    return clDSPOp::PeakLevel(dpSrc, lCount);
}


void dsp_mixf (float *fpDest, const float *fpSrc, long lCount)
{
    clDSPOp::Mix(fpDest, fpSrc, lCount);
}


void dsp_mix (double *dpDest, const double *dpSrc, long lCount)
{
    clDSPOp::Mix(dpDest, dpSrc, lCount);
}


void dsp_mix2f (float *fpDest, const float *fpSrc1, const float *fpSrc2,
    long lCount)
{
    clDSPOp::Mix(fpDest, fpSrc1, fpSrc2, lCount);
}


void dsp_mix2 (double *dpDest, const double *dpSrc1, const double *dpSrc2,
    long lCount)
{
    clDSPOp::Mix(dpDest, dpSrc1, dpSrc2, lCount);
}


void dsp_mixnf (float *fpDest, const float *fpSrc, long lChCount,
    long lCount)
{
    clDSPOp::Mix(fpDest, fpSrc, lChCount, lCount);
}


void dsp_mixn (double *dpDest, const double *dpSrc, long lChCount,
    long lCount)
{
    clDSPOp::Mix(dpDest, dpSrc, lChCount, lCount);
}


void dsp_extractf (float *fpDest, const float *fpSrc, long lChannel,
    long lChCount, long lCount)
{
    clDSPOp::Extract(fpDest, fpSrc, lChannel, lChCount, lCount);
}


void dsp_extract (double *dpDest, const double *dpSrc, long lChannel,
    long lChCount, long lCount)
{
    clDSPOp::Extract(dpDest, dpSrc, lChannel, lChCount, lCount);
}


void dsp_packf (float *fpDest, const float *fpSrc, long lChannel,
    long lChCount, long lCount)
{
    clDSPOp::Pack(fpDest, fpSrc, lChannel, lChCount, lCount);
}


void dsp_pack (double *dpDest, const double *dpSrc, long lChannel,
    long lChCount, long lCount)
{
    clDSPOp::Pack(dpDest, dpSrc, lChannel, lChCount, lCount);
}


void dsp_fftw_convertf2cf (stpSCplx spDest, const float *fpSrc, long lCount)
{
    clDSPOp::FFTWConvert(spDest, fpSrc, lCount);
}


void dsp_fftw_convertf2cd (stpDCplx spDest, const float *fpSrc, long lCount)
{
    clDSPOp::FFTWConvert(spDest, fpSrc, lCount);
}


void dsp_fftw_convertd2cf (stpSCplx spDest, const double *dpSrc, long lCount)
{
    clDSPOp::FFTWConvert(spDest, dpSrc, lCount);
}


void dsp_fftw_convertd2cd (stpDCplx spDest, const double *dpSrc, long lCount)
{
    clDSPOp::FFTWConvert(spDest, dpSrc, lCount);
}


void dsp_fftw_convertcf2f (float *fpDest, const stpSCplx spSrc, long lCount)
{
    clDSPOp::FFTWConvert(fpDest, spSrc, lCount);
}


void dsp_fftw_convertcd2f (float *fpDest, const stpDCplx spSrc, long lCount)
{
    clDSPOp::FFTWConvert(fpDest, spSrc, lCount);
}


void dsp_fftw_convertcf2d (double *dpDest, const stpSCplx spSrc, long lCount)
{
    clDSPOp::FFTWConvert(dpDest, spSrc, lCount);
}


void dsp_fftw_convertcd2d (double *dpDest, const stpDCplx spSrc, long lCount)
{
    clDSPOp::FFTWConvert(dpDest, spSrc, lCount);
}


dsp_t dsp_new ()
{
    return ((dsp_t) new clDSPOp);
}


void dsp_delete (dsp_t dspInst)
{
    delete ((clDSPOp *) dspInst);
}


void dsp_win_bartlettf (dsp_t dspInst, float *fpDest, long lCount)
{
    ((clDSPOp *) dspInst)->WinBartlett(fpDest, lCount);
}


void dsp_win_bartlett (dsp_t dspInst, double *dpDest, long lCount)
{
    ((clDSPOp *) dspInst)->WinBartlett(dpDest, lCount);
}


void dsp_win_blackmanf (dsp_t dspInst, float *fpDest, long lCount, 
    float fAlpha)
{
    ((clDSPOp *) dspInst)->WinBlackman(fpDest, lCount, fAlpha);
}


void dsp_win_blackman (dsp_t dspInst, double *dpDest, long lCount,
    double dAlpha)
{
    ((clDSPOp *) dspInst)->WinBlackman(dpDest, lCount, dAlpha);
}


void dsp_win_blackman_harrisf (dsp_t dspInst, float *fpDest, long lCount)
{
    ((clDSPOp *) dspInst)->WinBlackmanHarris(fpDest, lCount);
}


void dsp_win_blackman_harris (dsp_t dspInst, double *dpDest, long lCount)
{
    ((clDSPOp *) dspInst)->WinBlackmanHarris(dpDest, lCount);
}


void dsp_win_cos_taperedf (dsp_t dspInst, float *fpDest, long lCount)
{
    ((clDSPOp *) dspInst)->WinCosTapered(fpDest, lCount);
}


void dsp_win_cos_tapered (dsp_t dspInst, double *dpDest, long lCount)
{
    ((clDSPOp *) dspInst)->WinCosTapered(dpDest, lCount);
}


void dsp_win_exact_blackmanf (dsp_t dspInst, float *fpDest, long lCount)
{
    ((clDSPOp *) dspInst)->WinExactBlackman(fpDest, lCount);
}


void dsp_win_exact_blackman (dsp_t dspInst, double *dpDest, long lCount)
{
    ((clDSPOp *) dspInst)->WinExactBlackman(dpDest, lCount);
}


void dsp_win_flat_topf (dsp_t dspInst, float *fpDest, long lCount)
{
    ((clDSPOp *) dspInst)->WinFlatTop(fpDest, lCount);
}


void dsp_win_flat_top (dsp_t dspInst, double *dpDest, long lCount)
{
    ((clDSPOp *) dspInst)->WinFlatTop(dpDest, lCount);
}


void dsp_win_generic_cosf (dsp_t dspInst, float *fpDest, long lCount,
    const float *fpCoeff, long lCoeffCount)
{
    ((clDSPOp *) dspInst)->WinGenericCos(fpDest, lCount, fpCoeff,
        lCoeffCount);
}


void dsp_win_generic_cos (dsp_t dspInst, double *dpDest, long lCount,
    const double *dpCoeff, long lCoeffCount)
{
    ((clDSPOp *) dspInst)->WinGenericCos(dpDest, lCount, dpCoeff,
        lCoeffCount);
}


void dsp_win_hammingf (dsp_t dspInst, float *fpDest, long lCount)
{
    ((clDSPOp *) dspInst)->WinHamming(fpDest, lCount);
}


void dsp_win_hamming (dsp_t dspInst, double *dpDest, long lCount)
{
    ((clDSPOp *) dspInst)->WinHamming(dpDest, lCount);
}


void dsp_win_hanningf (dsp_t dspInst, float *fpDest, long lCount)
{
    ((clDSPOp *) dspInst)->WinHanning(fpDest, lCount);
}


void dsp_win_hanning (dsp_t dspInst, double *dpDest, long lCount)
{
    ((clDSPOp *) dspInst)->WinHanning(dpDest, lCount);
}


void dsp_win_kaiserf (dsp_t dspInst, float *fpDest, float fBeta, long lCount)
{
    ((clDSPOp *) dspInst)->WinKaiser(fpDest, fBeta, lCount);
}


void dsp_win_kaiser (dsp_t dspInst, double *dpDest, double dBeta, long lCount)
{
    ((clDSPOp *) dspInst)->WinKaiser(dpDest, dBeta, lCount);
}


void dsp_win_kaiser_besself (dsp_t dspInst, float *fpDest, float fAlpha,
    long lCount)
{
    ((clDSPOp *) dspInst)->WinKaiserBessel(fpDest, fAlpha, lCount);
}


void dsp_win_kaiser_bessel (dsp_t dspInst, double *dpDest, double dAlpha,
    long lCount)
{
    ((clDSPOp *) dspInst)->WinKaiserBessel(dpDest, dAlpha, lCount);
}


void dsp_win_tukeyf (dsp_t dspInst, float *fpDest, long lCount)
{
    ((clDSPOp *) dspInst)->WinTukey(fpDest, lCount);
}


void dsp_win_tukey (dsp_t dspInst, double *dpDest, long lCount)
{
    ((clDSPOp *) dspInst)->WinTukey(dpDest, lCount);
}


void dsp_win_dolph_chebyshevf (dsp_t dspInst, float *fpDest, float fGamma,
    long lCount)
{
    ((clDSPOp *) dspInst)->WinDolphChebyshev(fpDest, fGamma, lCount);
}


void dsp_win_dolph_chebyshev (dsp_t dspInst, double *dpDest, double dGamma,
    long lCount)
{
    ((clDSPOp *) dspInst)->WinDolphChebyshev(dpDest, dGamma, lCount);
}


long dsp_rebufferf (dsp_t dspInst, float *fpDest, const float *fpSrc,
    long lDestCount, long lSrcCount)
{
    return ((clDSPOp *) dspInst)->ReBuffer(fpDest, fpSrc, lDestCount,
        lSrcCount);
}


long dsp_rebuffer (dsp_t dspInst, double *dpDest, const double *dpSrc,
    long lDestCount, long lSrcCount)
{
    return ((clDSPOp *) dspInst)->ReBuffer(dpDest, dpSrc, lDestCount,
        lSrcCount);
}


float dsp_deg2radf (dsp_t dspInst, float fSource)
{
    return ((clDSPOp *) dspInst)->DegToRad(fSource);
}


double dsp_deg2rad (dsp_t dspInst, double dSource)
{
    return ((clDSPOp *) dspInst)->DegToRad(dSource);
}


float dsp_rad2degf (dsp_t dspInst, float fSource)
{
    return ((clDSPOp *) dspInst)->RadToDeg(fSource);
}


double dsp_rad2deg (dsp_t dspInst, double dSource)
{
    return ((clDSPOp *) dspInst)->RadToDeg(dSource);
}


void dsp_fir_allocatef (dsp_t dspInst, const float *fpCoeff, long lCount)
{
    ((clDSPOp *) dspInst)->FIRAllocate(fpCoeff, lCount);
}


void dsp_fir_allocate (dsp_t dspInst, const double *dpCoeff, long lCount)
{
    ((clDSPOp *) dspInst)->FIRAllocate(dpCoeff, lCount);
}


void dsp_fir_free (dsp_t dspInst)
{
    ((clDSPOp *) dspInst)->FIRFree();
}


void dsp_fir_filterf (dsp_t dspInst, float *fpVect, long lCount)
{
    ((clDSPOp *) dspInst)->FIRFilter(fpVect, lCount);
}


void dsp_fir_filter (dsp_t dspInst, double *dpVect, long lCount)
{
    ((clDSPOp *) dspInst)->FIRFilter(dpVect, lCount);
}


void dsp_fir_filterf_nip (dsp_t dspInst, float *fpDest, const float *fpSrc,
    long lCount)
{
    ((clDSPOp *) dspInst)->FIRFilter(fpDest, fpSrc, lCount);
}


void dsp_fir_filter_nip (dsp_t dspInst, double *dpDest, const double *dpSrc,
    long lCount)
{
    ((clDSPOp *) dspInst)->FIRFilter(dpDest, dpSrc, lCount);
}


void dsp_fir_filterf_fst (dsp_t dspInst, float *fpDest, float *fpSrc,
    long lCount)
{
    ((clDSPOp *) dspInst)->FIRFilterF(fpDest, fpSrc, lCount);
}


void dsp_fir_filter_fst (dsp_t dspInst, double *dpDest, double *dpSrc,
    long lCount)
{
    ((clDSPOp *) dspInst)->FIRFilterF(dpDest, dpSrc, lCount);
}


void dsp_iir_initf (dsp_t dspInst, const float *fpCoeffs)
{
    ((clDSPOp *) dspInst)->IIRInitialize(fpCoeffs);
}


void dsp_iir_init (dsp_t dspInst, const double *dpCoeffs)
{
    ((clDSPOp *) dspInst)->IIRInitialize(dpCoeffs);
}


void dsp_iir_filterf (dsp_t dspInst, float *fpVect, long lCount)
{
    ((clDSPOp *) dspInst)->IIRFilter(fpVect, lCount);
}


void dsp_iir_filter (dsp_t dspInst, double *dpVect, long lCount)
{
    ((clDSPOp *) dspInst)->IIRFilter(dpVect, lCount);
}


void dsp_iir_filterf_nip (dsp_t dspInst, float *fpDest, const float *fpSrc,
    long lCount)
{
    ((clDSPOp *) dspInst)->IIRFilter(fpDest, fpSrc, lCount);
}


void dsp_iir_filter_nip (dsp_t dspInst, double *dpDest, const double *dpSrc,
    long lCount)
{
    ((clDSPOp *) dspInst)->IIRFilter(dpDest, dpSrc, lCount);
}


void dsp_iir_clear (dsp_t dspInst)
{
    ((clDSPOp *) dspInst)->IIRClear();
}


void dsp_fft_init (dsp_t dspInst, long lSize, int iReal)
{
    ((clDSPOp *) dspInst)->FFTInitialize(lSize, (iReal) ? true : false);
}


void dsp_fft_uninit (dsp_t dspInst)
{
    ((clDSPOp *) dspInst)->FFTUninitialize();
}


void dsp_fftf (dsp_t dspInst, stpSCplx spDest, float *fpSrc)
{
    ((clDSPOp *) dspInst)->FFTi(spDest, fpSrc);
}


void dsp_fft (dsp_t dspInst, stpDCplx spDest, double *dpSrc)
{
    ((clDSPOp *) dspInst)->FFTi(spDest, dpSrc);
}


void dsp_fftf_nip (dsp_t dspInst, stpSCplx spDest, const float *fpSrc)
{
    ((clDSPOp *) dspInst)->FFTo(spDest, fpSrc);
}


void dsp_fft_nip (dsp_t dspInst, stpDCplx spDest, const double *dpSrc)
{
    ((clDSPOp *) dspInst)->FFTo(spDest, dpSrc);
}


void dsp_cfftf_nip (dsp_t dspInst, stpSCplx spDest, const stpSCplx spSrc)
{
    ((clDSPOp *) dspInst)->FFTo(spDest, spSrc);
}


void dsp_cfft_nip (dsp_t dspInst, stpDCplx spDest, const stpDCplx spSrc)
{
    ((clDSPOp *) dspInst)->FFTo(spDest, spSrc);
}


void dsp_ifftf_nip (dsp_t dspInst, float *fpDest, const stpSCplx spSrc)
{
    ((clDSPOp *) dspInst)->IFFTo(fpDest, spSrc);
}


void dsp_ifft_nip (dsp_t dspInst, double *dpDest, const stpDCplx spSrc)
{
    ((clDSPOp *) dspInst)->IFFTo(dpDest, spSrc);
}


void dsp_cifftf_nip (dsp_t dspInst, stpSCplx spDest, const stpSCplx spSrc)
{
    ((clDSPOp *) dspInst)->IFFTo(spDest, spSrc);
}


void dsp_cifft_nip (dsp_t dspInst, stpDCplx spDest, const stpDCplx spSrc)
{
    ((clDSPOp *) dspInst)->IFFTo(spDest, spSrc);
}


void dsp_dct_init (dsp_t dspInst, long lSize)
{
    ((clDSPOp *) dspInst)->DCTInitialize(lSize);
}


void dsp_dct_uninit (dsp_t dspInst)
{
    ((clDSPOp *) dspInst)->FFTUninitialize();
}


void dsp_dctf (dsp_t dspInst, float *fpVect)
{
    ((clDSPOp *) dspInst)->DCTi(fpVect);
}


void dsp_dct (dsp_t dspInst, double *dpVect)
{
    ((clDSPOp *) dspInst)->DCTi(dpVect);
}


void dsp_dctf_nip (dsp_t dspInst, float *fpDest, const float *fpSrc)
{
    ((clDSPOp *) dspInst)->DCTo(fpDest, fpSrc);
}


void dsp_dct_nip (dsp_t dspInst, double *dpDest, const double *dpSrc)
{
    ((clDSPOp *) dspInst)->DCTo(dpDest, dpSrc);
}


void dsp_idctf (dsp_t dspInst, float *fpVect)
{
    ((clDSPOp *) dspInst)->IDCTi(fpVect);
}


void dsp_idct (dsp_t dspInst, double *dpVect)
{
    ((clDSPOp *) dspInst)->IDCTi(dpVect);
}


void dsp_idctf_nip (dsp_t dspInst, float *fpDest, const float *fpSrc)
{
    ((clDSPOp *) dspInst)->IDCTo(fpDest, fpSrc);
}


void dsp_idct_nip (dsp_t dspInst, double *dpDest, const double *dpSrc)
{
    ((clDSPOp *) dspInst)->IDCTo(dpDest, dpSrc);
}


/* clIIRCascade */


dsp_iircas_t dsp_iir_cas_new ()
{
    return ((dsp_iircas_t) new clIIRCascade);
}


void dsp_iir_cas_delete (dsp_iircas_t dspiircasInst)
{
    delete ((clIIRCascade *) dspiircasInst);
}


int dsp_iir_cas_initf (dsp_iircas_t dspiircasInst, 
    const float fpCoeffs[][5], long lStageCount)
{
    if (!((clIIRCascade *) dspiircasInst)->Initialize(fpCoeffs, lStageCount))
        return 0;
    return 1;
}


int dsp_iir_cas_init (dsp_iircas_t dspiircasInst, 
    const double dpCoeffs[][5], long lStageCount)
{
    if (!((clIIRCascade *) dspiircasInst)->Initialize(dpCoeffs, lStageCount))
        return 0;
    return 1;
}


void dsp_iir_cas_uninit (dsp_iircas_t dspiircasInst)
{
    ((clIIRCascade *) dspiircasInst)->Uninitialize();
}


void dsp_iir_cas_processf (dsp_iircas_t dspiircasInst, 
    float *fpVect, long lCount)
{
    ((clIIRCascade *) dspiircasInst)->Process(fpVect, lCount);
}


void dsp_iir_cas_process (dsp_iircas_t dspiircasInst,
    double *dpVect, long lCount)
{
    ((clIIRCascade *) dspiircasInst)->Process(dpVect, lCount);
}


void dsp_iir_cas_processf_nip (dsp_iircas_t dspiircasInst,
    float *fpDest, const float *fpSrc, long lCount)
{
    ((clIIRCascade *) dspiircasInst)->Process(fpDest, fpSrc, lCount);
}


void dsp_iir_cas_process_nip (dsp_iircas_t dspiircasInst,
    double *dpDest, const double *dpSrc, long lCount)
{
    ((clIIRCascade *) dspiircasInst)->Process(dpDest, dpSrc, lCount);
}


void dsp_iir_cas_clear (dsp_iircas_t dspiircasInst)
{
    ((clIIRCascade *) dspiircasInst)->Clear();
}


/* clFFTDecimator */


dsp_decfft_t dsp_dec_fft_new ()
{
    return ((dsp_decfft_t) new clFFTDecimator);
}


void dsp_dec_fft_delete (dsp_decfft_t dspdecfftInst)
{
    delete ((clFFTDecimator *) dspdecfftInst);
}


int dsp_dec_fft_initf (dsp_decfft_t dspdecfftInst, long lDecFact,
    long lFiltSize, int iHighPass)
{
    float *fpNullPtr = NULL;
    
    if (!((clFFTDecimator *) dspdecfftInst)->Initialize(lDecFact, lFiltSize,
        fpNullPtr, (iHighPass) ? true : false)) return 0;
    return 1;
}


int dsp_dec_fft_init (dsp_decfft_t dspdecfftInst, long lDecFact,
    long lFiltSize, int iHighPass)
{
    double *dpNullPtr = NULL;
    
    if (!((clFFTDecimator *) dspdecfftInst)->Initialize(lDecFact, lFiltSize,
        dpNullPtr, (iHighPass) ? true : false)) return 0;
    return 1;
}


void dsp_dec_fft_uninit (dsp_decfft_t dspdecfftInst)
{
    ((clFFTDecimator *) dspdecfftInst)->Uninitialize();
}


void dsp_dec_fft_putf (dsp_decfft_t dspdecfftInst, const float *fpSrcData,
    long lSrcCount)
{
    ((clFFTDecimator *) dspdecfftInst)->Put(fpSrcData, lSrcCount);
}


void dsp_dec_fft_put (dsp_decfft_t dspdecfftInst, const double *dpSrcData,
    long lSrcCount)
{
    ((clFFTDecimator *) dspdecfftInst)->Put(dpSrcData, lSrcCount);
}


int dsp_dec_fft_getf (dsp_decfft_t dspdecfftInst, float *fpDestData,
    long lDestCount)
{
    if (!((clFFTDecimator *) dspdecfftInst)->Get(fpDestData, lDestCount))
        return 0;
    return 1;
}


int dsp_dec_fft_get (dsp_decfft_t dspdecfftInst, double *dpDestData,
    long lDestCount)
{
    if (!((clFFTDecimator *) dspdecfftInst)->Get(dpDestData, lDestCount))
        return 0;
    return 1;
}


/* clFIRDecimator */


dsp_decfir_t dsp_dec_fir_new ()
{
    return ((dsp_decfir_t) new clFIRDecimator);
}


void dsp_dec_fir_delete (dsp_decfir_t dspdecfirInst)
{
    delete ((clFIRDecimator *) dspdecfirInst);
}


int dsp_dec_fir_initf (dsp_decfir_t dspdecfirInst, long lDecFact, int iHighPass)
{
    float *fpNullPtr = NULL;
    
    if (!((clFIRDecimator *) dspdecfirInst)->Initialize(lDecFact, fpNullPtr,
        (iHighPass) ? true : false)) return 0;
    return 1;
}


int dsp_dec_fir_init (dsp_decfir_t dspdecfirInst, long lDecFact, int iHighPass)
{
    double *dpNullPtr = NULL;
    
    if (!((clFIRDecimator *) dspdecfirInst)->Initialize(lDecFact, dpNullPtr,
        (iHighPass) ? true : false)) return 0;
    return 1;
}


void dsp_dec_fir_uninit (dsp_decfir_t dspdecfirInst)
{
    ((clFIRDecimator *) dspdecfirInst)->Uninitialize();
}


void dsp_dec_fir_putf (dsp_decfir_t dspdecfirInst, const float *fpSrcData,
    long lSrcCount)
{
    ((clFIRDecimator *) dspdecfirInst)->Put(fpSrcData, lSrcCount);
}


void dsp_dec_fir_put (dsp_decfir_t dspdecfirInst, const double *dpSrcData,
    long lSrcCount)
{
    ((clFIRDecimator *) dspdecfirInst)->Put(dpSrcData, lSrcCount);
}


int dsp_dec_fir_getf (dsp_decfir_t dspdecfirInst, float *fpDestData,
    long lDestCount)
{
    if (!((clFIRDecimator *) dspdecfirInst)->Get(fpDestData, lDestCount))
        return 0;
    return 1;
}


int dsp_dec_fir_get (dsp_decfir_t dspdecfirInst, double *dpDestData,
    long lDestCount)
{
    if (!((clFIRDecimator *) dspdecfirInst)->Get(dpDestData, lDestCount))
        return 0;
    return 1;
}


/* clIIRDecimator */


dsp_deciir_t dsp_dec_iir_new ()
{
    return ((dsp_deciir_t) new clIIRDecimator);
}


void dsp_dec_iir_delete (dsp_deciir_t dspdeciirInst)
{
    delete ((clIIRDecimator *) dspdeciirInst);
}


int dsp_dec_iir_initf (dsp_deciir_t dspdeciirInst, long lDecFact, int iHighPass)
{
    float *fpNullPtr = NULL;
    
    if (!((clIIRDecimator *) dspdeciirInst)->Initialize(lDecFact, fpNullPtr,
        (iHighPass) ? true : false)) return 0;
    return 1;
}


int dsp_dec_iir_init (dsp_deciir_t dspdeciirInst, long lDecFact, int iHighPass)
{
    double *dpNullPtr = NULL;
    
    if (!((clIIRDecimator *) dspdeciirInst)->Initialize(lDecFact, dpNullPtr,
        (iHighPass) ? true : false)) return 0;
    return 1;
}


void dsp_dec_iir_uninit (dsp_deciir_t dspdeciirInst)
{
    ((clIIRDecimator *) dspdeciirInst)->Uninitialize();
}


void dsp_dec_iir_putf (dsp_deciir_t dspdeciirInst, const float *fpSrcData,
    long lSrcCount)
{
    ((clIIRDecimator *) dspdeciirInst)->Put(fpSrcData, lSrcCount);
}


void dsp_dec_iir_put (dsp_deciir_t dspdeciirInst, const double *dpSrcData,
    long lSrcCount)
{
    ((clIIRDecimator *) dspdeciirInst)->Put(dpSrcData, lSrcCount);
}


int dsp_dec_iir_getf (dsp_deciir_t dspdeciirInst, float *fpDestData,
    long lDestCount)
{
    if (!((clIIRDecimator *) dspdeciirInst)->Get(fpDestData, lDestCount))
        return 0;
    return 1;
}


int dsp_dec_iir_get (dsp_deciir_t dspdeciirInst, double *dpDestData,
    long lDestCount)
{
    if (!((clIIRDecimator *) dspdeciirInst)->Get(dpDestData, lDestCount))
        return 0;
    return 1;
}


/* clRecDecimator */


dsp_decrec_t dsp_dec_rec_new ()
{
    return ((dsp_decrec_t) new clRecDecimator);
}


void dsp_dec_rec_delete (dsp_decrec_t dspdecrecInst)
{
    delete ((clRecDecimator *) dspdecrecInst);
}


int dsp_dec_rec_initf (dsp_decrec_t dspdecrecInst, long lDecFact,
    long lFiltSize, float fBandCenter, int iFilterType)
{
    float *fpNullPtr = NULL;

    if (!((clRecDecimator *) dspdecrecInst)->Initialize(lDecFact,
        lFiltSize, fpNullPtr, fBandCenter, iFilterType))
        return 0;
    return 1;
}


int dsp_dec_rec_init (dsp_decrec_t dspdecrecInst, long lDecFact,
    long lFiltSize, double dBandCenter, int iFilterType)
{
    double *dpNullPtr = NULL;
    
    if (!((clRecDecimator *) dspdecrecInst)->Initialize(lDecFact,
        lFiltSize, dpNullPtr, dBandCenter, iFilterType))
        return 0;
    return 1;
}


void dsp_dec_rec_uninit (dsp_decrec_t dspdecrecInst)
{
    ((clRecDecimator *) dspdecrecInst)->Uninitialize();
}


void dsp_dec_rec_putf (dsp_decrec_t dspdecrecInst, const float *fpSrcData,
    long lSrcCount)
{
    ((clRecDecimator *) dspdecrecInst)->Put(fpSrcData, lSrcCount);
}


void dsp_dec_rec_put (dsp_decrec_t dspdecrecInst, const double *dpSrcData,
    long lSrcCount)
{
    ((clRecDecimator *) dspdecrecInst)->Put(dpSrcData, lSrcCount);
}


int dsp_dec_rec_getf (dsp_decrec_t dspdecrecInst, float *fpDestData,
    long lDestCount)
{
    if (!((clRecDecimator *) dspdecrecInst)->Get(fpDestData, lDestCount))
        return 0;
    return 1;
}


int dsp_dec_rec_get (dsp_decrec_t dspdecrecInst, double *dpDestData,
    long lDestCount)
{
    if (!((clRecDecimator *) dspdecrecInst)->Get(dpDestData, lDestCount))
        return 0;
    return 1;
}


/* clFFTInterpolator */


dsp_intfft_t dsp_int_fft_new ()
{
    return ((dsp_intfft_t) new clFFTInterpolator);
}


void dsp_int_fft_delete (dsp_intfft_t dspintfftInst)
{
    delete ((clFFTInterpolator *) dspintfftInst);
}


int dsp_int_fft_initf (dsp_intfft_t dspintfftInst, long lIntFact,
    long lFiltSize, int iHighPass)
{
    float *fpNullPtr = NULL;
    
    if (!((clFFTInterpolator *) dspintfftInst)->Initialize(lIntFact, lFiltSize,
        fpNullPtr, (iHighPass) ? true : false)) return 0;
    return 1;
}


int dsp_int_fft_init (dsp_intfft_t dspintfftInst, long lIntFact,
    long lFiltSize, int iHighPass)
{
    double *dpNullPtr = NULL;
    
    if (!((clFFTInterpolator *) dspintfftInst)->Initialize(lIntFact, lFiltSize,
        dpNullPtr, (iHighPass) ? true : false)) return 0;
    return 1;
}


void dsp_int_fft_uninit (dsp_intfft_t dspintfftInst)
{
    ((clFFTInterpolator *) dspintfftInst)->Uninitialize();
}


void dsp_int_fft_putf (dsp_intfft_t dspintfftInst, const float *fpSrcData,
    long lSrcCount)
{
    ((clFFTInterpolator *) dspintfftInst)->Put(fpSrcData, lSrcCount);
}


void dsp_int_fft_put (dsp_intfft_t dspintfftInst, const double *dpSrcData,
    long lSrcCount)
{
    ((clFFTInterpolator *) dspintfftInst)->Put(dpSrcData, lSrcCount);
}


int dsp_int_fft_getf (dsp_intfft_t dspintfftInst, float *fpDestData,
    long lDestCount)
{
    if (!((clFFTInterpolator *) dspintfftInst)->Get(fpDestData, lDestCount))
        return 0;
    return 1;
}


int dsp_int_fft_get (dsp_intfft_t dspintfftInst, double *dpDestData,
    long lDestCount)
{
    if (!((clFFTInterpolator *) dspintfftInst)->Get(dpDestData, lDestCount))
        return 0;
    return 1;
}


/* clFIRInterpolator */


dsp_intfir_t dsp_int_fir_new ()
{
    return ((dsp_intfir_t) new clFIRInterpolator);
}


void dsp_int_fir_delete (dsp_intfir_t dspintfirInst)
{
    delete ((clFIRInterpolator *) dspintfirInst);
}


int dsp_int_fir_initf (dsp_intfir_t dspintfirInst, long lIntFact, int iHighPass)
{
    float *fpNullPtr = NULL;
    
    if (!((clFIRInterpolator *) dspintfirInst)->Initialize(lIntFact, fpNullPtr,
        (iHighPass) ? true : false)) return 0;
    return 1;
}


int dsp_int_fir_init (dsp_intfir_t dspintfirInst, long lIntFact, int iHighPass)
{
    double *dpNullPtr = NULL;
    
    if (!((clFIRInterpolator *) dspintfirInst)->Initialize(lIntFact, dpNullPtr,
        (iHighPass) ? true : false)) return 0;
    return 1;
}


void dsp_int_fir_uninit (dsp_intfir_t dspintfirInst)
{
    ((clFIRInterpolator *) dspintfirInst)->Uninitialize();
}


void dsp_int_fir_putf (dsp_intfir_t dspintfirInst, const float *fpSrcData,
    long lSrcCount)
{
    ((clFIRInterpolator *) dspintfirInst)->Put(fpSrcData, lSrcCount);
}


void dsp_int_fir_put (dsp_intfir_t dspintfirInst, const double *dpSrcData,
    long lSrcCount)
{
    ((clFIRInterpolator *) dspintfirInst)->Put(dpSrcData, lSrcCount);
}


int dsp_int_fir_getf (dsp_intfir_t dspintfirInst, float *fpDestData,
    long lDestCount)
{
    if (!((clFIRInterpolator *) dspintfirInst)->Get(fpDestData, lDestCount))
        return 0;
    return 1;
}


int dsp_int_fir_get (dsp_intfir_t dspintfirInst, double *dpDestData,
    long lDestCount)
{
    if (!((clFIRInterpolator *) dspintfirInst)->Get(dpDestData, lDestCount))
        return 0;
    return 1;
}


/* clIIRInterpolator */


dsp_intiir_t dsp_int_iir_new ()
{
    return ((dsp_intiir_t) new clIIRInterpolator);
}


void dsp_int_iir_delete (dsp_intiir_t dspintiirInst)
{
    delete ((clIIRInterpolator *) dspintiirInst);
}


int dsp_int_iir_initf (dsp_intiir_t dspintiirInst, long lIntFact, int iHighPass)
{
    float *fpNullPtr = NULL;
    
    if (!((clIIRInterpolator *) dspintiirInst)->Initialize(lIntFact, fpNullPtr,
        (iHighPass) ? true : false)) return 0;
    return 1;
}


int dsp_int_iir_init (dsp_intiir_t dspintiirInst, long lIntFact, int iHighPass)
{
    double *dpNullPtr = NULL;
    
    if (!((clIIRInterpolator *) dspintiirInst)->Initialize(lIntFact, dpNullPtr,
        (iHighPass) ? true : false)) return 0;
    return 1;
}


void dsp_int_iir_uninit (dsp_intiir_t dspintiirInst)
{
    ((clIIRInterpolator *) dspintiirInst)->Uninitialize();
}


void dsp_int_iir_putf (dsp_intiir_t dspintiirInst, const float *fpSrcData,
    long lSrcCount)
{
    ((clIIRInterpolator *) dspintiirInst)->Put(fpSrcData, lSrcCount);
}


void dsp_int_iir_put (dsp_intiir_t dspintiirInst, const double *dpSrcData,
    long lSrcCount)
{
    ((clIIRInterpolator *) dspintiirInst)->Put(dpSrcData, lSrcCount);
}


int dsp_int_iir_getf (dsp_intiir_t dspintiirInst, float *fpDestData,
    long lDestCount)
{
    if (!((clIIRInterpolator *) dspintiirInst)->Get(fpDestData, lDestCount))
        return 0;
    return 1;
}


int dsp_int_iir_get (dsp_intiir_t dspintiirInst, double *dpDestData,
    long lDestCount)
{
    if (!((clIIRInterpolator *) dspintiirInst)->Get(dpDestData, lDestCount))
        return 0;
    return 1;
}


/* clRecInterpolator */


dsp_intrec_t dsp_int_rec_new ()
{
    return ((dsp_intrec_t) new clRecInterpolator);
}


void dsp_int_rec_delete (dsp_intrec_t dspintrecInst)
{
    delete ((clRecInterpolator *) dspintrecInst);
}


int dsp_int_rec_initf (dsp_intrec_t dspintrecInst, long lIntFact,
    long lFiltSize, float fBandCenter, int iFilterType)
{
    float *fpNullPtr = NULL;

    if (!((clRecInterpolator *) dspintrecInst)->Initialize(lIntFact,
        lFiltSize, fpNullPtr, fBandCenter, iFilterType))
        return 0;
    return 1;
}


int dsp_int_rec_init (dsp_intrec_t dspintrecInst, long lIntFact,
    long lFiltSize, double dBandCenter, int iFilterType)
{
    double *dpNullPtr = NULL;
    
    if (!((clRecInterpolator *) dspintrecInst)->Initialize(lIntFact,
        lFiltSize, dpNullPtr, dBandCenter, iFilterType))
        return 0;
    return 1;
}


void dsp_int_rec_uninit (dsp_intrec_t dspintrecInst)
{
    ((clRecInterpolator *) dspintrecInst)->Uninitialize();
}


void dsp_int_rec_putf (dsp_intrec_t dspintrecInst, const float *fpSrcData,
    long lSrcCount)
{
    ((clRecInterpolator *) dspintrecInst)->Put(fpSrcData, lSrcCount);
}


void dsp_int_rec_put (dsp_intrec_t dspintrecInst, const double *dpSrcData,
    long lSrcCount)
{
    ((clRecInterpolator *) dspintrecInst)->Put(dpSrcData, lSrcCount);
}


int dsp_int_rec_getf (dsp_intrec_t dspintrecInst, float *fpDestData,
    long lDestCount)
{
    if (!((clRecInterpolator *) dspintrecInst)->Get(fpDestData, lDestCount))
        return 0;
    return 1;
}


int dsp_int_rec_get (dsp_intrec_t dspintrecInst, double *dpDestData,
    long lDestCount)
{
    if (!((clRecInterpolator *) dspintrecInst)->Get(dpDestData, lDestCount))
        return 0;
    return 1;
}


/* clFilter */


dsp_filter_t dsp_filter_new ()
{
    return ((dsp_filter_t) new clFilter);
}


void dsp_filter_delete (dsp_filter_t dspfilterInst)
{
    delete ((clFilter *) dspfilterInst);
}


int dsp_filter_initf (dsp_filter_t dspfilterInst, long lWindowSize)
{
    float *fpNullPtr = NULL;

    if (!((clFilter *) dspfilterInst)->Initialize(lWindowSize, fpNullPtr))
        return 0;
    return 1;
}


int dsp_filter_init (dsp_filter_t dspfilterInst, long lWindowSize)
{
    double *dpNullPtr = NULL;
    
    if (!((clFilter *) dspfilterInst)->Initialize(lWindowSize, dpNullPtr))
        return 0;
    return 1;
}


int dsp_filter_init2f (dsp_filter_t dspfilterInst, long lWindowSize,
    const float *fpFiltCoeffs, float fOverlap, float fBeta)
{
    if (!((clFilter *) dspfilterInst)->Initialize(lWindowSize, fpFiltCoeffs,
        fOverlap, fBeta)) return 0;
    return 1;
}


int dsp_filter_init2 (dsp_filter_t dspfilterInst, long lWindowSize,
    const double *dpFiltCoeffs, double dOverlap, double dBeta)
{
    if (!((clFilter *) dspfilterInst)->Initialize(lWindowSize, dpFiltCoeffs,
        dOverlap, dBeta)) return 0;
    return 1;
}


int dsp_filter_init_lpf (dsp_filter_t dspfilterInst, float fPassBand,
    float fStopBand, float fRippleRatio, float fOverlap)
{
    if (!((clFilter *) dspfilterInst)->InitializeLP(fPassBand, fStopBand,
        fRippleRatio, fOverlap)) return 0;
    return 1;
}


int dsp_filter_init_lp (dsp_filter_t dspfilterInst, double dPassBand,
    double dStopBand, double dRippleRatio, double dOverlap)
{
    if (!((clFilter *) dspfilterInst)->InitializeLP(dPassBand, dStopBand,
        dRippleRatio, dOverlap)) return 0;
    return 1;
}


int dsp_filter_init_hpf (dsp_filter_t dspfilterInst, float fPassBand,
    float fStopBand, float fRippleRatio, float fOverlap)
{
    if (!((clFilter *) dspfilterInst)->InitializeHP(fPassBand, fStopBand,
        fRippleRatio, fOverlap)) return 0;
    return 1;
}


int dsp_filter_init_hp (dsp_filter_t dspfilterInst, double dPassBand,
    double dStopBand, double dRippleRatio, double dOverlap)
{
    if (!((clFilter *) dspfilterInst)->InitializeHP(dPassBand, dStopBand,
        dRippleRatio, dOverlap)) return 0;
    return 1;
}


void dsp_filter_set_coeffsf (dsp_filter_t dspfilterInst, 
    const float *fpFiltCoeffs)
{
    ((clFilter *) dspfilterInst)->SetCoeffs(fpFiltCoeffs);
}


void dsp_filter_set_coeffs (dsp_filter_t dspfilterInst,
    const double *dpFiltCoeffs)
{
    ((clFilter *) dspfilterInst)->SetCoeffs(dpFiltCoeffs);
}


void dsp_filter_set_ccoeffsf (dsp_filter_t dspfilterInst,
    const stpSCplx spFiltCoeffs, int iSmooth)
{
    ((clFilter *) dspfilterInst)->SetCoeffs(spFiltCoeffs,
        (iSmooth) ? true : false);
}


void dsp_filter_set_ccoeffs (dsp_filter_t dspfilterInst,
    const stpDCplx spFiltCoeffs, int iSmooth)
{
    ((clFilter *) dspfilterInst)->SetCoeffs(spFiltCoeffs,
        (iSmooth) ? true : false);
}


void dsp_filter_get_coeffsf (dsp_filter_t dspfilterInst,
    float *fpFiltCoeffs)
{
    ((clFilter *) dspfilterInst)->GetCoeffs(fpFiltCoeffs);
}


void dsp_filter_get_coeffs (dsp_filter_t dspfilterInst,
    double *dpFiltCoeffs)
{
    ((clFilter *) dspfilterInst)->GetCoeffs(dpFiltCoeffs);
}


void dsp_filter_get_ccoeffsf (dsp_filter_t dspfilterInst,
    stpSCplx spFiltCoeffs)
{
    ((clFilter *) dspfilterInst)->GetCoeffs(spFiltCoeffs);
}


void dsp_filter_get_ccoeffs (dsp_filter_t dspfilterInst,
    stpDCplx spFiltCoeffs)
{
    ((clFilter *) dspfilterInst)->GetCoeffs(spFiltCoeffs);
}


void dsp_filter_putf (dsp_filter_t dspfilterInst, const float *fpSrcData,
    long lSrcCount)
{
    ((clFilter *) dspfilterInst)->Put(fpSrcData, lSrcCount);
}


void dsp_filter_put (dsp_filter_t dspfilterInst, const double *dpSrcData,
    long lSrcCount)
{
    ((clFilter *) dspfilterInst)->Put(dpSrcData, lSrcCount);
}


void dsp_filter_put2f (dsp_filter_t dspfilterInst, const float *fpSrcData,
    long lSrcCount, const stpSCplx spCoeffs)
{
    ((clFilter *) dspfilterInst)->Put(fpSrcData, lSrcCount, spCoeffs);
}


void dsp_filter_put2 (dsp_filter_t dspfilterInst, const double *dpSrcData,
    long lSrcCount, const stpDCplx spCoeffs)
{
    ((clFilter *) dspfilterInst)->Put(dpSrcData, lSrcCount, spCoeffs);
}


void dsp_filter_getf (dsp_filter_t dspfilterInst, float *fpDestData,
    long lDestCount)
{
    ((clFilter *) dspfilterInst)->Get(fpDestData, lDestCount);
}


void dsp_filter_get (dsp_filter_t dspfilterInst, double *dpDestData,
    long lDestCount)
{
    ((clFilter *) dspfilterInst)->Get(dpDestData, lDestCount);
}


void dsp_filter_design_lpf (dsp_filter_t dspfilterInst, float *fpCorner,
    int iDCBlock)
{
    ((clFilter *) dspfilterInst)->DesignLP(fpCorner,
        (iDCBlock) ? true : false);
}


void dsp_filter_design_lp (dsp_filter_t dspfilterInst, double *dpCorner,
    int iDCBlock)
{
    ((clFilter *) dspfilterInst)->DesignLP(dpCorner,
        (iDCBlock) ? true : false);
}


void dsp_filter_design_lp2f (dsp_filter_t dspfilterInst, float *fpCorner,
    float fSampleRate, int iDCBlock)
{
    ((clFilter *) dspfilterInst)->DesignLP(fpCorner, fSampleRate,
        (iDCBlock) ? true : false);
}


void dsp_filter_design_lp2 (dsp_filter_t dspfilterInst, double *dpCorner,
    double dSampleRate, int iDCBlock)
{
    ((clFilter *) dspfilterInst)->DesignLP(dpCorner, dSampleRate,
        (iDCBlock) ? true : false);
}


void dsp_filter_design_hpf (dsp_filter_t dspfilterInst, float *fpCorner)
{
    ((clFilter *) dspfilterInst)->DesignHP(fpCorner);
}


void dsp_filter_design_hp (dsp_filter_t dspfilterInst, double *dpCorner)
{
    ((clFilter *) dspfilterInst)->DesignHP(dpCorner);
}


void dsp_filter_design_hp2f (dsp_filter_t dspfilterInst, float *fpCorner,
    float fSampleRate)
{
    ((clFilter *) dspfilterInst)->DesignHP(fpCorner, fSampleRate);
}


void dsp_filter_design_hp2 (dsp_filter_t dspfilterInst, double *dpCorner,
    double dSampleRate)
{
    ((clFilter *) dspfilterInst)->DesignHP(dpCorner, dSampleRate);
}


void dsp_filter_design_bpf (dsp_filter_t dspfilterInst, float *fpLowCorner,
    float *fpHighCorner)
{
    ((clFilter *) dspfilterInst)->DesignBP(fpLowCorner, fpHighCorner);
}


void dsp_filter_design_bp (dsp_filter_t dspfilterInst, double *dpLowCorner,
    double *dpHighCorner)
{
    ((clFilter *) dspfilterInst)->DesignBP(dpLowCorner, dpHighCorner);
}


void dsp_filter_design_bp2f (dsp_filter_t dspfilterInst, float *fpLowCorner,
    float *fpHighCorner, float fSampleRate)
{
    ((clFilter *) dspfilterInst)->DesignBP(fpLowCorner, fpHighCorner,
        fSampleRate);
}


void dsp_filter_design_bp2 (dsp_filter_t dspfilterInst, double *dpLowCorner,
    double *dpHighCorner, double dSampleRate)
{
    ((clFilter *) dspfilterInst)->DesignBP(dpLowCorner, dpHighCorner,
        dSampleRate);
}


void dsp_filter_design_brf (dsp_filter_t dspfilterInst, float *fpLowCorner,
    float *fpHighCorner)
{
    ((clFilter *) dspfilterInst)->DesignBR(fpLowCorner, fpHighCorner);
}


void dsp_filter_design_br (dsp_filter_t dspfilterInst, double *dpLowCorner,
    double *dpHighCorner)
{
    ((clFilter *) dspfilterInst)->DesignBR(dpLowCorner, dpHighCorner);
}


void dsp_filter_design_br2f (dsp_filter_t dspfilterInst, float *fpLowCorner,
    float *fpHighCorner, float fSampleRate)
{
    ((clFilter *) dspfilterInst)->DesignBR(fpLowCorner, fpHighCorner,
        fSampleRate);
}


void dsp_filter_design_br2 (dsp_filter_t dspfilterInst, double *dpLowCorner,
    double *dpHighCorner, double dSampleRate)
{
    ((clFilter *) dspfilterInst)->DesignBR(dpLowCorner, dpHighCorner,
        dSampleRate);
}


/* clFlipBand */


dsp_flip_t dsp_flip_new ()
{
    return ((dsp_flip_t) new clFlipBand);
}


void dsp_flip_delete (dsp_flip_t dspflipInst)
{
    delete ((clFlipBand *) dspflipInst);
}


void dsp_flip_initf (dsp_flip_t dspflipInst, long lSize)
{
    float *fpNullPtr = NULL;

    ((clFlipBand *) dspflipInst)->Initialize(lSize, fpNullPtr);
}


void dsp_flip_init (dsp_flip_t dspflipInst, long lSize)
{
    double *dpNullPtr = NULL;
    
    ((clFlipBand *) dspflipInst)->Initialize(lSize, dpNullPtr);
}


void dsp_flip_uninit (dsp_flip_t dspflipInst)
{
    ((clFlipBand *) dspflipInst)->Uninitialize();
}


void dsp_flip_putf (dsp_flip_t dspflipInst, const float *fpSrc, long lCount)
{
    ((clFlipBand *) dspflipInst)->Put(fpSrc, lCount);
}


void dsp_flip_put (dsp_flip_t dspflipInst, const double *dpSrc, long lCount)
{
    ((clFlipBand *) dspflipInst)->Put(dpSrc, lCount);
}


void dsp_flip_getf (dsp_flip_t dspflipInst, float *fpDest, long lCount)
{
    ((clFlipBand *) dspflipInst)->Get(fpDest, lCount);
}


void dsp_flip_get (dsp_flip_t dspflipInst, double *dpDest, long lCount)
{
    ((clFlipBand *) dspflipInst)->Get(dpDest, lCount);
}


void dsp_flip_clear (dsp_flip_t dspflipInst)
{
    ((clFlipBand *) dspflipInst)->Clear();
}


/* clReBuffer */


dsp_rebuf_t dsp_rebuf_new ()
{
    return ((dsp_rebuf_t) new clReBuffer);
}


void dsp_rebuf_delete (dsp_rebuf_t dsprebufInst)
{
    delete ((clReBuffer *) dsprebufInst);
}


void dsp_rebuf_putf (dsp_rebuf_t dsprebufInst, const float *fpSrcData,
    long lSrcCount)
{
    ((clReBuffer *) dsprebufInst)->Put(fpSrcData, lSrcCount);
}


void dsp_rebuf_put (dsp_rebuf_t dsprebufInst, const double *dpSrcData,
    long lSrcCount)
{
    ((clReBuffer *) dsprebufInst)->Put(dpSrcData, lSrcCount);
}


int dsp_rebuf_getf (dsp_rebuf_t dsprebufInst, float *fpDestData,
    long lDestCount)
{
    if (!((clReBuffer *) dsprebufInst)->Get(fpDestData, lDestCount))
        return 0;
    return 1;
}


int dsp_rebuf_get (dsp_rebuf_t dsprebufInst, double *dpDestData,
    long lDestCount)
{
    if (!((clReBuffer *) dsprebufInst)->Get(dpDestData, lDestCount))
        return 0;
    return 1;
}


long dsp_rebuf_size (dsp_rebuf_t dsprebufInst)
{
    return ((clReBuffer *) dsprebufInst)->GetCount();
}


void dsp_rebuf_clear (dsp_rebuf_t dsprebufInst)
{
    ((clReBuffer *) dsprebufInst)->Clear();
}


void dsp_rebuf_copy (dsp_rebuf_t dsprebufInst, dsp_rebuf_t dsprebufCopySrc)
{
    *((clReBuffer *) dsprebufInst) = *((clReBuffer *) dsprebufCopySrc);
    /* or */
    /*((clReBuffer *) dsprebufInst)->operator=(*((clReBuffer *) dsprebufCopySrc));*/
    /* if you like... :) */
}

}
