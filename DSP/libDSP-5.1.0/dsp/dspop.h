/*

    DSP operations C-wrapper, header
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


#include <dsp/dsptypes.h>

#ifndef DSPOP_H
    #define DSPOP_H

    /*
        Everything is documented in corresponding class::method's and HTML
        documentation of this is hyperlinked to corresponding methods.
    */


    /** DSP object type */
    typedef void * dsp_t;
    /** Cascaded IIR filter object type */
    typedef void * dsp_iircas_t;
    /** FFT decimator object type */
    typedef void * dsp_decfft_t;
    /** FIR decimator object type */
    typedef void * dsp_decfir_t;
    /** IIR decimator object type */
    typedef void * dsp_deciir_t;
    /** Recursive decimator object type */
    typedef void * dsp_decrec_t;
    /** FFT interpolator object type */
    typedef void * dsp_intfft_t;
    /** FIR interpolator object type */
    typedef void * dsp_intfir_t;
    /** IIR interpolator object type */
    typedef void * dsp_intiir_t;
    /** Recursive interpolator object type */
    typedef void * dsp_intrec_t;
    /** FFT filter object type */
    typedef void * dsp_filter_t;
    /** Band flipping object type */
    typedef void * dsp_flip_t;
    /** Rebuffering object type */
    typedef void * dsp_rebuf_t;


    #ifdef __cplusplus
    extern "C"
    {
    #endif

    /* clDSPOp */
    
    /**
        Initialize DSP library.
        
        This should be called first to enable possible hardware dependent
        optimizations.
    */
    void dsp_init (void);

    /** clDSPOp::Round */
    signed long dsp_roundf (float);
    signed long dsp_round (double);
    
    /** clDSPOp::Add */
    void dsp_addf (float *, float, long);
    void dsp_add (double *, double, long);
    void dsp_caddf (stpSCplx, stSCplx, long);
    void dsp_cadd (stpDCplx, stDCplx, long);
    void dsp_add2f (float *, const float *, long);
    void dsp_add2 (double *, const double *, long);
    void dsp_cadd2f (stpSCplx, const stpSCplx, long);
    void dsp_cadd2 (stpDCplx, const stpDCplx, long);
    void dsp_add3f (float *, const float *, const float *, long);
    void dsp_add3 (double *, const double *, const double *, long);
    void dsp_cadd3f (stpSCplx, const stpSCplx, const stpSCplx, long);
    void dsp_cadd3 (stpDCplx, const stpDCplx, const stpDCplx, long);
    
    /** clDSPOp::Sub */
    void dsp_subf (float *, float, long);
    void dsp_sub (double *, double, long);
    void dsp_csubf (stpSCplx, stSCplx, long);
    void dsp_csub (stpDCplx, stDCplx, long);
    void dsp_sub2f (float *, const float *, long);
    void dsp_sub2 (double *, const double *, long);
    void dsp_csub2f (stpSCplx, const stpSCplx, long);
    void dsp_csub2 (stpDCplx, const stpDCplx, long);
    void dsp_sub3f (float *, const float *, const float *, long);
    void dsp_sub3 (double *, const double *, const double *, long);
    void dsp_csub3f (stpSCplx, const stpSCplx, const stpSCplx, long);
    void dsp_csub3 (stpDCplx, const stpDCplx, const stpDCplx, long);
    
    /** clDSPOp::Mul */
    void dsp_mulf (float *, float, long);
    void dsp_mul (double *, double, long);
    void dsp_chmulf (stpSCplx, float, long);
    void dsp_chmul (stpDCplx, double, long);
    void dsp_cmulf (stpSCplx, stSCplx, long);
    void dsp_cmul (stpDCplx, stDCplx, long);
    void dsp_mulf_nip (float *, const float *, float, long);
    void dsp_mul_nip (double *, const double *, double, long);
    void dsp_mul2f (float *, const float *, long);
    void dsp_mul2 (double *, const double *, long);
    void dsp_chmul2f (stpSCplx, const float *, long);
    void dsp_chmul2 (stpDCplx, const double *, long);
    void dsp_cmul2f (stpSCplx, const stpSCplx, long);
    void dsp_cmul2 (stpDCplx, const stpDCplx, long);
    void dsp_mul3f (float *, const float *, const float *, long);
    void dsp_mul3 (double *, const double *, const double *, long);
    void dsp_cmul3f (stpSCplx, const stpSCplx, const stpSCplx, long);
    void dsp_cmul3 (stpDCplx, const stpDCplx, const stpDCplx, long);
    void dsp_ccmulf (stpSCplx, const stpSCplx, long);
    void dsp_ccmul (stpDCplx, const stpDCplx, long);
    void dsp_ccmulf_nip (stpSCplx, const stpSCplx, const stpSCplx, long);
    void dsp_ccmul_nip (stpDCplx, const stpDCplx, const stpDCplx, long);
    
    /** clDSPOp::Div */
    void dsp_divf (float *, float, long);
    void dsp_div (double *, double, long);
    void dsp_cdivf (stpSCplx, stSCplx, long);
    void dsp_cdiv (stpDCplx, stDCplx, long);
    void dsp_div2f (float *, const float *, long);
    void dsp_div2 (double *, const double *, long);
    void dsp_cdiv2f (stpSCplx, const stpSCplx, long);
    void dsp_cdiv2 (stpDCplx, const stpDCplx, long);
    void dsp_div3f (float *, const float *, const float *, long);
    void dsp_div3 (double *, const double *, const double *, long);
    void dsp_cdiv3f (stpSCplx, const stpSCplx, const stpSCplx, long);
    void dsp_cdiv3 (stpDCplx, const stpDCplx, const stpDCplx, long);
    
    /** clDSPOp::Div1x */
    void dsp_div1xf (float *, long);
    void dsp_div1x (double *, long);
    void dsp_div1xf_nip (float *, const float *, long);
    void dsp_div1x_nip (double *, const double *, long);

    /** clDSPOp::MulAdd */
    void dsp_muladdf (float *, float, float, long);
    void dsp_muladd (double *, double, double, long);
    void dsp_muladdf_nip (float *, const float *, float, float, long);
    void dsp_muladd_nip (double *, const double *, double, double, long);
    void dsp_cmuladdf (stpSCplx, const stpSCplx, const stpSCplx, long);
    void dsp_cmuladd (stpDCplx, const stpDCplx, const stpDCplx, long);
    
    /** clDSPOp::Abs */
    void dsp_absf (float *, long);
    void dsp_abs (double *, long);
    void dsp_absf_nip (float *, const float *, long);
    void dsp_abs_nip (double *, const double *, long);
    
    /** clDSPOp::Sqrt */
    void dsp_sqrtf (float *, long);
    void dsp_sqrt (double *, long);
    void dsp_sqrtf_nip (float *, const float *, long);
    void dsp_sqrt_nip (double *, const double *, long);
    
    /** clDSPOp::Zero */
    void dsp_zerof (float *, long);
    void dsp_zero (double *, long);
    void dsp_czerof (stpSCplx, long);
    void dsp_czero (stpDCplx, long);
    
    /** clDSPOp::Set */
    void dsp_setf (float *, float, long);
    void dsp_set (double *, double, long);
    void dsp_csetf (stpSCplx, stSCplx, long);
    void dsp_cset (stpDCplx, stDCplx, long);
    void dsp_set2f (float *, float, long, long, long);
    void dsp_set2 (double *, double, long, long, long);
    void dsp_cset2f (stpSCplx, stSCplx, long, long, long);
    void dsp_cset2 (stpDCplx, stDCplx, long, long, long);
    
    /** clDSPOp::Clip */
    void dsp_clipf (float *, float, long);
    void dsp_clip (double *, double, long);
    void dsp_clipf_nip (float *, const float *, float, long);
    void dsp_clip_nip (double *, const double *, double, long);
    void dsp_clip2f (float *, float, float, long);
    void dsp_clip2 (double *, double, double, long);
    void dsp_clip2f_nip (float *, const float *, float, float, long);
    void dsp_clip2_nip (double *, const double *, double, double, long);
    
    /** clDSPOp::ClipZero */
    void dsp_clipzerof (float *, long);
    void dsp_clipzero (double *, long);
    void dsp_clipzerof_nip (float *, const float *, long);
    void dsp_clipzero_nip (double *, const double *, long);
    
    /** clDSPOp::Copy */
    void dsp_copyf (float *, const float *, long);
    void dsp_copy (double *, const double *, long);
    /*void dsp_ccopyf (stpSCplx, const stpSCplx, long);
    void dsp_ccopy (stpDCplx, const stpDCplx, long);*/
    
    /** clDSPOp::Convolve */
    float dsp_convolvef (const float *, const float *, long);
    double dsp_convolve (const double *, const double *, long);
    void dsp_convolve2f (float *, const float *, const float *, long);
    void dsp_convolve2 (double *, const double *, const double *, long);
    
    /** clDSPOp::Correlate */
    float dsp_correlatef (const float *, const float *, long);
    double dsp_correlate (const double *, const double *, long);
    void dsp_correlate2f (float *, const float *, const float *, long);
    void dsp_correlate2 (double *, const double *, const double *, long);
    
    /** clDSPOp::AutoCorrelate */
    float dsp_autocorrf (const float *, long);
    double dsp_autocorr (const double *, long);
    void dsp_autocorr2f (float *, const float *, long);
    void dsp_autocorr2 (double *, const double *, long);
    
    /** clDSPOp::DotProduct */
    float dsp_dotproductf (const float *, const float *, long);
    double dsp_dotproduct (const double *, const double *, long);
    
    /** clDSPOp::MinMax */
    void dsp_minmaxf (float *, float *, const float *, long);
    void dsp_minmax (double *, double *, const double *, long);
    
    /** clDSPOp::Mean */
    float dsp_meanf (const float *, long);
    double dsp_mean (const double *, long);
    
    /** clDSPOp::Median */
    float dsp_meadianf (const float *, long);
    double dsp_median (const double *, long);
    
    /** clDSPOp::Negate */
    void dsp_negatef (float *, long);
    void dsp_negate (double *, long);
    void dsp_negatef_nip (float *, const float *, long);
    void dsp_negate_nip (double *, const double *, long);
    
    /** clDSPOp::Normalize */
    void dsp_normalizef (float *, long);
    void dsp_normalize (double *, long);
    void dsp_normalizef_nip (float *, const float *, long);
    void dsp_normalize_nip (double *, const double *, long);
    
    /** clDSPOp::Product */
    float dsp_productf (const float *, long);
    double dsp_product (const double *, long);
    
    /** clDSPOp::Reverse */
    void dsp_reversef (float *, long);
    void dsp_reverse (double *, long);
    void dsp_reversef_nip (float *, const float *, long);
    void dsp_reverse_nip (double *, const double *, long);
    
    /** clDSPOp::Scale */
    void dsp_scalef (float *, long);
    void dsp_scale (double *, long);
    void dsp_scalef_nip (float *, const float *, long);
    void dsp_scale_nip (double *, const double *, long);
    /** clDSPOp::Scale01 */
    void dsp_scale01f (float *, long);
    void dsp_scale01 (double *, long);
    void dsp_scale01f_nip (float *, const float *, long);
    void dsp_Scale01_nip (double *, const double *, long);
    
    /** clDSPOp::Sort */
    void dsp_sortf (float *, long);
    void dsp_sort (double *, long);
    void dsp_sortl (long *, long);
    
    /** clDSPOp::StdDev */
    void dsp_stddevf (float *, float *, const float *, long);
    void dsp_stddev (double *, double *, const double *, long);
    
    /** clDSPOp::Sum */
    float dsp_sumf (const float *, long);
    double dsp_sum (const double *, long);
    
    /** clDSPOp::Square */
    void dsp_squaref (float *, long);
    void dsp_square (double *, long);
    void dsp_squaref_nip (float *, const float *, long);
    void dsp_square_nip (double *, const double *, long);
    
    /** clDSPOp::Convert */
    void dsp_convertu8f (float *, const unsigned char *, long);
    void dsp_convertu8 (double *, const unsigned char *, long);
    void dsp_converts16f (float *, const signed short *, long, int);
    void dsp_converts16 (double *, const signed short *, long, int);
    void dsp_converts32f (float *, const signed int *, long, int);
    void dsp_converts32 (double *, const signed int *, long, int);
    void dsp_convertd64f (float *, const double *, long);
    void dsp_convertf32 (double *, const float *, long);
    
    void dsp_convertf32c (unsigned char *, const float *, long);
    void dsp_convertd64c (unsigned char *, const double *, long);
    void dsp_convertf32s (signed short *, const float *, long, int);
    void dsp_convertd64s (signed short *, const double *, long, int);
    void dsp_convertf32i (signed int *, const float *, long, int);
    void dsp_convertd64i (signed int *, const double *, long, int);
    
    /** clDSPOp::CartToPolar */
    void dsp_cart2polarf (float *, float *, const float *, const float *,
        long);
    void dsp_cart2polar (double *, double *, const double *, const double *,
        long);
    void dsp_cart2polar2f (float *, float *, const stpSCplx, long);
    void dsp_cart2polar2 (double *, double *, const stpDCplx, long);
    void dsp_cart2polar3f (stpSPolar, const stpSCplx, long);
    void dsp_cart2polar3 (stpDPolar, const stpDCplx, long);
    void dsp_cart2polar4f (utpSCoord, long);
    void dsp_cart2polar4 (utpDCoord, long);
    
    /** clDSPOp::PolarToCart */
    void dsp_polar2cartf (float *, float *, const float *, const float *,
        long);
    void dsp_polar2cart (double *, double *, const double *, const double *,
        long);
    void dsp_polar2cart2f (stpSCplx, const float *, const float *, long);
    void dsp_polar2cart2 (stpDCplx, const double *, const double *, long);
    void dsp_polar2cart3f (stpSCplx, const stpSPolar, long);
    void dsp_polar2cart3 (stpDCplx, const stpDPolar, long);
    void dsp_polar2cart4f (utpSCoord, long);
    void dsp_polar2cart4 (utpDCoord, long);
    
    /** clDSPOp::CrossCorr */
    float dsp_crosscorrf (const float *, const float *, long);
    double dsp_crosscorr (const double *, const double *, long);
    /** clDSPOp::DelCrossCorr */
    float dsp_crosscorr2f (const float *, const float *, long, long);
    double dsp_crosscorr2 (const double *, const double *, long, long);
    void dsp_crosscorr3f (float *, const float *, const float *, long,
        const long *, long);
    void dsp_crosscorr3 (double *, const double *, const double *, long,
        const long *, long);
    
    /** clDSPOp::Energy */
    float dsp_energyf (const float *, long);
    double dsp_energy (const double *, long);
    
    /** clDSPOp::Magnitude */
    void dsp_magnitudef (float *, const stpSCplx, long);
    void dsp_magnitude (double *, const stpDCplx, long);
    
    /** clDSPOp::Power */
    void dsp_powerf (float *, const stpSCplx, long);
    void dsp_power (double *, const stpDCplx, long);
    
    /** clDSPOp::Phase */
    void dsp_phasef (float *, const stpSCplx, long);
    void dsp_phase (double *, const stpDCplx, long);
    
    /** clDSPOp::PowerPhase */
    void dsp_powerphasef (float *, float *, const stpSCplx, long);
    void dsp_powerphase (double *, double *, const stpDCplx, long);
    
    /** clDSPOp::Decimate */
    void dsp_decimatef (float *, const float *, long, long);
    void dsp_decimate (double *, const double *, long, long);
    /** clDSPOp::DecimateAvg */
    void dsp_decimateavgf (float *, const float *, long, long);
    void dsp_decimateavg (double *, const double *, long, long);
    
    /** clDSPOp::Interpolate */
    void dsp_interpolatef (float *, const float *, long, long);
    void dsp_interpolate (double *, const double *, long, long);
    /** clDSPOp::InterpolateAvg */
    void dsp_interpolateavgf (float *, const float *, long, long);
    void dsp_interpolateavg (double *, const double *, long, long);

    /** clDSPOp::Resample */
    void dsp_resamplef (float *, long, const float *, long);
    void dsp_resample (double *, long, const double *, long);
    /** clDSPOp::ResampleAvg */
    void dsp_resampleavgf (float *, long, const float *, long);
    void dsp_resampleavg (double *, long, const double *, long);
    
    /** clDSPOp::RMS */
    float dsp_rmsf (const float *, long);
    double dsp_rms (const double *, long);
    
    /** clDSPOp::Variance */
    float dsp_variancef (float *, float *, const float *, long);
    double dsp_variance (double *, double *, const double *, long);
    
    /** clDSPOp::PeakLevel */
    float dsp_peaklevelf (const float *fpSrc, long);
    double dsp_peaklevel (const double *dpSrc, long);
    
    /** clDSPOp::Mix */
    void dsp_mixf (float *, const float *, long);
    void dsp_mix (double *, const double *, long);
    void dsp_mix2f (float *, const float *, const float *, long);
    void dsp_mix2 (double *, const double *, const double *, long);
    void dsp_mixnf (float *, const float *, long, long);
    void dsp_mixn (double *, const double *, long, long);
    
    /** clDSPOp::Extract */
    void dsp_extractf (float *, const float *, long, long, long);
    void dsp_extract (double *, const double *, long, long, long);
    
    /** clDSPOp::Pack */
    void dsp_packf (float *, const float *, long, long, long);
    void dsp_pack (double *, const double *, long, long, long);
    
    /** clDSPOp::FFTWConvert */
    void dsp_fftw_convertf2cf (stpSCplx, const float *, long);
    void dsp_fftw_convertf2cd (stpDCplx, const float *, long);
    void dsp_fftw_convertd2cf (stpSCplx, const double *, long);
    void dsp_fftw_convertd2cd (stpDCplx, const double *, long);
    void dsp_fftw_convertcf2f (float *, const stpSCplx, long);
    void dsp_fftw_convertcd2f (float *, const stpDCplx, long);
    void dsp_fftw_convertcf2d (double *, const stpSCplx, long);
    void dsp_fftw_convertcd2d (double *, const stpDCplx, long);
    
    /**
        Creates new instance of DSP object.
        
        \return DSP object instance
    */
    dsp_t dsp_new (void);
    /**
        Deletes DSP object instance.
        
        \param dspInst DSP object instance
    */
    void dsp_delete (dsp_t);

    /** clDSPOp::WinBartlett */
    void dsp_win_bartlettf (dsp_t, float *, long);
    void dsp_win_bartlett (dsp_t, double *, long);
    /** clDSPOp::WinBlackman */
    void dsp_win_blackmanf (dsp_t, float *, long, float);
    void dsp_win_blackman (dsp_t, double *, long, double);
    /** clDSPOp::WinBlackmanHarris */
    void dsp_win_blackman_harrisf (dsp_t, float *, long);
    void dsp_win_blackman_harris (dsp_t, double *, long);
    /** clDSPOp::WinCosTapered */
    void dsp_win_cos_taperedf (dsp_t, float *, long);
    void dsp_win_cos_tapered (dsp_t, double *, long);
    /** clDSPOp::WinExactBlackman */
    void dsp_win_exact_blackmanf (dsp_t, float *, long);
    void dsp_win_exact_blackman (dsp_t, double *, long);
    /** clDSPOp::WinExp */
    void dsp_win_expf (dsp_t, float *, float, long);
    void dsp_win_exp (dsp_t, double *, double, long);
    /** clDSPOp::WinFlatTop */
    void dsp_win_flat_topf (dsp_t, float *, long);
    void dsp_win_flat_top (dsp_t, double *, long);
    /** clDSPOp::WinGenericCos */
    void dsp_win_generic_cosf (dsp_t, float *, long, const float *, long);
    void dsp_win_generic_cos (dsp_t, double *, long, const double *, long);
    /** clDSPOp::WinHamming */
    void dsp_win_hammingf (dsp_t, float *, long);
    void dsp_win_hamming (dsp_t, double *, long);
    /** clDSPOp::WinHanning */
    void dsp_win_hanningf (dsp_t, float *, long);
    void dsp_win_hanning (dsp_t, double *, long);
    /** clDSPOp::WinKaiser */
    void dsp_win_kaiserf (dsp_t, float *, float, long);
    void dsp_win_kaiser (dsp_t, double *, double, long);
    /** clDSPOp::WinKaiserBessel */
    void dsp_win_kaiser_besself (dsp_t, float *, float, long);
    void dsp_win_kaiser_bessel (dsp_t, double *, double, long);
    /** clDSPOp::WinTukey */
    void dsp_win_tukeyf (dsp_t, float *, long);
    void dsp_win_tukey (dsp_t, double *, long);
    /** clDSPOp::WinDolphChebyshev */
    void dsp_win_dolph_chebyshevf (dsp_t, float *, float, long);
    void dsp_win_dolph_chebyshev (dsp_t, double *, double, long);
    
    /** clDSPOp::ReBuffer */
    long dsp_rebufferf (dsp_t, float *, const float *, long, long);
    long dsp_rebuffer (dsp_t, double *, const double *, long, long);

    /** clDSPOp::DegToRad */
    float dsp_deg2radf (dsp_t, float);
    double dsp_deg2rad (dsp_t, double);
    /** clDSPOp::RadToDeg */
    float dsp_rad2degf (dsp_t, float);
    double dsp_rad2deg (dsp_t, double);
    
    /** clDSPOp::FIRAllocate */
    void dsp_fir_allocatef (dsp_t, const float *, long);
    void dsp_fir_allocate (dsp_t, const double *, long);
    /** clDSPOp::FIRFree */
    void dsp_fir_free (dsp_t);
    /** clDSPOp::FIRFilter */
    void dsp_fir_filterf (dsp_t, float *, long);
    void dsp_fir_filter (dsp_t, double *, long);
    void dsp_fir_filterf_nip (dsp_t, float *, const float *, long);
    void dsp_fir_filter_nip (dsp_t, double *, const double *, long);
    /** clDSPOp::FIRFilterF */
    void dsp_fir_filterf_fst (dsp_t, float *, float *, long);
    void dsp_fir_filter_fst (dsp_t, double *, double *, long);
    
    /** clDSPOp::IIRInitialize */
    void dsp_iir_initf (dsp_t, const float *);
    void dsp_iir_init (dsp_t, const double *);
    /** clDSPOp::IIRFilter */
    void dsp_iir_filterf (dsp_t, float *, long);
    void dsp_iir_filter (dsp_t, double *, long);
    void dsp_iir_filterf_nip (dsp_t, float *, const float *, long);
    void dsp_iir_filter_nip (dsp_t, double *, const double *, long);
    /** clDSPOp::IIRClear */
    void dsp_iir_clear (dsp_t);
    
    /** clDSPOp::FFTInitialize */
    void dsp_fft_init (dsp_t, long, int);
    /** clDSPOp::FFTUninitialize */
    void dsp_fft_uninit (dsp_t);
    /** clDSPOp::FFTi */
    void dsp_fftf (dsp_t, stpSCplx, float *);
    void dsp_fft (dsp_t, stpDCplx, double *);
    /** clDSPOp::FFTo */
    void dsp_fftf_nip (dsp_t, stpSCplx, const float *);
    void dsp_fft_nip (dsp_t, stpDCplx, const double *);
    void dsp_cfftf_nip (dsp_t, stpSCplx, const stpSCplx);
    void dsp_cfft_nip (dsp_t, stpDCplx, const stpDCplx);
    /** clDSPOp::IFFTo */
    void dsp_ifftf_nip (dsp_t, float *, const stpSCplx);
    void dsp_ifft_nip (dsp_t, double *, const stpDCplx);
    void dsp_cifftf_nip (dsp_t, stpSCplx, const stpSCplx);
    void dsp_cifft_nip (dsp_t, stpDCplx, const stpDCplx);

    /** clDSPOp::DCTInitialize */
    void dsp_dct_init (dsp_t, long);
    /** clDSPOp::DCTUninitialize */
    void dsp_dct_uninit (dsp_t);
    /** clDSPOp::DCTi */
    void dsp_dctf (dsp_t, float *);
    void dsp_dct (dsp_t, double *);
    /** clDSPOp::DCTo */
    void dsp_dctf_nip (dsp_t, float *, const float *);
    void dsp_dct_nip (dsp_t, double *, const double *);
    /** clDSPOp::IDCTi */
    void dsp_idctf (dsp_t, float *);
    void dsp_idct (dsp_t, double *);
    /** clDSPOp::IDCTo */
    void dsp_idctf_nip (dsp_t, float *, const float *);
    void dsp_idct_nip (dsp_t, double *, const double *);


    /* clIIRCascade */

    /**
        Creates new instance of cascaded IIR object.
        
        \return Cascade IIR instance
    */
    dsp_iircas_t dsp_iir_cas_new (void);
    /**
        Deletes cascaded IIR object instance.
        
        \param dspiircasInst Cascaded IIR instance
    */
    void dsp_iir_cas_delete (dsp_iircas_t);
    /** clIIRCascade::Initialize */
    int dsp_iir_cas_initf (dsp_iircas_t, const float[][5], long);
    int dsp_iir_cas_init (dsp_iircas_t, const double[][5], long);
    void dsp_iir_cas_uninit (dsp_iircas_t);
    void dsp_iir_cas_processf (dsp_iircas_t, float *, long);
    void dsp_iir_cas_process (dsp_iircas_t, double *, long);
    void dsp_iir_cas_processf_nip (dsp_iircas_t, float *, const float *, long);
    void dsp_iir_cas_process_nip (dsp_iircas_t, double *, const double *, long);
    void dsp_iir_cas_clear (dsp_iircas_t);

    /* clFFTDecimator */

    /**
        Creates new instance of FFT decimator object.
        
        \return FFT decimator instance
    */
    dsp_decfft_t dsp_dec_fft_new (void);
    /**
        Deletes FFT decimator object instance.
        
        \param dspdecfftInst FFT decimator instance
    */
    void dsp_dec_fft_delete (dsp_decfft_t);
    /** clFFTDecimator::Initialize */
    int dsp_dec_fft_initf (dsp_decfft_t, long, long, int);
    int dsp_dec_fft_init (dsp_decfft_t, long, long, int);
    /** clFFTDecimator::Uninitialize */
    void dsp_dec_fft_uninit (dsp_decfft_t);
    /** clFFTDecimator::Put */
    void dsp_dec_fft_putf (dsp_decfft_t, const float *, long);
    void dsp_dec_fft_put (dsp_decfft_t, const double *, long);
    /** clFFTDecimator::Get */
    int dsp_dec_fft_getf (dsp_decfft_t, float *, long);
    int dsp_dec_fft_get (dsp_decfft_t, double *, long);

    /* clFIRDecimator */

    /**
        Creates new instance of FIR decimator object.
        
        \return FIR decimator instance
    */
    dsp_decfir_t dsp_dec_fir_new (void);
    /**
        Deletes FIR decimator object instance.
        
        \param dspdecfirInst FIR decimator instance
    */
    void dsp_dec_fir_delete (dsp_decfir_t);
    /** clFIRDecimator::Initialize */
    int dsp_dec_fir_initf (dsp_decfir_t, long, int);
    int dsp_dec_fir_init (dsp_decfir_t, long, int);
    /** clFIRDecimator::Uninitialize */
    void dsp_dec_fir_uninit (dsp_decfir_t);
    /** clFIRDecimator::Put */
    void dsp_dec_fir_putf (dsp_decfir_t, const float *, long);
    void dsp_dec_fir_put (dsp_decfir_t, const double *, long);
    /** clFIRDecimator::Get */
    int dsp_dec_fir_getf (dsp_decfir_t, float *, long);
    int dsp_dec_fir_get (dsp_decfir_t, double *, long);

    /* clIIRDecimator */

    /**
        Creates new instance of IIR decimator object.
        
        \return IIR decimator instance
    */
    dsp_deciir_t dsp_dec_iir_new (void);
    /**
        Deletes IIR decimator object instance.
        
        \param dspdeciirInst IIR decimator instance
    */
    void dsp_dec_iir_delete (dsp_deciir_t);
    /** clIIRDecimator::Initialize */
    int dsp_dec_iir_initf (dsp_deciir_t, long, int);
    int dsp_dec_iir_init (dsp_deciir_t, long, int);
    /** clIIRDecimator::Uninitialize */
    void dsp_dec_iir_uninit (dsp_deciir_t);
    /** clIIRDecimator::Put */
    void dsp_dec_iir_putf (dsp_deciir_t, const float *, long);
    void dsp_dec_iir_put (dsp_deciir_t, const double *, long);
    /** clIIRDecimator::Get */
    int dsp_dec_iir_getf (dsp_deciir_t, float *, long);
    int dsp_dec_iir_get (dsp_deciir_t, double *, long);

    /* clRecDecimator */

    /**
        Creates new instance of recursive decimator object.
        
        \return Recursive decimator instance
    */
    dsp_decrec_t dsp_dec_rec_new (void);
    /**
        Deletes recursive decimator object instance.
        
        \param dspdecrecInst Recursive decimator instance
    */
    void dsp_dec_rec_delete (dsp_decrec_t);
    /** clRecDecimator::Initialize */
    int dsp_dec_rec_initf (dsp_decrec_t, long, long, float, int);
    int dsp_dec_rec_init (dsp_decrec_t, long, long, double, int);
    /** clRecDecimator::Uninitialize */
    void dsp_dec_rec_uninit (dsp_decrec_t);
    /** clRecDecimator::Put */
    void dsp_dec_rec_putf (dsp_decrec_t, const float *, long);
    void dsp_dec_rec_put (dsp_decrec_t, const double *, long);
    /** clRecDecimator::Get */
    int dsp_dec_rec_getf (dsp_decrec_t, float *, long);
    int dsp_dec_rec_get (dsp_decrec_t, double *, long);

    /* clFFTInterpolator */

    /**
        Creates new instance of FFT interpolator object.
        
        \return FFT interpolator instance
    */
    dsp_intfft_t dsp_int_fft_new (void);
    /**
        Deletes FFT interpolator object instance.
        
        \param dspintfftInst FFT interpolator instance
    */
    void dsp_int_fft_delete (dsp_intfft_t);
    /** clFFTInterpolator::Initialize */
    int dsp_int_fft_initf (dsp_intfft_t, long, long, int);
    int dsp_int_fft_init (dsp_intfft_t, long, long, int);
    /** clFFTInterpolator::Uninitialize */
    void dsp_int_fft_uninit (dsp_intfft_t);
    /** clFFTInterpolator::Put */
    void dsp_int_fft_putf (dsp_intfft_t, const float *, long);
    void dsp_int_fft_put (dsp_intfft_t, const double *, long);
    /** clFFTInterpolator::Get */
    int dsp_int_fft_getf (dsp_intfft_t, float *, long);
    int dsp_int_fft_get (dsp_intfft_t, double *, long);

    /* clFIRInterpolator */

    /**
        Creates new instance of FIR interpolator object.
        
        \return FIR interpolator instance
    */
    dsp_intfir_t dsp_int_fir_new (void);
    /**
        Deletes FIR interpolator object instance.
        
        \param dspintfirInst FIR interpolator instance
    */
    void dsp_int_fir_delete (dsp_intfir_t);
    /** clFIRInterpolator::Initialize */
    int dsp_int_fir_initf (dsp_intfir_t, long, int);
    int dsp_int_fir_init (dsp_intfir_t, long, int);
    /** clFIRInterpolator::Uninitialize */
    void dsp_int_fir_uninit (dsp_intfir_t);
    /** clFIRInterpolator::Put */
    void dsp_int_fir_putf (dsp_intfir_t, const float *, long);
    void dsp_int_fir_put (dsp_intfir_t, const double *, long);
    /** clFIRInterpolator::Get */
    int dsp_int_fir_getf (dsp_intfir_t, float *, long);
    int dsp_int_fir_get (dsp_intfir_t, double *, long);

    /* clIIRInterpolator */

    /**
        Creates new instance of IIR interpolator object.
        
        \return IIR interpolator instance
    */
    dsp_intiir_t dsp_int_iir_new (void);
    /**
        Deletes IIR interpolator object instance.
        
        \param dspintiirInst IIR interpolator instance
    */
    void dsp_int_iir_delete (dsp_intiir_t);
    /** clIIRInterpolator::Initialize */
    int dsp_int_iir_initf (dsp_intiir_t, long, int);
    int dsp_int_iir_init (dsp_intiir_t, long, int);
    /** clIIRInterpolator::Uninitialize */
    void dsp_int_iir_uninit (dsp_intiir_t);
    /** clIIRInterpolator::Put */
    void dsp_int_iir_putf (dsp_intiir_t, const float *, long);
    void dsp_int_iir_put (dsp_intiir_t, const double *, long);
    /** clIIRInterpolator::Get */
    int dsp_int_iir_getf (dsp_intiir_t, float *, long);
    int dsp_int_iir_get (dsp_intiir_t, double *, long);

    /* clRecInterpolator */

    /**
        Creates new instance of recursive interpolator object.
        
        \return Recursive interpolator instance
    */
    dsp_intrec_t dsp_int_rec_new (void);
    /**
        Deletes recursive interpolator object instance.
        
        \param dspdecrecInst Recursive interpolator instance
    */
    void dsp_int_rec_delete (dsp_intrec_t);
    /** clRecInterpolator::Initialize */
    int dsp_int_rec_initf (dsp_intrec_t, long, long, float, int);
    int dsp_int_rec_init (dsp_intrec_t, long, long, double, int);
    /** clRecInterpolator::Uninitialize */
    void dsp_int_rec_uninit (dsp_intrec_t);
    /** clRecInterpolator::Put */
    void dsp_int_rec_putf (dsp_intrec_t, const float *, long);
    void dsp_int_rec_put (dsp_intrec_t, const double *, long);
    /** clRecInterpolator::Get */
    int dsp_int_rec_getf (dsp_intrec_t, float *, long);
    int dsp_int_rec_get (dsp_intrec_t, double *, long);

    /* clFilter */

    /**
        Creates new instance of FFT filter object.
        
        \return FFT filter instance
    */
    dsp_filter_t dsp_filter_new (void);
    /**
        Deletes instance of FFT filter object.
        
        \param dspfilterInst FFT filter instance
    */
    void dsp_filter_delete (dsp_filter_t);
    /** clFilter::Initialize */
    int dsp_filter_initf (dsp_filter_t, long);
    int dsp_filter_init (dsp_filter_t, long);
    int dsp_filter_init2f (dsp_filter_t, long, const float *, float, float);
    int dsp_filter_init2 (dsp_filter_t, long, const double *, double, double);
    /** clFilter::InitializeLP */
    int dsp_filter_init_lpf (dsp_filter_t, float, float, float, float);
    int dsp_filter_init_lp (dsp_filter_t, double, double, double, double);
    /** clFilter::InitializeHP */
    int dsp_filter_init_hpf (dsp_filter_t, float, float, float, float);
    int dsp_filter_init_hp (dsp_filter_t, double, double, double, double);
    /** clFilter::SetCoeffs */
    void dsp_filter_set_coeffsf (dsp_filter_t, const float *);
    void dsp_filter_set_ceoffs (dsp_filter_t, const double *);
    void dsp_filter_set_ccoeffsf (dsp_filter_t, const stpSCplx, int);
    void dsp_filter_set_ccoeffs (dsp_filter_t, const stpDCplx, int);
    /** clFilter::GetCoeffs */
    void dsp_filter_get_coeffsf (dsp_filter_t, float *);
    void dsp_filter_get_coeffs (dsp_filter_t, double *);
    void dsp_filter_get_ccoeffsf (dsp_filter_t, stpSCplx);
    void dsp_filter_get_ccoeffs (dsp_filter_t, stpDCplx);
    /** clFilter::Put */
    void dsp_filter_putf (dsp_filter_t, const float *, long);
    void dsp_filter_put (dsp_filter_t, const double *, long);
    void dsp_filter_put2f (dsp_filter_t, const float *, long, const stpSCplx);
    void dsp_filter_put2 (dsp_filter_t, const double *, long, const stpDCplx);
    /** clFilter::Get */
    void dsp_filter_getf (dsp_filter_t, float *, long);
    void dsp_filter_get (dsp_filter_t, double *, long);
    /** clFilter::DesignLP */
    void dsp_filter_design_lpf (dsp_filter_t, float *, int);
    void dsp_filter_design_lp (dsp_filter_t, double *, int);
    void dsp_filter_design_lp2f (dsp_filter_t, float *, float, int);
    void dsp_filter_design_lp2 (dsp_filter_t, double *, double, int);
    /** clFilter::DesignHP */
    void dsp_filter_design_hpf (dsp_filter_t, float *);
    void dsp_filter_design_hp (dsp_filter_t, double *);
    void dsp_filter_design_hp2f (dsp_filter_t, float *, float);
    void dsp_filter_design_hp2 (dsp_filter_t, double *, double);
    /** clFilter::DesignBP */
    void dsp_filter_design_bpf (dsp_filter_t, float *, float *);
    void dsp_filter_design_bp (dsp_filter_t, double *, double *);
    void dsp_filter_design_bp2f (dsp_filter_t, float *, float *, float);
    void dsp_filter_design_bp2 (dsp_filter_t, double *, double *, double);
    /** clFilter::DesignBR */
    void dsp_filter_design_brf (dsp_filter_t, float *, float *);
    void dsp_filter_design_br (dsp_filter_t, double *, double *);
    void dsp_filter_design_br2f (dsp_filter_t, float *, float *, float);
    void dsp_filter_design_br2 (dsp_filter_t, double *, double *, double);


    /* clFlipBand */
    
    /**
        Creates new instance of band flipping object.
        
        \return FlipBand instance
    */
    dsp_flip_t dsp_flip_new (void);
    /**
        Deletes instance of band flipping object.
        
        \param dspflipInst FlipBand instance
    */
    void dsp_flip_delete (dsp_flip_t);
    /** clFlipBand::Initialize */
    void dsp_flip_initf (dsp_flip_t, long);
    void dsp_flip_init (dsp_flip_t, long);
    /** clFlipBand::Uninitialize */
    void dsp_flip_uninit (dsp_flip_t);
    /** clFlipBand::Put */
    void dsp_flip_putf (dsp_flip_t, const float *, long);
    void dsp_flip_put (dsp_flip_t, const double *, long);
    /** clFlipBand::Get */
    void dsp_flip_getf (dsp_flip_t, float *, long);
    void dsp_flip_get (dsp_flip_t, double *, long);
    /** clFlipBand::Clear */
    void dsp_flip_clear (dsp_flip_t);


    /* clReBuffer */

    /**
        Creates new instance of rebuffering object.
        
        \return Rebuffer instance
    */
    dsp_rebuf_t dsp_rebuf_new (void);
    /**
        Deletes instance of rebuffering object.
        
        \param dsprebufInst Rebuffer instance
    */
    void dsp_rebuf_delete (dsp_rebuf_t);
    /** clReBuffer::Put */
    void dsp_rebuf_putf (dsp_rebuf_t, const float *, long);
    void dsp_rebuf_put (dsp_rebuf_t, const double *, long);
    /** clReBuffer::Get */
    int dsp_rebuf_getf (dsp_rebuf_t, float *, long);
    int dsp_rebuf_get (dsp_rebuf_t, double *, long);
    /** clReBuffer::GetCount */
    long dsp_rebuf_size (dsp_rebuf_t);
    /** clReBuffer::Clear */
    void dsp_rebuf_clear (dsp_rebuf_t);
    /** clReBuffer::operator= */
    void dsp_rebuf_copy (dsp_rebuf_t, dsp_rebuf_t);
    
    #ifdef __cplusplus
    }
    #endif

#endif
