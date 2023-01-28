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

#ifndef X86_H
    #define X86_H

    #if ((__GNUC__ == 3) && (__GNUC_MINOR__ >= 1))
        #define likely(x)           __builtin_expect(x, 1)
        #define unlikely(x)         __builtin_expect(x, 0)
        #define prefetch(x, w, l)   __builtin_prefetch(x, w, l)
    #else
        #define likely(x)           (x)
        #define unlikely(x)         (x)
        #define prefetch(x, w, l)
    #endif

    #define X86_ASM     asm volatile

    typedef float v2sf __attribute__((vector_size(8)));
    typedef float v4sf __attribute__((vector_size(16)));
    typedef double v2df __attribute__((vector_size(16)));

    typedef v2sf * pv2sf;
    typedef v4sf * pv4sf;
    typedef v2df * pv2df;

    #pragma pack(1)

    typedef struct _stm64
    {
        float f[2];
    } stm64, *stpm64;

    typedef struct _stm128
    {
        union
        {
            float f[4];
            double d[2];
        };
    } stm128, *stpm128;

    #pragma pack()


    #ifdef __cplusplus
    extern "C"
    {
    #endif

    #ifndef DSP_X86_64
    const char *dsp_x86_cpuid (void);
    void dsp_x86_features (unsigned int *, unsigned int *);
    unsigned int dsp_x86_amd_features (void);
    #endif
    int dsp_x86_have_e3dnow (void);
    int dsp_x86_have_sse (void);
    void dsp_x86_3dnow_copyf (float *, const float *, int);
    void dsp_x86_3dnow_copyd (double *, const double *, int);
    void dsp_x86_3dnow_addf (float *, float, int);
    void dsp_x86_sse_addf (float *, float, int);
    void dsp_x86_sse_add (double *, double, int);
    void dsp_x86_3dnow_mulf (float *, float, int);
    void dsp_x86_sse_mulf (float *, float, int);
    void dsp_x86_sse_mul (double *, double, int);
    void dsp_x86_3dnow_mulf_nip (float *, const float *, float, int);
    void dsp_x86_sse_mulf_nip (float *, const float *, float, int);
    void dsp_x86_sse_mul_nip (double *, const double *, double, int);
    void dsp_x86_3dnow_add2f (float *, const float *, int);
    void dsp_x86_sse_add2f (float *, const float *, int);
    void dsp_x86_sse_add2 (double *, const double *, int);
    void dsp_x86_3dnow_mul2f (float *, const float *, int);
    void dsp_x86_sse_mul2f (float *, const float *, int);
    void dsp_x86_sse_mul2 (double *, const double *, int);
    void dsp_x86_3dnow_add3f (float *, const float *, const float *, int);
    void dsp_x86_sse_add3f (float *, const float *, const float *, int);
    void dsp_x86_sse_add3 (double *, const double *, const double *, int);
    void dsp_x86_3dnow_mul3f (float *, const float *, const float *, int);
    void dsp_x86_sse_mul3f (float *, const float *, const float *, int);
    void dsp_x86_sse_mul3 (double *, const double *, const double *, int);
    void dsp_x86_3dnow_cmulf (float *, const float *, int);
    void dsp_x86_sse_cmulf (float *, const float *, int);
    void dsp_x86_sse3_cmulf (float *, const float *, int);
    void dsp_x86_sse_cmul (double *, const double *, int);
    void dsp_x86_sse3_cmul (double *, const double *, int);
    void dsp_x86_3dnow_cmul2f (float *, const float *, int);
    void dsp_x86_sse_cmul2f (float *, const float *, int);
    void dsp_x86_sse3_cmul2f (float *, const float *, int);
    void dsp_x86_sse_cmul2 (double *, const double *, int);
    void dsp_x86_sse3_cmul2 (double *, const double *, int);
    void dsp_x86_3dnow_cmul3f (float *, const float *, const float *, int);
    void dsp_x86_sse_cmul3f (float *, const float *, const float *, int);
    void dsp_x86_sse3_cmul3f (float *, const float *, const float *, int);
    void dsp_x86_sse_cmul3 (double *, const double *, const double *, int);
    void dsp_x86_sse3_cmul3 (double *, const double *, const double *, int);
    void dsp_x86_3dnow_maf (float *, float, float, int);
    void dsp_x86_sse_maf (float *, float, float, int);
    void dsp_x86_sse_ma (double *, double, double, int);
    void dsp_x86_3dnow_ma2f (float *, const float *, float, float, int);
    void dsp_x86_sse_ma2f (float *, const float *, float, float, int);
    void dsp_x86_sse_ma2 (double *, const double *, double, double, int);
    void dsp_x86_3dnow_cmaf (float *, const float *, const float *, int);
    void dsp_x86_sse_cmaf (float *, const float *, const float *, int);
    void dsp_x86_sse3_cmaf (float *, const float *, const float *, int);
    void dsp_x86_sse_cma (double *, const double *, const double *, int);
    void dsp_x86_sse3_cma (double *, const double *, const double *, int);
    void dsp_x86_3dnow_amf (float *, float, float, int);
    void dsp_x86_sse_amf (float *, float, float, int);
    void dsp_x86_sse_am (double *, double, double, int);
    float dsp_x86_3dnow_macf (const float *, const float *, int);
    float dsp_x86_sse_macf (const float *, const float *, int);
    double dsp_x86_sse_mac (const double *, const double *, int);
    void dsp_x86_3dnow_minmaxf (float *, float *, const float *, int);
    void dsp_x86_sse_minmaxf (float *, float *, const float *, int);
    void dsp_x86_sse_minmax (double *, double *, const double *, int);
    float dsp_x86_3dnow_crosscorrf (const float *, const float *, int);
    float dsp_x86_sse_crosscorrf (const float *, const float *, int);
    double dsp_x86_sse_crosscorr (const double *, const double *, int);
    void dsp_x86_3dnow_i16tof (float *, const short *, int, int);
    void dsp_x86_3dnow_i32tof (float *, const int *, int, int);
    void dsp_x86_3dnow_ftoi16 (short *, const float *, int);
    /* not in-place capable */
    void dsp_x86_3dnow_firf (float *, const float *, int, const float *, int);
    void dsp_x86_sse_firf (float *, const float *, int, const float *, int);
    void dsp_x86_sse_fir (double *, const double *, int, const double *, int);
    void dsp_x86_3dnow_iirf (float *, int, const float *, float *, float *);
    void dsp_x86_sse_iirf (float *, int, const float *, float *, float *);
    void dsp_x86_sse_iir (double *, int, const double *, double *, double *);
    void dsp_x86_3dnow_iirf_nip (float *, const float *, int, const float *, 
        float *, float *);
    void dsp_x86_sse_iirf_nip (float *, const float *, int, const float *, 
        float *, float *);
    void dsp_x86_sse_iir_nip (double *, const double *, int, const double *, 
        double *, double *);

    #ifdef __cplusplus
    }
    #endif

#endif

#endif
