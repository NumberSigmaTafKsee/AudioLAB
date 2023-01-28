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

#ifndef X86_64_H
    #define X86_64_H

    #define X86_64_ASM      asm volatile

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

    void dsp_x86_64_minmaxf (float *, float *, const float *, long);
    void dsp_x86_64_minmax (double *, double *, const double *, long);
    float dsp_x86_64_crosscorrf (const float *, const float *, long);
    double dsp_x86_64_crosscorr (const double *, const double *, long);
    /* not in-place capable */
    void dsp_x86_64_firf (float *, const float *, long, const float *, long);
    void dsp_x86_64_fir (double *, const double *, long, const double *, long);
    void dsp_x86_64_iirf (float *, long, const float *, float *, float *);
    void dsp_x86_64_iir (double *, long, const double *, double *, double *);
    void dsp_x86_64_iirf_nip (float *, const float *, long, const float *, 
        float *, float *);
    void dsp_x86_64_iir_nip (double *, const double *, long, const double *, 
        double *, double *);

    #ifdef __cplusplus
    }
    #endif

#endif

#endif
