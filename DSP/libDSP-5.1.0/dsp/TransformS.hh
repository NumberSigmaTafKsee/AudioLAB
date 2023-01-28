/*

    DFT, DCT and DST routines
    Copyright (C) 1996-2001 Takuya OOURA
    Copyright (C) 1999-2005 Jussi Laako

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

    Version 1.0 : Converted T.OOURA's functions to C++ class
    Version 1.1 : Updated to match T.OOURA's Oct/99 code and all versions
    Version 1.2 : Some conversion and precision updates
    Version 1.3 : Updated to match T.OOURA's Oct/01 code

*/


#include <dsp/DSPConfig.hh>


/* This is decimation-in-frequency split-radix transform */


#ifndef TRANSFORMS_HH
    #define TRANSFORMS_HH

    #ifdef TRANSFORM_INLINE
        #define TS_INLINE           inline
    #else
        #define TS_INLINE
    #endif

    #if (defined(_GNU_SOURCE) && (!defined(ALPHA)))
        #ifndef TRANSFORM_EXT_PREC
            #define TRANSFORM_EXT_PREC
        #endif
    #endif

    #if (defined(__BORLANDC__) || defined(SOLARIS) || defined(OSX))
        #define sinf(x)             sin(x)
        #define cosf(x)             cos(x)
        #define atanf(x)            atan(x)
    #endif

    // #define TRANSFORM_EXT_PREC for extended precision for table calculation
    // for single precision transforms

    // #define USE_CDFT_THREADS to enable parallelization
    
    #ifdef USE_CDFT_THREADS
        #include <pthread.h>

        #ifndef CDFT_THREADS_BEGIN_N
            #define CDFT_THREADS_BEGIN_N 8192
        #endif
        #ifndef CDFT_4THREADS_BEGIN_N
            #define CDFT_4THREADS_BEGIN_N 65536
        #endif

        #define cdft_thread_t pthread_t
        
        typedef struct cdft_arg_st {
            long n0;
            long n;
            void *a;  // float * / double *
            long nw;
            void *w;  // float * / double *
            void *klass;
        } cdft_arg_t;
    #endif

    
    /**
        Decimation-in-frequency split-radix transform

        \author Takuya OOURA
        \author Jussi Laako

        \note See clTransform4 for details.
    */
    class clTransformS
    {
            void makeipt(long, long *);
            void makewt(long, long *, float *);
            void makewt(long, long *, double *);
            void makect(long, long *, float *);
            void makect(long, long *, double *);
            void cftfsub(long, float *, long *, long, float *);
            void cftfsub(long, double *, long *, long, double *);
            void cftbsub(long, float *, long *, long, float *);
            void cftbsub(long, double *, long *, long, double *);
            void bitrv2(long, long *, float *);
            void bitrv2(long, long *, double *);
            void bitrv2conj(long, long *, float *);
            void bitrv2conj(long, long *, double *);
            void bitrv216(float *);
            void bitrv216(double *);
            void bitrv216neg(float *);
            void bitrv216neg(double *);
            void bitrv208(float *);
            void bitrv208(double *);
            void bitrv208neg(float *);
            void bitrv208neg(double *);
            void cftf1st(long, float *, float *);
            void cftf1st(long, double *, double *);
            void cftb1st(long, float *, float *);
            void cftb1st(long, double *, double *);
            void cftrec4(long, float *, long, float *);
            void cftrec4(long, double *, long, double *);
            long cfttree(long, long, long, float *, long, float *);
            long cfttree(long, long, long, double *, long, double *);
            void cftleaf(long, long, float *, long, float *);
            void cftleaf(long, long, double *, long, double *);
            void cftmdl1(long, float *, float *);
            void cftmdl1(long, double *, double *);
            void cftmdl2(long, float *, float *);
            void cftmdl2(long, double *, double *);
            void cftfx41(long, float *, long, float *);
            void cftfx41(long, double *, long, double *);
            void cftf161(float *, float *);
            void cftf161(double *, double *);
            void cftf162(float *, float *);
            void cftf162(double *, double *);
            void cftf081(float *, float *);
            void cftf081(double *, double *);
            void cftf082(float *, float *);
            void cftf082(double *, double *);
            void cftf040(float *);
            void cftf040(double *);
            void cftb040(float *);
            void cftb040(double *);
            void cftx020(float *);
            void cftx020(double *);
            void rftfsub(long, float *, long, float *);
            void rftfsub(long, double *, long, double *);
            void rftbsub(long, float *, long, float *);
            void rftbsub(long, double *, long, double *);
            void dctsub(long, float *, long, float *);
            void dctsub(long, double *, long, double *);
            void dstsub(long, float *, long, float *);
            void dstsub(long, double *, long, double *);
            #ifdef USE_CDFT_THREADS
            void cftrec4_th(long, float *, long, float *);
            void cftrec4_th(long, double *, long, double *);
            void cdft_thread_create(cdft_thread_t *, void * (*)(void *), cdft_arg_t *);
            void cdft_thread_wait(cdft_thread_t);
            #endif
        public:
            clTransformS()
                { }
            virtual ~clTransformS()
                { }
            // Complex DFT
            // (2 * size, <1|-1>, input/output, work area [0] = 0 init,
            // cos/sin table)
            void cdft (long, long, float *, long *, float *);
            void cdft (long, long, double *, long *, double *);
            // Real DFT
            // (size, <1|-1>, input/output, work area, cos/sin table)
            void rdft (long, long, float *, long *, float *);
            void rdft (long, long, double *, long *, double *);
            // DCT
            void ddct (long, long, float *, long *, float *);
            void ddct (long, long, double *, long *, double *);
            // DST
            void ddst (long, long, float *, long *, float *);
            void ddst (long, long, double *, long *, double *);
            // CT of RDFT
            void dfct (long, float *, float *, long *, float *);
            void dfct (long, double *, double *, long *, double *);
            // ST of RDFT
            void dfst (long, float *, float *, long *, float *);
            void dfst (long, double *, double *, long *, double *);
            // Internal use only
            #ifdef USE_CDFT_THREADS
            void * cftrec1f_th(void *);
            void * cftrec1_th(void *);
            void * cftrec2f_th(void *);
            void * cftrec2_th(void *);
            #endif
    };

#endif

