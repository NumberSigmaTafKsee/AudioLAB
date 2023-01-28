/*

    DFT, DCT and DST routines
    Copyright (C) 1996-1999 Takuya OOURA
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

*/


/* This is decimation-in-frequency radix-2/4/8 transform */


#ifndef TRANSFORM8_HH
    #define TRANSFORM8_HH

    #ifdef TRANSFORM_INLINE
        #define T8_INLINE           inline
    #else
        #define T8_INLINE
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


    /**
        Decimation-in-frequency radix-2/4/8 transform

        \author Takuya OOURA
        \author Jussi Laako

        \note See clTransform4 for details.
    */
    class clTransform8
    {
            void makewt(long, long *, float *);
            void makewt(long, long *, double *);
            void makect(long, long *, float *);
            void makect(long, long *, double *);
            void bitrv2(long, long *, float *);
            void bitrv2(long, long *, double *);
            void bitrv2conj(long, long *, float *);
            void bitrv2conj(long, long *, double *);
            void cftfsub(long, float *, float *);
            void cftfsub(long, double *, double *);
            void cftbsub(long, float *, float *);
            void cftbsub(long, double *, double *);
            void cft1st(long, float *, float *);
            void cft1st(long, double *, double *);
            void cftmdl(long, long, float *, float *);
            void cftmdl(long, long, double *, double *);
            void rftfsub(long, float *, long, float *);
            void rftfsub(long, double *, long, double *);
            void rftbsub(long, float *, long, float *);
            void rftbsub(long, double *, long, double *);
            void dctsub(long, float *, long, float *);
            void dctsub(long, double *, long, double *);
            void dstsub(long, float *, long, float *);
            void dstsub(long, double *, long, double *);
        public:
            clTransform8 ()
                { }
            virtual ~clTransform8 ()
                { }
            // Complex DFT
            // (2 * size, <1|-1>, input/output, work area [0] = 0 init,
            // cos/sin table)
            void cdft(long, long, float *, long *, float *);
            void cdft(long, long, double *, long *, double *);
            // Real DFT
            // (size, <1|-1>, input/output, work area, cos/sin table)
            void rdft(long, long, float *, long *, float *);
            void rdft(long, long, double *, long *, double *);
            // DCT
            void ddct(long, long, float *, long *, float *);
            void ddct(long, long, double *, long *, double *);
            // DST
            void ddst(long, long, float *, long *, float *);
            void ddst(long, long, double *, long *, double *);
            // CT of RDFT
            void dfct(long, float *, float *, long *, float *);
            void dfct(long, double *, double *, long *, double *);
            // ST of RDFT
            void dfst(long, float *, float *, long *, float *);
            void dfst(long, double *, double *, long *, double *);
    };

#endif

