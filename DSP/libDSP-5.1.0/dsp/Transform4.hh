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


/* This is decimation-in-frequency radix-2/4 transform */


#ifndef TRANSFORM4_HH
    #define TRANSFORM4_HH

    #ifdef TRANSFORM_INLINE
        #define T4_INLINE           inline
    #else
        #define T4_INLINE
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
        Decimation-in-frequency radix-2/4 transform

        \author Takuya OOURA
        \author Jussi Laako

        \note DFT parts of this class has been wrapped into clDSPOp with
        easier interface.
    */
    class clTransform4
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
            clTransform4 ()
                { }
            virtual ~clTransform4 ()
                { }
            /**
                Complex DFT

                Forward:
                \f[X(k)=\sum_{j=0}^{N-1}x(j)\exp(2\pi{ijk/N}), 0\leq{k}\leq{N-1}\f]
                Inverse:
                \f[X(k)=\sum_{j=0}^{N-1}x(j)\exp(-2\pi{ijk/N}), 0\leq{k}\leq{N-1}\f]

                Input data:
                \f[\left\vert\begin{array}{l}
                a(2j)=\Re_{x}(j)\\
                a(2j+1)=\Im_{x}(j)
                \end{array}\right., 0\leq{j}\leq{N-1}\f]
                Output data:
                \f[\left\vert\begin{array}{l}
                a(2k)=\Re_{X}(k)\\
                a(2k+1)=\Im_{X}(k)
                \end{array}\right., 0\leq{k}\leq{N-1}\f]

                Work area size:
                \f[2+\sqrt{N}\f]

                Scaling factor:
                \f[\frac{1}{N}\f]

                \param Size Window size, 2*N
                \param Direction Direction of transform; 1 = fwd, -1 = inv
                \param Data Source & destination data, size 2*N
                \param WorkArea Working area, [0] = 0 to initialize
                \param CosSin cos/sin table, size N/2
            */
            // (2 * size, <1|-1>, input/output, work area [0] = 0 init,
            // cos/sin table)
            void cdft(long, long, float *, long *, float *);
            /// \overload
            void cdft(long, long, double *, long *, double *);
            /**
                Real DFT

                Forward:
                \f[\left\vert\begin{array}{ll}
                \Re_{X}(k)=\sum_{j=0}^{N-1}a(j)\cos(2\pi{jk/N})&, 0\leq{k}\leq{N/2}\\
                \Im_{X}(k)=\sum_{j=0}^{N-1}a(j)\sin(2\pi{jk/N})&, 0<k<N/2
                \end{array}\right.\f]
                Inverse:
                \f[
                a(k)=\frac{\left(\Re_{X}(0)+\Re_{X}(\frac{N}{2})cos(\pi{k})\right)}{2}+
                \sum_{j=1}^{N/2-1}\Re_{X}(j)\cos(2\pi{jk/N})+
                \sum_{j=1}^{N/2-1}\Im_{X}(j)\sin(2\pi{jk/N}), 0\leq{k}\leq{N-1}
                \f]

                Output data:
                \f[\left\vert\begin{array}{ll}
                a(2k)=\Re_{X}(k)&, 0\leq{k}<N/2\\
                a(2k+1)=\Im_{X}(k)&, 0<k<N/2\\
                a(1)=\Re_{X}(N/2)
                \end{array}\right.\f]
                Input data:
                \f[\left\vert\begin{array}{ll}
                a(2j)=\Re_{X}(j)&, 0\leq{j}<N/2\\
                a(2j+1)=\Im_{X}(j)&, 0<j<N/2\\
                a(1)=\Re_{X}(N/2)
                \end{array}\right.\f]

                Work area size:
                \f[2+\sqrt{\frac{N}{2}}\f]

                Scaling factor:
                \f[\frac{2}{N}\f]

                \param Size Window size, N
                \param Direction Direction of transform; 1 = fwd, -1 = inv
                \param Data Source & destination data, size N
                \param WorkArea Working area, [0] = 0 to initialize
                \param CosSin cos/sin table, size N/2
            */
            // (size, <1|-1>, input/output, work area, cos/sin table)
            void rdft(long, long, float *, long *, float *);
            /// \overload
            void rdft(long, long, double *, long *, double *);
            /**
                DCT

                Forward:
                \f[C(k)=\sum_{j=0}^{N-1}a(j)\cos(\pi{j(k+1/2)/N}), 0\leq{k}\leq{N-1}\f]
                Inverse:
                \f[C(k)=\sum_{j=0}^{N-1}a(j)\cos(\pi(j+1/2)k/N), 0\leq{k}\leq{N-1}\f]

                Input/output data:
                \f[a(k)=C(k), 0\leq{k}\leq{N-1}\f]

                Work area size:
                \f[2+\sqrt{\frac{N}{2}}\f]

                Scaling factors:
                \f[(0)0.5, \frac{2}{N}\f]

                \param Size Window size, N
                \param Direction Direction of transform; 1 = fwd, -1 = inv
                \param Data Source & destination data, size N
                \param WorkArea Working area, [0] = 0 to initialize
                \param CosSin cos/sin table, size N*5/4
            */
            void ddct(long, long, float *, long *, float *);
            /// \overload
            void ddct(long, long, double *, long *, double *);
            /**
                DST

                Forward:
                \f[S(k)=\sum_{j=1}^{N}A(j)\sin(\pi{j(k+1/2)/N}), 0\leq{k}\leq{N-1}\f]
                Inverse:
                \f[S(k)=\sum_{j=0}^{N-1}a(j)\sin(\pi(j+1/2)k/N), 0<k\leq{N}\f]

                Input data (forward):
                \f[\left\vert\begin{array}{ll}
                a(j)=A(j)&, 0<j<N\\
                a(0)=A(N)
                \end{array}\right.\f]
                Output data (forward):
                \f[a(k)=S(k), 0\leq{k}\leq{N-1}\f]
                Output data (inverse):
                \f[\left\vert\begin{array}{ll}
                a(k)=S(k), 0<k<N\\
                a(0)=S(N)
                \end{array}\right.\f]

                Work area size:
                \f[2+\sqrt{\frac{N}{2}}\f]

                Scaling factors:
                \f[(0)0.5, \frac{2}{N}\f]

                \param Size Window size, N
                \param Direction Direction of transform; 1 = fwd, -1 = inv
                \param Data Source & destination data, size N
                \param WorkArea Working area, [0] = 0 to initialize
                \param CosSin cos/sin table, size N*5/4
            */
            void ddst(long, long, float *, long *, float *);
            /// \overload
            void ddst(long, long, double *, long *, double *);
            /**
                CT of RDFT

                Definition:
                \f[C(k)=\sum_{j=0}^{N}a(j)\cos(\pi{jk/N}), 0\leq{k}\leq{N}\f]

                Output data:
                \f[a(k)=C(k), 0\leq{k}\leq{N}\f]

                Work area size:
                \f[2+\sqrt{\frac{N}{4}}\f]

                Scaling factors:
                \f[(0)0.5, (N)0.5, \frac{2}{N}\f]

                \param Size Window size, N
                \param Data Source & destination data, size N+1
                \param Scratch Scratch pad, size N/2+1
                \param WorkArea Working area, [0] = 0 to initialize
                \param CosSin cos/sin table, size N*5/8
            */
            void dfct(long, float *, float *, long *, float *);
            /// \overload
            void dfct(long, double *, double *, long *, double *);
            /**
                ST of RDFT

                Definition:
                \f[S(k)=\sum_{j=1}^{N-1}a(j)\sin(\pi{jk/N}), 0<k<N\f]

                Output data:
                \f[a(k)=S(k), 0<k<N\f]

                Work area size:
                \f[2+\sqrt{\frac{N}{4}}\f]

                Scaling factor:
                \f[\frac{2}{N}\f]

                \param Size Window size, N
                \param Data Source & destination data, size N+1
                \param Scratch Scratch pad, size N/2
                \param WorkArea Working area, [0] = 0 to initialize
                \param CosSin cos/sin table, size N*5/8
            */
            void dfst(long, float *, float *, long *, float *);
            /// \overload
            void dfst(long, double *, double *, long *, double *);
    };

#endif

