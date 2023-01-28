/*

    DSP operations, header
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


#ifdef USE_INTEL_MATH
    #include <mathimf.h>
#else
    #include <math.h>
#endif
#include <float.h>

#include <commoncpp/Alloc.hh>

#include <dsp/dsptypes.h>
#include <dsp/DSPConfig.hh>
#include <dsp/Transform4.hh>
#include <dsp/Transform8.hh>
#include <dsp/TransformS.hh>
 
#ifdef DSP_USE_FFTW
    #include <fftw3.h>
#endif


#ifndef DSPOP_HH
    #define DSPOP_HH

    // To help compilers common subexpression elimination optimization
    #ifdef __GNUG__
        #define CONSTFUNC       __const__
    #else
        #define CONSTFUNC
    #endif
    
    // Enable/disable function inlining
    #ifndef __INTEL_COMPILER
        #define INLINE          inline
    #else
        #define INLINE
    #endif

    // function definitions for C-libraries lacking ISO C 9x standard
    #if (defined(__BCPLUSPLUS__) || defined(SOLARIS) || defined(OSX))
        #if (!defined(OSX))
            #define isgreater(x, y) (x > y)
            #define isless(x, y)    (x < y)
        #endif
        #define fmax(x, y)      ((x >= y) ? x : y)
        #define fmin(x, y)      ((x <= y) ? x : y)
        #define fmaxf(x, y)     fmax(x, y)
        #define fminf(x, y)     fmin(x, y)
        #define fma(x, y, z)    (z + x * y)
        #define fmaf(x, y, z)   fma(x, y, z)
        #define fabsf(x)        fabs(x)
        #define powf(x, y)      pow(x, y)
        #define sqrtf(x)        sqrt(x)
        #define expf(x)         exp(x)
        #define logf(x)         log(x)
        #define log10f(x)       log10(x)
        #define sinf(x)         sin(x)
        #define sinhf(x)        sinh(x)
        #define cosf(x)         cos(x)
        #define coshf(x)        cosh(x)
        #define atanf(x)        atan(x)
        #define atan2f(x, y)    atan2(x, y)
        #define acosf(x)        acos(x)
        #define acosh(x)        (log(x + sqrt(x * x - 1)))
        #define acoshf(x)       (logf(x + sqrtf(x * x - 1)))
        #define hypotf(x, y)    hypot(x, y)
    #endif

    // Microsoft VisualC++ partially conforms with ISO 9x standard
    #ifdef _MSC_VER
        #define fmax(x, y)      ((x >= y) ? x : y)
        #define fmin(x, y)      ((x <= y) ? x : y)
        #define fmaxf(x, y)     fmax(x, y)
        #define fminf(x, y)     fmin(x, y)
        #define fma(x, y, z)    (z + x * y)
        #define fmaf(x, y, z)   fma(x, y, z)
        #define acosh(x)        (log(x + sqrt(x * x - 1)))
        #define acoshf(x)       (logf(x + sqrtf(x * x - 1)))
        #define hypotf(x, y)    hypot(x, y)
        #define isgreater(x, y) (x > y)
        #define isless(x, y)    (x < y)
    #endif

    // Intel C++ partially conforms with ISO 9x standard
    #ifdef __INTEL_COMPILER
        #undef isgreater
        #undef isless
        #define isgreater(x, y) (x > y)
        #define isless(x, y)    (x < y)
    #endif

    // Maximum iterations for modified zero order Bessel function of first kind
    #define DSP_MAXBESSEL       32L

    // Filename for FFTW wisdom
    #define DSP_WISDOM_FILE     "fftw.wisdom"


    /**
        Class specialization to support automatic typecasts to cartesian/polar
        datatypes.
    */
    class clDSPAlloc : public clAlloc
    {
        public:
            clDSPAlloc () {}
            clDSPAlloc (const clDSPAlloc &CopySrc) : clAlloc(CopySrc) {}
            clDSPAlloc (long lAllocSize) : clAlloc(lAllocSize) {}
            operator stSCplx *()
                { return ((stpSCplx) vpPtr); }
            operator stDCplx *()
                { return ((stpDCplx) vpPtr); }
            operator stSPolar *()
                { return ((stpSPolar) vpPtr); }
            operator stDPolar *()
                { return ((stpDPolar) vpPtr); }
            #if defined(DSP_USE_FFTW)
            operator fftwf_complex *()
                { return ((fftwf_complex *) vpPtr); }
            operator fftw_complex *()
                { return ((fftw_complex *) vpPtr); }
            #endif
    };

    /**
        Class of DSP operations
    */
    class clDSPOp
    {
            long lPrevSrcCount;
            long lPrevDestCount;
            float fPI;
            double dPI;
            // --- Dynamically allocated arrays
            // FIR
            int iFIRDlyIdx;
            long lFIRLength;
            clDSPAlloc FIRCoeff;
            clDSPAlloc FIRBuf;
            clDSPAlloc FIRWork;
            // IIR
            float fpIIR_C[5];
            float fpIIR_X[3];
            float fpIIR_Y[2];
            double dpIIR_C[5];
            double dpIIR_X[3];
            double dpIIR_Y[2];
            // FFT (and other transforms)
            bool bFFTInitialized;
            bool bRealTransform;
            long lFFTLength;
            float fFFTScale;
            double dFFTScale;
            long *lpSBitRevWork;
            long *lpDBitRevWork;
            float *fpCosSinTable;
            double *dpCosSinTable;
            void *vpSTfrm;
            void *vpDTfrm;
            clDSPAlloc SBitRevWork;
            clDSPAlloc DBitRevWork;
            clDSPAlloc SCosSinTable;
            clDSPAlloc DCosSinTable;
            clDSPAlloc FFTBuf;
            #if defined(DSP_USE_FFTW)
            fftwf_plan fftwpSPlan;
            fftwf_plan fftwpSIPlan;
            fftw_plan fftwpDPlan;
            fftw_plan fftwpDIPlan;
            #endif
            #if defined(DSP_USE_RADIX4)
            clTransform4 Tfrm;
            #elif defined(DSP_USE_RADIX8)
            clTransform8 Tfrm;
            #else  // Split-radix
            clTransformS Tfrm;
            #endif
        protected:
            /**
                Cartesian to polar conversion.

                \f[V=\sqrt{\Re^2+\Im^2}\f]
                \f[\varphi=\arctan{\left(\frac{\Im}{\Re}\right)}\f]

                \param Magn Magnitude
                \param Phase Phase
                \param Real Real
                \param Imag Imaginary
            */
            static void Cart2Polar (float *, float *, float, float);
            /// \overload
            static void Cart2Polar (double *, double *, double, double);
            /** 
                \overload
                \param Magn Magnitude
                \param Phase Phase
                \param Cplx Cartesian
            */
            static void Cart2Polar (float *, float *, const stpSCplx);
            /// \overload
            static void Cart2Polar (double *, double *, const stpDCplx);
            /**
                \overload
                \param Polar Polar
                \param Cplx Cartesian
            */
            static void Cart2Polar (stpSPolar, const stpSCplx);
            /// \overload
            static void Cart2Polar (stpDPolar, const stpDCplx);
            /**
                \overload
                \param Coord Polar & Cartesian
            */
            static void Cart2Polar (utpSCoord);
            /// \overload
            static void Cart2Polar (utpDCoord);
            /**
                Polar to cartesian conversion.

                \f[\Re = V\cos(\varphi)\f]
                \f[\Im = V\sin(\varphi)\f]

                \note See Cart2Polar for details
            */
            static void Polar2Cart (float *, float *, float, float);
            /// \overload
            static void Polar2Cart (double *, double *, double, double);
            /// \overload
            static void Polar2Cart (stpSCplx, float, float);
            /// \overload
            static void Polar2Cart (stpDCplx, double, double);
            /// \overload
            static void Polar2Cart (stpSCplx, const stpSPolar);
            /// \overload
            static void Polar2Cart (stpDCplx, const stpDPolar);
            /// \overload
            static void Polar2Cart (utpSCoord);
            /// \overload
            static void Polar2Cart (utpDCoord);
            /**
                Complex addition

                \f[(z_{r},z_{i})=(x_{r}+y_{r},x_{i}+y_{i})\f]

                \param CplxDest Source & destination
                \param CplxSrc Source
            */
            static void CplxAdd (stpSCplx, const stpSCplx);
            /// \overload
            static void CplxAdd (stpDCplx, const stpDCplx);
            /**
                \overload
                \param CplxDest Destination
                \param CplxSrc1 Source 1
                \param CplxSrc2 Source 2
            */
            static void CplxAdd (stpSCplx, const stpSCplx, const stpSCplx);
            /// \overload
            static void CplxAdd (stpDCplx, const stpDCplx, const stpDCplx);
            /**
                Complex subtraction

                \f[(z_{r},z_{i})=(x_{r}-y_{r},x_{i}-y_{i})\f]

                \param CplxDest Source & destination
                \param CplxSrc Source
            */
            static void CplxSub (stpSCplx, const stpSCplx);
            /// \overload
            static void CplxSub (stpDCplx, const stpDCplx);
            /**
                \overload
                \param CplxDest Destination
                \param CplxSrc1 Source 1
                \param CplxSrc2 Source 2
            */
            static void CplxSub (stpSCplx, const stpSCplx, const stpSCplx);
            /// \overload
            static void CplxSub (stpDCplx, const stpDCplx, const stpDCplx);
            /**
                Complex multiply with real value

                \f[(z_r,z_i)=(x_{r}y,x_{i}y)\f]

                \param CplxDest Source & destination
                \param Src Source
            */
            static void CplxMul (stpSCplx, float);
            /// \overload
            static void CplxMul (stpDCplx, double);
            /**
                Complex multiply

                \f[\Re_{z}=\Re_{x}\Re_{y}-\Im_{x}\Im_{y}\f]
                \f[\Im_{z}=\Re_{x}\Im_{y}+\Re_{y}\Im_{x}\f]

                \param CplxDest Source & destination
                \param CplxSrc Source
            */
            static void CplxMul (stpSCplx, const stpSCplx);
            /// \overload
            static void CplxMul (stpDCplx, const stpDCplx);
            /**
                \overload
                \param CplxDest Destination
                \param CplxSrc Source 1
                \param CplxSrc Source 2
            */
            static void CplxMul (stpSCplx, const stpSCplx, const stpSCplx);
            /// \overload
            static void CplxMul (stpDCplx, const stpDCplx, const stpDCplx);
            /**
                Complex multiply with complex conjugate

                \f[\Re_{z}=\Re_{x}\Re_{y}-\Im_{x}(-\Im_{y})\f]
                \f[\Im_{z}=\Re_{x}(-\Im_{y})+\Re_{y}\Im_{x}\f]

                \param CplxDest Source & destination
                \param CplxSrc Source
            */
            static void CplxMulC (stpSCplx, const stpSCplx);
            /// \overload
            static void CplxMulC (stpDCplx, const stpDCplx);
            /**
                \overload
                \param CplxDest Destination
                \param CplxSrc1 Source 1
                \param CplxSrc2 Source 2
            */
            static void CplxMulC (stpSCplx, const stpSCplx, const stpSCplx);
            /// \overload
            static void CplxMulC (stpDCplx, const stpDCplx, const stpDCplx);
            /**
                Complex division

                \f[\Re_{z}=\frac{\Re_{x}\Re_{y}+\Im_{x}\Im_{y}}{\Re_{y}^2+\Im_{y}^2}\f]
                \f[\Im_{z}=\frac{\Im_{x}\Re_{y}-\Re_{x}\Im_{y}}{\Re_{y}^2+\Im_{y}^2}\f]

                \param CplxDest Source & destination
                \param CplxSrc Source
            */
            static void CplxDiv (stpSCplx, const stpSCplx);
            /// \overload
            static void CplxDiv (stpDCplx, const stpDCplx);
            /**
                \overload
                \param CplxDest Destination
                \param CplxSrc Source 1
                \param CplxSrc Source 2
            */
            static void CplxDiv (stpSCplx, const stpSCplx, const stpSCplx);
            /// \overload
            static void CplxDiv (stpDCplx, const stpDCplx, const stpDCplx);
            // In-place is allowed for following
            /**
                Complex exp function (e raised to y)

                \f[\Re_{z}=\exp(\Re_{y})\cos(\Im_{y})\f]
                \f[\Im_{z}=\exp(\Re_{y})\sin(\Im_{y})\f]

                \param CplxDest Destination
                \param CplxSrc Source
            */
            static void CplxExp (stpSCplx, const stpSCplx);
            static void CplxExp (stpDCplx, const stpDCplx);
            /**
                Complex natural logarithm

                \f[\Re_{z}=\ln(|y|)\f]
                \f[\Im_{z}=\arg(y)\f]
            */
            static void CplxLog (stpSCplx, const stpSCplx);
            static void CplxLog (stpDCplx, const stpDCplx);
            /**
                Complex 10-base logarithm

                \f[\Re_{z}=\log(|y|)\f]
                \f[\Im_{z}=\arg(y)\f]
            */
            static void CplxLog10 (stpSCplx, const stpSCplx);
            static void CplxLog10 (stpDCplx, const stpDCplx);
            /**
                Complex x raised to y

                \f[z=\exp(\ln(x)y)\f]
            */
            static void CplxPow (stpSCplx, const stpSCplx, const stpSCplx);
            static void CplxPow (stpDCplx, const stpDCplx, const stpDCplx);
            /**
                Complex root y of x

                \f[z=x^{\frac{1}{y}}\f]
            */
            static void CplxRoot (stpSCplx, const stpSCplx, const stpSCplx);
            static void CplxRoot (stpDCplx, const stpDCplx, const stpDCplx);
            static void CplxConj (stpSCplx spCplx)
                { spCplx->I = -(spCplx->I); }
            static void CplxConj (stpDCplx spCplx)
                { spCplx->I = -(spCplx->I); }
            static void CplxConj (stpSCplx, const stpSCplx);
            static void CplxConj (stpDCplx, const stpDCplx);
            /**
                Return multiple of n

                \f[!n\f]

                \param Value Value of n
            */
            static double Multiple (long);
            /**
                Zero-order modified Bessel function of the first kind.

                \f[I_{0}(x)=\sum_{k=0}^{K}\left[\frac{(x/2)^{k}}{k!}\right]^{2}\f]

                \param Value Value of x
            */
            static float ModZeroBessel (float);
            /// \overload
            static double ModZeroBessel (double);
            /**
                nth-order Chebyshev polynomial
                
                \f[T_{n}(x)=\left\{\begin{array}{l}\cos(n\cos^{-1}x), |x|\leq{1}\\
                \cosh(n\cosh^{-1}x), |x|>1\end{array}\right.\f]
                
                \param Order Order of polynomial
                \param Value Value of x
            */
            static float ChebyshevPolynom (float, float);
            static double ChebyshevPolynom (double, double);
        public:
            clDSPOp();
            virtual ~clDSPOp();
            /**
                Round floatingpoint number to integer
            */
            static signed long Round (float);
            static signed long Round (double);
            /**
                Add single value to vector

                \f[x(i)=x(i)+y, 0\leq{i}\leq{N-1}\f]

                \param Dest Source & destination, x()
                \param Src Source, y
                \param Count Vector length, N
            */
            static void Add (float *, float, long);
            /// \overload
            static void Add (double *, double, long);
            /// \overload
            static void Add (stpSCplx, stSCplx, long);
            /// \overload
            static void Add (stpDCplx, stDCplx, long);
            /**
                Add two vectors

                \f[x(i)=x(i)+y(i), 0\leq{i}\leq{N-1}\f]

                \param Dest Source & destination, x()
                \param Src Source, y()
                \param Count Vector length, N
            */
            static void Add (float *, const float *, long);
            /// \overload
            static void Add (double *, const double *, long);
            /// \overload
            static void Add (stpSCplx, const stpSCplx, long);
            /// \overload
            static void Add (stpDCplx, const stpDCplx, long);
            /**
                Add two vectors and return result in third

                \f[z(i)=x(i)+y(i), 0\leq{i}\leq{N-1}\f]

                \param Dest Destination, z()
                \param Src1 Source 1, x()
                \param Src2 Source 2, y()
                \param Count Vector length, N
            */
            static void Add (float *, const float *, const float *, long);
            /// \overload
            static void Add (double *, const double *, const double *, long);
            /// \overload
            static void Add (stpSCplx, const stpSCplx, const stpSCplx, long);
            /// \overload
            static void Add (stpDCplx, const stpDCplx, const stpDCplx, long);
            /**
                Subtract single value from vector

                \f[x(i)=x(i)-y, 0\leq{i}\leq{N-1}\f]

                \param Dest Source & destination, x()
                \param Src Source, y
                \param Count Vector length, N
            */
            static void Sub (float *, float, long);
            /// \overload
            static void Sub (double *, double, long);
            /// \overload
            static void Sub (stpSCplx, stSCplx, long);
            /// \overload
            static void Sub (stpDCplx, stDCplx, long);
            /**
                Subtract two vectors

                \f[x(i)=x(i)-y(i), 0\leq{i}\leq{N-1}\f]

                \param Dest Source & destination, x()
                \param Src Source, y()
                \param Count Vector length, N
            */
            static void Sub (float *, const float *, long);
            /// \overload
            static void Sub (double *, const double *, long);
            /// \overload
            static void Sub (stpSCplx, const stpSCplx, long);
            /// \overload
            static void Sub (stpDCplx, const stpDCplx, long);
            /**
                Subtract two vectors and return result in third

                \f[z(i)=x(i)-y(i), 0\leq{i}\leq{N-1}\f]

                \param Dest Destination, z()
                \param Src1 Source 1, x()
                \param Src2 Source 2, y()
                \param Count Vector length, N
            */
            static void Sub (float *, const float *, const float *, long);
            /// \overload
            static void Sub (double *, const double *, const double *, long);
            /// \overload
            static void Sub (stpSCplx, const stpSCplx, const stpSCplx, long);
            /// \overload
            static void Sub (stpDCplx, const stpDCplx, const stpDCplx, long);
            /**
                Multiply vector with single value in-place

                \f[x(i)=x(i)y, 0\leq{i}\leq{N-1}\f]

                \param Dest Source & destination, x()
                \param Src Source, y
                \param Count Vector length, N
            */
            static void Mul (float *, float, long);
            /// \overload
            static void Mul (double *, double, long);
            /// \overload
            static void Mul (stpSCplx, float, long);
            /// \overload
            static void Mul (stpDCplx, double, long);
            /// \overload
            static void Mul (stpSCplx, stSCplx, long);
            /// \overload
            static void Mul (stpDCplx, stDCplx, long);
            /**
                Multiply vector with single value out-of-place

                \f[z(i)=x(i)y, 0\leq{i}\leq{N-1}\f]

                \param Dest Destination, z()
                \param Src1 Source 1, x()
                \param Src2 Source 2, y
                \param Count Vector length, N
            */
            static void Mul (float *, const float *, float, long);
            /// \overload
            static void Mul (double *, const double *, double, long);
            /**
                Multiply two vectors

                \f[x(i)=x(i)y(i), 0\leq{i}\leq{N-1}\f]

                \param Dest Source & destination, x()
                \param Src Source, y()
                \param Count Vector length, N
            */
            static void Mul (float *, const float *, long);
            /// \overload
            static void Mul (double *, const double *, long);
            /// \overload
            static void Mul (stpSCplx, const float *, long);
            /// \overload
            static void Mul (stpDCplx, const double *, long);
            /// \overload
            static void Mul (stpSCplx, const stpSCplx, long);
            /// \overload
            static void Mul (stpDCplx, const stpDCplx, long);
            /**
                Multiply two vectors and return result in third

                \f[z(i)=x(i)y(i), 0\leq{i}\leq{N-1}\f]

                \param Dest Destination, z()
                \param Src1 Source 1, x()
                \param Src2 Source 2, y()
                \param Count Vector length, N
            */
            static void Mul (float *, const float *, const float *, long);
            /// \overload
            static void Mul (double *, const double *, const double *, long);
            /// \overload
            static void Mul (stpSCplx, const stpSCplx, const stpSCplx, long);
            /// \overload
            static void Mul (stpDCplx, const stpDCplx, const stpDCplx, long);
            /**
                Multiply vector with others complex conjugate

                \f[X(i)=X(i)Y(i)^{*}, 0\leq{i}\leq{N-1}\f]
                \f[Z(i)=X(i)Y(i)^{*}, 0\leq{i}\leq{N-1}\f]

                \param Dest Source & destination, X()
                \param Src Source, Y()
                \param Count Vector length, N
            */
            static void MulC (stpSCplx, const stpSCplx, long);
            /// \overload
            static void MulC (stpDCplx, const stpDCplx, long);
            /// \overload
            static void MulC (stpSCplx, const stpSCplx, const stpSCplx, long);
            /// \overload
            static void MulC (stpDCplx, const stpDCplx, const stpDCplx, long);
            /**
                Multiply two vectors in-place with third vector

                \f[\left\vert\begin{array}{l}x(i)=x(i)z(i)\\
                y(i)=y(i)z(i)\end{array}\right., 0\leq{i}\leq{N-1}\f]

                \param Dest1 Source & destination 1, x()
                \param Dest2 Source & destination 2, y()
                \param Src Source, z()
                \param Count Vector length, N
            */
            static void Mul2 (float *, float *, const float *, long);
            /// \overload
            static void Mul2 (double *, double *, const double *, long);
            /**
                Multiply two vectors out-of-place with third vector and
                return results in separate two vectors

                \param Dest1 Destination 1
                \param Dest2 Destination 2
                \param Src1 Source 1
                \param Src2 Source 2
                \param Count Vector length
            */
            static void Mul2 (float *, float *, const float *, 
                const float *, const float *, long);
            /// \overload
            static void Mul2 (double *, double *, const double *, 
                const double *, const double *, long);
            /**
                Divide vector with single value

                \f[x(i)=\frac{x(i)}{y}, 0\leq{i}\leq{N-1}\f]

                \param Dest Source & destination, x()
                \param Src Source, y
                \param Count Vector length, N
            */
            static void Div (float *, float, long);
            /// \overload
            static void Div (double *, double, long);
            /// \overload
            static void Div (stpSCplx, stSCplx, long);
            /// \overload
            static void Div (stpDCplx, stDCplx, long);
            /**
                Divide two vectors

                \f[x(i)=\frac{x(i)}{y(i)}, 0\leq{i}\leq{N-1}\f]

                \param Dest Source & destination, x()
                \param Src Source, y()
                \param Count Vector length, N
            */
            static void Div (float *, const float *, long);
            /// \overload
            static void Div (double *, const double *, long);
            /// \overload
            static void Div (stpSCplx, const stpSCplx, long);
            /// \overload
            static void Div (stpDCplx, const stpDCplx, long);
            /**
                Divide two vectors and return result in third

                \f[z(i)=\frac{x(i)}{y(i)}, 0\leq{i}\leq{N-1}\f]

                \param Dest Destination, z()
                \param Src1 Source 1, x()
                \param Src2 Source 2, y()
                \param Count Vector length, N
            */
            static void Div (float *, const float *, const float *, long);
            /// \overload
            static void Div (double *, const double *, const double *, long);
            /// \overload
            static void Div (stpSCplx, const stpSCplx, const stpSCplx, long);
            /// \overload
            static void Div (stpDCplx, const stpDCplx, const stpDCplx, long);
            /**
                Calculate 1 / value in-place

                \f[x(i)=\frac{1}{x(i)}, 0\leq{i}\leq{N-1}\f]

                \param Vect Source & destination
                \param Count Vector length, N
            */
            static void Div1x (float *, long);
            /// \overload
            static void Div1x (double *, long);
            /**
                Calculate 1 / value out-of-place

                \f[y(i)=\frac{1}{x(i)}, 0\leq{i}\leq{N-1}\f]

                \param Dest Destination, y()
                \param Src Source, x()
                \param Count Vector length, N
            */
            static void Div1x (float *, const float *, long);
            /// \overload
            static void Div1x (double *, const double *, long);
            /**
                Multiply-add vector in-place
                
                \f[x(i)=x(i)s+o, 0\leq{i}\leq{N-1}\f]
                
                \param Vect Source & destination, x()
                \param Mul Multiply value, s
                \param Add Addition value, o
                \param Count Vector length, N
            */
            static void MulAdd (float *, float, float, long);
            static void MulAdd (double *, double, double, long);
            /**
                Multiply-add vector out-of-place
                
                \f[y(i)=x(i)s+o, 0\leq{i}\leq{N-1}\f]
                
                \param Dest Destination, y()
                \param Src Source, x()
                \param Mul Multiply value, s
                \param Add Addition value, o
                \param Count Vector length, N
            */
            static void MulAdd (float *, const float *, float, float, long);
            static void MulAdd (double *, const double *, double, double, long);
            /**
                Multiply-add complex vectors out-of-place
                
                \f[Z(i)=X(i)Y(i), 0\leq{i}\leq{N-1}\f]
                
                \param Dest Destination, Z()
                \param Src1 Source 1, X()
                \param Src2 Source 2, Y()
                \param Count Vector length, N
            */
            static void MulAddC (stpSCplx, const stpSCplx, const stpSCplx, 
                long);
            static void MulAddC (stpDCplx, const stpDCplx, const stpDCplx,
                long);
            /**
                Get absolute value in-place

                \f[x(i)=|x(i)|, 0\leq{i}\leq{N-1}\f]

                \param Vect Source & destination, x()
                \param Count Vector length, N
            */
            static void Abs (float *, long);
            /// \overload
            static void Abs (double *, long);
            /**
                Get absolute value out-of-place

                \f[y(i)=|x(i)|, 0\leq{i}\leq{N-1}\f]

                \param Dest Destination, y()
                \param Src Source, x()
                \param Count Vector length, N
            */
            static void Abs (float *, const float *, long);
            /// \overload
            static void Abs (double *, const double *, long);
            /**
                Calculate square-roots in-place

                \f[x(i)=\sqrt{x(i)}, 0\leq{i}\leq{N-1}\f]

                \param Vect Source & destination, x()
                \param Count Vector length, N
            */
            static void Sqrt (float *, long);
            /// \overload
            static void Sqrt (double *, long);
            /**
                Calculate squre-roots out-of-place

                \f[y(i)=\sqrt{x(i)}, 0\leq{i}\leq{N-1}\f]

                \param Dest Destination, y()
                \param Src Source, x()
                \param Count Vector length, N
            */
            static void Sqrt (float *, const float *, long);
            /// \overload
            static void Sqrt (double *, const double *, long);
            /**
                Set all vector elements to zero

                \param Vect Vector
                \param Count Vector length
            */
            static void Zero (float *, long);
            /// \overload
            static void Zero (double *, long);
            /// \overload
            static void Zero (stpSCplx, long);
            /// \overload
            static void Zero (stpDCplx, long);
            /**
                Set all vector elements to specified value

                \param Dest Destination
                \param Src Source value
                \param Count Vector length
            */
            static void Set (float *, float, long);
            /// \overload
            static void Set (double *, double, long);
            /// \overload
            static void Set (stpSCplx, stSCplx, long);
            /// \overload
            static void Set (stpDCplx, stDCplx, long);
            /**
                Set specified vector elements to specified value with
                boundary check

                \param Dest Destination
                \param Src Source value
                \param Start Starting index
                \param Count Number of elements to set
                \param Length Vector length
            */
            static void Set (float *, float, long, long, long);
            /// \overload
            static void Set (double *, double, long, long, long);
            /// \overload
            static void Set (stpSCplx, stSCplx, long, long, long);
            /// \overload
            static void Set (stpDCplx, stDCplx, long, long, long);
            /**
                Clip vector in-place to specified value

                \f[x(i)=\left\{\begin{array}{ll}x(i)&,x(i)\leq{y}\\
                y&,x(i)>y\end{array}\right., 0\leq{i}\leq{N-1}\f]

                \param Vect Source & destination, x()
                \param Value Clipping value, y
                \param Count Vector length, N
            */
            static void Clip (float *, float, long);
            /// \overload
            static void Clip (double *, double, long);
            /**
                Clip vector out-of-place to specified value

                \f[z(i)=\left\{\begin{array}{ll}x(i)&,x(i)\leq{y}\\
                y&,x(i)>y\end{array}\right., 0\leq{i}\leq{N-1}\f]

                \param Dest Destination, z()
                \param Src Source, x()
                \param Value Clipping value, y
                \param Count Vector length, N
            */
            static void Clip (float *, const float *, float, long);
            /// \overload
            static void Clip (double *, const double *, double, long);
            /**
                Clip vector in-place to fit between specified values

                \f[x(i)=\left\{\begin{array}{ll}a&,x(i)<a\\
                x(i)&,a\leq{x(i)}\leq{b}\\
                b&,x(i)>b\end{array}\right., 0\leq{i}\leq{N-1}\f]

                \param Vect Source & destination, x()
                \param Min Minimum value, a
                \param Max Maximum value, b
                \param Count Vector length, N
            */
            static void Clip (float *, float, float, long);
            /// \overload
            static void Clip (double *, double, double, long);
            /**
                Clip vector out-of-place to fit between specified values

                \f[y(i)=\left\{\begin{array}{ll}a&,x(i)<a\\
                x(i)&,a\leq{x(i)}\leq{b}\\
                b&,x(i)>b\end{array}\right., 0\leq{i}\leq{N-1}\f]

                \param Dest Destination, y()
                \param Src Source, x()
                \param Min Minimum value, a
                \param Max Maximum value, b
                \param Count Vector length, N
            */
            static void Clip (float *, const float *, float, float, long);
            /// \overload
            static void Clip (double *, const double *, double, double, long);
            /**
                Clip smallest values in vector to zero (in-place)

                \f[x(i)=\left\{\begin{array}{ll}0&,x(i)<0\\
                x(i)&,x(i)\leq{0}\end{array}\right., 0\leq{i}\leq{N-1}\f]

                \param Vect Source & destination, x()
                \param Count Vector length, N
            */
            static void ClipZero (float *, long);
            /// \overload
            static void ClipZero (double *, long);
            /**
                Clip smallest values in vector to zero (out-of-place)

                \f[y(i)=\left\{\begin{array}{ll}0&,x(i)<0\\
                x(i)&,x(i)\leq{0}\end{array}\right., 0\leq{i}\leq{N-1}\f]

                \param Dest Destination, y()
                \param Src Source, x()
                \param Count Vector length, N
            */
            static void ClipZero (float *, const float *, long);
            /// \overload
            static void ClipZero (double *, const double *, long);
            /**
                Copy source vector to destination vector(s)

                \param Dest Destination vector
                \param Src Source vector
                \param Count Vector length
            */
            static void Copy (float *, const float *, long);
            /// \overload
            static void Copy (double *, const double *, long);
            /// \overload
            static void Copy (stpSCplx, const stpSCplx, long);
            /// \overload
            static void Copy (stpDCplx, const stpDCplx, long);
            /// \overload
            static void Copy (float *, float *, const float *, long);
            /// \overload
            static void Copy (double *, double *, const double *, long);
            /**
                Convolve two vectors (same length, finite) without lag
                \note Circular convolution, result is not scaled

                \f[z=\sum_{i=0}^{K}x(i)y(K-i), K=N-1\f]

                \param Src1 Source 1, x()
                \param Src2 Source 2, y()
                \param Count Vector length, N
                \return Convolution result, z
            */
            static float Convolve (const float *, const float *, long);
            /// \overload
            static double Convolve (const double *, const double *, long);
            /**
                Convolve two vectors (same length, finite) with lag
                \note Circular convolution, result is not scaled

                \f[z(i)=\sum_{k=0}^{N-1}x(i-k)y(k)\f]

                \param Dest Destination, z()
                \param Src1 Source 1, x()
                \param Src2 Source 2, y()
                \param Count Vector length, N
            */
            static void Convolve (float *, const float *, const float *, 
                long);
            /// \overload
            static void Convolve (double *, const double *, const double *, 
                long);
            /**
                Correlate two vectors (same length, finite) without lag
                \note Circular correlation, result is scaled

                \f[z=\frac{1}{N}\sum_{i=0}^{N-1}x(i)y(i)\f]

                \param Src1 Source1, x()
                \param Src2 Source2, y()
                \param Count Vector length, N
                \return Correlation result, z
            */
            static float Correlate (const float *, const float *, long);
            /// \overload
            static double Correlate (const double *, const double *, long);
            /**
                Correlate two vectors (same length, finite) with lag
                \note Circular correlation, result is scaled

                \f[z(i)=\frac{1}{N}\sum_{k=0}^{N-1}x(k)y(k+i)\f]

                \param Dest Destination, z()
                \param Src1 Source 1, x()
                \param Src2 Source 2, y()
                \param Count Vector length, N
            */
            static void Correlate (float *, const float *, const float *, 
                long);
            /// \overload
            static void Correlate (double *, const double *, const double *, 
                long);
            /**
                Autocorrelate vector
                \note Circular, result is scaled

                \f[z=\frac{1}{N}\sum_{i=0}^{N-1}x(i)^2\f]

                \param Src Source, x()
                \param Count Vector length, N
                \return Autocorrelation (energy)
            */
            static float AutoCorrelate (const float *, long);
            /// \overload
            static double AutoCorrelate (const double *, long);
            /**
                Autocorrelate vector
                \note Circular, result is scaled

                \f[y(i)=\frac{1}{N}\sum_{k=0}^{N-1}x(k)x(k+i)\f]

                \param Dest Destination, y()
                \param Src Source, x()
                \param Count Vector length, N
            */
            static void AutoCorrelate (float *, const float *, long);
            /// \overload
            static void AutoCorrelate (double *, const double *, long);
            /**
                Dot product of two vectors

                \f[z=\sum_{i=0}^{N-1}x(i)y(i)\f]

                \param Src1 Source 1, x()
                \param Src2 Source 2, y()
                \param Count Vector length, N
                \return Dot product
            */
            static float DotProduct (const float *, const float *, long);
            /// \overload
            static double DotProduct (const double *, const double *, long);
            /**
                Find minimum and maximum values of vector

                \param Min Minimum
                \param Max Maximum
                \param Src Source
                \param Count Vector length
            */
            static void MinMax (float *, float *, const float *, long);
            /// \overload
            static void MinMax (double *, double *, const double *, long);
            /**
                Calculate mean of vector

                \f[a=\frac{1}{N}\sum_{i=0}^{N-1}x(i)\f]

                \param Src Source
                \param Count Vector length
                \return Mean
            */
            static float Mean (const float *, long);
            /// \overload
            static double Mean (const double *, long);
            /**
                Calculate median of vector

                \f[m=\left\{\begin{array}{ll}x_{sorted}\left(\frac{N}{2}\right)&,N\bmod{2}\neq{0}\\
                0.5\times\left(x_{sorted}\left(\frac{N}{2}-1\right)+x_{sorted}\left(\frac{N}{2}\right)\right)&,N\bmod{2}=0
                \end{array}\right.\f]

                \param Src Source, x()
                \param Count Vector length, N
                \return Median, m
            */
            static float Median (const float *, long);
            /// \overload
            static double Median (const double *, long);
            /**
                Negate vector elements

                \f[x(i)=-x(i), 0\leq{i}\leq{N-1}\f]

                \param Vect Source & destination
                \param Count Vector length
            */
            static void Negate (float *, long);
            /// \overload
            static void Negate (double *, long);
            /**
                Negate vector elements

                \f[y(i)=-x(i), 0\leq{i}\leq{N-1}\f]

                \param Dest Destination, y()
                \param Src Source, x()
                \param Count Vector length, N
            */
            static void Negate (float *, const float *, long);
            /// \overload
            static void Negate (double *, const double *, long);
            /**
                Normalize vector elements

                \f[x(i)=\frac{x(i)-\mu}{\sigma}, 0\leq{i}\leq{N-1}\f]

                \param Vect Source & destination, x()
                \param Count Vector length, N
            */
            static void Normalize (float *, long);
            /// \overload
            static void Normalize (double *, long);
            /**
                Normalize vector elements

                \f[y(i)=\frac{x(i)-\mu}{\sigma}, 0\leq{i}\leq{N-1}\f]

                \param Dest Destination, y()
                \param Src Source, x()
                \param Count Vector length, N
            */
            static void Normalize (float *, const float *, long);
            /// \overload
            static void Normalize (double *, const double *, long);
            /**
                Product of vector elements

                \f[p=\prod{x(i)}, 0\leq{i}\leq{N-1}\f]

                \param Src Source, x()
                \param Count Vector length, N
                \return Product, p
            */
            static float Product (const float *, long);
            /// \overload
            static double Product (const double *, long);
            /**
                Reverse vector (in-place)
                
                \f[y(i)=x(K-i), K=N-1, 0\leq{i}\leq{N-1}\f]
                
                \param Vect Source & destination x()
                \param Count Vector length, N
            */
            static void Reverse (float *, long);
            /// \overload
            static void Reverse (double *, long);
            /// \overload
            static void Reverse (stpSCplx, long);
            /// \overload
            static void Reverse (stpDCplx, long);
            /**
                Reverse vector (out-of-place)

                \f[y(i)=x(K-i), K=N-1, 0\leq{i}\leq{N-1}\f]

                \param Dest Destination, y()
                \param Src Source, x()
                \param Count Vector length
            */
            static void Reverse (float *, const float *, long);
            /// \overload
            static void Reverse (double *, const double *, long);
            /// \overload
            static void Reverse (stpSCplx, const stpSCplx, long);
            /// \overload
            static void Reverse (stpDCplx, const stpDCplx, long);
            /**
                Scale (normalize) vector to range [-1:1]

                \f[x(i)=x(i)\frac{2}{\max-\min}+\left(1-\max\frac{2}{\max-\min}\right), 0\leq{i}\leq{N-1}\f]

                \param Vect Source & destination, x()
                \param Count Vector length, N
            */
            static void Scale (float *, long);
            /// \overload
            static void Scale (double *, long);
            /**
                \overload
                \param Dest Destination
                \param Src Source
                \param Count Vector length
            */
            static void Scale (float *, const float *, long);
            /// \overload
            static void Scale (double *, const double *, long);
            /**
                Scale (normalize) vector to range [0:1]

                \f[x(i)=x(i)\frac{1}{\max-\min}+\left(-\min\frac{1}{\max-\min}\right), 0\leq{i}\leq{N-1}\f]

                \param Vect Source & destination, x()
                \param Count Vector length, N
            */
            static void Scale01 (float *, long);
            /// \overload
            static void Scale01 (double *, long);
            /**
                \overload
                \param Dest Destination
                \param Src Source
                \param Count Vector length
            */
            static void Scale01 (float *, const float *, long);
            /// \overload
            static void Scale01 (double *, const double *, long);
            /**
                Sort vector elements (in-place).
                Vector is sorted using quick-sort algorithm.

                \param Vect Source & destination
                \param Count Vector length
            */                
            static void Sort (float *, long);
            /// \overload
            static void Sort (double *, long);
            /// \overload
            static void Sort (long *, long);
            /**
                Calculate standard deviation and mean of vector

                \f[\mu=\frac{1}{N}\sum_{i=0}^{N-1}x(i)\f]
                \f[\sigma=\sqrt{\frac{1}{N}\sum_{i=0}^{N-1}(x(i)-\mu)^2}\f]

                \param StdDev Standard deviation, \f$\sigma\f$
                \param Mean Mean, \f$\mu\f$
                \param Src Source, x()
                \param Count Vector length, N
            */
            static void StdDev (float *, float *, const float *, long);
            /// \overload
            static void StdDev (double *, double *, const double *, long);
            /**
                Calculate sum of vector elements

                \f[y=\sum_{i=0}^{N-1}x(i)\f]

                \param Src Source, x()
                \param Count Vector length, N
                \return Sum, y
            */
            static float Sum (const float *, long);
            /// \overload
            static double Sum (const double *, long);
            /**
                Square vector in-place

                \f[x(i)=x(i)^2, 0\leq{i}\leq{N-1}\f]

                \param Vect Source & destination, x()
                \param Count Vector length, N
            */
            static void Square (float *, long);
            /// \overload
            static void Square (double *, long);
            /**
                Square vector out-of-place

                \f[y(i)=x(i)^2, 0\leq{i}\leq{N-1}\f]

                \param Dest Destination, y()
                \param Src Source, x()
                \param Count Vector length, N
            */
            static void Square (float *, const float *, long);
            /// \overload
            static void Square (double *, const double *, long);
            /**
                Vector datatype conversion

                \note 12-bit data in LSBs, 24-bit data in MSBs

                \param Dest Destination
                \param Src Source
                \param Count Vector length
            */
            static void Convert (float *, const unsigned char *, long);
            /**
                \overload
                \param Dest Destination
                \param Src Source
                \param Count Vector length
                \param 12bit 12-bit data
            */
            static void Convert (float *, const signed short *, long, bool);
            /**
                \overload
                \param Dest Destination
                \param Src Source
                \param Count Vector length
                \param 24bit 24-bit data
            */
            static void Convert (float *, const signed int *, long, bool);
            /// \overload
            static void Convert (float *, const double *, long);
            /// \overload
            static void Convert (double *, const unsigned char *, long);
            /// \overload
            static void Convert (double *, const signed short *, long, bool);
            /// \overload
            static void Convert (double *, const signed int *, long, bool);
            /// \overload
            static void Convert (double *, const float *, long);
            /// \overload
            static void Convert (unsigned char *, const float *, long);
            /// \overload
            static void Convert (unsigned char *, const double *, long);
            /// \overload
            static void Convert (signed short *, const float *, long, bool);
            /// \overload
            static void Convert (signed short *, const double *, long, bool);
            /// \overload
            static void Convert (signed int *, const float *, long, bool);
            /// \overload
            static void Convert (signed int *, const double *, long, bool);
            /**
                Convert cartesian to polar vectors (out-of-place).
                See ::Cart2Polar for formula

                \param Magn Magnitude
                \param Phase Phase
                \param Real Real
                \param Imag Imaginary
                \param Count Vector length
            */
            static void CartToPolar (float *, float *, const float *, 
                const float *, long);
            /// \overload
            static void CartToPolar (double *, double *, const double *, 
                const double *, long);
            /**
                \overload
                \param Magn Magnitude
                \param Phase Phase
                \param Cart Cartesian
                \param Count Vector length
            */
            static void CartToPolar (float *, float *, const stpSCplx, long);
            /// \overload
            static void CartToPolar (double *, double *, const stpDCplx, long);
            /**
                \overload
                \param Polar Polart
                \param Cart Cartesian
                \param Count Vector length
            */
            static void CartToPolar (stpSPolar, const stpSCplx, long);
            /// \overload
            static void CartToPolar (stpDPolar, const stpDCplx, long);
            /**
                Convert cartesian to polar vectors (in-place)
            */
            static void CartToPolar (utpSCoord, long);
            /// \overload
            static void CartToPolar (utpDCoord, long);
            /**
                Convert polar to cartesian vectors (out-of-place).
                See ::Polar2Cart for formula

                \param Real Real
                \param Imag Imaginary
                \param Magn Magnitude
                \param Phase Phase
                \param Count Vector length
            */
            static void PolarToCart (float *, float *, const float *,
                const float *, long);
            /// \overload
            static void PolarToCart (double *, double *, const double *,
                const double *, long);
            /**
                \overload
                \param Cart Cartesian
                \param Magn Magnitude
                \param Phase Phase
                \param Count Vector length
            */
            static void PolarToCart (stpSCplx, const float *, const float *, 
                long);
            /// \overload
            static void PolarToCart (stpDCplx, const double *, const double *, 
                long);
            /**
                \overload
                \param Polar Polar
                \param Cart Cartesian
                \param Count Vector length
            */
            static void PolarToCart (stpSCplx, const stpSPolar, long);
            /// \overload
            static void PolarToCart (stpDCplx, const stpDPolar, long);
            /**
                Convert cartesian to polar vectors (in-place)
            */
            static void PolarToCart (utpSCoord, long);
            /// \overload
            static void PolarToCart (utpDCoord, long);
            /**
                Calculate normalized cross correlation of two vectors without
                delay

                \f[z=\frac{\frac{1}{N}\sum_{i=0}^{N-1}x(i)y(i)}
                {\frac{1}{N}\sqrt{\sum_{i=0}^{N-1}x(i)^2\sum_{i=0}^{N-1}y(i)^2}}\f]

                \param Src1 Source 1, x()
                \param Src2 Source 2, y()
                \param Count Vector length, N
                \return Normalized cross-correlation
            */
            static float CrossCorr (const float *, const float *, long);
            /// \overload
            static double CrossCorr (const double *, const double *, long);
            /**
                Calculate normalized cross correlation of two vectors with
                delay in second

                \f[z=\frac{\frac{1}{K}\sum_{i=0}^{K-1}x(i)y(i+k)}
                {\frac{1}{K}\sqrt{\sum_{i=0}^{K-1}x(i)^2\sum_{i=0}^{K-1}y(i+k)^2}}
                ,K=N-k\f]

                \param Src1 Source 1, x()
                \param Src2 Source 2, y()
                \param Delay Delay, k
                \param Count Vector length, N
                \return Normalized cross-correlation
            */
            static float DelCrossCorr (const float *, const float *, long, 
                long);
            /// \overload
            static double DelCrossCorr (const double *, const double *, long, 
                long);
            /**
                Calculate normalized cross correlation of two vectors
                with delays specified in vector 

                \f[z(j)=\frac{\frac{1}{K}\sum_{i=0}^{K-1}x(i)y(i+k(j))}
                {\frac{1}{K}\sqrt{\sum_{i=0}^{K-1}x(i)^2\sum_{i=0}^{K-1}y(i+k(j))^2}}
                ,K=N-k(j),0\leq{j}\leq{M-1}\f]

                \param Dest Destination vector, z()
                \param Src1 Source 1, x()
                \param Src2 Source 2, y()
                \param Count Vector length, N
                \param Delay Delays, k()
                \param DelayCount Delay vector length, M
            */
            static void DelCrossCorr (float *, const float *, const float *, 
                long, const long *, long);
            /// \overload
            static void DelCrossCorr (double *, const double *, const double *, 
                long, const long *, long);
            /**
                Calculates energy (square) of vector

                \f[w=\sum_{i=0}^{N-1}x(i)^2\f]

                \param Src Source, x()
                \param Count Vector length, N
                \return Energy, w
            */
            static float Energy (const float *, long);
            /// \overload
            static double Energy (const double *, long);
            /**
                Calculates magnitudes (linear) of vector

                \f[V(i)=\sqrt{\Re_{X}(i)^2+\Im_{X}(i)^2}, 0\leq{i}\leq{N-1}\f]

                \param Dest Destination, V()
                \param Src Source, X()
                \param Count Vector length
            */
            static void Magnitude (float *, const stpSCplx, long);
            /// \overload
            static void Magnitude (double *, const stpDCplx, long);
            /**
                Calculates powers (in dB) of vector

                \f[P(i)=20\log\sqrt{\Re_{X}(i)^2+\Im_{X}(i)^2}, 0\leq{i}\leq{N-1}\f]

                \param Dest Destination, P()
                \param Src Source, X()
                \param Count Vector length
            */
            static void Power (float *, const stpSCplx, long);
            /// \overload
            static void Power (double *, const stpDCplx, long);
            /**
                Calculated phases (in rad) of vector

                \f[\varphi(i)=\arctan{\frac{\Im_{X}(i)}{\Re_{X}(i)}}, 0\leq{i}\leq{N-1}\f]

                \param Dest Destination, \f$\varphi(i)\f$
                \param Src Source, X()
                \param Count Vector length, N
            */
            static void Phase (float *, const stpSCplx, long);
            /// \overload
            static void Phase (double *, const stpDCplx, long);
            /**
                Calculates powers (in dB) and phases (in rad) of vector
                See ::Power() and ::Phase() for formulas

                \param Power Powers
                \param Phase Phases
                \param Src Source
                \param Count Vector length
            */
            static void PowerPhase (float *, float *, const stpSCplx, long);
            /// \overload
            static void PowerPhase (double *, double *, const stpDCplx, long);
            /**
                Decimate vector without average
                \note This can be used in-place also

                \param Dest Destination
                \param Src Source
                \param Factor Decimation factor
                \param Count Vector length (source)
            */
            static void Decimate (float *, const float *, long, long);
            /// \overload
            static void Decimate (double *, const double *, long, long);
            /**
                Decimate vector with average.
                Linear (arithmetic) mean is used to evaluate new values.

                \param Dest Destination
                \param Src Source
                \param Factor Decimation factor
                \param Count Vector length (source)
            */
            static void DecimateAvg (float *, const float *, long, long);
            /// \overload
            static void DecimateAvg (double *, const double *, long, long);
            /**
                Iterpolate vector without average (nulling)

                \param Dest Destination
                \param Src Source
                \param Factor Interpolation factor
                \param Count Vector length (source)
            */
            static void Interpolate (float *, const float *, long, long);
            /// \overload
            static void Interpolate (double *, const double *, long, long);
            /**
                Interpolate vector with average (linear).
                Lagrange interpolation is used to estimate new values.

                \f[y(k)=\frac{k-(i+1)}{i-(i+1)}\times{x(i)}+\frac{k-i}{(i+1)-i}\times{x(i+1)}, 0\leq{k}\leq{M-1}, 0\leq{i}\leq{N-1}\f]

                Where M is destination length and N is source length.

                \param Dest Destination
                \param Src Source
                \param Factor Interpolation factor
                \param Count Vector length (source)
            */
            static void InterpolateAvg (float *, const float *, long, long);
            static void InterpolateAvg (double *, const double *, long, long);
            /**
                Resample vector to different length.

                \param Dest Destination
                \param DestCount Destionation vector length
                \param Src Source
                \param SrcCount Source vector length
            */
            static void Resample (float *, long, const float *, long);
            static void Resample (double *, long, const double *, long);
            /**
                Resample vector to different length with average.

                \param Dest Destination
                \param DestCount Destination vector length
                \param Src Source
                \param SrcCount Source vector length
            */
            static void ResampleAvg (float *, long, const float *, long);
            static void ResampleAvg (double *, long, const double *, long);
            /**
                Calculate RMS (root mean square) of vector

                \f[I=\sqrt{\frac{\sum_{i=0}^{N-1}x(i)^2}{N}}\f]

                \param Src Source, x()
                \param Count Vector length, N
                \return RMS, I
            */
            static float RMS (const float *, long);
            /// \overload
            static double RMS (const double *, long);
            /**
                Calculate variance and mean of vector
                \note variance and mean pointers can be NULL

                \f[\mu=\frac{1}{N}\sum_{i=0}^{N-1}x(i)\f]
                \f[\sigma^2=\frac{\sum_{i=0}^{N}(x(i)-\mu)^2}{N}\f]

                \param Variance Variance, \f$\sigma^2\f$
                \param Mean Mean, \f$\mu\f$
                \param Src Source, x()
                \param Count Vector length, N
                \return Variance \f$\sigma^2\f$
            */
            static float Variance (float *, float *, const float *, long);
            /// \overload
            static double Variance (double *, double *, const double *, long);
            /**
                Find peak level of vector and return result in dB

                \param Src Source
                \param Count Vector length
                \return Peak level (dB)
            */
            static float PeakLevel (const float *, long);
            /// \overload
            static double PeakLevel (const double *, long);
            /**
                Calculate Bartlett (triangle) window coefficients

                \param Dest Destination
                \param Count Vector length
            */
            void WinBartlett (float *, long);
            /// \overload
            void WinBartlett (double *, long);
            /**
                Calculate Blackman window coefficients

                \f[y(n)=0.42-0.5\cos\left(\frac{2\pi{n}}{N}\right)+0.08\cos\left(\frac{4\pi{n}}{N}\right), 0\leq{n}\leq{N-1}\f]

                \param Dest Destination, y()
                \param Count Vector length, N
            */
            // (destination vector, length)
            void WinBlackman (float *, long);
            /// \overload
            void WinBlackman (double *, long);
            /**
                Calculate blackman window coefficients with specified alpha
                
                \f[y(n)=\frac{\alpha+1}{2}-0.5\cos\left(\frac{2\pi{n}}{N-1}\right)-\frac{\alpha}{2}\cos\left(\frac{4\pi{n}}{N-1}\right), 0\leq{n}\leq{N-1}\f]
                
                or with optimal alpha when alpha is specified as 0
                
                \f[\alpha=\frac{0.5}{1+\cos\left(\frac{2\pi}{N-1}\right)}\f]
                
                (around -0.25 for large windows)
            */
            void WinBlackman (float *, long, float);
            void WinBlackman (double *, long, double);
            /**
                Calculate Blackman-Harris window coefficients

                \f[y(n)=0.42323-0.49855\cos\left(\frac{2\pi{n}}{N}\right)+0.07922\cos\left(\frac{4\pi{n}}{N}\right), 0\leq{n}\leq{N-1}\f]

                \param Dest Destination, y()
                \param Count Vector length, N
            */
            void WinBlackmanHarris (float *, long);
            /// \overload
            void WinBlackmanHarris (double *, long);
            /**
                Calculate cosine tapered window coefficients

                \f[y(n)=0.5\left(1-\cos\left(\frac{2\pi{n}}{N}\right)\right), 0\leq{n}\leq{N-1}\f]

                \param Dest Destination, y()
                \param Count Vector length, N
            */
            void WinCosTapered (float *, long);
            /// \overload
            void WinCosTapered (double *, long);
            /**
                \overload
                Apply cosine tapered window to vector (in-place)
            */
            void WinCosTaperedA (float *, long);
            /// \overload
            void WinCosTaperedA (double *, long);
            /**
                \overload
                Apply cosine tapered window to vector (out-of-place)
            */
            void WinCosTaperedA (float *, const float *, long);
            /// \overload
            void WinCosTaperedA (double *, const double *, long);
            /**
                Calculate exact Blackman window coefficients

                \f[y(n)=\frac{7938}{18608}-\frac{9240}{18608}\cos\left(\frac{2\pi{n}}{N}\right)+\frac{1430}{18608}\cos\left(\frac{4\pi{n}}{N}\right), 0\leq{n}\leq{N-1}\f]

                \param Dest Destination, y()
                \param Count Vector length, N
            */
            void WinExactBlackman (float *, long);
            /// \overload
            void WinExactBlackman (double *, long);
            /**
                Calculate exponential window coefficients

                \f[y(n)=y(N-1-n)=\exp\left(\frac{\ln(z+1)}{N/2}n\right)-1, 0\leq{n}\leq{N/2}\f]

                \param Dest Destination, y()
                \param Final Final value, z
                \param Count Vector length, N
            */
            void WinExp (float *, float, long);
            /// \overload
            void WinExp (double *, double, long);
            /**
                Calculate flat top window coefficients

                \f[y(n)=0.2810639-0.5208972\cos\left(\frac{2\pi{n}}{N}\right)+0.1980399\cos\left(\frac{4\pi{n}}{N}\right), 0\leq{n}\leq{N-1}\f]

                \param Dest Destination, y()
                \param Count Vector length, N
            */
            void WinFlatTop (float *, long);
            /// \overload
            void WinFlatTop (double *, long);
            /**
                Calculate generic cosine window coefficients

                \f[y(n)=\sum_{i=0}^{M-1}-1^{i}x(i)\cos\left(\frac{2\pi{in}}{N}\right), 0\leq{n}\leq{N-1}\f]

                \param Dest Destination, y()
                \param Count Vector length, N
                \param Coeff Coefficients, x()
                \param CoeffCount Coefficient vector length, M
            */
            void WinGenericCos (float *, long, const float *, long);
            /// \overload
            void WinGenericCos (double *, long, const double *, long);
            /**
                Calculate Hamming window coefficients

                \f[y(n)=0.54-0.46\cos\left(\frac{2\pi{n}}{N}\right), 0\leq{n}\leq{N-1}\f]

                \param Dest Destination vector, y()
                \param Count Vector length, N
            */
            void WinHamming (float *, long);
            /// \overload
            void WinHamming (double *, long);
            /**
                Calculate Hanning (Hann?) window coefficients

                \f[y(n)=0.5-0.5\cos\left(\frac{2\pi{n}}{N}\right), 0\leq{n}\leq{N-1}\f]

                \param Dest Destination vector, y()
                \param Count Vector length, N
            */
            void WinHanning (float *, long);
            /// \overload
            void WinHanning (double *, long);
            /**
                Calculate Kaiser window coefficients
                See ::ModZeroBessel() for \f$I_{0}\f$

                \f[y(n)=\frac{I_{0}\left\{\beta\sqrt{1-\left|1-\frac{2n}{N}\right|^2}\right\}}
                {I_{0}\left\{\beta\right\}}, 0\leq{n}\leq{N-1}\f]

                \param Dest Destination, y()
                \param Beta Beta, \f$\beta\f$
                \param Count Vector length, N
            */
            void WinKaiser (float *, float, long);
            /// \overload
            void WinKaiser (double *, double, long);
            /**
                Calculate Kaiser-Bessel window coefficients
                See ::ModZeroBessel() for \f$I_{0}\f$

                \f[w(n_{DFT})=I_{0}\left\{\pi\alpha\left[1.0-\left(\frac{n_{DFT}-N/2}{N/2}\right)^{2}\right]^{1/2}\right\}/I_{0}(\pi\alpha), 0\leq{n_{DFT}}\leq{N-1}\f]
                
                \param Dest Destination, w()
                \param Alpha Alpha, \f$\alpha\f$
                \param Count Vector length, N
            */
            void WinKaiserBessel (float *, float, long);
            /// \overload
            void WinKaiserBessel (double *, double, long);
            /**
                Calculate Tukey window coefficients

                \f[y(n)=0.5\left[1+\cos\left(\frac{(n-N/2)\pi}{N/2}\right)\right], 0\leq{n}\leq{N-1}\f]

                \param Dest Destination, y()
                \param Count Vector length, N
            */
            void WinTukey (float *, long);
            /// \overload
            void WinTukey (double *, long);
            /**
                Calculate Dolph-Chebyshev window coefficients
                
                \f[w(n)=\frac{1}{N}\left[1+2r\sum_{m=1}^{M}T_{2M}\left(x_{0}\cos{\frac{\theta_{m}}{2}}\right)\cos{m\theta_{n}}\right], {-M}\leq{n}\leq{M}\f]
                \f[x_{0}=\cosh\left(\frac{1}{2M}\cosh^{-1}\frac{1}{r}\right)\f]

                Gamma (r) is sidelobe / mainlobe ratio.
                
                \param Dest Destination w()
                \param Gamma Gamma (r)
                \param Count Vector length, N
            */
            void WinDolphChebyshev (float *, float, long);
            void WinDolphChebyshev (double *, double, long);
            /**
                Mix two channels interleaved in vector

                \param Dest Destination
                \param Src Source
                \param Count Vector length (destination)
            */
            static void Mix (float *, const float *, long);
            /// \overload
            static void Mix (double *, const double *, long);
            /**
                \overload
                Mix two channels in separate vectors
            */
            static void Mix (float *, const float *, const float *, long);
            /// \overload
            static void Mix (double *, const double *, const double *, long);
            /**
                \overload
                Mix n channels interleaved in vector
            */
            static void Mix (float *, const float *, long, long);
            /// \overload
            static void Mix (double *, const double *, long, long);
            /**
                Spatialize one channel to two channels by inversing others 
                phase. For playing mono sound on stereo headphones.

                \param Dest1 Destination 1
                \param Dest2 Destination 2
                \param Src Source
                \param Count Vector length
            */
            void Spatialize (float *, float *, const float *, long);
            /// \overload
            void Spatialize (double *, double *, const double *, long);
            /// \overload
            void Spatialize (float *, const float *, long);
            /// \overload
            void Spatialize (double *, const double *, long);
            /**
                Extract channel n from N channel interleaved input

                \param Dest Destination
                \param Src Source
                \param Channel Channel to extract
                \param ChCount Channel count
                \param Count Vector length (source)
            */
            static void Extract (float *, const float *, long, long, long);
            /// \overload
            static void Extract (double *, const double *, long, long, long);
            /**
                Pack channel n to N channel interleaved output

                \param Dest Destination
                \param Src Source
                \param Channel Channel to pack
                \param ChCount Channel count
                \param Count Vector length (source)
            */
            static void Pack (float *, const float *, long, long, long);
            /// \overload
            static void Pack (double *, const double *, long, long, long);
            /**
                Rebuffer to different sized buffer

                Returns 0 if there was no complete result block,
                recall with new source buffer and same destination buffer.

                Returns 1 if there was complete result block,
                recall with new source and new destination buffer next time.

                Returns 2 if source buffer is not yet empty, but one result
                buffer available, recall with same source buffer and new
                destination buffer next time.
                
                \param Dest Destination
                \param Src Source
                \param DestCount Destination vector length
                \param SrcCount Source vector length
                \return Rebuffering result
            */
            long ReBuffer (float *, const float *, long, long);
            /// \overload
            long ReBuffer (double *, const double *, long, long);
            /**
                Convert degrees to radians

                \f[y=\frac{\pi}{180}x\f]

                \param Source Degrees, x
                \return Radians, y
            */
            double DegToRad (double dSource)
                { return ((dPI / 180.0) * dSource); }
            /// \overload
            float DegToRad (float fSource)
                { return ((fPI / 180.0F) * fSource); }
            /**
                Convert radians to degrees

                \f[y=\frac{180}{\pi}x\f]

                \param Source Radians, x
                \param Degrees, y
            */
            double RadToDeg (double dSource)
                { return ((180.0 / dPI) * dSource); }
            /// \overload
            float RadToDeg (float fSource)
                { return ((180.0F / fPI) * fSource); }
            /**
                Convert Real-FFTW complex output to our Cplx
                \note Destination length is source length / 2 + 1
                \note This is only for FFTW < 3.x
            */
            static void FFTWConvert (stpSCplx, const float *, long);
            /// \overload
            static void FFTWConvert (stpDCplx, const float *, long);
            /// \overload
            static void FFTWConvert (stpSCplx, const double *, long);
            /// \overload
            static void FFTWConvert (stpDCplx, const double *, long);
            /**
                Convert our Cplx to Real-FFTW complex input
            */
            static void FFTWConvert (float *, const stpSCplx, long);
            /// \overload
            static void FFTWConvert (float *, const stpDCplx, long);
            /// \overload
            static void FFTWConvert (double *, const stpSCplx, long);
            /// \overload
            static void FFTWConvert (double *, const stpDCplx, long);
            /**
                Allocate FIR filter

                \param Coeff Coefficients
                \param Count Vector length
            */
            void FIRAllocate (const float *, long);
            /// \overload
            void FIRAllocate (const double *, long);
            /**
                Filter using FIR filter

                \param Vect Source & destination
                \param Count Vector length
            */
            void FIRFilter (float *, long);
            /// \overload
            void FIRFilter (double *, long);
            /**
                \overload
                \param Dest Destination
                \param Src Source
                \param Count Vector length
            */
            void FIRFilter (float *, const float *, long);
            /// \overload
            void FIRFilter (double *, const double *, long);
            /**
                Filter using FIR filter (fast version).
                
                \note New data at source[FIRlength], must have
                FIRlength size scratch pad at start of source vector.

                \param Dest Destination
                \param Src Source
                \param Count Vector length
            */
            void FIRFilterF (float *, float *, long);
            /// \overload
            void FIRFilterF (double *, double *, long);
            /**
                Free FIR filter
            */
            void FIRFree ();
            /**
                Initialize IIR filter.

                Data format is:
                [0..2] = b[0..2]
                [3..4] = a[0..1]

                \param Coeffs Coefficient vector of length 5
            */
            void IIRInitialize (const float *);
            /// \overload
            void IIRInitialize (const double *);
            /**
                Filter using IIR filter.

                \f[H(z)=\frac{b_{0}+b_{1}z^{-1}+b_{2}z^{-2}}{1-a_{0}z^{-1}-a_{1}z^{-2}}\f]
                \f[y(n)=\sum_{i=1}^{M}a_{i}y(n-i)+\sum_{j=0}^{N}b_{j}x(n-j)\f]

                \param Vect Source & destination
                \param Count Vector length
            */
            void IIRFilter (float *, long);
            /// \overload
            void IIRFilter (double *, long);
            /**
                \overload
                \param Dest Destination
                \param Src Source
                \param Count Vector length
            */
            void IIRFilter (float *, const float *, long);
            /// \overload
            void IIRFilter (double *, const double *, long);
            /**
                Clear IIR filter feedback chain.
            */
            void IIRClear ();
            /**
                Initialize FFT

                \param Size Transform size
                \param Real Real transform?
            */
            void FFTInitialize (long, bool);
            /**
                Uninitialize FFT
            */
            void FFTUninitialize ();
            /**
                FFT half in-place, source vector is modified.
                \note Output is FFT size / 2 + 1

                \param Dest Destination
                \param Src Source
            */
            void FFTi (stpSCplx, float *);
            /// \overload
            void FFTi (stpDCplx, double *);
            /**
                FFT out-of-place.
                \note Output is FFT size / 2 + 1 for real input

                \param Dest Destination
                \param Src Source
            */
            void FFTo (stpSCplx, const float *);
            /// \overload
            void FFTo (stpDCplx, const double *);
            /// \overload
            void FFTo (stpSCplx, const stpSCplx);
            /// \overload
            void FFTo (stpDCplx, const stpDCplx);
            /**
                IFFT out-of-place.
                \note Input is FFT size / 2 + 1 for real output

                \param Dest Destination
                \param Src Source
            */
            void IFFTo (float *, const stpSCplx);
            /// \overload
            void IFFTo (double *, const stpDCplx);
            /// \overload
            void IFFTo (stpSCplx, const stpSCplx);
            /// \overload
            void IFFTo (stpDCplx, const stpDCplx);
            /**
                DCT initialize
                
                \note Internal state is shared with FFT

                \param Size Transform size
            */
            void DCTInitialize (long);
            /**
                Uninitialize DCT
            */
            void DCTUninitialize ();
            /**
                DCT in-place
                
                \param Vect Source & destination vector
            */
            void DCTi (float *);
            /// \overload
            void DCTi (double *);
            /**
                DCT out-of-place
                
                \param Dest Destination vector
                \param Src Source vector
            */
            void DCTo (float *, const float *);
            /// \overload
            void DCTo (double *, const double *);
            /**
                Inverse DCT in-place
                
                \param Vect Source & destination vector
            */
            void IDCTi (float *);
            /// \overload
            void IDCTi (double *);
            /**
                Inverse DCT out-of-place
                
                \param Dest Destination vector
                \param Src Source vector
            */
            void IDCTo (float *, const float *);
            /// \overload
            void IDCTo (double *, const double *);
    };

#endif

