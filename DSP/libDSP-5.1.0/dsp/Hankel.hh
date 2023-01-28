/*

    Hankel transform (Abel transform followed by FFT)
    Copyright (C) 1995 Colorado School of Mines. All rights reserved.
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

    Source converted from cwplib (by Dave Hale and Lydia Deng / Colorado
    School of Mines).

    References:
    Hansen, E. W., 1985, Fast Hankel transform algorithm:  IEEE Trans. on
    Acoustics, Speech and Signal Processing, v. ASSP-33, n. 3, p. 666-671.
    (Beware of several errors in the equations in this paper!)

*/


#include <dsp/DSPOp.hh>
//#include <dsp/PFAFFT.hh>


#ifndef HANKEL_HH
    #define HANKEL_HH

    /**
        Class implementation of (modified) Hankel-transform.

        \author Colorado School of Mines
        \author Jussi Laako

        \par
        Source converted from cwplib (by Dave Hale and Lydia Deng / Colorado
        School of Mines).

        \par
        Transform is significantly sped-up by doing Abel-transform followed
        by inverse-FFT.

        \par References
        Hansen, E. W., 1985, Fast Hankel transform algorithm:  IEEE Trans. on
        Acoustics, Speech and Signal Processing, v. ASSP-33, n. 3, p. 666-671.
        (Beware of several errors in the equations in this paper!)
    */
    class clHankel
    {
            long lSize;
            long lFFTSize;
            float fOutScale0;
            float fOutScale1;
            double dOutScale0;
            double dOutScale1;
            clDSPAlloc A;
            clDSPAlloc B0;
            clDSPAlloc B1;
            clDSPAlloc GX;
            clDSPAlloc GK;
            clDSPOp DSP;
            //clPFAFFT FFT;
            /**
                Initialize Abel transformer

                \param NullPtr NULL pointer
            */
            void InitAbel (const float *);
            /// \overload
            void InitAbel (const double *);
            /**
                Do Abel transform

                \param Dest Destination
                \param Src Source
            */
            void DoAbel (float *, const float *);
            /// \overload
            void DoAbel (double *, const double *);
            /**
                Unintialize Abel transformer
            */
            void UninitAbel ();
        public:
            clHankel ();
            virtual ~clHankel ();
            /**
                Initialize Hankel transform

                The NULL pointer is used make C++ compiler select correct
                datatype for the transform from overloaded method set.

                \param Size Window size
                \param NullPtr NULL pointer
            */
            void Initialize (long, const float *);
            /// \overload
            void Initialize (long, const double *);
            /**
                Uninitialize, destructor also does this
            */
            void Uninitialize ();
            /**
                Process data (0th order)

                Definition:
                \f[h_{0}(k)=\int_{0}^{\infty}rj_{0}(kr)f(r)\mathrm{d}r\f]

                \note Input and output is window size / 2 + 1!

                \param Dest Destination
                \param Src Source
            */
            void Process0 (float *, const float *);
            /// \overload
            void Process0 (double *, const double *);
            /**
                Process data (1st order)

                Definition:
                \f[h_{1}(k)=\int_{0}^{\infty}rj_{1}(kr)f(r)\mathrm{d}r\f]

                \note Input and output is window size / 2 + 1!

                \param Dest Destination
                \param Src Source
            */
            void Process1 (float *, const float *);
            /// \overload
            void Process1 (double *, const double *);
    };

#endif
