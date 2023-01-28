/*

    Recursive interpolation class, header
    Copyright (C) 2000-2005 Jussi Laako

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


#include <dsp/DSPConfig.hh>
#include <dsp/DSPOp.hh>
#include <dsp/FFTInterpolator.hh>
#include <dsp/FIRInterpolator.hh>
#include <dsp/IIRInterpolator.hh>


#ifndef RECINTERPOLATOR_HH
    #define RECINTERPOLATOR_HH

    #define RECINT_MAX_SUB_ROUNDS       32  ///< Maximum number of subsequent rounds


    /**
        Recursive interpolation filter class implementation.

        Data is interpolated and filtered in recursive manner in rounds of
        factor 2.
    */
    class clRecInterpolator
    {
            bool bInitialized;
            int iType;  ///< Type filter; 0 = FFT, 1= FIR, 2 = IIR
            long lFactor;  ///< Interpolation factor
            long lFilterSize;  ///< Size of input buffer
            long lSubRounds;  ///< Number of recursive rounds
            long lIntSize;  ///< Size of recursion buffer
            bool bpHalves[RECINT_MAX_SUB_ROUNDS];
            clDSPAlloc IntBuf;
            clFFTInterpolator FFTIntBank[RECINT_MAX_SUB_ROUNDS];
            clFIRInterpolator FIRIntBank[RECINT_MAX_SUB_ROUNDS];
            clIIRInterpolator IIRIntBank[RECINT_MAX_SUB_ROUNDS];
            void InitHalves (double);
        public:
            enum eFilterType
            {
                FILTER_TYPE_FFT = 0,
                FILTER_TYPE_FIR = 1,
                FILTER_TYPE_IIR = 2
            };
            clRecInterpolator ();
            virtual ~clRecInterpolator ();
            /**
                Initialize interpolator, interpolation factor must be
                powers of two.
                You can re-initialize interpolator without uninitializing first.

                The NULL pointer is used select correct overloaded function
                matching input datatype.

                \param IntFact Interpolation factor
                \param FiltSize Filter size (for FFT)
                \param NullPtr NULL pointer
                \param BandCenter Normalized center frequency of filter band (destination)
                \param FilterType Filter type; 0 = FFT, 1 = FIR, 2 = IIR
            */
            bool Initialize (long, long, const float *, float = 0.0f, 
                int = 0);
            /// \overload
            bool Initialize (long, long, const double *, double = 0.0, 
                int = 0);
            /// \overload
            bool Initialize (long, long, const float *, float = 0.0f, 
                bool = false);
            /// \overload
            bool Initialize (long, long, const double *, double = 0.0, 
                bool = false);
            /**
                Uninitialize interpolator.
            */
            void Uninitialize ();
            /**
                Feed data into interpolator.

                \param SrcData Source data
                \param SrcCount Source data count
            */
            void Put (const float *, long);
            /// \overload
            void Put (const double *, long);
            /**
                Get data from interpolator.

                Returns false if there's not enough data feeded into
                the interpolator.

                \param DestData Destination buffer
                \param DestCount Number of samples
                \return Data was available?
            */
            bool Get (float *, long);
            /// \overload
            bool Get (double *, long);
    };

#endif
