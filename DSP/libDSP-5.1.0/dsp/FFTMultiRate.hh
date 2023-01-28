/*

    Base class for FFT based multirate filters, header
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


#include <dsp/Filter.hh>


#ifndef FFTMULTIRATE_HH
    #define FFTMULTIRATE_HH

    #define FFTMULTIRATE_RIPPLERATIO        144.49f
    #define FFTMULTIRATE_OVERLAP            0.75f
    #define FFTMULTIRATE_DELTAOMEGA         0.1f


    /**
        Base class for FFT based multirate filters.

        Data is is filtered using FFT filter and then decimated or interpolated
        by derived classes.
    */
    class clFFTMultiRate
    {
        protected:
            bool bInitialized;
            long lFactor;  ///< Rate change factor
            long lFilterSize;  ///< Size of filter FFT
            clFilter Filter;
        public:
            clFFTMultiRate ();
            virtual ~clFFTMultiRate ();
            /**
                Initialize filter.
                You can re-initialize without uninitializing first.

                The NULL pointer is used to select correct overloaded function
                matching input datatype.

                \param FactorP Rate change factor
                \param FiltSize Filter size, specify <= 0 for automatic selection
                \param NullPtr NULL pointer
                \param HighPass High-pass filter?
            */
            bool Initialize (long, long, const float *, bool = false);
            /// \overload
            bool Initialize (long, long, const double *, bool = false);
            /**
                Uninitialize filter.
            */
            void Uninitialize ();
    };

#endif

