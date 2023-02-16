/*

    IIR Interpolation class, header
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


#include <dsp/DSPConfig.hh>
#include <dsp/DSPOp.hh>
#include <dsp/ReBuffer.hh>
#include <dsp/IIRMultiRate.hh>


#ifndef IIRINTERPOLATOR_HH
    #define IIRINTERPOLATOR_HH


    /**
        IIR interpolation filter class implementation.

        Data is interpolated by factor 2 or 3 and then filtered using IIR
        filter.
    */
    class clIIRInterpolator : public clIIRMultiRate
    {
            clDSPAlloc IntBuf;
            clDSPOp DSP;
            clReBuffer OutBuf;
        public:
            clIIRInterpolator ();
            virtual ~clIIRInterpolator ();
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

                \param DestData Destination buffer
                \param DestCount Number of samples to fetch
                \return Success
            */
            bool Get (float *, long);
            /// \overload
            bool Get (double *, long);
    };

#endif
