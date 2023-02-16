/*

    IIR Decimation class, header
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


#ifndef IIRDECIMATOR_HH
    #define IIRDECIMATOR_HH


    /**
        IIR decimation filter class implementation.

        Data is filtered using IIR filter and then decimated by factor 
        2 or 3.
    */
    class clIIRDecimator : public clIIRMultiRate
    {
            clDSPAlloc DecBuf;
            clDSPOp DSP;
            clReBuffer InBuf;
        public:
            clIIRDecimator ();
            virtual ~clIIRDecimator ();
            /**
                Uninitialize decimator.
            */
            void Uninitialize ();
            /**
                Feed data into decimator.

                \param SrcData Source data
                \param SrcCount Source data count
            */
            void Put (const float *, long);
            /// \overload
            void Put (const double *, long);
            /**
                Get data from decimator.

                \param DestData Destination buffer
                \param DestCount Number of samples to fetch
                \return Success
            */
            bool Get (float *, long);
            /// \overload
            bool Get (double *, long);
    };

#endif
