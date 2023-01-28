/*

    Class for flipping odd frequency bands, header
    Copyright (C) 2003-2005 Jussi Laako

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


#include <dsp/DSPOp.hh>
#include <dsp/ReBuffer.hh>


#ifndef FLIPBAND_HH
    #define FLIPBAND_HH


    class clFlipBand
    {
            bool bInitialized;
            long lBlockSize;
            long lCBlockSize;
            clDSPAlloc Proc;
            clDSPAlloc CProc;
            clDSPOp DSP;
            clReBuffer InBuf;
            clReBuffer OutBuf;
        public:
            clFlipBand ();
            virtual ~clFlipBand ();
            /**
                Initialize band flipping
                
                \param TransformSize Size of transform to use for flipping
                \param NullPtr NULL pointer
            */
            void Initialize (long, const float *);
            /// \overload
            void Initialize (long, const double *);
            /**
                Uninitialize
            */
            void Uninitialize ();
            /**
                Put data into flipping buffer
                
                \param Src Pointer to source data
                \pararm Count Number of samples to put
            */
            void Put (const float *, long);
            /// \overload
            void Put (const double *, long);
            /**
                Get data from flipping buffer
                
                \param Dest Pointer to destination buffer
                \param Count Number of samples to get
                \return Success?
            */
            bool Get (float *, long);
            /// \overload
            bool Get (double *, long);
            /**
                Clear buffer contents
            */
            void Clear ();
    };

#endif
