/*

    Rebuffering class
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

*/


#include <typeinfo>

#include <commoncpp/Alloc.hh>


#ifndef REBUFFER_HH
    #define REBUFFER_HH


    /**
        Class for splitting data into buffers of different sizes.
        It's implemented using dynamically growing circular FIFO buffer.
    */
    class clReBuffer
    {
            long lSize;
            long lPutIndex;
            long lGetIndex;
            long lCount;
            clAlloc Buffer;
            void CheckSize (long, long);
        protected:
            void *Index (const std::type_info &, long);
            void *GetPtr (const std::type_info &);
        public:
            clReBuffer ();
            clReBuffer (const clReBuffer &CopySrc)
                {
                    *this = CopySrc;
                }
            virtual ~clReBuffer ();
            /**
                Put data into FIFO.

                \param fpSrcData Source buffer
                \param lSrcCount Number of samples in source.
            */
            void Put (const float *, long);
            /// \overload
            void Put (const double *, long);
            /**
                Get data from FIFO.
                Returns false if there's not enough samples in FIFO to
                fill requested buffer.

                \param fpDestData Destination buffer
                \param lDestCount Number of samples to fetch
                \return Data available
            */
            bool Get (float *, long);
            /// \overload
            bool Get (double *, long);
            /**
                Get number of samples in FIFO.

                \return Number of samples
            */
            long GetCount () const
                { 
                    return lCount; 
                }
            /**
                Clear FIFO contents and reset buffer size.
            */
            void Clear ()
                {
                    lSize = 0;
                    lPutIndex = 0;
                    lGetIndex = 0;
                    lCount = 0;
                    Buffer.Free();
                }
            clReBuffer & operator= (const clReBuffer &);
    };


#endif
