/*

    Template rebuffering class
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


#include <cmath>
#include <cfloat>
#include <cstring>
#include <typeinfo>

#include <commoncpp/Exception.hh>
#include <commoncpp/Compilers.hh>

#include <dsp/DSPOp.hh>


#ifndef REBUFFERT_HH
    #define REBUFFERT_HH


    /**
        Exception class for clReBufferT.
    */
    class clXReBufferT : public clException
    {
        public:
            clXReBufferT (const char *cpErrorMsg) : 
                clException(cpErrorMsg)
                { }
    };


    /**
        Template class for splitting data into buffers of different sizes.
        It's implemented using dynamically growing circular FIFO buffer.
        The ring buffer is linearized when linear access is requested.
    */
    template <class TReBuffer_t> class clReBufferT
    {
            long lSize;
            long lPutIndex;
            long lGetIndex;
            long lCount;
            clDSPAlloc Buffer;

            inline void CheckSize (long lDataCount)
                {
                    long lTotalSize = (lCount + lDataCount) * 
                        sizeof(TReBuffer_t);
                    if (lTotalSize > Buffer.GetSize())
                    {
                        #ifndef _ISOC9X_SOURCE
                        long lTwosExp = (long) std::ceil(
                            std::log((double) lTotalSize) / 
                            std::log((double) 2));
                        long lNewBufferSize = (long) 
                            std::pow(2.0, (double) lTwosExp);
                        #else
                        long lTwosExp = (long) ceil(log2(lTotalSize));
                        long lNewBufferSize = (long) exp2(lTwosExp);
                        #endif
                        Buffer.Resize(lNewBufferSize);

                        long lTailSize = lSize - lGetIndex;
                        if (lCount > lTailSize)
                        {
                            long lHeadSize = lCount - lTailSize;
                            long lNewSize = 
                                lNewBufferSize / sizeof(TReBuffer_t);
                            unsigned char *cpBuffer = Buffer;
                            unsigned char *cpDest = 
                                cpBuffer + lSize * sizeof(TReBuffer_t);
                            memcpy(cpDest, cpBuffer, 
                                lHeadSize * sizeof(TReBuffer_t));
                            lPutIndex = lSize + lHeadSize;
                            if (lPutIndex >= lNewSize)
                                lPutIndex -= lNewSize;
                        }

                        lSize = lNewBufferSize / sizeof(TReBuffer_t);
                    }
                }

        protected:

            /*long lSize;
            long lPutIndex;
            long lGetIndex;
            long lCount;
            clDSPAlloc Buffer;*/

            TReBuffer_t * GetPtr ()
                {
                    if (lGetIndex == 0 || lCount == 0)
                        return ((TReBuffer_t *) Buffer);
                    long lTempCount = lCount;
                    clDSPAlloc Temp(lTempCount * sizeof(TReBuffer_t));
                    Get((TReBuffer_t *) Temp, lTempCount);
                    Clear();
                    Put((TReBuffer_t *) Temp, lTempCount);
                    return ((TReBuffer_t *) Buffer);
                }
            TReBuffer_t * GetRawPtr ()
                {
                    return ((TReBuffer_t *) Buffer.GetPtr());
                }
            void CopyGet (TReBuffer_t *fpDestData, long lDestCount) const
                {
                    long lTailSize = lSize - lGetIndex;
                    const TReBuffer_t *fpBuffer = 
                        (const TReBuffer_t *) Buffer.GetPtr();
                    if (lDestCount > lTailSize)
                    {
                        long lHeadSize = lDestCount - lTailSize;
                        memcpy(fpDestData, &fpBuffer[lGetIndex],
                            lTailSize * sizeof(TReBuffer_t));
                        memcpy(&fpDestData[lTailSize], fpBuffer,
                            lHeadSize * sizeof(TReBuffer_t));
                    }
                    else
                    {
                        memcpy(fpDestData, &fpBuffer[lGetIndex],
                            lDestCount * sizeof(TReBuffer_t));
                    }
                }

        public:

            clReBufferT ()
                {
                    lSize = 0;
                    lPutIndex = 0;
                    lGetIndex = 0;
                    lCount = 0;
                }
            clReBufferT (const clReBufferT &CopySrc)
                {
                    lSize = 0;
                    lPutIndex = 0;
                    lGetIndex = 0;
                    lCount = 0;
                    *this = CopySrc;
                }
            clReBufferT (long lNewSize)
                {
                    lSize = 0;
                    lPutIndex = 0;
                    lGetIndex = 0;
                    lCount = 0;
                    SetSize(lNewSize);
                }
            clReBufferT (const TReBuffer_t *fpSrcData, long lSrcCount)
                {
                    lSize = 0;
                    lPutIndex = 0;
                    lGetIndex = 0;
                    lCount = 0;
                    Put(fpSrcData, lSrcCount);
                }
            virtual ~clReBufferT ()
                { }
            /**
                Put data into FIFO.

                \param fpSrcData Source buffer
                \param lSrcCount Number of samples in source.
            */
            void Put (const TReBuffer_t *fpSrcData, long lSrcCount)
                {
                    CheckSize(lSrcCount);
                    if (lPutIndex >= lSize) lPutIndex = 0;
                    long lTailSpace = lSize - lPutIndex;
                    TReBuffer_t *fpBuffer = Buffer;
                    if (lSrcCount > lTailSpace)
                    {
                        long lHeadSize = lSrcCount - lTailSpace;
                        memcpy(&fpBuffer[lPutIndex], fpSrcData, 
                            lTailSpace * sizeof(TReBuffer_t));
                        memcpy(fpBuffer, &fpSrcData[lTailSpace], 
                            lHeadSize * sizeof(TReBuffer_t));
                        lPutIndex = lHeadSize;
                    }
                    else
                    {
                        memcpy(&fpBuffer[lPutIndex], fpSrcData, 
                            lSrcCount * sizeof(TReBuffer_t));
                        lPutIndex += lSrcCount;
                    }
                    lCount += lSrcCount;
                }
            /**
                Put data into fifo from another FIFO.

                \param Src Source FIFO
            */
            void Put (clReBufferT &Src)
                {
                    Put(Src.GetPtr(), Src.Size());
                }
            /**
                Get data from FIFO.
                Returns false if there's not enough samples in FIFO to
                fill requested buffer.

                \param fpDestData Destination buffer
                \param lDestCount Number of samples to fetch
                \return Data available
            */
            bool Get (TReBuffer_t *fpDestData, long lDestCount)
                {
                    if (lCount < lDestCount)
                        return false;
                    long lTailSize = lSize - lGetIndex;
                    TReBuffer_t *fpBuffer = Buffer;
                    if (lDestCount > lTailSize)
                    {
                        long lHeadSize = lDestCount - lTailSize;
                        memcpy(fpDestData, &fpBuffer[lGetIndex], 
                            lTailSize * sizeof(TReBuffer_t));
                        memcpy(&fpDestData[lTailSize], fpBuffer, 
                            lHeadSize * sizeof(TReBuffer_t));
                        lGetIndex = lHeadSize;
                    }
                    else
                    {
                        memcpy(fpDestData, &fpBuffer[lGetIndex], 
                            lDestCount * sizeof(TReBuffer_t));
                        lGetIndex += lDestCount;
                    }
                    lCount -= lDestCount;
                    return true;
                }
            /**
                Get number of samples in FIFO.

                \return Number of samples
            */
            long GetCount () const
                { 
                    return lCount; 
                }
            long Size () const
                {
                    return lCount;
                }
            /**
                Set number of samples in FIFO.

                \note This is destructive resize.
            */
            void SetSize (long lNewCount)
                {
                    Clear();
                    CheckSize(lNewCount);
                    lCount = lNewCount;
                    lPutIndex = lCount;
                }
            /**
                Set number of samples in FIFO.

                \note This is non-destructive resize.
            */
            void Resize (long lNewCount)
                {
                    if (lNewCount > lSize)
                        CheckSize(lNewCount - lCount);
                    lCount = lNewCount;
                    lPutIndex = lCount;
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
            clReBufferT & operator= (const clReBufferT &Src)
                {
                    SetSize(Src.Size());
                    Src.CopyGet((TReBuffer_t *) Buffer, Src.Size());
                    lPutIndex = Src.Size();
                    lCount = Src.Size();
                    return (*this);
                }
            TReBuffer_t & operator[] (long lIndex)
                {
                    if (unlikely(lIndex >= lCount))
                    {
                        throw clXReBufferT(
                            "clReBufferT<>::operator[]: Index >= Size()");
                    }
                    long lBufIdx = lGetIndex + lIndex;
                    if (lBufIdx >= lSize) lBufIdx -= lSize;
                    TReBuffer_t *fpBuffer = Buffer;
                    return fpBuffer[lIndex];
                }
    };


#endif
