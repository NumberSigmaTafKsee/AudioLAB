/*

    Lockless ring buffer
    Copyright (C) 2004-2005 Jussi Laako

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


#include <cstdio>
#include <cstddef>

#include <commoncpp/Exception.hh>
#include <commoncpp/Compilers.hh>

#include <dsp/DSPOp.hh>


#ifndef RINGBUFFER_HH
    #define RINGBUFFER_HH


    /**
        Exception class for clRingBuffer.
    */
    class clXRingBuffer : public clException
    {
        public:
            clXRingBuffer (const char *cpErrorMsg) : 
                clException(cpErrorMsg)
                { }
    };


    /**
        Class for predefined size lockless ring buffer.
    */
    template <class TRingBuffer_t> class clRingBuffer
    {
            volatile size_t sizeBuf;
            volatile size_t sizeCount;
            volatile size_t sizePutPos;
            volatile size_t sizeGetPos;
            clDSPAlloc Buffer;
        public:
            clRingBuffer ()
                {
                    sizeBuf = 0;
                    sizeCount = 0;
                    sizePutPos = 0;
                    sizeGetPos = 0;
                }
            clRingBuffer (const clRingBuffer<TRingBuffer_t> &CopySrc)
                {
                    sizeBuf = CopySrc.sizeBuf;
                    sizeCount = CopySrc.sizeCount;
                    sizePutPos = CopySrc.sizePutPos;
                    sizeGetPos = CopySrc.sizeGetPos;
                    Buffer = CopySrc.Buffer;
                }
            clRingBuffer (size_t sizeNew)
                { SetSize(sizeNew); }
            virtual ~clRingBuffer ()
                {
                    try
                    {
                        Free();
                    }
                    catch (...)
                    { }
                }
            void SetSize (size_t sizeNew)
                {
                    Buffer.Size(sizeNew * sizeof(TRingBuffer_t));
                    sizeBuf = sizeNew;
                    sizeCount = 0;
                    sizePutPos = 0;
                    sizeGetPos = 0;
                }
            void Free ()
                {
                    sizeBuf = 0;
                    sizeCount = 0;
                    sizePutPos = 0;
                    sizeGetPos = 0;
                    Buffer.Free();
                }
            void Lock ()
                { Buffer.Lock(); }
            void UnLock ()
                { Buffer.UnLock(); }
            size_t Size () const
                {
                    return sizeBuf;
                }
            size_t GetCount () const
                {
                    return sizeCount;
                }
            void Put (const TRingBuffer_t *pSrc, size_t sizeSrc)
                {
                    register size_t sizeHead = 0;
                    register size_t sizeTail;
                    TRingBuffer_t *pBuffer = Buffer;

                    if (unlikely(sizeSrc > sizeBuf))
                        throw clXRingBuffer(
                            "clRingBuffer::Put(): sizeSrc > sizeBuf");
                    sizeTail = sizeSrc;
                    if ((sizePutPos + sizeTail) > sizeBuf)
                    {
                        sizeHead = (sizePutPos + sizeTail) % sizeBuf;
                        sizeTail -= sizeHead;
                    }
                    clDSPOp::Copy(pBuffer + sizePutPos, pSrc, sizeTail);
                    sizePutPos += sizeTail;
                    sizePutPos %= sizeBuf;
                    if (sizeHead > 0)
                    {
                        clDSPOp::Copy(pBuffer + sizePutPos, pSrc + sizeTail,
                            sizeHead);
                        sizePutPos += sizeHead;
                        sizePutPos %= sizeBuf;
                    }
                    sizeCount += sizeSrc;
                    if (sizeCount > sizeBuf)
                        sizeCount = sizeBuf;
                    return true;
                }
            bool Put2 (const TRingBuffer_t *pSrc, size_t sizeSrc)
                {
                    register size_t sizeHead = 0;
                    register size_t sizeTail;
                    TRingBuffer_t *pBuffer = Buffer;

                    if (unlikely(sizeSrc > sizeBuf))
                        throw clXRingBuffer(
                            "clRingBuffer::Put(): sizeSrc > sizeBuf");
                    if (sizeSrc > (sizeBuf - sizeCount))
                        return false;
                    sizeTail = sizeSrc;
                    if ((sizePutPos + sizeTail) > sizeBuf)
                    {
                        sizeHead = (sizePutPos + sizeTail) % sizeBuf;
                        sizeTail -= sizeHead;
                    }
                    clDSPOp::Copy(pBuffer + sizePutPos, pSrc, sizeTail);
                    sizePutPos += sizeTail;
                    sizePutPos %= sizeBuf;
                    if (sizeHead > 0)
                    {
                        clDSPOp::Copy(pBuffer + sizePutPos, pSrc + sizeTail,
                            sizeHead);
                        sizePutPos += sizeHead;
                        sizePutPos %= sizeBuf;
                    }
                    sizeCount += sizeSrc;
                    return true;
                }
            bool Get (TRingBuffer_t *pDst, size_t sizeDst)
                {
                    register size_t sizeHead = 0;
                    register size_t sizeTail;
                    TRingBuffer_t *pBuffer = Buffer;

                    if (unlikely(sizeDst > sizeBuf))
                        throw clXRingBuffer("clRingBuffer::Get(): sizeDst > sizeBuf");
                    if (sizeDst > sizeCount)
                        return false;
                    sizeTail = sizeDst;
                    if ((sizeGetPos + sizeTail) > sizeBuf)
                    {
                        sizeHead = (sizeGetPos + sizeTail) % sizeBuf;
                        sizeTail -= sizeHead;
                    }
                    clDSPOp::Copy(pDst, pBuffer + sizeGetPos, sizeTail);
                    sizeGetPos += sizeTail;
                    sizeGetPos %= sizeBuf;
                    if (sizeHead > 0)
                    {
                        clDSPOp::Copy(pDst + sizeTail, pBuffer + sizeGetPos,
                            sizeHead);
                        sizeGetPos += sizeHead;
                        sizeGetPos %= sizeBuf;
                    }
                    sizeCount -= sizeDst;
                    return true;
                }
    };

#endif
