/*

    Rebuffering class
    Copyright (C) 2001-2003 Jussi Laako

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
#include <cstring>
#include <cstdio>

#include <commoncpp/Exception.hh>

#include "dsp/ReBuffer.hh"


inline void clReBuffer::CheckSize (long lDataCount, long lTypeSize)
{
    long lTotalSize = (lCount + lDataCount) * lTypeSize;
    if (lTotalSize > Buffer.GetSize())
    {
        #ifndef _ISOC9X_SOURCE
        long lTwosExp = (long) ceil(log(lTotalSize) / log(2));
        long lNewBufferSize = (long) pow(2, lTwosExp);
        #else
        long lTwosExp = (long) ceil(log2(lTotalSize));
        long lNewBufferSize = (long) exp2(lTwosExp);
        #endif
        Buffer.Resize(lNewBufferSize);

        long lTailSize = lSize - lGetIndex;
        if (lCount > lTailSize)
        {
            long lHeadSize = lCount - lTailSize;
            long lNewSize = lNewBufferSize / lTypeSize;
            unsigned char *cpBuffer = Buffer;
            unsigned char *cpDest = cpBuffer + lSize * lTypeSize;
            memcpy(cpDest, cpBuffer, lHeadSize * lTypeSize);
            lPutIndex = lSize + lHeadSize;
            if (lPutIndex >= lNewSize)
                lPutIndex -= lNewSize;
        }

        lSize = lNewBufferSize / lTypeSize;
    }
}


void *clReBuffer::Index (const std::type_info &TypeInfo, long lIndex)
{
    if (lIndex >= lSize)
        throw clException("clReBuffer::Index(): lIndex >= lSize");
    long lBufIdx = lGetIndex + lIndex;
    if (lBufIdx >= lSize) lBufIdx -= lSize;
    if (TypeInfo == typeid(float))
    {
        float *fpBuffer = Buffer;
        return (&fpBuffer[lBufIdx]);
    }
    else if (TypeInfo == typeid(double))
    {
        double *dpBuffer = Buffer;
        return (&dpBuffer[lBufIdx]);
    }
    else
        throw clException("clReBuffer::Index(): typeid()");
}


void *clReBuffer::GetPtr (const std::type_info &TypeInfo)
{
    if (lGetIndex == 0 || lCount == 0) 
        return ((void *) Buffer);
    long lTempCount = lCount;
    clAlloc Temp;
    if (TypeInfo == typeid(float))
    {
        Temp.Size(lTempCount * sizeof(float));
        Get((float *) Temp, lTempCount);
        Clear();
        Put((float *) Temp, lTempCount);
    }
    else if (TypeInfo == typeid(double))
    {
        Temp.Size(lTempCount * sizeof(double));
        Get((double *) Temp, lTempCount);
        Clear();
        Put((double *) Temp, lTempCount);
    }
    else
        throw clException("clReBuffer::GetPtr(): typeid()");
    return ((void *) Buffer);
}


clReBuffer::clReBuffer ()
{
    lSize = 0;
    lPutIndex = 0;
    lGetIndex = 0;
    lCount = 0;
}


clReBuffer::~clReBuffer ()
{
}


void clReBuffer::Put (const float *fpSrcData, long lSrcCount)
{
    CheckSize(lSrcCount, sizeof(float));
    if (lPutIndex >= lSize) lPutIndex = 0;
    long lTailSpace = lSize - lPutIndex;
    float *fpBuffer = Buffer;
    if (lSrcCount > lTailSpace)
    {
        long lHeadSize = lSrcCount - lTailSpace;
        memcpy(&fpBuffer[lPutIndex], fpSrcData, lTailSpace * sizeof(float));
        memcpy(fpBuffer, &fpSrcData[lTailSpace], lHeadSize * sizeof(float));
        lPutIndex = lHeadSize;
    }
    else
    {
        memcpy(&fpBuffer[lPutIndex], fpSrcData, lSrcCount * sizeof(float));
        lPutIndex += lSrcCount;
    }
    lCount += lSrcCount;
}


void clReBuffer::Put (const double *dpSrcData, long lSrcCount)
{
    CheckSize(lSrcCount, sizeof(double));
    if (lPutIndex >= lSize) lPutIndex = 0;
    long lTailSpace = lSize - lPutIndex;
    double *dpBuffer = Buffer;
    if (lSrcCount > lTailSpace)
    {
        long lHeadSize = lSrcCount - lTailSpace;
        memcpy(&dpBuffer[lPutIndex], dpSrcData, lTailSpace * sizeof(double));
        memcpy(dpBuffer, &dpSrcData[lTailSpace], lHeadSize * sizeof(double));
        lPutIndex = lHeadSize;
    }
    else
    {
        memcpy(&dpBuffer[lPutIndex], dpSrcData, lSrcCount * sizeof(double));
        lPutIndex += lSrcCount;
    }
    lCount += lSrcCount;
}


bool clReBuffer::Get (float *fpDestData, long lDestCount)
{
    if (lCount < lDestCount)
        return false;
    long lTailSize = lSize - lGetIndex;
    float *fpBuffer = Buffer;
    if (lDestCount > lTailSize)
    {
        long lHeadSize = lDestCount - lTailSize;
        memcpy(fpDestData, &fpBuffer[lGetIndex], lTailSize * sizeof(float));
        memcpy(&fpDestData[lTailSize], fpBuffer, lHeadSize * sizeof(float));
        lGetIndex = lHeadSize;
    }
    else
    {
        memcpy(fpDestData, &fpBuffer[lGetIndex], lDestCount * sizeof(float));
        lGetIndex += lDestCount;
    }
    lCount -= lDestCount;
    return true;
}


bool clReBuffer::Get (double *dpDestData, long lDestCount)
{
    if (lCount < lDestCount)
        return false;
    long lTailSize = lSize - lGetIndex;
    double *dpBuffer = Buffer;
    if (lDestCount > lTailSize)
    {
        long lHeadSize = lDestCount - lTailSize;
        memcpy(dpDestData, &dpBuffer[lGetIndex], lTailSize * sizeof(double));
        memcpy(&dpDestData[lTailSize], dpBuffer, lHeadSize * sizeof(double));
        lGetIndex = lHeadSize;
    }
    else
    {
        memcpy(dpDestData, &dpBuffer[lGetIndex], lDestCount * sizeof(double));
        lGetIndex += lDestCount;
    }
    lCount -= lDestCount;
    return true;
}


clReBuffer & clReBuffer::operator= (const clReBuffer &Src)
{
    lSize = Src.lSize;
    lPutIndex = Src.lPutIndex;
    lGetIndex = Src.lGetIndex;
    lCount = Src.lCount;
    Buffer.Copy(Src.Buffer);
    return (*this);
}
