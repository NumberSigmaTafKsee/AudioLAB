/*

    Class for flipping odd frequency bands
    Copyright (C) 2003 Jussi Laako

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


#include "dsp/FlipBand.hh"


clFlipBand::clFlipBand ()
{
    bInitialized = false;
    lBlockSize = 0;
}


clFlipBand::~clFlipBand ()
{
    if (bInitialized)
        Uninitialize();
}


void clFlipBand::Initialize (long lTransformSize, const float *fpNullPtr)
{
    lBlockSize = lTransformSize;
    lCBlockSize = lTransformSize / 2 + 1;
    DSP.FFTInitialize(lTransformSize, true);
    Proc.Size(lBlockSize * sizeof(float));
    CProc.Size(lCBlockSize * sizeof(stSCplx));
    bInitialized = true;
}


void clFlipBand::Initialize (long lTransformSize, const double *dpNullPtr)
{
    lBlockSize = lTransformSize;
    lCBlockSize = lTransformSize / 2 + 1;
    DSP.FFTInitialize(lTransformSize, true);
    Proc.Size(lBlockSize * sizeof(double));
    CProc.Size(lCBlockSize * sizeof(stDCplx));
    bInitialized = true;
}


void clFlipBand::Uninitialize ()
{
    if (bInitialized)
    {
        DSP.FFTUninitialize();
        Clear();
        Proc.Free();
        CProc.Free();
        bInitialized = false;
    }
}


void clFlipBand::Put (const float *fpSrc, long lCount)
{
    float *fpProc;
    stpSCplx spCProc;

    fpProc = Proc;
    spCProc = CProc;
    InBuf.Put(fpSrc, lCount);
    while (InBuf.Get(fpProc, lBlockSize))
    {
        DSP.FFTi(spCProc, fpProc);
        DSP.Reverse(spCProc, lCBlockSize);
        DSP.IFFTo(fpProc, spCProc);
        OutBuf.Put(fpProc, lBlockSize);
    }
}


void clFlipBand::Put (const double *dpSrc, long lCount)
{
    double *dpProc;
    stpDCplx spCProc;

    dpProc = Proc;
    spCProc = CProc;
    InBuf.Put(dpSrc, lCount);
    while (InBuf.Get(dpProc, lBlockSize))
    {
        DSP.FFTi(spCProc, dpProc);
        DSP.Reverse(spCProc, lCBlockSize);
        DSP.IFFTo(dpProc, spCProc);
        OutBuf.Put(dpProc, lBlockSize);
    }
}


bool clFlipBand::Get (float *fpDest, long lCount)
{
    return OutBuf.Get(fpDest, lCount);
}


bool clFlipBand::Get (double *dpDest, long lCount)
{
    return OutBuf.Get(dpDest, lCount);
}


void clFlipBand::Clear ()
{
    InBuf.Clear();
    OutBuf.Clear();
}
