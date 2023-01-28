/*

    Cascaded IIR stages
    Copyright (C) 2002 Jussi Laako

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


#include "dsp/IIRCascade.hh"


clIIRCascade::clIIRCascade ()
{
    bInitialized = false;
    lStages = 0;
    IIR = NULL;
}


clIIRCascade::~clIIRCascade ()
{
    if (bInitialized)
        Uninitialize();
}


bool clIIRCascade::Initialize (const float fppCoeffs[][5], long lStageCount)
{
    long lStageCntr;

    if (bInitialized) Uninitialize();
    lStages = lStageCount;
    IIR = new clDSPOp[lStages];
    if (IIR == NULL)
        return false;
    for (lStageCntr = 0; lStageCntr < lStages; lStageCntr++)
        IIR[lStageCntr].IIRInitialize(fppCoeffs[lStageCntr]);
    return true;
}


bool clIIRCascade::Initialize (const double dppCoeffs[][5], long lStageCount)
{
    long lStageCntr;

    if (bInitialized) Uninitialize();
    lStages = lStageCount;
    IIR = new clDSPOp[lStages];
    if (IIR == NULL)
        return false;
    for (lStageCntr = 0; lStageCntr < lStages; lStageCntr++)
        IIR[lStageCntr].IIRInitialize(dppCoeffs[lStageCntr]);
    return true;
}


void clIIRCascade::Uninitialize ()
{
    if (IIR != NULL)
    {
        delete [] IIR;
        IIR = NULL;
    }
    lStages = 0;
    bInitialized = false;
}


void clIIRCascade::Process (float *fpVect, long lCount)
{
    long lStageCntr;
    
    for (lStageCntr = 0; lStageCntr < lStages; lStageCntr++)
        IIR[lStageCntr].IIRFilter(fpVect, lCount);
}


void clIIRCascade::Process (double *dpVect, long lCount)
{
    long lStageCntr;
    
    for (lStageCntr = 0; lStageCntr < lStages; lStageCntr++)
        IIR[lStageCntr].IIRFilter(dpVect, lCount);
}


void clIIRCascade::Process (float *fpDest, const float *fpSrc, long lCount)
{
    long lStageCntr;
    
    IIR[0].IIRFilter(fpDest, fpSrc, lCount);
    for (lStageCntr = 1; lStageCntr < lStages; lStageCntr++)
        IIR[lStageCntr].IIRFilter(fpDest, lCount);
}


void clIIRCascade::Process (double *dpDest, const double *dpSrc, long lCount)
{
    long lStageCntr;
    
    IIR[0].IIRFilter(dpDest, dpSrc, lCount);
    for (lStageCntr = 1; lStageCntr < lStages; lStageCntr++)
        IIR[lStageCntr].IIRFilter(dpDest, lCount);
}


void clIIRCascade::Clear ()
{
    long lStageCntr;
    
    for (lStageCntr = 0; lStageCntr < lStages; lStageCntr++)
        IIR[lStageCntr].IIRClear();
}
