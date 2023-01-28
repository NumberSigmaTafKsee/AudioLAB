/*

    Class/template system for threading using pthreads
    Copyright (C) 2002-2005 Jussi Laako

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with this library; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

*/


#include <cstdio>
#include <climits>
#ifndef USE_NPTL
#include <pthread.h>
#else
#include <nptl/pthread.h>
#endif

#include "commoncpp/Threads.hh"


extern "C"
{

void *WrapThreadsBase (void *vpParam)
{
    clThreadsBase::stpParams spParams =
        (clThreadsBase::stpParams) vpParam;
    clThreadsBase *Klass = spParams->Klass;
    void *vpICParam = spParams->vpParam;

    pthread_setcancelstate(PTHREAD_CANCEL_ENABLE, NULL);
    pthread_setcanceltype(PTHREAD_CANCEL_ASYNCHRONOUS, NULL);
    
    delete spParams;
    return Klass->InternalCaller(vpICParam);
}

}


clThreadsBase::clThreadsBase ()
{
    iThreadCount = 0;
}


clThreadsBase::~clThreadsBase ()
{
    try
    {
        MtxBase.Wait();

        ThreadMap_t::iterator iterThreads;
    
        iterThreads = mapThreads.begin();
        while (iterThreads != mapThreads.end())
        {
            pthread_cancel((*iterThreads).second);
            pthread_join((*iterThreads).second, NULL);
            iterThreads++;
        }

        MtxBase.Release();
    }
    catch (...)
    {
    }
}


int clThreadsBase::Create (void *vpParam, bool bDetached)
{
    stpParams spParams;
    pthread_t ptidThread;
    int iThreadHandle = -1;
    
    spParams = new stParams;
    spParams->Klass = this;
    spParams->vpParam = vpParam;

    pthread_create(&ptidThread, NULL, WrapThreadsBase, (void *) spParams);
    if (!bDetached)
    {
        MtxBase.Wait();

        while (mapThreads.find(iThreadCount) != mapThreads.end())
        {
            iThreadCount++;
            if (iThreadCount >= INT_MAX)
                iThreadCount = 0;
        }
        iThreadHandle = iThreadCount;
        mapThreads[iThreadCount++] = ptidThread;
        if (iThreadCount >= INT_MAX)
            iThreadCount = 0;

        MtxBase.Release();
    }
    else
    {
        pthread_detach(ptidThread);
    }
    return iThreadHandle;
}


void *clThreadsBase::Wait (int iThreadHandle)
{
    MtxBase.Wait();

    ThreadMap_t::iterator iterThread;
    void *vpReturn = NULL;

    iterThread = mapThreads.find(iThreadHandle);
    if (iterThread != mapThreads.end())
    {
        MtxBase.Release();

        pthread_join((*iterThread).second, &vpReturn);

        MtxBase.Wait();

        mapThreads.erase(iterThread);
    }

    MtxBase.Release();

    return vpReturn;
}


void clThreadsBase::Cancel (int iThreadHandle)
{
    MtxBase.Wait();

    ThreadMap_t::iterator iterThread;

    iterThread = mapThreads.find(iThreadHandle);
    if (iterThread != mapThreads.end())
    {
        MtxBase.Release();

        pthread_cancel((*iterThread).second);

        MtxBase.Wait();

        mapThreads.erase(iterThread);
    }

    MtxBase.Release();
}


bool clThreadsBase::SetSched (pthread_t ptidThread, int iPolicy, int iDelta)
{
    int iEC;
    uid_t uidCurrent;
    struct sched_param sSchedParam;

    // Serialize uid operations
    MtxBase.Wait();
    uidCurrent = getuid();
    setuid(0);
#   ifndef BSDSYS
    sSchedParam.sched_priority = sched_get_priority_min(iPolicy) + iDelta;
#   else
    sSchedParam.sched_priority = iDelta;
#   endif
    iEC = pthread_setschedparam(pthread_self(), iPolicy, &sSchedParam);
    setuid(uidCurrent);
    MtxBase.Release();
    if (iEC != 0) return false;
    return true;
}
