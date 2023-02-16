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


#include <unistd.h>
#include <sys/types.h>
#ifndef USE_NPTL
#include <pthread.h>
#else
#include <nptl/pthread.h>
#endif
#include <map>

#include <commoncpp/Mutex.hh>


#ifndef THREADS_HH
#define THREADS_HH


    /// Thread ID container type.
    typedef std::map<int, pthread_t> ThreadMap_t;


    /**
        Base class for threads.
        
        \note Not to be used directly.
    */
    class clThreadsBase
    {
            int iThreadCount;
            ThreadMap_t mapThreads;
            clMutex MtxBase;
        public:
            typedef struct _stParams
            {
                clThreadsBase *Klass;
                void *vpParam;
            } stParams, *stpParams;

            /**
                Constructor.
            */
            clThreadsBase ();
            /**
                Destructor.
            */
            virtual ~clThreadsBase ();
            /**
                Create thread.
                
                \param vpParam Parameter to the thread (as a stParams instance)
                \param bDetached Start in detached mode
                \return Thread handle
            */
            int Create (void *, bool);
            /**
                Wait for thread to terminate (join thread).
                
                \param iThreadHandle Handle to the thread to wait for
                \return Thread return value
            */
            void *Wait (int);
            /**
                Cancel thread.
                
                Thread will terminate at cancellation point. See IEEE 1003.1
                standard for details.
                
                \param iThreadHandle Handle to the thread to terminate
            */
            void Cancel (int);
            /**
                Return POSIX thread id to the current thread.
                
                \return pthread_t type thread id to the current thread
            */
            pthread_t Self ()
                { return pthread_self(); }
            /**
                Set thread's scheduling parameters.
                
                \note You have to call this unless you start thread in detached mode.

                \param ptidThread pthread_t thread id
                \param iPolicy Scheduling policy
                \param iDelta Delta to the minimum priority
                \return Success?
            */
            bool SetSched (pthread_t, int, int);
            /**
                Abstract (pure) virtual function for thread's function.
            */
            virtual void *InternalCaller (void *) = 0;
    };
    
    
    /**
        Template thread class.
        
        This template class allows individual functions (methods) of class
        to be started as threads.
    */
    template <class TThreads> class clThreads :
        public clThreadsBase
    {
        public:
            typedef void *(TThreads::*Method_t)(void *);
            typedef struct _stParams2
            {
                Method_t Method;
                void *vpParam;
            } stParams2, *stpParams2;

            TThreads *Klass;    ///< Pointer template type class instance

            /**
                Constructor.
                
                \param KlassInst Pointer to template type class instance
            */
            clThreads (TThreads &KlassInst)
                { Klass = &KlassInst; }
            /**
                Create thread.
                
                \param Method Pointer to function (method) in template type class instance to start as a thread
                \param vpParam Parameter to the thread
                \param bDetached Start in detached mode?
                \return Handle to the thread
            */
            int Create (Method_t Method, void *vpParam, bool bDetached = false)
                {
                    stpParams2 spParams;
                    
                    spParams = new stParams2;
                    spParams->Method = Method;
                    spParams->vpParam = vpParam;
                    
                    return clThreadsBase::Create((void *) spParams, 
                        bDetached);
                }
            /**
                Internal thread function (method) caller.
                
                \note This is for purely internal use.
            */
            virtual void *InternalCaller (void *vpParam)
                {
                    stpParams2 spParams = (stpParams2) vpParam;
                    Method_t Method = spParams->Method;
                    void *vpThreadParam = spParams->vpParam;
                    
                    delete spParams;
                    try
                    {
                        return (Klass->*Method)(vpThreadParam);
                    }
                    catch (...)
                    {
                        pthread_exit(NULL);
                    }
                    return NULL;  // dummy for Intel C++
                }
    };

#endif
