/*

    Mutex class, header
    Copyright (C) 1999-2005 Jussi Laako

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


#ifndef USE_NPTL
#include <pthread.h>
#else
#include <nptl/pthread.h>
#endif


#ifndef MUTEX_HH
    #define MUTEX_HH


    /**
        Class implementation of POSIX mutex semaphore
    */
    class clMutex
    {
            pthread_mutex_t pthmMutex;
        public:
            /**
                Constructor; creates and initializes the mutex
            */
            clMutex () { pthread_mutex_init(&pthmMutex, NULL); }
            /**
                Desctructor; destroys the mutex
            */
            virtual ~clMutex () { pthread_mutex_destroy(&pthmMutex); }
            /**
                Lock mutex.
                Process is blocked until mutex lock becomes available.

                \return Success
            */
            bool Wait ()
                { 
                    if (pthread_mutex_lock(&pthmMutex) != 0)
                        return false;
                    return true;
                }
            /**
                Unlock mutex

                \return Success
            */
            bool Release ()
                { 
                    if (pthread_mutex_unlock(&pthmMutex) != 0)
                        return false;
                    return true;
                }
            /**
                Try locking mutex without blocking the calling process.

                \return Success
            */
            bool TryLock ()
                {
                    if (pthread_mutex_trylock(&pthmMutex) != 0)
                        return false;
                    return true;
                }
            /**
                Return pointer to mutex variable. For use with condition
                variables.

                \return Pointer to pthread mutex variable
            */
            pthread_mutex_t *GetPtr () { return (&pthmMutex); }
    };


    /**
        Exception safe mutex acquirer class.
    */
    class clSafeMutex
    {
            bool bAcquired;
            clMutex *Mtx;
        public:
            clSafeMutex (clMutex &MtxP)
                {
                    bAcquired = false;
                    Mtx = &MtxP;
                }
            virtual ~clSafeMutex ()
                {
                    if (bAcquired)
                        Mtx->Release();
                }
            /**
                Lock mutex.
                Process is blocked until mutex lock becomes available.

                \return Success
            */
            bool Wait ()
                {
                    if (!Mtx->Wait())
                        return false;
                    bAcquired = true;
                    return true;
                }
            /**
                Unlock mutex

                \return Success
            */
            bool Release ()
                {
                    bAcquired = false;
                    if (!Mtx->Release())
                        return false;
                    return true;
                }
    };

#endif
