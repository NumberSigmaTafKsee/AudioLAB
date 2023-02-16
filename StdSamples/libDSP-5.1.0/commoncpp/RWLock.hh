/*

    Read-write locks class, header
    Copyright (C) 2001-2005 Jussi Laako

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


#ifdef _XOPEN_SOURCE
#   if (_XOPEN_SOURCE < 600)
#       undef _XOPEN_SOURCE
#       define _XOPEN_SOURCE 600
#   endif
#else
#   define _XOPEN_SOURCE 600
#endif


#ifndef USE_NPTL
#include <pthread.h>
#else
#include <nptl/pthread.h>
#endif


#ifndef RWLOCK_HH
#define RWLOCK_HH


    /**
        Class implementation of SUSv2 read-write locks
    */
    class clRWLock
    {
            pthread_rwlock_t ptrwlLock;
        public:
            /**
                Constructor; creates and initializes the read-write lock
            */
            clRWLock ()
                { pthread_rwlock_init(&ptrwlLock, NULL); }
            /**
                Desctructor; destroys the read-write lock
            */
            virtual ~clRWLock ()
                { pthread_rwlock_destroy(&ptrwlLock); }
            /**
                Obtain read-lock.
                There can be multiple read-locks at the same time.

                \return Success
            */
            bool WaitRead ()
                {
                    if (pthread_rwlock_rdlock(&ptrwlLock) != 0)
                        return false;
                    return true;
                }
            /**
                Obtain exclusive write-lock.
                There can be only one write-lock and no read-locks at the
                same time.

                \return Success
            */
            bool WaitWrite ()
                {
                    if (pthread_rwlock_wrlock(&ptrwlLock) != 0)
                        return false;
                    return true;
                }
            /**
                Release lock

                \return Success
            */
            bool Release ()
                {
                    if (pthread_rwlock_unlock(&ptrwlLock) != 0)
                        return false;
                    return true;
                }
            /**
                Try to obtain read-lock (non-blocking)

                \return Success
            */
            bool TryRead ()
                {
                    if (pthread_rwlock_tryrdlock(&ptrwlLock) != 0)
                        return false;
                    return true;
                }
            /**
                Try to obtain exclusive write-lock (non-blocking)

                \return Success
            */
            bool TryWrite ()
                {
                    if (pthread_rwlock_trywrlock(&ptrwlLock) != 0)
                        return false;
                    return true;
                }
    };

#endif
