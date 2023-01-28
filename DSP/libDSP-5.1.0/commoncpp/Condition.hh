/*

    Condition class, header
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
#include <ctime>
#include <sys/time.h>
#include <unistd.h>


#ifndef CONDITION_HH
    #define CONDITION_HH


    /**
        Class implementation of POSIX condition variable
    */
    class clCondition
    {
            pthread_cond_t pthcCond;
        public:
            /**
                Constructor; creates and initializes the condition variable
            */
            clCondition() { pthread_cond_init(&pthcCond, NULL); }
            /**
                Destructor; destroys the condition variable
            */
            virtual ~clCondition() { pthread_cond_destroy(&pthcCond); }
            /**
                Wait for condition to be signalled. 

                \note Mutex must be locked before calling this. Mutex is 
                freed for the waiting time and relocked on return.

                \param pthmCond Pointer to mutex variable
                \return Success
            */
            bool Wait(pthread_mutex_t *pthmCond)
                {
                    if (pthread_cond_wait(&pthcCond, pthmCond) == 0)
                    {
                        return true;
                    }
                    else
                    {
                        return false;
                    }
                }
            /**
                \overload

                \param pthmCond Pointer to mutex varibale
                \param iCondTO Timeout in milliseconds
                \return Success
            */
            bool Wait(pthread_mutex_t *pthmCond, int iCondTO)
                {
                    struct timeval sTimeVal;
                    struct timespec sTimeSpec;
                    gettimeofday(&sTimeVal, NULL);
                    sTimeSpec.tv_sec = sTimeVal.tv_sec + iCondTO / 1000;
                    sTimeSpec.tv_nsec = (sTimeVal.tv_usec + iCondTO * 1000) %
                        1000000 * 1000;
                    if (pthread_cond_timedwait(&pthcCond, pthmCond, 
                        &sTimeSpec) == 0)
                    {
                        return true;
                    }
                    else
                    {
                        return false;
                    }
                }
            /**
                Signal condition variable.
                Mutex can be locked to achieve desired scheduling behaviour.
                Only single waiting thread is woken up.
            */
            void Notify()
                { pthread_cond_signal(&pthcCond); }
            /**
                Signal condition variable.
                Mutex can be locked to achieve desired scheduling behaviour.
                All waiting threads are woken up.
            */
            void NotifyAll() 
                { pthread_cond_broadcast(&pthcCond); }
    };

#endif

