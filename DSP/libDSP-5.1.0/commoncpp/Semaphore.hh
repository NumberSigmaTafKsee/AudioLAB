/*

    Semaphore wrapper class
    Copyright (C) 2000-2005 Jussi Laako

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


#ifndef BSDSYS
#include <semaphore.h>
#else
#include <sched.h>
#endif


#ifndef SEMAPHORE_HH
    #define SEMAPHORE_HH


    /**
        Class implementation of POSIX counting semaphores
    */
    class clSemaphore
    {
            #ifndef BSDSYS
            sem_t semSemaphore;
            #else
            volatile unsigned int uiSemaphore;
            #endif
        public:
            /**
                Constructor; create and initialize the semaphore.
                Semaphore value is initialized to zero.
            */
            clSemaphore ()
                #ifndef BSDSYS
                { sem_init(&semSemaphore, 0, 0); }
                #else
                { uiSemaphore = 0; }
                #endif
            /**
                Constructor; create and initialize the semaphore.
                Semaphore is initialized to caller specified value.

                \param uiSemValue Initial semaphore value
            */
            clSemaphore (unsigned int uiSemValue)
                #ifndef BSDSYS
                { sem_init(&semSemaphore, 0, uiSemValue); }
                #else
                { uiSemaphore = uiSemValue; }
                #endif
            /**
                Destructor; destroy the semaphore
            */
            virtual ~clSemaphore ()
                #ifndef BSDSYS
                { sem_destroy(&semSemaphore); }
                #else
                { }
                #endif
            /**
                Initialize the semaphore to specified value.
                This can be used to explicitly set the semaphore value.

                \param uiSemValue New value of semaphore
                \return Success
            */
            bool Initialize (unsigned int uiSemValue)
                #ifndef BSDSYS
                { 
                    if (sem_init(&semSemaphore, 0, uiSemValue) < 0) 
                        return false;
                    return true;
                }
                #else
                {
                    uiSemaphore = uiSemValue;
                    return true;
                }
                #endif
            /**
                Wait until semaphore count becomes non-zero and then
                decrement the count.
            */
            void Wait ()
                #ifndef BSDSYS
                { sem_wait(&semSemaphore); }
                #else
                {
                    while (uiSemaphore == 0) sched_yield();
                    uiSemaphore--;
                }
                #endif
            /**
                Test if semaphore count is non-zero, if it is, then count
                is decremented.
            */
            bool TryWait ()
                #ifndef BSDSYS
                { 
                    if (sem_trywait(&semSemaphore) < 0) return false;
                    return true;
                }
                #else
                {
                    if (uiSemaphore > 0) uiSemaphore--;
                    else return false;
                    return true;
                }
                #endif
            /**
                Post (increment) semaphore

                \return Success
            */
            bool Post ()
                #ifndef BSDSYS
                {
                    if (sem_post(&semSemaphore) < 0) return false;
                    return true;
                }
                #else
                {
                    uiSemaphore++;
                    return true;
                }
                #endif
            /**
                Get value of semaphore

                \return Semaphore count
            */
            int GetValue ()
                #ifndef BSDSYS
                {
                    int iSemValue;
                    sem_getvalue(&semSemaphore, &iSemValue);
                    return iSemValue;
                }
                #else
                { return ((int) uiSemaphore); }
                #endif
    };
    
#endif

