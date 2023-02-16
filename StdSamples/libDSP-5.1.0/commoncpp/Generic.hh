/*

    Portable generic functions
    Copyright (C) 2004-2005 Jussi Laako

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
#include <sys/time.h>
#include <sys/select.h>


#ifndef GENERIC_HH
#define GENERIC_HH


    inline void uSleep (long lUSecs)
        {
            struct timeval sTimeVal;
            
            sTimeVal.tv_sec = lUSecs / 1000000L;
            sTimeVal.tv_usec = lUSecs % 1000000L;
            select(0, NULL, NULL, NULL, &sTimeVal);
        }

#ifdef __i686__
    inline unsigned long long read_tsc ()
        {
            unsigned long long ullTSC = 0;
            
            asm volatile (
                "rdtsc\n"
                : "=A" (ullTSC)
                : "0" (ullTSC));
            return ullTSC;
        }
#endif


#endif
