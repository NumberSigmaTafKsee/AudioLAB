/*

    Helper macros for special compiler features
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


#ifndef COMPILERS_HH
    #define COMPILERS_HH

    #if ((__GNUC__ == 3) && (__GNUC_MINOR__ >= 1))
        #define likely(x)           __builtin_expect(x, 1)
        #define unlikely(x)         __builtin_expect(x, 0)
        #define prefetch(x, w, l)   __builtin_prefetch(x, w, l)
    #else
        #define likely(x)           (x)
        #define unlikely(x)         (x)
        #define prefetch(x, w, l)
    #endif

#endif
