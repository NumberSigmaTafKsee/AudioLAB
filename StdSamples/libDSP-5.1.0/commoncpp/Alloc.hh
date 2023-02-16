/*

    Memory allocation wrapper class
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


#include <cstdlib>
#include <sys/types.h>
#include <sys/mman.h>
#include <string>
#include <stdexcept>
#include <cstring>

#ifndef ALLOC_HH
    #define ALLOC_HH

    #if (!defined(DSP_IPP) && !defined(ALLOC_DEBUG))
        /*#define ALLOC_USE_IPP
        #include <ipp.h>*/
        #define ALLOC_USE_ALIGN
        #define ALLOC_ALIGNMENT     32
    #elif (defined(DSP_IPP))
        #define ALLOC_USE_IPP
        #include <ipp.h>
    #endif


    /**
        Class for memory allocation operations. Mainly to avoid memory leaks
        and to simplify pointer typecasts.
    */
    class clAlloc
    {
        protected:
            bool bLocked;
            long lSize;
            void *vpPtr;
        public:
            /** 
                Constructor; initializes empty allocation
            */
            clAlloc ()
                {
                    bLocked = false;
                    lSize = 0l;
                    vpPtr = NULL;
                }
            /**
                Copy constructor; copies instance
            */
            clAlloc (const clAlloc &CopySrc)
                {
                    Copy(CopySrc);
                }
            /** 
                Constructor; allocates specified amount of memory

                \param lAllocSize Size, in bytes, of allocation
                \throw runtime_error
            */
            clAlloc (long lAllocSize)
                {
                    bLocked = false;
                    lSize = 0l;
                    vpPtr = NULL;
                    Size(lAllocSize);
                }
            /** 
                Destructor; frees allocated memory block
            */
            virtual ~clAlloc ()
                { Free(); }
            /** 
                Allocate specified amount of memory, previously allocated
                block is freed before allocating new one.

                \param lAllocSize Size, in bytes, of new allocation
                \return Pointer to memory block
                \throw runtime_error
            */
            void *SetSize (long lAllocSize)
                {
                    if (lSize == lAllocSize) return vpPtr;
                    Free();
                    if (lAllocSize <= 0l) return vpPtr;
                    #if defined(ALLOC_USE_ALIGN)
                    if (posix_memalign(&vpPtr, ALLOC_ALIGNMENT, lAllocSize))
                        throw std::runtime_error(std::string("Out of memory!"));
                    #elif defined(ALLOC_USE_IPP)
                    vpPtr = ippMalloc(lAllocSize);
                    if (vpPtr == NULL)
                        throw std::runtime_error(std::string("Out of memory!"));
                    #else
                    vpPtr = (void *) malloc(lAllocSize);
                    if (vpPtr == NULL) 
                        throw std::runtime_error(std::string("Out of memory!"));
                    #endif
                    lSize = lAllocSize;
                    return vpPtr;
                }
            /// \overload
            void *Size (size_t sizeNew)
                { return SetSize((long) sizeNew); }
            /**
                Resize memory block to specified size

                \param lAllocSize New size, in bytes, of memory block
                \return Pointer to memory block
                \throw runtime_error
            */
            void *Resize (long lAllocSize)
                {
                    if (bLocked) UnLock();
                    if (lSize == lAllocSize) return vpPtr;
                    if (lAllocSize <= 0l)
                    {
                        Free();
                        return vpPtr;
                    }
                    #if defined(ALLOC_USE_ALIGN)
                    void *vpOldPtr = vpPtr;
                    vpPtr = NULL;
                    if (posix_memalign(&vpPtr, ALLOC_ALIGNMENT, lAllocSize))
                        throw std::runtime_error(std::string("Out of memory!"));
                    memcpy(vpPtr, vpOldPtr, 
                        (lSize <= lAllocSize) ? lSize : lAllocSize);
                    free(vpOldPtr);
                    #elif defined(ALLOC_USE_IPP)
                    void *vpOldPtr = vpPtr;
                    vpPtr = ippMalloc(lAllocSize);
                    if (vpPtr == NULL)
                        throw std::runtime_error(std::string("Out of memory!"));
                    ippsCopy_8u((Ipp8u *) vpOldPtr, (Ipp8u *) vpPtr, 
                        (lSize <= lAllocSize) ? lSize : lAllocSize);
                    ippFree(vpOldPtr);
                    #else
                    vpPtr = (void *) realloc(vpPtr, lAllocSize);
                    if (vpPtr == NULL)
                        throw std::runtime_error(std::string("Out of memory!"));
                    #endif
                    lSize = lAllocSize;
                    return vpPtr;
                }
            /**
                Free allocated memory block. It is not an error to free
                already freed or non-allocated block.
            */
            void Free ()
                {
                    if (bLocked) UnLock();
                    if (vpPtr != NULL) 
                    {
                        #if defined(ALLOC_USE_IPP)
                        ippFree(vpPtr);
                        #else
                        free(vpPtr);
                        #endif
                        lSize = 0l;
                        vpPtr = NULL;
                    }
                }
            /**
                Get pointer to allocated block

                \return Pointer to memory block
            */
            void *GetPtr () const
                { return vpPtr; }
            /**
                Get size of allocation

                \return Size, in bytes, of allocation
            */
            long GetSize () const
                { return lSize; }
            /// \overload
            size_t Size () const
                { return ((size_t) lSize); }
            /**
                Lock memory block (prevent paging out)
            */
            void Lock ()
                {
                    bLocked = true;
                    #if (defined(SOLARIS))
                    mlock((char *) vpPtr, lSize);
                    #else
                    mlock(vpPtr, lSize);
                    #endif
                }
            /**
                Unlock memory block (allow paging out)
            */
            void UnLock ()
                {
                    bLocked = false;
                    #if (defined(SOLARIS))
                    munlock((char *) vpPtr, lSize);
                    #else
                    munlock(vpPtr, lSize);
                    #endif
                }
            /**
                Copy specified memory block to this instance

                \param Src Source of copy
                \throw runtime_error
            */
            void Copy (const clAlloc &Src)
                {
                    Free();
                    Size(Src.lSize);
                    memcpy(vpPtr, Src.vpPtr, lSize);
                }
            /**
                Copy contents of this instance to specified memory block

                \param Dest Destination of copy
                \throw runtime_error
            */
            void CopyTo (clAlloc &Dest)
                {
                    Dest.Copy(*this);
                }
            /**
                Return copy of this instance

                \return Copy of this memory block
            */
            clAlloc GetCopy ()
                {
                    clAlloc Temp;

                    Temp.Copy(*this);
                    return Temp;
                }
            /**
                Return pointer to memory block
            */
            operator char *() const
                { return ((char *) vpPtr); }
            /// \overload
            operator unsigned char *() const
                { return ((unsigned char *) vpPtr); }
            /// \overload
            operator short *() const
                { return ((short *) vpPtr); }
            /// \overload
            operator unsigned short *() const
                { return ((unsigned short *) vpPtr); }
            /// \overload
            operator int *() const
                { return ((int *) vpPtr); }
            /// \overload
            operator unsigned int *() const
                { return ((unsigned int *) vpPtr); }
            /// \overload
            operator long *() const
                { return ((long *) vpPtr); }
            /// \overload
            operator unsigned long *() const
                { return ((unsigned long *) vpPtr); }
            /// \overload
            operator float *() const
                { return ((float *) vpPtr); }
            /// \overload
            operator double *() const
                { return ((double *) vpPtr); }
            /// \overload
            operator long double *() const
                { return ((long double *) vpPtr); }
            /// \overload
            operator void *() const
                { return vpPtr; }
            operator void *()
                { return vpPtr; }
            /**
                Self explanatory

                \param Src Source of copy
                \return Reference to this
                \throw runtime_error
            */
            clAlloc & operator = (const clAlloc &Src)
                { 
                    Copy(Src);
                    return (*this);
                }
    };

#endif

