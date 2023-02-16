#pragma once

#include <cstddef>
#include <cstdlib>
#include <cstring>

#include "fxdsp.hpp"

namespace FXDSP
{
    template<typename T>
    struct CircularBuffer
    {
        unsigned    length;
        unsigned    wrap;
        T*          buffer;
        unsigned    read_index;
        unsigned    write_index;
        unsigned    count;

        CircularBuffer(size_t len)
        {
             // use next power of two so we can do a bitwise wrap
            length = next_pow2(length);
            buffer = new T[length];

            ClearBuffer(buffer, length);

            cb->length = length;
            cb->wrap = length - 1;
            cb->buffer = buffer;
            cb->read_index = 0;
            cb->write_index = 0;
            cb->count = 0;
        }
        ~CircularBuffer() {
            if(buffer) delete [] buffer;
        }

        size_t size() const { return length; }
        
        /*******************************************************************************
        CircularBufferWrite */
        void write(const T src, unsigned n_samples)
        {
            #pragma omp simd
            for (unsigned i=0; i < n_samples; ++i)
            {
                buffer[++write_index & cb->wrap] = *src++;
            }
            count += n_samples;

            if (count > length)
            {
                throw std::runtime_error("Value error");
            }            
        }

        void read(T* dest, unsigned n_samples)
        {
            #pragma omp simd
            for (unsigned i=0; i < n_samples; ++i)
            {
                *dest++ = buffer[++read_index & wrap];
            }
            count -= n_samples;

            if (cb->count > cb->length)
            {
                throw std::runtime_error("Value error");
            }            
        }

        void flush()
        {
            ClearBuffer(buffer, length);
            count = 0;            
        }


        void rewind(unsigned samples)
        {
            read_index = ((read_index + length) - samples) % length;
            count += samples;

            if (cb->count > cb->length)
            {
                throw std::runtime_error("Value error");
            }            
        }

    };

}