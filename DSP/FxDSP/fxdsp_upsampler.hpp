#pragma once

#include <cmath>
#include <cstddef>
#include <cstdlib>
#include "fxdsp.hpp"
#include "fxdsp_fir.hpp"
#include "fxdsp_polyphase.hpp"

namespace FXDSP
{


    template<typename T, class FIRFilter>    
    struct Upsampler
    {
        unsigned factor;
        FIRFilter** polyphase;


        Upsampler(ResampleFactor_t factor)
        {
            unsigned n_filters = 1;
            switch(factor)
            {
                case X2:
                    n_filters = 2;
                    break;
                case X4:
                    n_filters = 4;
                    break;
                case X8:
                    n_filters = 8;
                    break;
                /*
                case X16:
                    n_filters = 16;
                    break;
                */
                default:
                    return NULL;
            }

            
            // Allocate memory for the polyphase array
            polyphase = new FIRFilter* [n_filters];
            assert(polyphase != nullptr);
            polyphase = polyphase;

            // Create polyphase filters
            unsigned idx;
            for(idx = 0; idx < n_filters; ++idx)
            {
                polyphase[idx] = new FIRFilter(PolyphaseCoeffs[factor][idx], 64, DIRECT);
            }

            // Add factor
            factor = n_filters;                        
        }
        ~Upsampler()
        {
            if (polyphase)
            {
                for (unsigned i = 0; i < factor; ++i)
                {
                    delete polyphase[i];
                }
                delete [] polyphase;
            }            
        }
        void flush() {
            unsigned idx;
            for (idx = 0; idx < factor; ++idx)
            {
                polyphase[idx]->flush();;
            }
        }
        void ProcessBlock(size_t n, T * inBuffer, T * outBuffer)
        {
            T tempbuf[n];
            #pragma omp simd aligned(inBuffer,outBuffer)
            for (unsigned filt = 0; filt < factor; ++filt)
            {
                polyphase[filt]->ProcessBlock(n, inBuffer, tempbuf);
                CopyBufferStride(outBuffer+filt, factor, tempbuf, 1, n_samples);
            }
            VectorScalarMultiply(outBuffer, (const T*)outBuffer,
                                factor, n_samples * factor);
        }            
    };
}
