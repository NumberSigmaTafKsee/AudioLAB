#pragma once

#include <cstddef>
#include <cstdlibh>

#include "fxdsp.hpp"
#include "fxdsp_fir.hpp"

namespace FXDSP
{
    
    template<typename T, class FIRFilter>
    struct DecimatorT
    {
        unsigned        factor;
        FIRFilter  t ** polyphase;
        unsigned        n_filters=0;

        Decimator(ResampleFactor_t factor)
        {
            n_filters = 1;
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
            /*  case X16:
                    n_filters = 16;
                    break; */
                default:
                    return NULL;
            }

            
            // Allocate memory for the polyphase array
            polyphase = new FIRFilter* [n_filters];
            assert(polyphase != nullptr);
            
            // Create polyphase filters
            unsigned idx;
            for(idx = 0; idx < n_filters; ++idx)
            {
                polyphase[idx] = new FIRFilter(PolyphaseCoeffs[factor][idx], 64, DIRECT);
            }

            // Add factor
            factor = n_filters;            
        }
        ~Decimator() {
            if(polyphase) {
                for(size_t i = 0; i < n_filters; i++)
                    if(polyphase[i]) delete polyphase[i];
                delete [] polyphase;
            }
        }
        void flush() {
            unsigned idx;
            for (idx = 0; idx < factor; ++idx)
            {
               polyphase[idx]->flush();
            }
        }
        void ProcessBlock(size_t n, const T * inBuffer, T * outBuffer)
        {
            if (outBuffer)
            {
                unsigned declen = n_samples / factor;                
                ClearBuffer(outBuffer, declen);

                T temp[n];
                for (unsigned filt = 0; filt < factor; ++filt)
                {
                    CopyBufferStride(temp, 1, inBuffer, factor, declen);
                    polyphase[filt]->ProcessBlock(declen,temp,temp);
                    VectorVectorAdd(outBuffer, (const T*)outBuffer, temp, declen);
                }                
            }
            else
            {
                throw std::runtime_error("Null pointer");
            }
        }
    };

    
}
