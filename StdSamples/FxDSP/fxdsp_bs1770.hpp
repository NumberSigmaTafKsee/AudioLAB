#pragma once

#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#include "fxdsp.hpp"
#include "fxdsp_circular_buffer.hpp"
#include "fxdsp_biquad.hpp"
#include "fxdsp_upsampler.hpp"

#define PREFILTER_FC    (1500.12162162162167078350)
#define PREFILTER_GAIN  (3.99976976976976983380)
#define PREFILTER_Q     (sqrt(2.0)/2.0) // Close enough....
#define RLBFILTER_FC    (40.2802802803)
#define RLBFILTER_Q     (0.92792792793)
#define GATE_LENGTH_S   (0.4)
#define GATE_OVERLAP    (0.75)


namespace FXDSP
{
    /* Channel id numbers */
    enum
    {
        LEFT = 0,
        RIGHT,
        CENTER,
        LEFT_SURROUND,
        RIGHT_SURROUND,
        N_CHANNELS
    };


    double CHANNEL_GAIN[N_CHANNELS] =
    {
        1.0,    /* LEFT */
        1.0,    /* RIGHT */
        1.0,    /* CENTER */
        1.41,   /* LEFT_SURROUND */
        1.41    /* RIGHT_SURROUND */
    };


    template<typename T>
    struct KWeightingFilter
    {
        BiquadFilter<T>*   pre_filter;
        BiquadFilter<T>*   rlb_filter;

        KWeightingFilter(T sample_rate)
        {
            T b[3] = {0.};
            T a[2] = {0.};         
            calc_prefilter(b, a, sample_rate);
            filter->pre_filter = new BiquadFilter<T>(b, a);
            calc_rlbfilter(b, a, sample_rate);
            filter->rlb_filter = new BiquadFilter<T>(b, a);
        }
        ~KWeightingFilter()
        {
            if(pre_filter) delete pre_filter;
            if(rlb_filter) delete rlb_filter;
        }
        void ProcessBlock(size_t n, T * src, T * dst)
        {
            T scratch[length];
            filter->pre_filter->ProcessBlock(n,src,scratch);
            filter->rlb_filter->ProcessBlock(n,scratch,dst);
        }
    };
    };

    
    tempalte<typename T, class Upsampler>
    struct BS1770Meter
    {
        KWeightingFilter<T>**  filters;
        Upsampler**         upsamplers;
        CircularBuffer<T>**    buffers;
        unsigned            n_channels;
        unsigned            sample_count;
        unsigned            gate_len;
        unsigned            overlap_len;

        BS1770MeterT(unsigned n_channels, T sample_rate)
        {
            
            filters = new KWeightingFilter* [n_channels]; 
            upsampler =  new Upsampler*[n_channels];(Upsampler**)malloc(n_channels * sizeof(Upsampler*));
            buffers = new CircularBuffer*[n_channels]; 
            for (unsigned i = 0; i < n_channels; ++i)
            {
                filters[i] = KWeightingFilterInit(sample_rate);
                upsamplers[i] = UpsamplerInit(X4);
                buffers[i] = CircularBufferInit((unsigned)(2 * GATE_LENGTH_S * sample_rate));
            }

            sample_count = 0;
            n_channels = n_channels;
            gate_len = (unsigned)(GATE_LENGTH_S * sample_rate);
            overlap_len = (unsigned)(GATE_OVERLAP * gate_len);
            filters= filters;
            upsamplers = upsamplers;
            buffers = buffers;
        }
        void ProcessBlock(T * loudness, T ** peaks, const T ** samples, size_t n_samples)
        {
            unsigned os_length = 4 * n_samples;
            T filtered[n_samples];
            T os_sig[os_length];
            T gate[gate_len];
            T sum = 0.0;

            if (meter)
            {
                *loudness = 0.0;
				#pragma omp simd aligned(peaks,samples)
                for (unsigned i = 0; i < n_channels; ++i)
                {
                    // Calculate peak for each channel
                    upsamplers[i]->ProcessBlock(n_samples,samples[i], os_sig);
                    VectorAbs(os_sig, (const T*)os_sig, os_length);
                    *peaks[i] = AmpToDb(VectorMax(os_sig, os_length));

                    filters[i]->ProcessBlock(n_samples,samples[i], filtered);
                    CircularBufferWrite(buffers[i], (const T*)filtered, n_samples);

                    if (CircularBufferCount(buffers[i]) >= gate_len)
                    {
                        CircularBufferRead(buffers[i], gate, gate_len);
                        CircularBufferRewind(buffers[i], overlap_len);
                        sum += CHANNEL_GAIN[i] * MeanSquare(gate, gate_len);
                    }
                }
                *loudness = -0.691 + 10 * std::log10(sum);                
            }
        }
    };
              
    
}
