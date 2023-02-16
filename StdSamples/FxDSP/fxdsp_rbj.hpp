#pragma once

#include "fxdsp.hpp"
#include "fxdsp_biquad.hpp"

namespace FXDSP
{
    /** Filter types */
    typedef enum Filter_t
    {
        /** Lowpass */
        LOWPASS,

        /** Highpass */
        HIGHPASS,

        /** Bandpass */
        BANDPASS,

        /** Allpass */
        ALLPASS,

        /** Notch */
        NOTCH,

        /** Peaking */
        PEAK,

        /** Low Shelf */
        LOW_SHELF,

        /** High Shelf */
        HIGH_SHELF,

        /** Number of Filter types */
        N_FILTER_TYPES
    }Filter_t;


    /* RBJFilter ***********************************************************/
    template<typename T>
    struct RBJFilter
    {
        BiquadFilter<T>* biquad;
        Filter_t type;
        T omega;
        T Q;
        T cosOmega;
        T sinOmega;
        T alpha;
        T A;
        T dbGain;
        T b[3];
        T a[3];
        T sampleRate;


        void update()
        {
            filter->cosOmega = cos(filter->omega);
            filter->sinOmega = sin(filter->omega);

            switch (filter->type)
            {
            case LOWPASS:
                filter->alpha = filter->sinOmega / (2.0 * filter->Q);
                filter->b[0] = (1 - filter->cosOmega) / 2;
                filter->b[1] = 1 - filter->cosOmega;
                filter->b[2] = filter->b[0];
                filter->a[0] = 1 + filter->alpha;
                filter->a[1] = -2 * filter->cosOmega;
                filter->a[2] = 1 - filter->alpha;
                break;

            case HIGHPASS:
                filter->alpha = filter->sinOmega / (2.0 * filter->Q);
                filter->b[0] = (1 + filter->cosOmega) / 2;
                filter->b[1] = -(1 + filter->cosOmega);
                filter->b[2] = filter->b[0];
                filter->a[0] = 1 + filter->alpha;
                filter->a[1] = -2 * filter->cosOmega;
                filter->a[2] = 1 - filter->alpha;
                break;

            case BANDPASS:
                filter->alpha = filter->sinOmega * sinhf(logf(2.0) / 2.0 * \
                    filter->Q * filter->omega/filter->sinOmega);
                filter->b[0] = filter->sinOmega / 2;
                filter->b[1] = 0;
                filter->b[2] = -filter->b[0];
                filter->a[0] = 1 + filter->alpha;
                filter->a[1] = -2 * filter->cosOmega;
                filter->a[2] = 1 - filter->alpha;
                break;

            case ALLPASS:
                filter->alpha = filter->sinOmega / (2.0 * filter->Q);
                filter->b[0] = 1 - filter->alpha;
                filter->b[1] = -2 * filter->cosOmega;
                filter->b[2] = 1 + filter->alpha;
                filter->a[0] = filter->b[2];
                filter->a[1] = filter->b[1];
                filter->a[2] = filter->b[0];
                break;

            case NOTCH:
                filter->alpha = filter->sinOmega * sinhf(logf(2.0) / 2.0 * \
                    filter->Q * filter->omega/filter->sinOmega);
                filter->b[0] = 1;
                filter->b[1] = -2 * filter->cosOmega;
                filter->b[2] = 1;
                filter->a[0] = 1 + filter->alpha;
                filter->a[1] = filter->b[1];
                filter->a[2] = 1 - filter->alpha;
                break;

            case PEAK:
                filter->alpha = filter->sinOmega * sinhf(logf(2.0) / 2.0 * \
                    filter->Q * filter->omega/filter->sinOmega);
                filter->b[0] = 1 + (filter->alpha * filter->A);
                filter->b[1] = -2 * filter->cosOmega;
                filter->b[2] = 1 - (filter->alpha * filter->A);
                filter->a[0] = 1 + (filter->alpha / filter->A);
                filter->a[1] = filter->b[1];
                filter->a[2] = 1 - (filter->alpha / filter->A);
                break;

            case LOW_SHELF:
                filter->alpha = filter->sinOmega / 2.0 * sqrt( (filter->A + 1.0 / \
                    filter->A) * (1.0 / filter->Q - 1.0) + 2.0);
                filter->b[0] = filter->A * ((filter->A + 1) - ((filter->A - 1) *       \
                    filter->cosOmega) + (2 * sqrtf(filter->A) * filter->alpha));
                filter->b[1] = 2 * filter->A * ((filter->A - 1) - ((filter->A + 1) *   \
                    filter->cosOmega));
                filter->b[2] = filter->A * ((filter->A + 1) - ((filter->A - 1) *       \
                    filter->cosOmega) - (2 * sqrtf(filter->A) * filter->alpha));
                filter->a[0] = ((filter->A + 1) + ((filter->A - 1) *                   \
                    filter->cosOmega) + (2 * sqrtf(filter->A) * filter->alpha));
                filter->a[1] = -2 * ((filter->A - 1) + ((filter->A + 1) *              \
                    filter->cosOmega));
                filter->a[2] = ((filter->A + 1) + ((filter->A - 1) *                   \
                    filter->cosOmega) - (2 * sqrtf(filter->A) * filter->alpha));
                break;

            case HIGH_SHELF:
                filter->alpha = filter->sinOmega / 2.0 * sqrt( (filter->A + 1.0 / \
                    filter->A) * (1.0 / filter->Q - 1.0) + 2.0);
                filter->b[0] = filter->A * ((filter->A + 1) + ((filter->A - 1) *       \
                    filter->cosOmega) + (2 * sqrtf(filter->A) * filter->alpha));
                filter->b[1] = -2 * filter->A * ((filter->A - 1) + ((filter->A + 1) *  \
                    filter->cosOmega));
                filter->b[2] = filter->A * ((filter->A + 1) + ((filter->A - 1) *       \
                    filter->cosOmega) - (2 * sqrtf(filter->A) * filter->alpha));
                filter->a[0] = ((filter->A + 1) - ((filter->A - 1) *                   \
                    filter->cosOmega) + (2 * sqrtf(filter->A) * filter->alpha));
                filter->a[1] = 2 * ((filter->A - 1) - ((filter->A + 1) *               \
                    filter->cosOmega));
                filter->a[2] = ((filter->A + 1) - ((filter->A - 1) *                   \
                    filter->cosOmega) - (2 * sqrtf(filter->A) * filter->alpha));
                break;

            default:
                return ERROR;
                break;
            }

            // Normalize filter coefficients
            T factor = 1.0 / filter->a[0];
            T norm_a[2];
            T norm_b[3];
            VectorScalarMultiply(norm_a, &filter->a[1], factor, 2);
            VectorScalarMultiply(norm_b, filter->b, factor, 3);
            BiquadCoefficients c(norm_b,norm_a);
            biquad->setCoefficients(c);
        }

        RBJFilter(Filter_t type, T cutoff, T sampleRate)
        {
            // Initialization
            this->type = type;
            filter->omega =  HZ_TO_RAD(cutoff) / sampleRate; //hzToRadians(cutoff, sampleRate);
            filter->Q = 1;
            filter->A = 1;
            filter->dbGain = 0;
            this->sampleRate = sampleRate;

            // Initialize biquad
            float b[3] = {0, 0, 0};
            float a[2] = {0, 0};            
            biquad = new BiquadFilter<T>();

            // Calculate coefficients
            update();
        }  
        ~RBJFilter
        {
            if(biquad) delete biquad;
        }

        void setType(Filter_t type)
        {
            filter->type = type;
            update();
        }
        void setCutoff(T cutoff)
        {
            filter->omega = HZ_TO_RAD(cutoff) / filter->sampleRate;
            update();
        }
        void setQ(T q) {
            Q = q;
            update();
        }
        void setParams(Filter_t type, T c, T q)
        {
            this->type = type;
            omega = HZ_TO_RAD(cutoff) / filter->sampleRate;
            this->Q = q;
            update();
        }
        void ProcessBlock(size_t n, T * inBuffer, T * outBuffer)
        {
            biquad->ProcessBlock(n,inBuffer,outBuffer);
        }
        void flush() {
            biquad->flush();
        }
        T Tick(T I, T A=1, T X=1, T Y=1)
        {
            return 0;
        }
    };

}