#pragma once

#include <cfloat>
#include <cmath>
#include <cstdlib>

namespace FXDSP
{
    typedef enum _Opto_t
    {
        /** Light-Dependent-Resistor output. Based
        on Vactrol VTL series. datasheet:
        http://pdf.datasheetcatalog.com/datasheet/perkinelmer/VT500.pdf
        Midpoint Delay values:
        Turn-on delay:   ~10ms
        Turn-off delay:  ~133ms
        */
        OPTO_LDR,

        /** TODO: Add Phototransistor/voltage output opto model*/
        OPTO_PHOTOTRANSISTOR
    } Opto_t;


    /* Scale a sample to the output curve given for the given optocoupler type
    */
    static inline double
    scale_sample(double sample, Opto_t opto_type)
    {

        double out = 0.0;

        switch (opto_type)
        {
            case OPTO_LDR:
                out = 3.0e-6 / pow(sample+DBL_MIN, 0.7);
                out = out > 1.0? 1.0 : out;
                break;
            case OPTO_PHOTOTRANSISTOR:
                out = sample;
                break;
            default:
                break;
        }
        return out;
    }

    /* Calculate the turn-on times [in seconds] for the given optocoupler type with
    the specified delay value
    */
    static inline double
    calculate_on_time(double delay, Opto_t opto_type)
    {
        /* Prevent Denormals */
        double time = DBL_MIN;

        double delay_sq = delay*delay;

        switch (opto_type)
        {
            case OPTO_LDR:
                time = 0.01595 * delay_sq + 0.02795 * delay + 1e-5;
                break;

            case OPTO_PHOTOTRANSISTOR:
                time = 0.01595 * delay_sq + 0.02795 * delay + 1e-5;
                break;

            default:
                break;
        }
        return time;
    }

    /* Calculate the turn-off times [in seconds] for the given optocoupler type with
    the specified delay value
    */
    static inline double
    calculate_off_time(double delay, Opto_t opto_type)
    {
        /* Prevent Denormals */
        double time = DBL_MIN;

        switch (opto_type)
        {
            case OPTO_LDR:
                time = 1.5*powf(delay+FLT_MIN,3.5);
                break;

            case OPTO_PHOTOTRANSISTOR:
                time = 1.5*powf(delay+FLT_MIN,3.5);
                break;
            default:
                break;
        }
        return time;
    }


    /* Opto *******************************************************************/
    template<typename T>
    struct Opto
    {
        Opto_t      type;           //model type
        T       sample_rate;
        T       previous;
        T       delay;
        T       on_cutoff;
        T       off_cutoff;
        T        delta_sign;     // sign of signal dv/dt
        OnePole<T>*  lp;

        Opto(T delay, T sample_rate)
        {
            // Initialization
            this->type = opto_type;
            this->sample_rate = sample_rate;
            this->previous = 0;
            this->delta_sign = 1;
            setDelay(delay);
            lp = new OnePole<T>(opto->on_cutoff, opto->sample_rate, LOWPASS);
        }
        ~Opto()
        {
            if(lp) delete lp;
        }

        void setDelay(T delay)
        {
            this->delay = delay;
            on_cutoff = 1.0/(float)calculate_on_time((double)delay, type);
            off_cutoff = 1.0/(float)calculate_off_time((double)delay, type);
        }
        void ProcessBlock(int n, T * inBuffer, T * outBuffer)
        {
            #pragma omp simd aligned(inBuffer,outBuffer)
            for (unsigned i = 0; i < n; ++i)
            {
                out_buffer[i] = Tick(in_buffer[i]);
            }
        }
        T Tick(T I, T A=1, T X=1, T Y=1)
        {
            T out;
            char prev_delta;

            /* Check sign of dv/dt */
            prev_delta = opto->delta_sign;
            opto->delta_sign = (I - opto->previous) >= 0 ? 1 : -1;

            /* Update lopwass if sign changed */
            if (opto->delta_sign != prev_delta)
            {
                if (opto->delta_sign == 1)
                {
                    lp->setCutoff(opto->on_cutoff);
                }
                else
                {
                    lp->setCutoff(opto->off_cutoff);
                }
            }

            /* Do Delay model */
            out = lp->Tick(I);
            opto->previous = out;
            out = (float)scale_sample((double)out, opto->type);

            /* spit out sample */
            return out;
        }
    };

}
