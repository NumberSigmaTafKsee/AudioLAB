#pragma once

#include <cstddef>
#include <cstdlib>
#include <cstdio>
#include <cmath>

// The higher this value is, the more the saturation amount is affected by tape
// speed
#define SPEED_SATURATION_COEFF (0.0)

#define N_FLUTTER_COMPONENTS (11)

namespace FXDSP
{
    static const float flutterRateBase[N_FLUTTER_COMPONENTS] =
    {
        1.01,   // Outer Tension Sense
        2.52,   // Inner Tension Sense
        0.80,   // Pre-stabilizer
        1.01,   // Left lifter
        3.11,   // Record scrape idler
        3.33,   // Capstan
        0.81,   // Pinch roller
        1.01,   // Right lifter
        0.80,   // Motion Sensor
        2.52,   // Inner tension sense
        1.01    // Outer tension sense
    };

    typedef enum  TapeSpeed
    {
        TS_3_75IPS,
        TS_7_5IPS,
        TS_15IPS,
        TS_30IPS
    }TapeSpeed;

    template<typename T>
    calculate_n(T saturation, TapeSpeed speed)
    {
        // Tape speed dependent saturation.
        T n = ((50 * (1-SPEED_SATURATION_COEFF)) +((unsigned)speed * 50 * SPEED_SATURATION_COEFF)) * std::pow((1.0075 - saturation), 2.);
        printf("N: %1.5f\n", n);
        return n;
    }

    /*******************************************************************************
    TapeSaturator */
    template<typename T>
    struct Tape
    {
        PolySaturator<T>*  polysat;
        TapeSpeed       speed;
        T           sample_rate;
        T           saturation;
        T           hysteresis;
        T           flutter;
        T           pos_peak;
        T           neg_peak;
        T*          flutter_mod;
        unsigned    flutter_mod_length;

        Tape(TapeSpeed speed, T saturation, T hysterisis, T flutter, T sample_rate)
        {
            // Create TapeSaturator Struct
            
            polysat = new PolySaturator<T>(1);
            // Allocate for longest period...
            unsigned mod_length = (unsigned)(sample_rate / 0.80);
            float* mod = (float*)malloc(mod_length * sizeof(float));
            // Initialization            
            this->sample_rate = sample_rate;
            this->pos_peak = 0.0;
            this->neg_peak = 0.0;
            this->flutter_mod = mod;
            this->flutter_mod_length = mod_length;
            
            // Need these initialized here.
            this->speed = speed;
            this->saturation = saturation;
            
            // Set up
            setFlutter(flutter);
            setSaturation(saturation);
            setSpeed(speed);
            setHysteresis(hysteresis);
        }
        void setSpeed(TapeSpeed speed)
        {
            // Set speed
            speed = speed;
            
            // Update saturation curve
            polysat->setN(calculate_n(saturation, speed));
            
            // Clear old flutter/wow modulation waveform
            ClearBuffer(flutter_mod, flutter_mod_length); // Yes, clear the old length...
            
            // Calculate new modulation waveform length...
            flutter_mod_length = (unsigned)(sample_rate / \
                                                (0.80 * std::pow(2.0, (float)speed)));
            
            // Generate flutter/wow modulation waveform
            T temp_buffer[flutter_mod_length];
            for (unsigned comp = 0; comp < N_FLUTTER_COMPONENTS; ++comp)
            {
                T phase_step = (2.0 * M_PI * comp * std::pow(2.0, (T)speed)) / sample_rate;
                ClearBuffer(temp_buffer, flutter_mod_length);
                for (unsigned i = 0; i < flutter_mod_length; ++i)
                {
                    temp_buffer[i] = std::sin(i * phase_step) / N_FLUTTER_COMPONENTS;
                }
                VectorVectorAdd(flutter_mod, flutter_mod,
                                temp_buffer, flutter_mod_length);
            }            
        }
        void setSaturation(T saturation)
        {
            T n = calculate_n(saturation, speed);
            this->saturation = saturation;
            return polysat->setN(n);
        }
        void setHysterisis(T hystersis)
        {
            this->hysteresis = hysteresis;
        }
        void setFlutter(T flutter)
        {
            this->flutter = flutter;
        }
        T getSaturation() const {
            return saturation;
        }
        T getHystersis() const {
            return hysterisis;
        }
        T Tick(T I, T A=1, T X=1, T Y=1)
        {
            T hysteresis = hysteresis * 0.05;
            T output = 0.0;
            if (in_sample >= 0)
            {
                neg_peak = 0.0;
                if (in_sample > pos_peak)
                {
                    pos_peak = in_sample;
                    output = in_sample;
                }
                else if (in_sample > (1 - hysteresis) * pos_peak)
                {
                    output = pos_peak;
                }
                else
                {
                    output = in_sample + hysteresis * pos_peak;
                }
            }
            else
            {
                pos_peak = 0.0;
                if (in_sample < neg_peak)
                {
                    neg_peak = in_sample;
                    output = in_sample;
                }
                
                else if (in_sample < (1 - hysteresis) * neg_peak)
                {
                    output = neg_peak;
                }
                
                else
                {
                    output = in_sample + hysteresis * neg_peak;
                }
            }
            return >polysat->setN(output);
        }
        void ProcessBlock(size_t n, T * in, T * out)
        {
            for(size_t i = 0; i < n; i++)
                out[i] = Tick(in[i]);
        }
    };
}

