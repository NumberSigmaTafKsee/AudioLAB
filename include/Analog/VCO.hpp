#pragma once
#include <cmath>
#include "SoundObject.hpp"


namespace Analog
{
    enum {
        POLYBLEP,
        BLIT,
        MINBLEP,
        DPW,
    };

    const DspFloatType TWO_PI = 2 * M_PI;

    template<typename T>
    inline T square_number(const T &x) {
        return x * x;
    }

    // Adapted from "Phaseshaping Oscillator Algorithms for Musical Sound
    // Synthesis" by Jari Kleimola, Victor Lazzarini, Joseph Timoney, and Vesa
    // Valimaki.
    // http://www.acoustics.hut.fi/publications/papers/smc2010-phaseshaping/
    inline DspFloatType blep(DspFloatType t, DspFloatType dt) {
        if (t < dt) {
            return -square_number(t / dt - 1);
        } else if (t > 1 - dt) {
            return square_number((t - 1) / dt + 1);
        } else {
            return 0;
        }
    }

    // Derived from blep().
    inline DspFloatType blamp(DspFloatType t, DspFloatType dt) {
        if (t < dt) {
            t = t / dt - 1;
            return -1 / 3.0 * square_number(t) * t;
        } else if (t > 1 - dt) {
            t = (t - 1) / dt + 1;
            return 1 / 3.0 * square_number(t) * t;
        } else {
            return 0;
        }
    }

    template<typename T>
    inline int64_t bitwiseOrZero(const T &t) {
        return static_cast<int64_t>(t) | 0;
    }

    // this is vectorized
    struct VCOPolyBLEP : public OscillatorProcessor
    {
        std::vector<OscillatorProcessor*> slaves;
    
        int    m_waveform = 0;
        DspFloatType m_morph = 0;
        DspFloatType m_freq  = 440.0f;
        DspFloatType m_phase = 0;
        DspFloatType m_index = 1;
        DspFloatType m_gain = 1;
        DspFloatType m_fm = 0;
        DspFloatType m_pm = 0;
        DspFloatType m_fenv = 1;
        DspFloatType m_penv = 1;  
        DspFloatType m_drift = 0;  
        DspFloatType m_mod = 1;
        DspFloatType m_cmod = 1;
        DspFloatType m_env = 1;
        DspFloatType m_lfo = 1;
        DspFloatType m_pwm = 0.5;

        
        enum Waveform {
            SINE=0,
            COSINE,
            TRIANGLE,
            SQUARE,
            RECTANGLE,
            SAWTOOTH,
            RAMP,
            MODIFIED_TRIANGLE,
            MODIFIED_SQUARE,
            HALF_WAVE_RECTIFIED_SINE,
            FULL_WAVE_RECTIFIED_SINE,
            TRIANGULAR_PULSE,
            TRAPEZOID_FIXED,
            TRAPEZOID_VARIABLE
        };

        Waveform waveform;
        DspFloatType sampleRate;
        DspFloatType freqInSecondsPerSample;
        DspFloatType amplitude; // Frequency dependent gain [0.0..1.0]
        DspFloatType pulseWidth; // [0.0..1.0]
        DspFloatType t; // The current phase [0.0..1.0) of the oscillator.


        VCOPolyBLEP(DspFloatType sampleRate, Waveform wave = SINE)
        : sampleRate(sampleRate), amplitude(1.0), t(0.0), OscillatorProcessor()
        {     
            setSampleRate(sampleRate);
            setFrequency(440.0);        
            setWaveform(wave);
            setPulseWidth(0.5);
        }
        VCOPolyBLEP() = default;
        ~VCOPolyBLEP() = default;
        
        void init(DspFloatType sr)
        {
            sampleRate = sr;
            amplitude  = 1.0;
            t = 0.0;
            setSampleRate(sampleRate);
            setFrequency(440.0);        
            setWaveform(SAWTOOTH);
            setPulseWidth(0.5);
        }
        void setFrequency(DspFloatType freqInHz) {
            setdt(freqInHz / sampleRate);
            m_freq = freqInHz;
        }
        void setdt(DspFloatType time) {
            freqInSecondsPerSample = time;
            m_phase = t;
        }
        void setSampleRate(DspFloatType sampleRate)  {
            const DspFloatType freqInHz = getFreqInHz();
            this->sampleRate = sampleRate;
            setFrequency(freqInHz);
        }
        void setWaveform(Waveform waveform)
        {
            this->waveform = waveform;
        }
        void setPhase(DspFloatType p) {
            t = p;
            m_phase = p;
        }
        DspFloatType getPhase() {
            return t;
        }
        void setPulseWidth(DspFloatType pw)  {
            this->pulseWidth = pw;
            m_pwm = pw;
        }

        DspFloatType get()  {
            DspFloatType r = 0;
            DspFloatType tfreq  = m_freq;
            DspFloatType tphase = t;
            setPulseWidth(m_pwm);        
            setFrequency(m_index * m_fenv*(m_freq + m_freq*m_fm));
            t += m_pm;
            t *= m_penv;
            if(getFreqInHz() >= sampleRate / 2) {
                return sin();
            } else switch (waveform) {
                case SINE:                      
                    r = sin();
                    break;
                case COSINE:
                    r = cos();
                    break;
                case TRIANGLE:
                    r = tri();
                    break;
                case SQUARE:
                    r = sqr();
                    break;
                case RECTANGLE:
                    r = rect();
                    break;
                case SAWTOOTH:                                
                    r = saw();
                    break;
                case RAMP:
                    r = ramp();
                    break;
                case MODIFIED_TRIANGLE:
                    r = tri2();
                    break;
                case MODIFIED_SQUARE:
                    r = sqr2();
                    break;
                case HALF_WAVE_RECTIFIED_SINE:
                    r = half();
                    break;
                case FULL_WAVE_RECTIFIED_SINE:
                    r = full();
                    break;
                case TRIANGULAR_PULSE:
                    r = trip();
                    break;
                case TRAPEZOID_FIXED:
                    r = trap();
                    break;
                case TRAPEZOID_VARIABLE:
                    r = trap2();
                    break;
                default:
                    r = 0.0;
            }
            setFrequency(tfreq);
            t = tphase;
            m_phase = t;
            r *= m_env;
            r *= m_lfo;
            r *= m_mod;
            r  = fmod(r,m_cmod);
            t += m_drift;
            return m_gain*r;
        }
        

        DspFloatType Tick(DspFloatType I=0, DspFloatType A=1, DspFloatType X=0, DspFloatType Y=0)
        {
            return getAndInc();
        }

        void ProcessSIMD(size_t n, DspFloatType * out) {
            DspFloatType phase[n];            
            // generate the phase vector            
            for(size_t i = 0; i < n; i++)
            {
                phase[i] = t;
                inc();
            }
            // now blow this bitch
            if(getFreqInHz() >= sampleRate / 2) {
                sin(n,phase,out);
            } else switch (waveform) {
                case SINE:                      
                    sin(n,phase,out);
                    break;
                case COSINE:
                    cos(n,phase,out);
                    break;
                case TRIANGLE:
                    tri(n,phase,out);
                    break;
                case SQUARE:
                    sqr(n,phase,out);
                    break;
                case RECTANGLE:
                    rect(n,phase,out);
                    break;
                case SAWTOOTH:                                
                    saw(n,phase,out);
                    break;
                case RAMP:
                    ramp(n,phase,out);
                    break;
                case MODIFIED_TRIANGLE:
                    tri2(n,phase,out);
                    break;
                case MODIFIED_SQUARE:
                    sqr2(n,phase,out);
                    break;
                case HALF_WAVE_RECTIFIED_SINE:
                    half(n,phase,out);
                    break;
                case FULL_WAVE_RECTIFIED_SINE:
                    full(n,phase,out);
                    break;
                case TRIANGULAR_PULSE:
                    trip(n,phase,out);
                    break;
                case TRAPEZOID_FIXED:
                    trap(n,phase,out);
                    break;
                case TRAPEZOID_VARIABLE:
                    trap2(n,phase,out);
                    break;
                default:
                    memset(out,0x00,n*sizeof(DspFloatType));
            }
        }

        void inc() {
            t += freqInSecondsPerSample;
            t -= bitwiseOrZero(t);
        }

        DspFloatType getAndInc() {
            const DspFloatType sample = get();
            inc();
            return sample;
        }  

        DspFloatType getFreqInHz() 
        {
            return freqInSecondsPerSample * sampleRate;
        }

        void sync(DspFloatType phase)
        {
            t = phase;
            if (t >= 0) {
                t -= bitwiseOrZero(t);
            } else {
                t += 1 - bitwiseOrZero(t);
            }
        }    
        DspFloatType sin() {
            return amplitude * std::sin(TWO_PI * t);
        }
        void sin(size_t n, DspFloatType *phase, DspFloatType * out) {
            #pragma omp simd
            for(size_t i = 0; i < n; i++) {
                out[i] = amplitude * std::sin(TWO_PI * phase[i]);
            }
        }

        DspFloatType cos() {
            return amplitude * std::cos(TWO_PI * t);
        }
        void cos(size_t n, DspFloatType *phase, DspFloatType * out) {
            #pragma omp simd
            for(size_t i = 0; i < n; i++) {
                out[i] = amplitude * std::cos(TWO_PI * phase[i]);
            }
        }

        DspFloatType half() {
            DspFloatType t2 = t + 0.5;
            t2 -= bitwiseOrZero(t2);

            DspFloatType y = (t < 0.5 ? 2 * std::sin(TWO_PI * t) - 2 / M_PI : -2 / M_PI);
            y += TWO_PI * freqInSecondsPerSample * (blamp(t, freqInSecondsPerSample) + blamp(t2, freqInSecondsPerSample));

            return amplitude * y;
        }
        void half(size_t n, DspFloatType *phase, DspFloatType * out) {
            #pragma omp simd
            for(size_t i = 0; i < n; i++) {
                DspFloatType t  = phase[i];
                DspFloatType t2 = t + 0.5;
                t2 -= bitwiseOrZero(t2);

                DspFloatType y = (t < 0.5 ? 2 * std::sin(TWO_PI * t) - 2 / M_PI : -2 / M_PI);
                y += TWO_PI * freqInSecondsPerSample * (blamp(t, freqInSecondsPerSample) + blamp(t2, freqInSecondsPerSample));
                out[i] = amplitude * y;
            }
        }

        DspFloatType full() {
            DspFloatType _t = this->t + 0.25;
            _t -= bitwiseOrZero(_t);

            DspFloatType y = 2 * std::sin(M_PI * _t) - 4 / M_PI;
            y += TWO_PI * freqInSecondsPerSample * blamp(_t, freqInSecondsPerSample);

            return amplitude * y;
        }
        void full(size_t n, DspFloatType *phase, DspFloatType * out) {
            #pragma omp simd
            for(size_t i = 0; i < n; i++) {
                DspFloatType _t = phase[i] + 0.25;
                _t -= bitwiseOrZero(_t);

                DspFloatType y = 2 * std::sin(M_PI * _t) - 4 / M_PI;
                y += TWO_PI * freqInSecondsPerSample * blamp(_t, freqInSecondsPerSample);

                out[i] = amplitude * y;
            }
        }

        DspFloatType tri() {
            DspFloatType t1 = t + 0.25;
            t1 -= bitwiseOrZero(t1);

            DspFloatType t2 = t + 0.75;
            t2 -= bitwiseOrZero(t2);

            DspFloatType y = t * 4;

            if (y >= 3) {
                y -= 4;
            } else if (y > 1) {
                y = 2 - y;
            }

            y += 4 * freqInSecondsPerSample * (blamp(t1, freqInSecondsPerSample) - blamp(t2, freqInSecondsPerSample));

            return amplitude * y;
        }
        void tri(size_t n, DspFloatType *phase, DspFloatType * out) {
            #pragma omp simd
            for(size_t i = 0; i < n; i++) {
                DspFloatType t1 = phase[i] + 0.25;
                t1 -= bitwiseOrZero(t1);

                DspFloatType t2 = phase[i] + 0.75;
                t2 -= bitwiseOrZero(t2);

                DspFloatType y = phase[i] * 4;

                if (y >= 3) {
                    y -= 4;
                } else if (y > 1) {
                    y = 2 - y;
                }

                y += 4 * freqInSecondsPerSample * (blamp(t1, freqInSecondsPerSample) - blamp(t2, freqInSecondsPerSample));

                out[i] = amplitude * y;
            }
        }

        DspFloatType tri2() {
            
            DspFloatType pulseWidth = std::fmax(0.0001, std::fmin(0.9999, this->pulseWidth));

            DspFloatType t1 = t + 0.5 * pulseWidth;
            t1 -= bitwiseOrZero(t1);

            DspFloatType t2 = t + 1 - 0.5 * pulseWidth;
            t2 -= bitwiseOrZero(t2);

            DspFloatType y = t * 2;

            if (y >= 2 - pulseWidth) {
                y = (y - 2) / pulseWidth;
            } else if (y >= pulseWidth) {
                y = 1 - (y - pulseWidth) / (1 - pulseWidth);
            } else {
                y /= pulseWidth;
            }

            y += freqInSecondsPerSample / (pulseWidth - pulseWidth * pulseWidth) * (blamp(t1, freqInSecondsPerSample) - blamp(t2, freqInSecondsPerSample));

            return amplitude * y;
        }

        void tri2(size_t n, DspFloatType *phase, DspFloatType * out) {
            #pragma omp simd
            for(size_t i = 0; i < n; i++) {
                 DspFloatType pulseWidth = std::fmax(0.0001, std::fmin(0.9999, this->pulseWidth));

                DspFloatType t1 = phase[i] + 0.5 * pulseWidth;
                t1 -= bitwiseOrZero(t1);

                DspFloatType t2 = phase[i] + 1 - 0.5 * pulseWidth;
                t2 -= bitwiseOrZero(t2);

                DspFloatType y = phase[i] * 2;

                if (y >= 2 - pulseWidth) {
                    y = (y - 2) / pulseWidth;
                } else if (y >= pulseWidth) {
                    y = 1 - (y - pulseWidth) / (1 - pulseWidth);
                } else {
                    y /= pulseWidth;
                }

                y += freqInSecondsPerSample / (pulseWidth - pulseWidth * pulseWidth) * (blamp(t1, freqInSecondsPerSample) - blamp(t2, freqInSecondsPerSample));

                out[i] = amplitude * y;
            }
        }
        DspFloatType trip() {
            DspFloatType t1 = t + 0.75 + 0.5 * pulseWidth;
            t1 -= bitwiseOrZero(t1);

            DspFloatType y;
            if (t1 >= pulseWidth) {
                y = -pulseWidth;
            } else {
                y = 4 * t1;
                y = (y >= 2 * pulseWidth ? 4 - y / pulseWidth - pulseWidth : y / pulseWidth - pulseWidth);
            }

            if (pulseWidth > 0) {
                DspFloatType t2 = t1 + 1 - 0.5 * pulseWidth;
                t2 -= bitwiseOrZero(t2);

                DspFloatType t3 = t1 + 1 - pulseWidth;
                t3 -= bitwiseOrZero(t3);
                y += 2 * freqInSecondsPerSample / pulseWidth * (blamp(t1, freqInSecondsPerSample) - 2 * blamp(t2, freqInSecondsPerSample) + blamp(t3, freqInSecondsPerSample));
            }
            return amplitude * y;
        }
        void trip(size_t n, DspFloatType *phase, DspFloatType * out) {
            #pragma omp simd
            for(size_t i = 0; i < n; i++) {
                DspFloatType t1 = phase[i] + 0.75 + 0.5 * pulseWidth;
                t1 -= bitwiseOrZero(t1);

                DspFloatType y;
                if (t1 >= pulseWidth) {
                    y = -pulseWidth;
                } else {
                    y = 4 * t1;
                    y = (y >= 2 * pulseWidth ? 4 - y / pulseWidth - pulseWidth : y / pulseWidth - pulseWidth);
                }

                if (pulseWidth > 0) {
                    DspFloatType t2 = t1 + 1 - 0.5 * pulseWidth;
                    t2 -= bitwiseOrZero(t2);

                    DspFloatType t3 = t1 + 1 - pulseWidth;
                    t3 -= bitwiseOrZero(t3);
                    y += 2 * freqInSecondsPerSample / pulseWidth * (blamp(t1, freqInSecondsPerSample) - 2 * blamp(t2, freqInSecondsPerSample) + blamp(t3, freqInSecondsPerSample));
                }
                out[i] = amplitude * y;
            }
        }

        DspFloatType trap() {
            DspFloatType y = 4 * t;
            if (y >= 3) {
                y -= 4;
            } else if (y > 1) {
                y = 2 - y;
            }
            y = std::fmax(-1, std::fmin(1, 2 * y));

            DspFloatType t1 = t + 0.125;
            t1 -= bitwiseOrZero(t1);

            DspFloatType t2 = t1 + 0.5;
            t2 -= bitwiseOrZero(t2);

            // Triangle #1
            y += 4 * freqInSecondsPerSample * (blamp(t1, freqInSecondsPerSample) - blamp(t2, freqInSecondsPerSample));

            t1 = t + 0.375;
            t1 -= bitwiseOrZero(t1);

            t2 = t1 + 0.5;
            t2 -= bitwiseOrZero(t2);

            // Triangle #2
            y += 4 * freqInSecondsPerSample * (blamp(t1, freqInSecondsPerSample) - blamp(t2, freqInSecondsPerSample));

            return amplitude * y;
        }
        void trap(size_t n, DspFloatType *phase, DspFloatType * out) {
            #pragma omp simd
            for(size_t i = 0; i < n; i++) {
                DspFloatType y = 4 * phase[i];
                if (y >= 3) {
                    y -= 4;
                } else if (y > 1) {
                    y = 2 - y;
                }
                y = std::fmax(-1, std::fmin(1, 2 * y));

                DspFloatType t1 = phase[i] + 0.125;
                t1 -= bitwiseOrZero(t1);

                DspFloatType t2 = t1 + 0.5;
                t2 -= bitwiseOrZero(t2);

                // Triangle #1
                y += 4 * freqInSecondsPerSample * (blamp(t1, freqInSecondsPerSample) - blamp(t2, freqInSecondsPerSample));

                t1 = phase[i] + 0.375;
                t1 -= bitwiseOrZero(t1);

                t2 = t1 + 0.5;
                t2 -= bitwiseOrZero(t2);

                // Triangle #2
                y += 4 * freqInSecondsPerSample * (blamp(t1, freqInSecondsPerSample) - blamp(t2, freqInSecondsPerSample));

                out[i] = amplitude * y;
            }
        }

        DspFloatType trap2() {
            DspFloatType pulseWidth = std::fmin(0.9999, this->pulseWidth);
            DspFloatType scale = 1 / (1 - pulseWidth);

            DspFloatType y = 4 * t;
            if (y >= 3) {
                y -= 4;
            } else if (y > 1) {
                y = 2 - y;
            }
            y = std::fmax(-1, std::fmin(1, scale * y));

            DspFloatType t1 = t + 0.25 - 0.25 * pulseWidth;
            t1 -= bitwiseOrZero(t1);

            DspFloatType t2 = t1 + 0.5;
            t2 -= bitwiseOrZero(t2);

            // Triangle #1
            y += scale * 2 * freqInSecondsPerSample * (blamp(t1, freqInSecondsPerSample) - blamp(t2, freqInSecondsPerSample));

            t1 = t + 0.25 + 0.25 * pulseWidth;
            t1 -= bitwiseOrZero(t1);

            t2 = t1 + 0.5;
            t2 -= bitwiseOrZero(t2);

            // Triangle #2
            y += scale * 2 * freqInSecondsPerSample * (blamp(t1, freqInSecondsPerSample) - blamp(t2, freqInSecondsPerSample));

            return amplitude * y;
        }
        void trap2(size_t n, DspFloatType *phase, DspFloatType * out) {
            #pragma omp simd
            for(size_t i = 0; i < n; i++) {
                DspFloatType pulseWidth = std::fmin(0.9999, this->pulseWidth);
                DspFloatType scale = 1 / (1 - pulseWidth);

                DspFloatType y = 4 * phase[i];
                if (y >= 3) {
                    y -= 4;
                } else if (y > 1) {
                    y = 2 - y;
                }
                y = std::fmax(-1, std::fmin(1, scale * y));

                DspFloatType t1 = t + 0.25 - 0.25 * pulseWidth;
                t1 -= bitwiseOrZero(t1);

                DspFloatType t2 = t1 + 0.5;
                t2 -= bitwiseOrZero(t2);

                // Triangle #1
                y += scale * 2 * freqInSecondsPerSample * (blamp(t1, freqInSecondsPerSample) - blamp(t2, freqInSecondsPerSample));

                t1 = phase[i] + 0.25 + 0.25 * pulseWidth;
                t1 -= bitwiseOrZero(t1);

                t2 = t1 + 0.5;
                t2 -= bitwiseOrZero(t2);

                // Triangle #2
                y += scale * 2 * freqInSecondsPerSample * (blamp(t1, freqInSecondsPerSample) - blamp(t2, freqInSecondsPerSample));

                out[i] = amplitude * y;
            }
        }

        DspFloatType sqr() const {
            DspFloatType t2 = t + 0.5;
            t2 -= bitwiseOrZero(t2);

            DspFloatType y = t < 0.5 ? 1 : -1;
            y += blep(t, freqInSecondsPerSample) - blep(t2, freqInSecondsPerSample);

            return amplitude * y;
        }
        void sqr(size_t n, DspFloatType *phase, DspFloatType * out) {
            #pragma omp simd
            for(size_t i = 0; i < n; i++) {
                DspFloatType t2 = phase[i] + 0.5;
                t2 -= bitwiseOrZero(t2);

                DspFloatType y = phase[i] < 0.5 ? 1 : -1;
                y += blep(t, freqInSecondsPerSample) - blep(t2, freqInSecondsPerSample);

                out[i] = amplitude * y;
            }
        }

        DspFloatType sqr2() {
            DspFloatType t1 = t + 0.875 + 0.25 * (pulseWidth - 0.5);
            t1 -= bitwiseOrZero(t1);

            DspFloatType t2 = t + 0.375 + 0.25 * (pulseWidth - 0.5);
            t2 -= bitwiseOrZero(t2);

            // Square #1
            DspFloatType y = t1 < 0.5 ? 1 : -1;

            y += blep(t1, freqInSecondsPerSample) - blep(t2, freqInSecondsPerSample);

            t1 += 0.5 * (1 - pulseWidth);
            t1 -= bitwiseOrZero(t1);

            t2 += 0.5 * (1 - pulseWidth);
            t2 -= bitwiseOrZero(t2);

            // Square #2
            y += t1 < 0.5 ? 1 : -1;

            y += blep(t1, freqInSecondsPerSample) - blep(t2, freqInSecondsPerSample);

            return amplitude * 0.5 * y;
        }   
        void sqr2(size_t n, DspFloatType *phase, DspFloatType * out) {
            
            for(size_t i = 0; i < n; i++) {
                DspFloatType t1 = phase[i] + 0.875 + 0.25 * (pulseWidth - 0.5);
                t1 -= bitwiseOrZero(t1);

                DspFloatType t2 = phase[i] + 0.375 + 0.25 * (pulseWidth - 0.5);
                t2 -= bitwiseOrZero(t2);

                // Square #1
                DspFloatType y = t1 < 0.5 ? 1 : -1;

                y += blep(t1, freqInSecondsPerSample) - blep(t2, freqInSecondsPerSample);

                t1 += 0.5 * (1 - pulseWidth);
                t1 -= bitwiseOrZero(t1);

                t2 += 0.5 * (1 - pulseWidth);
                t2 -= bitwiseOrZero(t2);

                // Square #2
                y += t1 < 0.5 ? 1 : -1;

                y += blep(t1, freqInSecondsPerSample) - blep(t2, freqInSecondsPerSample);

                out[i] =  amplitude * 0.5 * y;
            }
        }

        DspFloatType rect() {
            DspFloatType t2 = t + 1 - pulseWidth;
            t2 -= bitwiseOrZero(t2);

            DspFloatType y = -2 * pulseWidth;
            if (t < pulseWidth) {
                y += 2;
            }

            y += blep(t, freqInSecondsPerSample) - blep(t2, freqInSecondsPerSample);

            return amplitude * y;
        }
        void rect(size_t n, DspFloatType *phase, DspFloatType * out) {
            #pragma omp simd
            for(size_t i = 0; i < n; i++) {
                DspFloatType t2 = phase[i] + 1 - pulseWidth;
                t2 -= bitwiseOrZero(t2);

                DspFloatType y = -2 * pulseWidth;
                if (phase[i] < pulseWidth) {
                    y += 2;
                }

                y += blep(phase[i], freqInSecondsPerSample) - blep(t2, freqInSecondsPerSample);

                out[i] =  amplitude * y;
            }
        }

        DspFloatType saw() {
            DspFloatType _t = t + 0.5;
            _t -= bitwiseOrZero(_t);

            DspFloatType y = 2 * _t - 1;
            y -= blep(_t, freqInSecondsPerSample);

            return (amplitude * y);
        }
        void saw(size_t n, DspFloatType *phase, DspFloatType * out) {
            #pragma omp simd
            for(size_t i = 0; i < n; i++) {        
                DspFloatType _t = phase[i] + 0.5;
                _t -= bitwiseOrZero(_t);

                DspFloatType y = 2 * _t - 1;
                y -= blep(_t, freqInSecondsPerSample);

                out[i]= (amplitude * y);
            }
        }
        DspFloatType ramp() {
            DspFloatType _t = t;
            _t -= bitwiseOrZero(_t);

            DspFloatType y = 1 - 2 * _t;
            y += blep(_t, freqInSecondsPerSample);

            return amplitude * y;
        }
        void ramp(size_t n, DspFloatType *phase, DspFloatType * out) {
            #pragma omp simd
            for(size_t i = 0; i < n; i++) {        
                DspFloatType _t = phase[i];
                _t -= bitwiseOrZero(_t);
                DspFloatType y = 1 - 2 * _t;
                y += blep(_t, freqInSecondsPerSample);
                out[i] = amplitude * y;
            }
        }

        void ProcessBlock(size_t n, DspFloatType * out) {
            ProcessSIMD(n,out);
        }
    };

    enum {
        SAWTOOTH,
        SQUARE,
        TRIANGLE,
        SINE,
        // other waveforms vary to the oscillator 
        // use the oscillator to change it
    };

    // PolyBLEP
    // Blip
    // minBLEP
    // DPW
    struct VCO
    {
        VCOPolyBLEP polyblep;

        VCO(DspFloatType sr) {
            polyblep.init(sr);
        }
        void setWaveForm(int type) {            
            switch(type) {
                case SAWTOOTH:
                    polyblep.setWaveform(VCOPolyBLEP::Waveform::SAWTOOTH); 
                    break;
                case SQUARE:
                    polyblep.setWaveform(VCOPolyBLEP::Waveform::SQUARE); 
                    break;
                case TRIANGLE:
                    polyblep.setWaveform(VCOPolyBLEP::Waveform::TRIANGLE); 
                    break;
                case SINE:
                    polyblep.setWaveform(VCOPolyBLEP::Waveform::SINE); 
                    break;        
            }
        }
        void setFrequency(DspFloatType f) {
            polyblep.setFrequency(f);
        }
        void setSampleRate(DspFloatType sr) {
            polyblep.setSampleRate(sr);
        }
        void setDuty(DspFloatType d) {
            polyblep.setPulseWidth(d);
        }

        DspFloatType Tick(DspFloatType I=1, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1)
        {
            return polyblep.Tick(I,A,X,Y);
        }
        void ProcessSIMD(size_t n, DspFloatType * out) {
            polyblep.ProcessSIMD(n,out);
        }
    };
}
