#pragma once

namespace SoundAlchemy::LFO
{

    template<typename T>
    struct TPhasor : TObject<T>
    {
        T phase,inc,f,fs;

        TPhasor(T Fs) : TObject<T>()
        {
            fs = Fs;
            f  = 0.2f;
            inc = f/fs;
            phase = 0.0f;
        }
        void setFrequency(T freq)
        {
            f = freq;
            inc = f/fs;
        }
        T Tick() {
            T out = phase;
            phase = fmod(phase + inc,1);
            return out;
        }
    };

    template<typename T>
    T sinwave(TPhasor<T> & phasor) {
        return std::sin(2*M_PI*phasor.Tick());
    }
    template<typename T>
    T coswave(TPhasor<T> & phasor) {
        return std::cos(2*M_PI*phasor.Tick());
    }
    template<typename T>
    T sawwave(TPhasor<T> & phasor) {
        return phasor.Tick();
    }
    template<typename T>
    T rampwave(TPhasor<T> & phasor) {
        return 1.0-phasor.Tick();
    }
    template<typename T>
    T squarewave(TPhasor<T> & phasor, T duty = (T)0.5) {
        T r = phasor.Tick();
        if(r < duty) return -1;
        return 1;
    }
    template<typename T>
    T trianglewave(TPhasor<T> & phasor) {
        return std::asin(std::sin(2*M_PI*phasor.Tick()))/1.5;
    }

    template<typename T>
    struct TLFO
    {
        TPhasor<T> phasor;
        T last;
        T duty = T(0.5);
    
        enum {
            POSITIVE,
            NEGATIVE,
            BIPOLAR,
            HALF,
            FULL,
        }
        polarity = POSITIVE;

        enum {
            SIN,
            COS,
            TAN,
            SAW,
            RAMP,
            SQUARE,
            TRIANGLE,
            NOISE,
        }
        type = SIN;

        TLFO(T Fs) : phasor(Fs)
        {
            last = T(0);
        }
        T Tick()
        {
            last = phasor.Tick();
            switch(type) {
                case SIN: return sin();
                case COS: return cos();
                case TAN: return tan();
                case SAW: return saw();
                case RAMP: return ramp();
                case SQUARE: return square();
                case TRIANGLE: return triangle();
            }
        }
        T sin() {
            T x = std::sin(2*M_PI*last);
            if(polarity == POSITIVE) return 0.5+0.5*x;
            if(polarity == NEGATIVE) return -(0.5+0.5*x);
            if(polarity == BIPOLAR)  return x;
            if(polarity == HALF) 
            {                
                if(x < 0) return 0;
                return x;
            }
            if(polarity == FULL) {
                return std::abs(x);
            }
        }
        T cos() {
            T x = std::cos(2*M_PI*last);
            if(polarity == POSITIVE) return 0.5+0.5*x;
            if(polarity == NEGATIVE) return -(0.5+0.5*x);
            if(polarity == BIPOLAR)  return x;
            if(polarity == HALF) 
            {                
                if(x < 0) return 0;
                return x;
            }
            if(polarity == FULL) {
                return std::abs(x);
            }
        }
        T tan() {
            T x =std::tan(M_PI/2*0.996*last);
            if(polarity == POSITIVE) return 0.5+0.5*x;
            if(polarity == NEGATIVE) return -(0.5+0.5*x);
            if(polarity == BIPOLAR)  return x;
            if(polarity == HALF) 
            {                
                if(x < 0) return 0;
                return x;
            }
            if(polarity == FULL) {
                return std::abs(x);
            }
            return x;
        }
        T saw() {
            T x = last;
            if(polarity == POSITIVE) return 0.5+0.5*x;
            if(polarity == NEGATIVE) return -(0.5+0.5*x);
            if(polarity == BIPOLAR)  return x;
            if(polarity == HALF) 
            {                
                if(x < 0) return 0;
                return x;
            }
            if(polarity == FULL) {
                return std::abs(x);
            }
            return x;
        }
        T ramp() {
            T x = 1-last;
            if(polarity == POSITIVE) return 0.5+0.5*x;
            if(polarity == NEGATIVE) return -(0.5+0.5*x);
            if(polarity == BIPOLAR)  return x;
            if(polarity == HALF) 
            {                
                if(x < 0) return 0;
                return x;
            }
            if(polarity == FULL) {
                return std::abs(x);
            }
            return x;
        }
        T square() {
            T x = last < duty? -1:1;
            if(polarity == POSITIVE) return 0.5+0.5*x;
            if(polarity == NEGATIVE) return -(0.5+0.5*x);
            if(polarity == BIPOLAR)  return x;
            if(polarity == HALF) 
            {                
                if(x < 0) return 0;
                return x;
            }
            if(polarity == FULL) {
                return std::abs(x);
            }
            return x;
        }
        T triangle()
        {
            T x = asin(sin(2*M_PI*last))/1.5;
            if(polarity == POSITIVE) return 0.5+0.5*x;
            if(polarity == NEGATIVE) return -(0.5+0.5*x);
            if(polarity == BIPOLAR)  return x;
            if(polarity == HALF) 
            {                
                if(x < 0) return 0;
                return x;
            }
            if(polarity == FULL) {
                return std::abs(x);
            }
            return x;
        }
        T noise()
        {
            return last * ((T)rand() / (T)RAND_MAX);
        }
    };
}

using LFO32 = SoundAlchemy::LFO::TLFO<float>;
using LFO64 = SoundAlchemy::LFO::TLFO<double>;
