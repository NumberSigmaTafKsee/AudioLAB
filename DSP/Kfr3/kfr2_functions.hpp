#pragma once

///////////////////////////////////////////////////////////////
// Functions
///////////////////////////////////////////////////////////////

namespace KfrDSP1
{
    template<typename T>
    kfr::univector<T> sinewave(size_t n, T freq, T sample_rate, T phase=(T)0) {        
        return DSP::sinewave(n,freq,sample_rate,phase);        
    }
    template<typename T>
    kfr::univector<T> squarewave(size_t n, T freq, T sample_rate, T phase=(T)0) {        
        return DSP::squarewave(n,freq,sample_rate,phase);        
    }
    template<typename T>
    kfr::univector<T> trianglewave(size_t n, T freq, T sample_rate, T phase=(T)0) {        
        return DSP::trianglewave(n,freq,sample_rate,phase);        
    }
    template<typename T>
    kfr::univector<T> sawtoothewave(size_t n, T freq, T sample_rate, T phase=(T)0) {        
        return DSP::squarewave(n,freq,sample_rate,phase);        
    }    

    template<typename T>
    struct FunctionGenerator
    {
        T phase,inc,f,sr;

        FunctionGenerator(T Fs = 44100.0f)
        {
            phase = 0;
            f     = 440.0f;
            sr    = Fs;
            inc   = f/sr;            
        }
        void Inc()
        {
            phase = fmod(phase + inc, 2*M_PI);
        }
        void rawsine()
        {
            return kfr::rawsine(phase);
        }
        void sine() {
            return kfr::sine(phase);
        }
        void sinenorm() {
            return kfr::sinenorm(phase);
        }
        void rawsquare() {
            return kfr::rawsquare(phase);
        }
        void square() {
            return kfr::square(phase);
        }
        void squarenorm() {
            return kfr::squarenorm(phase);
        }
        void rawtriangle() {
            return kfr::rawtriangle(phase);
        }
        void triangle() {
            return kfr::triangle(phase);
        }
        void trianglenorm() {
            return kfr::trianglenorm(phase);
        }
        void rawsawtooth() {
            return kfr::rawsawtooth(phase);
        }
        void sawtooth() {
            return kfr::sawtooth(phase);
        }
        void sawtoothnorm() {
            return kfr::sawtoothnorm(phase);
        }
        void isawtooth() {
            return kfr::isawtooth(phase);
        }
        void isawtoothnorm() {
            return kfr::isawtoothnorm(phase);
        }
    };
}