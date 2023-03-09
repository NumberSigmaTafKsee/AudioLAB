#pragma once

#include "SoundObject.hpp"
#include "stdsamples_function_generator.hpp"
#include "stdsamples_functions.hpp"

namespace AudioDSP
{
    template<typename T>
    sample_vector<T> generate_noise(size_t n)
    {
        sample_vector<T> r(n);
        std::generate(r.begin(),r.end(),[]() { return Oscillators::Generators::function_noise(); });
        return r;
    }
    template<typename T>
    sample_vector<T> generate_sin(DspFloatType freq, DspFloatType sampleRate, size_t n)
    {
        sample_vector<T> r(n);
        Oscillators::Generators::SineGenerator sine(freq,sampleRate);
        std::generate(r.begin(),r.end(),[&sine]() { return sine(); });
        return r;
    }
    template<typename T>
    sample_vector<T> generate_cos(DspFloatType freq, DspFloatType sampleRate, size_t n)
    {
        sample_vector<T> r(n);
        Oscillators::Generators::CosGenerator cose(freq,sampleRate);
        std::generate(r.begin(),r.end(),[&cose]() { return cose(); });
        return r;
    }
    template<typename T>
    sample_vector<T> generate_tan(DspFloatType freq, DspFloatType sampleRate, size_t n)
    {
        sample_vector<T> r(n);
        Oscillators::Generators::TanGenerator tane(freq,sampleRate);
        std::generate(r.begin(),r.end(),[&tane]() { return tane(); });
        return r;
    }
    template<typename T>
    sample_vector<T> generate_phasor(DspFloatType freq, DspFloatType sampleRate, size_t n)
    {
        sample_vector<T> r(n);
        Oscillators::Generators::PhasorGenerator phasore(freq,sampleRate);
        std::generate(r.begin(),r.end(),[&phasore]() { return phasore(); });
        return r;
    }
    template<typename T>
    sample_vector<T> generate_square(DspFloatType freq, DspFloatType sampleRate, size_t n)
    {
        sample_vector<T> r(n);
        Oscillators::Generators::SineGenerator squaree(freq,sampleRate);
        std::generate(r.begin(),r.end(),[&squaree]() { return squaree(); });
        return r;
    }
    template<typename T>
    sample_vector<T> generate_saw(DspFloatType freq, DspFloatType sampleRate, size_t n)
    {
        sample_vector<T> r(n);
        Oscillators::Generators::SawGenerator sawe(freq,sampleRate);
        std::generate(r.begin(),r.end(),[&sawe]() { return sawe(); });
        return r;
    }
    template<typename T>
    sample_vector<T> generate_triangle(DspFloatType freq, DspFloatType sampleRate, size_t n)
    {
        sample_vector<T> r(n);
        Oscillators::Generators::TriangleGenerator trie(freq,sampleRate);
        std::generate(r.begin(),r.end(),[&trie]() { return trie(); });
        return r;
    }
    /*
    template<typename T>
    sample_vector<T> generate_function(Oscillators::Generators::FunctionGenerator::Type type, DspFloatType freq, DspFloatType sampleRate, size_t n)
    {
        sample_vector<T> r(n);
        Oscillators::Generators::FunctionGenerator func(type,freq,sampleRate);
        std::generate(r.begin(),r.end(),[&func]() { return func(); });
    }
    */
    template<typename T>
    sample_vector<T> oscillator(OscillatorProcessor & osc, size_t n)
    {
        sample_vector<T> r(n);        
        std::generate(r.begin(),r.end(),[&osc]() { return osc.Tick(); });
        return r;
    }
    template<typename T>
    sample_vector<T> generator(GeneratorProcessor & osc, size_t n)
    {
        sample_vector<T> r(n);        
        std::generate(r.begin(),r.end(),[&osc]() { return osc.Tick(); });
        return r;
    }
    template<typename T>
    sample_vector<T> filter(FilterProcessor & filt, size_t n, const T * in, T * out)
    {
        sample_vector<T> r(n);
        memcpy(r.data(),in,n*sizeof(T));        
        std::for_each(r.begin(),r.end(),[&filt](T & x) { x = filt.Tick(x); });
        return r;
    }
    template<typename T>
    sample_vector<T> function(FunctionProcessor & func, size_t n, const T * in, T * out)
    {
        sample_vector<T> r(n);
        memcpy(r.data(),in,n*sizeof(T));        
        std::for_each(r.begin(),r.end(),[&func](T & x) { x = func.Tick(x); });
        return r;
    }
}
