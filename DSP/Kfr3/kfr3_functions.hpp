#pragma once

namespace kfr3
{
    template<typename T>
    sample_vector<T> sinewave(size_t n, T freq, T sample_rate, T phase=(T)0) {
        kfr::univector<T> r(n);                
        for(size_t i = 0; i < n; i++)
        {
            r[i] = kfr::sine(2*M_PI*phase);
            phase += freq/sample_rate;
            if(phase > (T)1.0) phase-=(T)1.0;
        }
        return make_univev(r);
    }
    template<typename T>
    sample_vector<T> squarewave(size_t n, T freq, T sample_rate, T phase=(T)0) {
        kfr::univector<T> r(n);
        for(size_t i = 0; i < n; i++)
        {
            r[i] = kfr::square(2*M_PI*phase);
            phase += freq/sample_rate;
            if(phase > (T)1.0) phase-=(T)1.0;
        }
        return make_univec(r);
    }
    template<typename T>
    kfr::univector<T> trianglewave(size_t n, T freq, T sample_rate, T phase=(T)0) {
        kfr::univector<T> r(n);
        for(size_t i = 0; i < n; i++)
        {
            r[i] = kfr::triangle(2*M_PI*phase);
            phase += freq/sample_rate;
            if(phase > (T)1.0) phase-=(T)1.0;
        }
        return r;
    }
    template<typename T>
    sample_vector<T> sawtoothwave(size_t n, T freq, T sample_rate, T phase=(T)0) {
        kfr::univector<T> r(n);
        for(size_t i = 0; i < n; i++)
        {
            r[i] = kfr::sawtooth(2*M_PI*phase);
            phase += freq/sample_rate;
            if(phase > (T)1.0) phase-=(T)1.0;
        }
        return make_univec(r);
    }
    template<typename T>
    sample_vector<T> generate_sin(size_t n, T start, T step) {
        kfr::univector<T> r(n);
        r = kfr::gen_sin(start,step);
        return make_univec(r);
    }

    template<typename T>
    sample_vector<T> generate_linear(size_t n, T start, T step) {
        kfr::univector<T> r(n);
        r = kfr::gen_linear(start,step);
        return make_univec(r);
    }

    template<typename T>
    sample_vector<T> generate_exp(size_t n, T start, T step) {
        kfr::univector<T> r(n);
        r = kfr::gen_exp(start,step);
        return make_univec(r);
    }

    template<typename T>
    sample_vector<T> generate_exp2(size_t n, T start, T step) {
        kfr::univector<T> r(n);
        r = kfr::gen_exp2(start,step);
        return make_univec(r);        
    }

    template<typename T>
    sample_vector<T> generate_cossin(size_t n, T start, T step) {
        kfr::univector<T> r(n);
        r = kfr::gen_cossin(start,step);
        return make_univec(r);
    }

    template<typename T>
    void deinterleave(sample_matrix<T> & out, const T * ptr, size_t size, size_t channels) {                        
        kfr::deinterleave(p.data(),ptr,channels,size);
    }
    
    template<typename T>
    sample_vector<T> interleave(const sample_matrix<T> & input) {        
        return make_univec(kfr::interleave(input));
    }

    template<typename T>
    sample_vector<T> generate_sequence(int n, int start, int step) {
        sample_vector<T> r(n);
        for(size_t i = 0; i < n; i++)
            r[i] = start + step*n;
        return r;
    }

    std::string dB_to_string(const double value, double min=-140.0f) { return kfr::dB_to_string(value,min); }
    std::string dB_to_utf8string(const double value, double min=-140.0f) { return kfr::dB_to_utf8string(value,min); }
}    