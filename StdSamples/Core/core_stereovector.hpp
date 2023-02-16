#pragma once

namespace Casino
{
    template<typename T>
    struct StereoVector
    {
        sample_vector<T> samples[2];

        StereoVector(size_t n) {
            for(size_t i = 0; i < n; i++) samples[i].resize(n);
        }

        sample_vector<T>& operator[](size_t channel) { return samples[channel]; }

        sample_vector<T> get_left_channel() { return samples[0]; }
        sample_vector<T> get_right_channel() { return samples[1]; }

        void set_left_channel(sample_vector<T> & v) { samples[0] = v; }
        void set_right_channel(sample_vector<T> & v) { samples[1] = v; }
    };
    template<typename T>
    StereoVector<T> stereo(const sample_vector<T> & left, const sample_vector<T> & right) {
        StereoVector<T> r(left.size());
        #pragma omp simd
        for(size_t i = 0; i < left.size(); i++)
        {
            r[0][i] = left[i];
            r[1][i] = right[i];
        }
    }
    template<typename T>
    StereoVector<T> pan(const sample_vector<T> & left, const sample_vector<T> & right, T amt) {
        StereoVector<T> r(left.size());
        T pan_map = ((amt+1)/2.0) * (M_PI/2.0);
        #pragma omp simd
        for(size_t i = 0; i < left.size(); i++)
        {
            r[0][i] = left[i] * sin(pan_map);
            r[1][i] = right[i] * cos(pan_map);
        }
        return r;
    }
    template<typename T>
    StereoVector<T> constant_power_pan(const sample_vector<T> & left, const sample_vector<T> & right, T pos) {
        StereoVector<T> r(left.size());        
        const T piover2 = 4.0*std::atan(1.0)*0.5;
        const T root2over2 = std::sqrt(2.0)*0.5;
        T thispos = pos * piover2;
        T angle   = thispos * 0.5;
        T pleft   = root2over2 * (std::cos(angle) - std::sin(angle));
        T pright  = root2over2 * (std::cos(angle) + std::sin(angle));
        #pragma omp simd
        for(size_t i = 0; i < left.size(); i++)
        {
            r[0][i] = left[i] * pleft;
            r[1][i] = right[i] * pright;
        }
        return r;
    }
    template<typename T>
    sample_vector<T> merge(const sample_vector<T> & left, const sample_vector<T> & right) {
        sample_vector<T> r(left.size()*2);
        size_t x = 0;
        #pragma omp simd
        for(size_t i = 0; i < left.size(); i++)
        {
            r[x++] = left[i];
            r[x++] = right[i];
        }
    }
}