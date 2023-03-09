#pragma once

namespace AudioDSP
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
}