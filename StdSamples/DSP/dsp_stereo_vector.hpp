#pragma once

namespace AudioDSP::Std
{
        template<typename T>
    StereoVector<T> stereo(const sample_vector<T> & left, const sample_vector<T> & right) {
        StereoVector<T> r(left.size());
        for(size_t i = 0; i < left.size(); i++)
        {
            r[0][i] = left[i];
            r[1][i] = right[i];
        }
    }

    template<typename T>
    sample_vector<T> merge(const sample_vector<T> & left, const sample_vector<T> & right) {
        sample_vector<T> r(left.size()*2);
        size_t x = 0;
        for(size_t i = 0; i < left.size(); i++)
        {
            r[x++] = left[i];
            r[x++] = right[i];
        }
        return r;
    }

    
    template<typename T>
    StereoVector<T> pan(const sample_vector<T> & left, const sample_vector<T> & right, T amt) {
        StereoVector<T> r(left.size());
        T pan_map = ((amt+1)/2.0) * (M_PI/2.0);
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
        for(size_t i = 0; i < left.size(); i++)
        {
            r[0][i] = left[i] * pleft;
            r[1][i] = right[i] * pright;
        }
        return r;
    }
}