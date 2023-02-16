#include <algorithm>

#include "effects.hpp"

Delay::Delay(
        const std::vector<size_t> &delays,
        const std::vector<double> &gains
):delays_(delays), gains_(gains){
    size_t max_delay = 0;
    for(auto &&i: delays){
        max_delay = std::max(i,max_delay);
    }
    older_data_.resize(max_delay,0);
    newer_data_.resize(max_delay,0);
}

void Delay::operator()(FFTwrapper &fft){
    fft.c2r_transform();
    operator()(fft.get_real_ptr(), fft.get_real_ptr(), fft.get_real_c());
}

//Only adds the delayed samples to output, does not copy original signal
//Input and output can be the same
void Delay::operator()(double *input, double* output, size_t size){
    //Input and output can be the same pointer, for simplicity, copy the samples
    std::vector<double> input_copy(input, input+size);

    //Add the delay from the second previous window, and zero older_data_
    size_t min = std::min(older_data_.size(), size);
    for(size_t i=0; i < min; ++i){
        output[i] += older_data_[i];
        older_data_[i] = 0;
    }
    if(size < older_data_.size())
        std::rotate(older_data_.begin(), older_data_.begin()+size, older_data_.end());
    
    std::swap(older_data_, newer_data_);

    size_t delaysSize = delays_.size();
    for(size_t i=0; i<size; ++i){
        for(size_t j=0; j< delaysSize; ++j){
            if(i<delays_[j])
                continue;
            output[i] += input_copy[i-delays_[j]] * gains_[j]; 
        }
    }
    
    size_t n_s = newer_data_.size();
    //Save the delayed samples that would lie out of the window
    for(size_t i=size; i<size + n_s; ++i){
        for(size_t j=0; j<delaysSize; ++j){
            if(i-delays_[j]>=size)
                continue;
            newer_data_[i-size] += input_copy[i-delays_[j]] * gains_[j];
        }
    }
}

Delay Delay::create_linear(size_t delay, size_t count){
    std::vector<size_t> delays;
    std::vector<double> gains;
    
    for(size_t i=0; i<count; ++i){
        delays.push_back((i+1)*delay);
    }
    for(size_t i=0; i<count; ++i){
        gains.push_back(1-((i+1.0)/(count+1.0)));
    }
    return Delay(delays, gains);
}
Delay Delay::create_exp(size_t delay, size_t count, double base){
    std::vector<size_t> delays;
    std::vector<double> gains;
    double gain;
    gain = base;
    for(size_t i=0; i<count; ++i){
        delays.push_back((i+1)*delay);
    }
    for(size_t i=0; i<count; ++i){
        gains.push_back(gain);
        gain *= base;
    }
    return Delay(delays, gains);
}