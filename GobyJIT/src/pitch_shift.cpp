#include <cmath>
#include <cstdint>
#include <vector>
#include <tuple>

#include "effects.hpp"
//in.size() has to be equal to out_fft.get_complex_c()
void PitchShift::shift_complex(const std::vector<std::tuple<double,double>> &in, FFTwrapper &out_fft, double ammount){
    out_fft.is_complex_state_ = true;
    //Zero out the complex part
    for(size_t i=0; i<out_fft.get_complex_c(); ++i){
        out_fft.get_complex_ptr()[i][0]=0;
        out_fft.get_complex_ptr()[i][1]=0;
    }
    //Do the pitch shift
    for(size_t i=0; i<out_fft.get_complex_c(); ++i){
        double x = i*ammount;
        size_t new_id = (size_t)std::floor(x+0.5);

        if(x<0 || new_id >= out_fft.get_complex_c())
            break;

        out_fft.get_complex_ptr()[new_id][0] += std::get<0>(in[i]);
        out_fft.get_complex_ptr()[new_id][1] += std::get<1>(in[i]);
    }
}

void PitchShift::operator()(FFTwrapper &fft){
    std::vector<std::tuple<double,double>> spectrum;
    fft.r2c_transform();
    for(size_t j=0; j<fft.get_complex_c(); ++j){
        spectrum.push_back(std::make_tuple(fft.get_complex_ptr()[j][0],fft.get_complex_ptr()[j][1]));
    }
    shift_complex(spectrum, fft, ammount_);
}

