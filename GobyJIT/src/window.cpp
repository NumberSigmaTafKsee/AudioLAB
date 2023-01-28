#include <vector>
#include <cmath>

#include "io.hpp"
#include "window.hpp"

const double PI = 3.141592653589793238463;

Window::Window(size_t size): size_(size){
    older_data_.resize(size / 2, 0);
    newer_data_.resize(size / 2 + size % 2, 0);
}
size_t Window::operator()(double *output, const IO &io){  
    if(next_end_)
        return 0;
    //Move data from new to old (and from old to new, but we are rewriting those)
    std::swap(older_data_, newer_data_);

    if(first_){
        //Initialize data to old on first read
        io.read_n(older_data_.data(), older_data_.size());
        first_ = false;
    }
    //Rewrite new
    size_t read_n = io.read_n(newer_data_.data(), newer_data_.size());
    
    for(size_t i=0; i<size_; ++i){
        if(i<older_data_.size()){
            output[i] = hann_val(i) * older_data_[i];
        }
        else{
            output[i] = hann_val(i) * newer_data_[i-older_data_.size()];
        }
    }
    
    if(read_n < newer_data_.size()){
        next_end_ = true;
        //Signalize to write everything to file
        return older_data_.size() + read_n;
    }
    return older_data_.size();
}
//Hann window with a bit of +-1 fiddling so we dont multiply the border values by 0
double Window::hann_val(size_t pos) const{
    double val = std::cos((PI/(size_+1))*(pos-((double)(size_-1)/2)));
    return val*val;
}