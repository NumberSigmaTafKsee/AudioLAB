#include <cstdio>
#include <cstring>
#include <sndfile.hh>
#include <string>
#include <algorithm>
#include <stdexcept>

#include "io.hpp"

IO::IO(const std::string &input_path, const std::string &output_path){
    info.format = 0;
    i_file_ = sf_open(input_path.c_str(), SFM_READ, &info);
    o_file_ = sf_open(output_path.c_str(), SFM_WRITE, &info);
}

IO::IO(IO &&io){
    i_file_ = io.i_file_;
    o_file_ = io.o_file_;
    info = io.info;
    io.i_file_ = nullptr;
    io.o_file_ = nullptr;
}
IO &IO::operator=(IO &&io){
    info = io.info;
    std::swap(i_file_, io.i_file_);
    std::swap(o_file_, io.o_file_);
    return *this;
}

IO::~IO(){
    if(i_file_)
        sf_close(i_file_);
    if(o_file_)
        sf_close(o_file_);
}
sf_count_t IO::read_n(double *dest, size_t n) const{
   sf_count_t count = sf_read_double(i_file_, dest, n);
   //Zero the data that the read didnt overwrite
   for(size_t i = count; i<n; ++i)
      dest[i] = 0;
   return count;
}
sf_count_t IO::write_n(double *data, size_t n) const{
    return sf_write_double(o_file_, data, n);
}

WriteOverlap::WriteOverlap(size_t size): size_(size){
    cache_.resize(size_,0);
}
//data buffer must be at least size_ long
void WriteOverlap::operator()(size_t shift, double *data, const IO &io){
    for(size_t i=0; i<size_; ++i){
        cache_[i] += data[i];
    }
    io.write_n(cache_.data(),shift);

    std::rotate(cache_.begin(),cache_.begin()+shift, cache_.end());
    for(size_t i=cache_.size()-shift; i<cache_.size(); ++i)
        cache_[i] = 0;
}