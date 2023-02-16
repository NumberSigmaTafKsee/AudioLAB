#ifndef IO_HPP_
#define IO_HPP_

#include <cstdio>
#include <cstring>
#include <sndfile.hh>
#include <string>
#include <vector>

//Supports only one channel
class IO{
public:
    IO(const std::string &input_path, const std::string &output_path);
    
    //Delete copying
    IO(const IO &io) = delete;
    IO &operator=(const IO &io) = delete;

    //Move
    IO(IO &&io);
    IO &operator=(IO &&io);

    ~IO();
    sf_count_t read_n(double *dest, size_t n) const;
    sf_count_t write_n(double *data, size_t n) const;
    SF_INFO get_info() const{return info;}
private:
    SNDFILE *i_file_, *o_file_;
    SF_INFO info;
};

///Writes to file with overlap
class WriteOverlap{
public:
    ///Default constructor
    /*
     * @param size the size of the window to be overlaped
     * @param max_count the maximum count of datapoints to be written at once
     */
    WriteOverlap(size_t size);
    ///Adds the data to cache and then shifts the cache by shift datapoints. The shifted data is written to a file via io
    void operator()(size_t shift, double *data, const IO &io);
private:
    size_t size_;
    std::vector<double> cache_;
};

#endif