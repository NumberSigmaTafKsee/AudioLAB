#ifndef WINDOW_HPP_
#define WINDOW_HPP_
#include <vector>

#include "io.hpp"

///Basic implementation of the Hann window
/*
 * Supposed to be used for reading
 */
class Window{
public:
    Window(size_t size);
    ///Reads data and returns the windowed result, includint the 50% overlap
    /*
     * Returns the number of samples that can be written to file after processing.
     * Output buffer for write should be used in conjuncion.
     */
    size_t operator()(double *output, const IO &io);
    ///Writes the buffered data to file
private:
    size_t size_;
    bool next_end_ = false, first_ = true;
    std::vector<double> older_data_, newer_data_;
    double hann_val(size_t pos) const;
};

#endif