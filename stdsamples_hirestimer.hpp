#pragma once

#include <chrono>


class HiResTimer
{
private:
    std::chrono::time_point<std::chrono::high_resolution_clock> start,stop;
public:    
    HiResTimer() {
    }

    void Start() {
        start = std::chrono::high_resolution_clock::now();
    }
    double Stop() {
        stop  = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> diff = stop - start;
        return diff.count();
    }
};
    