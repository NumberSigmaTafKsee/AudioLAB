//
//  sinusoid.hpp
//  Math
//
//  Copyright © 2015-2016 Dsperados (info@dsperados.com). All rights reserved.
//  Licensed under the BSD 3-clause license.
//

#ifndef DSPERADOS_MATH_SINUSOID_HPP
#define DSPERADOS_MATH_SINUSOID_HPP

#include <cmath>
#include <cstddef>
#include <vector>

#include "constants.hpp"

namespace math
{
    //! Generate a sine wave in a buffer
    template <typename T>
    std::vector<T> generateSineBuffer(std::size_t size, float order = 1, float amplitude = 1)
    {
        std::vector<T> buffer(size);
        for (auto i = 0; i < size; ++i)
            buffer[i] = sin(i / static_cast<double>(size) * TAU<double> * order) * amplitude;
        
        return buffer;
    }
    
    //! Generate a cosine wave in a buffer
    template <typename T>
    std::vector<T> generateCosineBuffer(std::size_t size, float order = 1, float amplitude = 1)
    {
        std::vector<T> buffer(size);
        for (auto i = 0; i < size; ++i)
            buffer[i] = cos((i / static_cast<double>(size)) * TAU<double> * order) * amplitude;
        
        return buffer;
    }
}

#endif
