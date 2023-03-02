//
//  Sigmoid.hpp
//  Math
//
//  Copyright © 2015-2016 Dsperados (info@dsperados.com). All rights reserved.
//  Licensed under the BSD 3-clause license.
//

#ifndef DSPERADOS_MATH_SIGMOID_HPP
#define DSPERADOS_MATH_SIGMOID_HPP

#include <cmath>
#include <stdexcept>

namespace math
{
    //! Normalized sigmoid function
    template <typename T>
    T sigmoid(const T& x, double negativeFactor, double positiveFactor)
    {
        if (negativeFactor <= 0 || positiveFactor <=0)
            throw std::runtime_error("Factor <= 0");
        
        auto yAtX1Negative = negativeFactor / (1 + std::abs(negativeFactor));
        auto yAtX1Positive = positiveFactor / (1 + std::abs(positiveFactor));
        
        if (x > 0)
            return ((x * positiveFactor) / (1 + std::abs(x * positiveFactor))) / yAtX1Positive;
        else if (x < 0)
            return ((x * negativeFactor) / (1 + std::abs(x * negativeFactor))) / yAtX1Negative;
        else
            return 0;
    }
    
    //! Normalized sigmoid function using tan
    template <typename T>
    T sigmoidTan(const T& x, double negativeFactor, double positiveFactor)
    {
        if (negativeFactor <= 0 || positiveFactor <=0)
            throw std::runtime_error("Factor <= 0");
        
        if (x > 0)
            return std::atan(x * positiveFactor) * (1 / std::atan(positiveFactor));
        else if (x < 0)
            return atan(x * negativeFactor) * (1 / atan(negativeFactor));
        else
            return 0;
    }
    
    //! Normalized sigmoid function using exp
    template <typename T>
    T sigmoidExp(const T& x, double negativeFactor, double positiveFactor)
    {
        if (negativeFactor <= 0 || positiveFactor <=0)
            throw std::runtime_error("Factor <= 0");
        
        if (x > 0)
            return (1 - std::exp(-x * positiveFactor)) / (1 - std::exp(-positiveFactor)) ;
        else if (x < 0)
            return (-1 + std::exp(x * negativeFactor)) / (1 - std::exp(-negativeFactor));
        else
            return 0;
    }
}

#endif /* DSPERADOS_MATH_SIGMOID_HPP */
