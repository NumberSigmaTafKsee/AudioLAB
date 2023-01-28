//
//  statistics.hpp
//  Math
//
//  Copyright © 2015-2016 Dsperados (info@dsperados.com). All rights reserved.
//  Licensed under the BSD 3-clause license.
//

#ifndef DSPERADOS_MATH_STATISTICS_HPP
#define DSPERADOS_MATH_STATISTICS_HPP

#include <cmath>
#include <cstddef>
#include <numeric>
#include <vector>

namespace math
{
    //! Calculate the mean
    template <typename T, typename Iterator>
    T mean(Iterator begin, Iterator end)
    {
        return std::accumulate(begin, end, T(0)) / static_cast<T>(std::distance(begin, end));
    }
    
    //! Calculate the mean square
    template <typename T, typename Iterator>
    T meanSquare(Iterator begin, Iterator end)
    {
        return std::accumulate(begin, end, T(0), [](const T& accumulatedValue, const T& currentValue) { return accumulatedValue + currentValue * currentValue; }) / static_cast<T>(std::distance(begin, end));
    }
    
    //! Calculate the root mean square
    template <typename T, typename Iterator>
    T rootMeanSquare(Iterator begin, Iterator end)
    {
        return std::sqrt(meanSquare<T>(begin, end));
    }
}

#endif
