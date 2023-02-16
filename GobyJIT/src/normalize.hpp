//
//  Normalize.hpp
//  Math
//
//  Copyright © 2015-2016 Dsperados (info@dsperados.com). All rights reserved.
//  Licensed under the BSD 3-clause license.
//

#ifndef MATH_NORMALIZE_HPP
#define MATH_NORMALIZE_HPP

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>

#include "analysis.hpp"

namespace math
{
    //! Normalize an area so the integral of the signal equals one
    template <typename InputIterator, typename OutputIterator>
    void normalizeArea(InputIterator inBegin, InputIterator inEnd, OutputIterator outBegin)
    {
        auto integral = std::accumulate(inBegin, inEnd, typename InputIterator::value_type{0});
        
        if (!integral)
            throw std::runtime_error("area equals zero");
            
        std::transform(inBegin, inEnd, outBegin, [&](const auto& x){ return x / integral; });
    }
    
    //! Normalize a range
    template <typename InputIterator, typename OutputIterator>
    void normalize(InputIterator inBegin, InputIterator inEnd, OutputIterator outBegin)
    {
        const auto absoluteExtrema = std::abs(*findExtrema(inBegin, inEnd));
        const auto factor = 1.0 / absoluteExtrema;
        std::transform(inBegin, inEnd, outBegin, [&](const auto& x){ return x * factor; });
    }
}

#endif
