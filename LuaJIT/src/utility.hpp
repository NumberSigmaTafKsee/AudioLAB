//
//  utility.hpp
//  Math
//
//  Copyright © 2015-2016 Dsperados (info@dsperados.com). All rights reserved.
//  Licensed under the BSD 3-clause license.
//

#ifndef DSPERADOS_MATH_UTILITY_HPP
#define DSPERADOS_MATH_UTILITY_HPP

#include <algorithm>
#include <cassert>
#include <cmath>
#include <stdexcept>

#include "constants.hpp"

namespace math
{
    //! Clamp a number within a given range
    /*! @throw std::invalid_argument if max < min */
    template <typename T>
    std::enable_if_t<!std::is_integral<T>::value || std::is_signed<T>::value, T> clamp(const T& value, const T& min, const T& max)
    {
        if (max < min)
            throw std::invalid_argument("max < min");
        
        return value < min ? min : max < value ? max : value;
    }

    //! Clamp a number within a given range
    /*! @throw std::invalid_argument if max < min */
    template <typename T>
    std::enable_if_t<std::is_unsigned<T>::value, T> clamp(T value, const T& min, const T& max)
    {
        using T2 = typename std::make_signed<T>::type;
        return clamp<T2>(static_cast<T2>(value), static_cast<T2>(min), static_cast<T2>(max));
    }
    
    //! Wrap a number within a given range
    /*! @throw std::invalid_argument if max <= min */
    template <typename T>
    std::enable_if_t<!std::is_integral<T>::value || std::is_signed<T>::value, T> wrap(T value, const T& min, const T& max)
    {
        if (max <= min)
            throw std::invalid_argument("max <= min");
        
        while (max <= value)
            value -= (max - min);
        
        while (value < min)
            value += (max - min);
        
        return value;
    }

    //! Wrap a number within a given range
    /*! @throw std::invalid_argument if max <= min */
    template <typename T>
    std::enable_if_t<std::is_unsigned<T>::value, T> wrap(T value, const T& min, const T& max)
    {
        using T2 = typename std::make_signed<T>::type;
        return wrap<T2>(static_cast<T2>(value), static_cast<T2>(min), static_cast<T2>(max));
    }
    
    //! Wrap a number within a given range, where min = 0
    /*! @throw std::invalid_argument if max <= 0 */
    template <typename T, typename U>
    auto wrap(T value, const U& max)
    {
        if (max <= 0)
            throw std::invalid_argument("max <= 0");
        
        while (max <= value)
            value -= max;
        
        while (value < 0)
            value += max;
        
        return value;
    }
    
    //! Returns whether an integral is a power of two
    /*! @return false for x == 0 */
    constexpr bool isPowerOf2(size_t x)
    {
        return (x == 0) ? false : (x == 1 || (x & (x - 1)) == 0);
    }
    
    //! Find the smallest power of 2 bigger than a given number
    /*! @return 1 for x == 0 */
    constexpr size_t ceilToPowerOf2(size_t x)
    {
        size_t y = 1;
        while (y < x)
            y <<= 1;
        
        return y;
    }
}

#endif
