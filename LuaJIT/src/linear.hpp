//
//  linear.hpp
//  Math
//
//  Copyright © 2015-2016 Dsperados (info@dsperados.com). All rights reserved.
//  Licensed under the BSD 3-clause license.
//

#ifndef DSPERADOS_MATH_LINEAR_HPP
#define DSPERADOS_MATH_LINEAR_HPP

#include <cmath>
#include <cstddef>
#include <vector>


namespace math
{
    //! Take the dot product of two containers
    template <class InputIterator1, class InputIterator2>
    auto dot(InputIterator1 begin1, std::size_t stride1, InputIterator2 begin2, std::size_t stride2, std::size_t size)
    {
        std::common_type_t<decltype(*begin1), decltype(*begin2)> out = {0};
        
        for (auto i = 0; i < size; ++i)
        {
            out += *begin1 * *begin2;
            begin1 += stride1;
            begin2 += stride2;
        }
        
        return out;
    }
    
    //! Take the dot product of two containers
    template <class InputIterator1, class InputIterator2>
    auto dot(InputIterator1 begin1, InputIterator2 begin2, std::size_t size)
    {
        return dot(begin1, 1, begin2, 1, size);
    }
}

#endif
