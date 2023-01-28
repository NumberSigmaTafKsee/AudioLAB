//
//  interleave.hpp
//  Math
//
//  Copyright © 2015-2016 Dsperados (info@dsperados.com). All rights reserved.
//  Licensed under the BSD 3-clause license.
//

#ifndef DSPERADOS_MATH_INTERLEAVE_HPP
#define DSPERADOS_MATH_INTERLEAVE_HPP

#include <cmath>
#include <complex>
#include <cstddef>
#include <iterator>
#include <vector>

namespace math
{
    //! Interleave two ranges
    template <class InputIt1, class InputIt2, class OutputIt>
    void interleave(InputIt1 inBegin1, InputIt1 inEnd, InputIt2 inBegin2, OutputIt outBegin)
    {
        for (auto it = inBegin1; it != inEnd; ++it)
        {
            *outBegin++ = *it;
            *outBegin++ = *inBegin2++;
        }
    }
    
    //! Deinterleave a range into two ranges
    template <class InputIt, class OutputIt1, class OutputIt2>
    void deinterleave(InputIt inBegin, InputIt inEnd, OutputIt1 outBegin1, OutputIt2 outBegin2)
    {
        for (auto it = inBegin; it != inEnd; ++it)
        {
            *outBegin1++ = *it++;
            *outBegin2++ = *it;
        }
    }
}

#endif
