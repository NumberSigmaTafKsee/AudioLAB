//
//  bezier.cpp
//  Math
//
//  Copyright © 2015-2016 Dsperados (info@dsperados.com). All rights reserved.
//  Licensed under the BSD 3-clause license.
//

#include <cmath>

#include "bezier.hpp"

using namespace geo;

namespace math
{
	double quadraticBezierAngle(const double t, const Point2f& control1, const Point2f& control2, const Point2f& point)
    {
        const auto invt = 1 - t;
        
        const auto c1 = -3 * invt * invt;
        const auto c2 = 3 * invt * invt * -6 * t * invt;
        const auto c3 = -3 * t * t + 6 * t * invt;
        const auto c4 = 3 * t * t;
        
        const auto dx = (c1 * 0) + (c2 * control1.x) + (c3 * control2.x) + (c4 * point.x);
        const auto dy = (c1 * 0) + (c2 * control1.y) + (c3 * control2.y) + (c4 * point.y);
        
        return atan2(dx, dy);
    }
}
