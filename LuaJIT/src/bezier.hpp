//
//  bezier.hpp
//  Math
//
//  Copyright © 2015-2016 Dsperados (info@dsperados.com). All rights reserved.
//  Licensed under the BSD 3-clause license.
//

#ifndef DSPERADOS_MATH_BEZIER_HPP
#define DSPERADOS_MATH_BEZIER_HPP

#include <dsperados/geometry/point.hpp>

namespace math
{
	//! The angle at any point of a quadratic bezier curve
	double quadraticBezierAngle(const double t, const geo::Point2f& control1, const geo::Point2f& control2, const geo::Point2f& point);
}

#endif
