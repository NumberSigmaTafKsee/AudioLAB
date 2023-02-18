#pragma once

#include <Eigen/Core>
#include <unsupported/Eigen/Polynomials>


namespace iir_filters
{
	struct RootSolver
	{
		// will solve any polynomial
		// H(s) = Y(s) / X(s)
		// zeros = Y(s)
		// poles = X(s)
		// if no zeros, then zeros=1
	};
}
