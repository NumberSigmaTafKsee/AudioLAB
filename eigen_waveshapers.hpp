#pragma once

#include <Eigen/Core>

typedef Eigen::Matrix<DspFloatType,1,Eigen::Dynamic,Eigen::RowMajor> DspVector;

inline void amplify_vector(size_t n, DspFloatType * buffer, DspFloatType gain) {
	Eigen::Map<DspVector> v(buffer,n);
	v *= gain;	
}
