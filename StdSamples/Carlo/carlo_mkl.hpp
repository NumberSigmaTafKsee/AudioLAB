#pragma once

#include <iostream>
#include <ccomplex>
#include <complex>
#include <vector>
#include <algorithm>
#include <cstring>
#include <chrono>
#include <random>
#include <functional>

#include "StdNoise.hpp"

#include "cppmkl/cppmkl_allocator.h"
#include "cppmkl/cppmkl_vml.h"
#include "cppmkl/cppmkl_cblas.h"
#include "cppmkl/matrix.h"

template<typename T> using vector_base = std::vector<T, cppmkl::cppmkl_allocator<T> >;

#include "carlo_vector.hpp"
#include "carlo_complex_vector.hpp"
#include "carlo_matrix.hpp"
#include "carlo_complex_matrix.hpp"
#include "carlo_mklfft.hpp"
