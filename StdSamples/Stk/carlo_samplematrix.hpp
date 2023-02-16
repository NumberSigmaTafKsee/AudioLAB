
#include "Samoles.hpp"
// samples.i
// mkl.i
// casino.i

// eigen
// #include "SimpleEigen.hpp"
// eigen.i
// simpleeigen.i
// samplevector.i

// armadillo
// #include "SimpleArmadillo.hpp"
// simplearmadillo.i

// octave
// #include "Octopus.hpp"
// octopus.i

// Arrayfire
// #include "GPUFFT/af.h"
// arrayfire.i

// cufft
// #include "GPUFFT/cufft.h"

#include "SampleVector.h"

namespace Casino::eigen
{
    template<typename T> using sample_vector = Casino::eigen::SampleVector<T>;
    template<typename T> using sample_matrix = Casino::eigen::SampleMatrix<T>;

    template<typename T>
    using complex_vector = sample_vector<std::complex<T>>;

    template<typename T>
    using complex_matrix = sample_matrix<std::complex<T>>;
};

namespace Casino::Octave
{

};

/*
namespace Casino::kfr
{
    template<typename T> using sample_vector = Casino::kfr::univector<T>;
    template<typename T> using sample_matrix = Casino::kfr::unimatrix<T>;
};
namespace Casino::armadillo
{
    template<typename T> using sample_vector = Casino::armadillo::vec;
    template<typename T> using sample_matrix = Casino::armadillo::mat;
};
namespace Casino::vectorclass
{
    Vector2d;
    Vector4d;
    Vector4f;
    Vector8f;
    Vector8d;
    Vector16f;
};
*/

// MKL
// Eigen + Armadillo Interface
// Octopus/Octave
// ArrayFire