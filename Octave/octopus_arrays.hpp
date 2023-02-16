#pragma once

namespace Octopus
{
    struct OctopusArrayXf : public Array<float>
    {
        OctopusArrayXf(const Array<float> & a) : Array<float>(a) {}
        OctopusArrayXf(size_t i) : Array<float>(i) {}
    };
    struct OctopusArrayXd : public Array<double>
    {
        OctopusArrayXd(const Array<double> & a) : Array<double>(a) {}
        OctopusArrayXd(size_t i) : Array<double>(i) {}
    };
    struct OctopusArrayXcf : public Array<std::complex<float>>
    {
        OctopusArrayXcf(const Array<std::complex<float>> & a) : Array<std::complex<float>>(a) {}
        OctopusArrayXcf(size_t i) : : Array<std::complex<float>>(i) {}
    };
    struct OctopusArrayXcd : public Array<std::complex<double>>
    {
        OctopusArrayXcd(const Array<std::complex<double>> & a) : Array<std::complex<double>>(a) {}
        OctopusArrayXcd(size_t i) : : Array<std::complex<double>>(i) {}
    };
}    