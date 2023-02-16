#include <vector>

typedef float DspFloatType;
#include "Plot.hpp"
#include "samples/sample.hpp"
#include "Spectrum/Spectrum.hpp"
#include "Spectrum/FFTProcessors.hpp"
int main()
{
    sample_vector<DspFloatType> v;
    Plot<DspFloatType> plotter;
    v = Spectrum::FFT::FFTWaveTableGenerator::sawtooth(440.0,44100.0f);
    plotter.plot_x(v.data(),v.size(),"Wave");
    char foo;
    std::cin >> foo;
}
