#include "audiodsp_samples_ptr.hpp"
#include "audiodsp_fvec.hpp"
#include "audiodsp_vtk.hpp"

int main()
{
    std::vector<float> v(10);
    AudioDSP::generate_sin(1000,44100,10,v.data());
    for(size_t i = 0; i < 10; i++) std::cout << v[i] << ",";
    std::cout << std::endl;
}