#include "Samples.hpp"
#include "SamplesDSP.hpp"

using namespace AudioDSP;
using namespace AudioDSP::Samples;


int main()
{
    sample_vector<float> v(10);
    v.fill(1);
    v.print();
}