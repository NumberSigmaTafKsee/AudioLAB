typedef float DspFloatType;
#include "StdSamples/stdsamples.hpp"
#include "StdSamples/stdsamples_fvec.hpp"
#include "StdSamples/stdsamples_vectorize.hpp"
#include "StdSamples/stdsamples_hirestimer.hpp"
#include <iomanip>

using namespace AudioDSP;

int main()
{
		sample_vector<float> a(256),b(256);
		fill(a,10.0f);
		HiResTimer timer;
		timer.Start();
		for(size_t i = 0; i < 100; i++)
			add(a.size(),a.data(),a.data(),b.data());
		std::cout << std::setprecision(2) << timer.Stop()/100.0 << std::endl;
		std::cout << sum(a.size(),a.data()) << std::endl;
}			
