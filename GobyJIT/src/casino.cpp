#include "Casino.hpp"
#include <omp.h>

using namespace Casino;

#define PRINT(foo) std::cout << foo << std::endl;

void DisplayIPP()
{
		const IppLibraryVersion *vers;
	int cache_size;	
	int mhz;
	int bCache;
	int numThreads;
	
	Ipp64u cpu_clocks;
	IppCache * pCache;
	vers = ippsGetLibVersion();
	ippGetL2CacheSize(&cache_size);
	ippGetCacheParams(&pCache);	
	cpu_clocks = ippGetCpuClocks();
	ippGetCpuFreqMhz(&mhz);
	ippGetMaxCacheSizeB(&bCache);
	// doesn't seem to do anything
	ippSetNumThreads(8);
	omp_set_num_threads(8);
	ippGetNumThreads(&numThreads);
	

	PRINT("--------------------------------------")
	PRINT("IPP Library Version");
	PRINT("Minor Version: " << vers->minor)
	PRINT("Major Version: " << vers->major)
	PRINT("Major Build  : " << vers->majorBuild)
	PRINT("Build        : " << vers->build);
	PRINT("TargetCPU[0] : " << vers->targetCpu[0]);
	PRINT("TargetCPU[1] : " << vers->targetCpu[1]);
	PRINT("TargetCPU[2] : " << vers->targetCpu[2]);
	PRINT("TargetCPU[3] : " << vers->targetCpu[3]);
	PRINT("Name         : " << vers->Name);
	PRINT("Version      : " << vers->Version);
	PRINT("Build Date   : " << vers->BuildDate);
	PRINT("--------------------------------------")

	
	PRINT("L2 Cache Size\t: " << cache_size);
	PRINT("Cache Type   \t: " << pCache->type);
	PRINT("Cache Level  \t: " << pCache->level);
	PRINT("Cache Size   \t: " << pCache->size);
	PRINT("CPU Clocks   \t: " << cpu_clocks);
	PRINT("CPU Mhz      \t: " << mhz);
	PRINT("B-Cache Size \t: " << bCache);
	PRINT("Num-Threads  \t: " << numThreads);
}

void acorr_test()
{
	FloatArray a(10),b(10),c(19);
	for(size_t i = 0; i < 10; i++) { a[i] = i+1; b[i] = a[i]; }
	CrossCorrFloat xcorr(10,10,19,-9);
	xcorr.Process(a.array,b.array,c.array);
	c.print();
}
void xcorr_test()
{
	FloatArray a(10),b(10),c(19);
	for(size_t i = 0; i < 10; i++) { a[i] = i+1; b[i] = a[i]; }
	CrossCorrFloat xcorr(10,10,19,-9);
	xcorr.Process(a.array,b.array,c.array);
	c.print();
}
void conv_test()
{
	FloatArray a(32),b(32),c(63);
	for(size_t i = 0; i < 32; i++) { a[i] = i+1; b[i] = a[i]; }
	ConvolveFloat conv(32,32);
	conv.Process(a.array,b.array,c.array);
	c.print();
}
void cfft_test()
{
	Ipp32fc * a = ippsMalloc_32fc(128);
	Ipp32fc * b = ippsMalloc_32fc(128);
	for(size_t i = 0; i < 128; i++) { a[i].re = i+1; a[i].im = i+1, b[i].re = 0; b[i].im = 0; }
	CFFTFC fft(128);
	
	fft.Forward(a,b);
	for(size_t i = 0; i < 128; i++) { 
		std::cout << "(" << b[i].re << "," << b[i].im << "),";
 	}
	std::cout << std::endl;
	
	memset(a,0,sizeof(Ipp32fc)*128);
	fft.Inverse(b,a);
	for(size_t i = 0; i < 128; i++) { 
		std::cout << "(" << a[i].re << "," << a[i].im << "),";
 	}
	std::cout << std::endl;
	
}
void rfft_test() {
	Ipp32f * a = ippsMalloc_32f(128);
	Ipp32f * b = ippsMalloc_32f(128);
	for(size_t i = 0; i < 128; i++) { a[i] = i+1; b[i] = 0; }
	RFFTF fft(128);
	
	fft.Forward(a,b);
	for(size_t i = 0; i < 128; i++) { 
		std::cout << b[i] << ",";
 	}
	std::cout << std::endl;
	
	memset(a,0,sizeof(Ipp32f)*128);
	fft.Inverse(b,a);
	for(size_t i = 0; i < 128; i++) { 
		std::cout << a[i] << ",";
 	}
	std::cout << std::endl;		
}
void cdft_test()
{
	Ipp32fc * a = ippsMalloc_32fc(128);
	Ipp32fc * b = ippsMalloc_32fc(128);
	for(size_t i = 0; i < 128; i++) { a[i].re = i+1; a[i].im = i+1, b[i].re = 0; b[i].im = 0; }
	CDFTF fft(128);
	
	fft.Forward(a,b);
	for(size_t i = 0; i < 128; i++) { 
		std::cout << "(" << b[i].re << "," << b[i].im << "),";
 	}
	std::cout << std::endl;
	
	memset(a,0,sizeof(Ipp32fc)*128);
	fft.Inverse(b,a);
	for(size_t i = 0; i < 128; i++) { 
		std::cout << "(" << a[i].re << "," << a[i].im << "),";
 	}
	std::cout << std::endl;
	
}
void rdft_test() {
	Ipp32f * a = ippsMalloc_32f(128);
	Ipp32f * b = ippsMalloc_32f(128);
	for(size_t i = 0; i < 128; i++) { a[i] = i+1; b[i] = 0; }
	RDFTF fft(128);
	
	fft.Forward(a,b);
	for(size_t i = 0; i < 128; i++) { 
		std::cout << b[i] << ",";
 	}
	std::cout << std::endl;
	
	memset(a,0,sizeof(Ipp32f)*128);
	fft.Inverse(b,a);
	for(size_t i = 0; i < 128; i++) { 
		std::cout << a[i] << ",";
 	}
	std::cout << std::endl;		
}
int main()
{
	cdft_test();
	rdft_test();	
}