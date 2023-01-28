#include "Std/Std.h"
#include "Std/StdFile.h"
#include "Std/StdPosix.h"
#include "Std/StdVector.h"
#include <ccomplex>

using namespace Std;

void print(const char * format, std::vector<std::string> & values) {
	std::string buffer;
	char * p = (char*)format;
	size_t x = 0;
	while(*p != 0)
	{
		if(*p != '{') buffer += *p;
		else {
			buffer += values[x++];
			while(*p != 0) {
				if(*p == '}') break;
				p++;
			}		
		}
		p++;
	}
	std::cout << buffer << std::endl;
}

template<typename T>
void println(const Std::ValVector<T> & v)
{
	std::cout << "Vector[" << v.size() << "]=";
	for(size_t i = 0; i < v.size()-1; i++)
		std::cout << v[i] << ",";
	std::cout << v.back() << std::endl;
}
int main()
{
	srand(time(NULL));
	ValVector<float> v(1024);
	v.random(-1.0f,1.0f);
	println(v);
	ValVector<std::complex<float>> v2(1024);
	ValVector<std::complex<float>> r(512);

	RealFFT1D<float> fft(1024);
	
	fft.Forward(v,v2);
	r = v2.slice(0,511);
	println(r);
	fft.Backward(v2,v);
	println(v);

}
