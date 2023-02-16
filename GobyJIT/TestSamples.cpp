#include "Samples.hpp"
#include "SamplesDSP.hpp"
using namespace Casino;
using namespace Casino::Samples;

using Complex = std::complex<double>;

template<typename T>
void sequence(complex_vector<T> & a)
{
    for(size_t i = 0; i < a.size(); i++) a[i] = Complex(i,-1);
}
template<typename T>
void sequence(complex_matrix<T> & a)
{
    for(size_t i = 0; i < a.rows(); i++) 
        for(size_t j=0; j < a.cols(); j++)
            a(i,j) = Complex(i,-1);
}
template<typename T>
void sequence(sample_vector<T> & a)
{
    for(size_t i = 0; i < a.size(); i++) a[i] = i;
}
template<typename T>
void sequence(sample_matrix<T> & a)
{
    for(size_t i = 0; i < a.rows(); i++) 
        for(size_t j=0; j < a.cols(); j++)
            a(i,j) = i;
}

void sequence(std::vector<Complex> & a)
{
    for(size_t i = 0; i < a.size(); i++) a[i] = Complex(i,-1);
}

template<typename T>
std::ostream& operator << (std::ostream& o, std::vector<T> & v)
{
    for(auto x : v) o << x << ",";
    o << std::endl;
    return o;
}

template<typename T>
std::ostream& operator << (std::ostream& o, complex_vector<T> & v)
{
    for(auto x : v) o << x << ",";
    o << std::endl;
    return o;
}



int main()
{
    complex_matrix<float>  v(3,16);
    complex_matrix<float>  r(3,16);
    sequence(v);
    std::cout << v << std::endl;
    C2CF2D fft(3,16,FORWARD);
    C2CF2D ifft(3,16,BACKWARD);
    fft.set_input(v);
    fft.Execute();
    fft.normalize();
    fft.get_output(r);    
    //std::reverse(r.begin(),r.end());
    //std::copy(r.begin(),r.begin()+r.size()/2+1,t.begin());
    ifft.set_input(r);
    ifft.Execute();    
    ifft.get_output(v);    
    std::cout << v << std::endl;
}