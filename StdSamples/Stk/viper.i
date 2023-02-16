// todo: table <=> array, table <=> vector, table <=> matrix etc
// csv read/write


%module viper
%{
#include "viper.hpp"
using namespace Viper;
%}
%include "viper.hpp"


//%include "stdint.i"
typedef signed char int8_t;
typedef unsigned char uint8_t;
typedef short int16_t;
typedef unsigned short uint16_t;
typedef int int32_t;
typedef unsigned int uint32_t;
typedef long int int64_t;
typedef unsigned long int uint64_t;

%include "std_math.i"
%include "std_vector.i"
%include "viper_cublas.i"
%include "viper.hpp"

%template(float_vector) std::vector<float>;
%template(double_vector) std::vector<double>;
%template(complex_vector) std::vector<std::complex<float>>;
%template(double_complex_vector) std::vector<std::complex<double>>;



%inline %{ 

Cublas  _cublas;
Cublas *cublas = &_cublas;

void CreateCublas() { 
    cublas = new Cublas();
    assert(cublas != NULL);
}
void DeleteCublas() { 
    if(cublas) delete cublas;
    cublas = NULL;
}
void synchronize() {
    cudaDeviceSynchronize();
}
unsigned seed = 0;

void set_seed()
{
    typedef std::chrono::high_resolution_clock myclock;
    myclock::time_point beginning = myclock::now();
    myclock::duration d = myclock::now() - beginning;    
    seed = d.count();
}


template<typename T>
float cast_float(T val) { return (float)val; }

template<typename T>
double cast_double(T val) { return (double)val; }

template<typename T>
int8_t cast_int8(T val) { return (int8_t)val; }

template<typename T>
uint8_t cast_uint8(T val) { return (uint8_t)val; }

template<typename T>
int16_t cast_int16(T val) { return (int16_t)val; }

template<typename T>
uint16_t cast_uint16(T val) { return (uint16_t)val; }

template<typename T>
int32_t cast_int32(T val) { return (int32_t)val; }

template<typename T>
uint32_t cast_uint32(T val) { return (uint32_t)val; }

template<typename T>
int64_t cast_int64(T val) { return (int64_t)val; }

template<typename T>
uint64_t cast_uint64(T val) { return (uint64_t)val; }

std::vector<float> vector_range(int start, int end, int inc=1) {
    std::vector<float> r;    
    for(int i = start; i <= end; i+=inc) {
        r.push_back((float)i);
    }
    return r;
}

%}

// lua only has double
%template(cast_double_float) cast_double<float>;

%template(VectorXf) Viper::Vector<float>;
%template(VectorXd) Viper::Vector<double>;
%template(MatrixXf) Viper::Matrix<float>;
%template(MatrixXd) Viper::Matrix<double>;
