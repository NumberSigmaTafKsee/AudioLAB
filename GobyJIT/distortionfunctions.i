%module distortionfunctions
%{
typedef float DspFloatType;
#include "FX/DistortionFunctions.hpp"
#include "FX/ClipFunctions.hpp"
%}

%include "std_vector.i"

typedef float DspFloatType;

%template(float_vector) std::vector<float>;
%template(double_vector) std::vector<double>;

%include "FX/DistortionFunctions.hpp"
%include "FX/ClipFunctions.hpp"