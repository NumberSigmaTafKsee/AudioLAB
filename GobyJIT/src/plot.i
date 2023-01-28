%module plot
%{
#include "Std/StdPlot.h"
using namespace std;
%}
%include "stdint.i"
%include "std_string.i"
%include "std_vector.i"

%template(double_vector) std::vector<double>;

%include "Std/StdPlot.h"
