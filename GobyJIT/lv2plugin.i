%module lv2plugin
%{
#include "FX/LV2Plugin.hpp"
#include <vector>
%}

%include "std_vector.i"
%include "std_string.i"
%include "FX/LV2Plugin.hpp"

%template(float_vector) std::vector<float>;

