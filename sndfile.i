%module sndfile
%{

#include <vector>
#include <cassert>
#include "SndFile.hpp"

using namespace std;
%}

%include "stdint.i"
%include "std_vector.i"

%template (byte_vector) std::vector<unsigned char>;
%template (short_vector) std::vector<short>;
%template (int_vector) std::vector<int>;
%template (float_vector) std::vector<float>;
%template (double_vector) std::vector<double>;

%include "SndFile.hpp"

