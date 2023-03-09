%module kfr_audiofft
%{
#include "KfrDSP/KfrDsp.hpp"
#include "KfrDSP/KfrAudio.hpp"
%}

%include "std_math.i"
%include "std_vector.i"
%include "std_list.i"
%include "std_map.i"
%include "lua_fnptr.i"

%template(float_vector) std::vector<float>;
%template(double_vector) std::vector<double>;

%include "samples/kfr_sample.hpp"
%include "samples/kfr_sample_dsp.hpp"

%template(sample_vector)  sample_vector<float>;
%template(complex_vector) complex_vector<float>;


%include "KfrDSP/KfrAudio.hpp"


%inline %{
float sampleRate=44100.0;
float invSampleRate=1.0/sampleRate;
Std::RandomMersenne noise;
%}