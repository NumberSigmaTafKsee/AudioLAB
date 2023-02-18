%module ccaubio 
%{

#include "ccaubio.h"
#include <cassert>
#include <vector>

using namespace std;

%}

%include "stdint.i"
%include "std_common.i"
%include "std_vector.i"

%include "types.h"

%template (float_vector) std::vector<float>;
%template (float_matrix) std::vector<std::vector<float>>;


%constant const char * sdesc_complex = "complex";
%constant const char * sdesc_energy = "energy";
%constant const char * sdesc_hfc = "hfc";
%constant const char * sdesc_kl = "kl";
%constant const char * sdesc_mkl = "mkl";
%constant const char * sdesc_phase = "phase";
%constant const char * sdesc_specdiff = "specdiff";
%constant const char * sdesc_wphase = "wphase";
%constant const char * sdesc_centroid = "centroid";
%constant const char * sdesc_decrease = "decrease";
%constant const char * sdesc_kurtosis = "kurtosis";
%constant const char * sdesc_rolloff = "rolloff";
%constant const char * sdesc_skewness = "skewness";
%constant const char * sdesc_slope = "slope";
%constant const char * sdesc_spread = "spread";

%constant const char * rectangle_window = "rectangle";
%constant const char * hamming_window = "hamming";
%constant const char * hanning_window = "hanning";
%constant const char * hanningz_window = "hanningz";
%constant const char * blackman_window = "blackman";
%constant const char * blackman_harris_window = "blackman_harris";
%constant const char * gaussian_window = "gaussian";
%constant const char * welch_window = "welch";
%constant const char * parzen_window = "parzen";
%constant const char * default_window = "default";

%constant const char* pitch_mcomb = "mcomb";
%constant const char* pitch_yinfast = "yinfast";
%constant const char* pitch_yinfft = "yinfft";
%constant const char* pitch_yin = "yin";
%constant const char* pitch_schmitt = "scmitt";
%constant const char* pitch_fcomb = "fcomb";
%constant const char* pitch_specacf = "specacf";
%constant const char* pitch_default = "default";

%include "ccaubio.h"

