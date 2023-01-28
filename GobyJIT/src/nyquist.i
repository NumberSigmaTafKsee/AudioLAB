%module nyquist
%{
#include "nyquist.h"
%}
%include "std_string.i"
%include "std_vector.i"

%template(float_vector) std::vector<float>;

namespace nqr
{
    enum PCMFormat
    {
        PCM_U8,
        PCM_S8,
        PCM_16,
        PCM_24,
        PCM_32,
        PCM_64,
        PCM_FLT,
        PCM_DBL,
        PCM_END
    };
}
%include "nyquist.h"
