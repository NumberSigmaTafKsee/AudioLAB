%module world
%{

#include "world/common.h"
#include "world/cheaptrick.h"
#include "world/constantnumbers.h"
#include "world/d4c.h"
#include "world/dio.h"
#include "world/fft.h"
#include "world/harvest.h"
#include "world/matlabfunctions.h"
#include "world/stonemask.h"
#include "world/synthesis.h"
#include "world/synthesisrealtime.h"
%}
#define WORLD_BEGIN_C_DECLS
#define WORLD_END_C_DECLS
%include "stdint.i"
%include "std_vector.i"
%include "std_string.i"
%include "std_complex.i"

%include "world/common.h"
%include "world/cheaptrick.h"
%include "world/constantnumbers.h"
%include "world/d4c.h"
%include "world/dio.h"
%include "world/fft.h"
%include "world/harvest.h"
%include "world/matlabfunctions.h"
%include "world/stonemask.h"
%include "world/synthesis.h"
%include "world/synthesisrealtime.h"
