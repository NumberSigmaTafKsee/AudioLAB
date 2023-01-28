%module audio_dc_block_filter
%{
#include "SoundObject.hpp"

#include <cassert>
#include <random>
#include <chrono>
#include <complex>
#include <cstdio>
#include <cstring>
#include <cstdlib>

#include "audio_dc_block_filter.hpp"
%}

%include "SoundObject.hpp"
%include "audio_dc_block_filter.hpp"