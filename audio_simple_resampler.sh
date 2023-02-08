swig -lua -c++ -Iinclude audio_simple_resampler.i
gcc -fmax-errors=1 -std=c++17 -Iinclude \
-O2 -fPIC -mavx2 -mfma -march=native -shared \
-o audio_simple_resampler.so audio_simple_resampler_wrap.cxx \
-lstdc++ -lm -lluajit
