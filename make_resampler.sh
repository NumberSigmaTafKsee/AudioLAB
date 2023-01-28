swig -lua -c++ -Iinclude Resampler.i
gcc -std=c++17 -O2 -fPIC -march=native -mavx2 -shared -o resampler.so Resampler_wrap.cxx -lstdc++ -lm -lsamplerate
