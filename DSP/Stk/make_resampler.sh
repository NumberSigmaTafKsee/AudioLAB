swig -lua -c++ -Iinclude src/Resampler.i
gcc -std=c++17 -O2 -fPIC -march=native -mavx2 -shared -o resampler.so src/Resampler_wrap.cxx -lstdc++ -lm -lsamplerate
