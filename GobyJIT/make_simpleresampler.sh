swig -lua -c++ simpleresampler.i
gcc -std=c++17 -I. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o simpleresampler.so simpleresampler_wrap.cxx -lstdc++ -lm -lluajit
