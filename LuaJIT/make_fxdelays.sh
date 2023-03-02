swig -lua -c++ -ICAnalog fxdelays.i
gcc -fmax-errors=1 -ICAnalog -O2 -fPIC -march=native -mavx2 -shared -o fxdelays.so fxdelays_wrap.cxx -lstdc++ -lm -lluajit -lfftw3 -lfftw3f
