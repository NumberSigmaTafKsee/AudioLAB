swig -lua -c++ -ICAnalog fxdistortion.i
gcc -fmax-errors=1 -ICAnalog -O2 -fPIC -march=native -mavx2 -shared -o fxdistortion.so fxdistortion_wrap.cxx -lstdc++ -lm -lluajit -lfftw3 -lfftw3f
