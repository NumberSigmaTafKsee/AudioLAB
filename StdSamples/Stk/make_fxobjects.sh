swig -lua -c++ -ICAnalog fxobjects.i
gcc -fmax-errors=1 -ICAnalog -O2 -fPIC -march=native -mavx2 -shared -o fxobjects.so fxobjects_wrap.cxx -lstdc++ -lm -lluajit -lfftw3 -lfftw3f
