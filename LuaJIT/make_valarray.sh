swig -lua -c++ src/valarray.i
gcc -fmax-errors=1  -O2 -fPIC -march=native -mavx2 -shared -o valarray.so src/valarray_wrap.cxx -lstdc++ -lm -lluajit
