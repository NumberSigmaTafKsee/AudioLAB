swig -lua -c++ valarray.i
gcc -fmax-errors=1  -O2 -fPIC -march=native -mavx2 -shared -o valarray.so valarray_wrap.cxx -lstdc++ -lm -lluajit
