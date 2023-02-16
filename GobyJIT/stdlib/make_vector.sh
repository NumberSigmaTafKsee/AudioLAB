swig -lua -c++ vector.i
gcc -fmax-errors=1  -O2 -fPIC -shared -o vector.so vector_wrap.cxx -lstdc++ -lm -lluajit
