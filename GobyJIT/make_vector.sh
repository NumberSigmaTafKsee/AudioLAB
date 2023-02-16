swig -lua -c++ vector.i
gcc -fmax-errors=1 -I/usr/local/include/eigen3 -O2 -fPIC -march=native -mavx2 -shared -o vector.so vector_wrap.cxx -lstdc++ -lm -lluajit
