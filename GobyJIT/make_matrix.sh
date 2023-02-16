swig -lua -c++ matrix.i
gcc -fmax-errors=1 -I/usr/local/include/eigen3 -O2 -fPIC -march=native -mavx2 -shared -o matrix.so  matrix_wrap.cxx -lstdc++ -lm -lluajit
