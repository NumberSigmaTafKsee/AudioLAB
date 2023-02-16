swig -lua -c++ -Iinclude/Eigen eigen-matrix.i
gcc -fmax-errors=1 -Iinclude/Eigen -O2 -fPIC -march=native -mavx2 -shared -o matrix.so eigen-matrix_wrap.cxx -lstdc++ -lm -lluajit
