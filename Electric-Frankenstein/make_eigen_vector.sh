swig -lua -c++ -Iinclude/Eigen eigen-vector.i
gcc -fmax-errors=1 -Iinclude/Eigen -O2 -fPIC -march=native -mavx2 -shared -o vector.so eigen-vector_wrap.cxx -lstdc++ -lm -lluajit
