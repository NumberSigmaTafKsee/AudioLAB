swig -lua -c++ -Iinclude/Eigen eigen-array.i
gcc -fmax-errors=1 -Iinclude/Eigen -O2 -fPIC -march=native -mavx2 -shared -o eigen_array.so eigen-array_wrap.cxx -lstdc++ -lm -lluajit
