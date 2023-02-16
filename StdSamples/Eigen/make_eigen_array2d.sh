swig -lua -c++ -Iinclude/Eigen eigen-array2d.i
gcc -std=c++17 -fmax-errors=1 -Iinclude/Eigen -O2 -fPIC -march=native -mavx2 -shared -o array2d.so eigen-array2d_wrap.cxx -lstdc++ -lm -lluajit
