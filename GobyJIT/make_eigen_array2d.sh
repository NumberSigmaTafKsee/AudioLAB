swig -lua -c++ -I/usr/local/include/eigen3 eigen-array2d.i
gcc -std=c++17 -fmax-errors=1 -I/usr/local/include/eigen3 -O2 -fPIC -march=native -mavx2 -shared -o eigen_array2d.so eigen-array2d_wrap.cxx -lstdc++ -lm -lluajit
