swig -python -c++ -Iinclude/Eigen eigen-array2d.i
gcc -std=c++17 -fmax-errors=1 -I/usr/local/include/python3.9 -Iinclude/Eigen -O2 -fPIC -march=native -mavx2 -shared -o _array2d.so eigen-array2d_wrap.cxx -lstdc++ -lm -lpython3.9
