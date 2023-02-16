swig -python -c++ -Iinclude/Eigen array.i
gcc -fmax-errors=1 -I/usr/local/include/python3.9 -Iinclude/Eigen -O2 -fPIC -march=native -mavx2 -shared -o _array.so array_wrap.cxx -lstdc++ -lm -lpython3.9
