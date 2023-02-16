swig -pyton -c++ matrix.i
gcc -fmax-errors=1 -I/usr/local/include/python3.9 -I/include -O2 -fPIC -march=native -mavx2 -shared -o _matrix.so  matrix_wrap.cxx -lstdc++ -lm -lpython3.9
