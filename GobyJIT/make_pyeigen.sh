swig -python -c++ -I/usr/local/include/python3.9 -Iinclude -Iinclude/SimpleEigen eigen.i
gcc -fmax-errors=1 -I/usr/local/include/python3.9 -Iinclude -Iinclude/SimpleEigen -O2 -fPIC -march=native -shared -o_se.so eigen_wrap.cxx -lstdc++ -lpython3.9
