swig -lua -c++ -I/usr/local/include/python3.9 -Iinclude -Iinclude/SimpleEigen src/eigen.i
gcc -fmax-errors=1 -I/usr/local/include/python3.9 -Iinclude -Iinclude/SimpleEigen -O2 -fPIC -march=native -shared -o se.so src/eigen_wrap.cxx -lstdc++ -lpython3.9
