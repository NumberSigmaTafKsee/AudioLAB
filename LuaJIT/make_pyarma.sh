swig -python -c++ -Iinclude src/armadillo.i
gcc -I/usr/local/include/python3.9 -fmax-errors=1 -O2 -fPIC -march=native -mavx2 -shared -o _arma.so src/armadillo_wrap.cxx -lstdc++ -lm -lpython3.9 -larmadillo
