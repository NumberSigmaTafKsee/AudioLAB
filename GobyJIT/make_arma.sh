swig -lua -c++ -Iinclude src/armadillo.i
gcc -Iinclude -fmax-errors=1 -O2 -fPIC -march=native -mavx2 -shared -o arma.so src/armadillo_wrap.cxx -lstdc++ -lm -lluajit -larmadillo
