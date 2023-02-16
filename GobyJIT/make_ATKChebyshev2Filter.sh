    swig -lua -c++ -Iinclude ATKChebyshev2Filter.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o ATKChebyshev2Filter.so ATKChebyshev2Filter_wrap.cxx -lstdc++ -lm -lluajit    
    