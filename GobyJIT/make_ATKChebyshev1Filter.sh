    swig -lua -c++ -Iinclude ATKChebyshev1Filter.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o ATKChebyshev1Filter.so ATKChebyshev1Filter_wrap.cxx -lstdc++ -lm -lluajit    
    