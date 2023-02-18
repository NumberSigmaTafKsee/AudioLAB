    swig -lua -c++ -I../include -I../include ATKChebyshev2Filter.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o ATKChebyshev2Filter.so ATKChebyshev2Filter_wrap.cxx -lstdc++ -lm -lluajit 
    