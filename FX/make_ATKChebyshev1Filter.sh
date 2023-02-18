    swig -lua -c++ -I../include -I../include ATKChebyshev1Filter.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o ATKChebyshev1Filter.so ATKChebyshev1Filter_wrap.cxx -lstdc++ -lm -lluajit 
    