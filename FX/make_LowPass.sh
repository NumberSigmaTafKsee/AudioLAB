    swig -lua -c++ -I../include -I../include LowPass.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o LowPass.so LowPass_wrap.cxx -lstdc++ -lm -lluajit 
    