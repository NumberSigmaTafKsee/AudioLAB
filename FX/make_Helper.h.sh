    swig -lua -c++ -I../include -I../include Helper.h.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o Helper.h.so Helper.h_wrap.cxx -lstdc++ -lm -lluajit 
    