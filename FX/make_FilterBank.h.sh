    swig -lua -c++ -I../include -I../include FilterBank.h.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o FilterBank.h.so FilterBank.h_wrap.cxx -lstdc++ -lm -lluajit 
    