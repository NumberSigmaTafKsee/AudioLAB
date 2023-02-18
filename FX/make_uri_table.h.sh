    swig -lua -c++ -I../include -I../include uri_table.h.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o uri_table.h.so uri_table.h_wrap.cxx -lstdc++ -lm -lluajit 
    