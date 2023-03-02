    swig -lua -c++ -Iinclude uri_table.h.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o uri_table.h.so uri_table.h_wrap.cxx -lstdc++ -lm -lluajit    
    