    swig -lua -c++ -I../include -I../include Filters.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o Filters.so Filters_wrap.cxx -lstdc++ -lm -lluajit 
    