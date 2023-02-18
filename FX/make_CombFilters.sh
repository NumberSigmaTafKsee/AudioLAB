    swig -lua -c++ -I../include -I../include CombFilters.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o CombFilters.so CombFilters_wrap.cxx -lstdc++ -lm -lluajit 
    