    swig -lua -c++ -I../include -I../include ATKRBJFilters.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o ATKRBJFilters.so ATKRBJFilters_wrap.cxx -lstdc++ -lm -lluajit 
    