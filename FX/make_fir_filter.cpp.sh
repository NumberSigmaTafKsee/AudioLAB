    swig -lua -c++ -I../include -I../include fir_filter.cpp.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o fir_filter.cpp.so fir_filter.cpp_wrap.cxx -lstdc++ -lm -lluajit 
    