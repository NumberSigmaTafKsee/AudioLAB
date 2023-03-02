    swig -lua -c++ -Iinclude fir_filter.cpp.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o fir_filter.cpp.so fir_filter.cpp_wrap.cxx -lstdc++ -lm -lluajit    
    