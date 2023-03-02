    swig -lua -c++ -Iinclude DCFilter.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o DCFilter.so DCFilter_wrap.cxx -lstdc++ -lm -lluajit    
    