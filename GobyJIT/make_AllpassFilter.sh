    swig -lua -c++ -Iinclude AllpassFilter.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o AllpassFilter.so AllpassFilter_wrap.cxx -lstdc++ -lm -lluajit    
    