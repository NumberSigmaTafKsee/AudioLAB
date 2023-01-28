    swig -lua -c++ -Iinclude VARCFilter.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o VARCFilter.so VARCFilter_wrap.cxx -lstdc++ -lm -lluajit    
    