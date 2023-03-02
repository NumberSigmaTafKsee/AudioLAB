    swig -lua -c++ -Iinclude VoltageControlledFilter.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o VoltageControlledFilter.so VoltageControlledFilter_wrap.cxx -lstdc++ -lm -lluajit    
    