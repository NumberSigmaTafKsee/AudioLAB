    swig -lua -c++ -Iinclude VCO.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o VCO.so VCO_wrap.cxx -lstdc++ -lm -lluajit    
    