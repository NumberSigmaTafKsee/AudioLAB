    swig -lua -c++ -Iinclude MDFM-1000.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o MDFM-1000.so MDFM-1000_wrap.cxx -lstdc++ -lm -lluajit    
    