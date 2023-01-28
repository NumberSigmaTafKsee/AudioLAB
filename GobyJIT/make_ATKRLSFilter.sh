    swig -lua -c++ -Iinclude ATKRLSFilter.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o ATKRLSFilter.so ATKRLSFilter_wrap.cxx -lstdc++ -lm -lluajit    
    